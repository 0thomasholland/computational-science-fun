# sde_push.jl — Fokker-Planck PDE solver (forward push of probability density).
# Uses the same SDEParams and potential as run.jl.

using Pkg;
Pkg.activate(@__DIR__);
using DualLimitCycle
using CairoMakie
using FFMPEG_jll



# t_ends = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 5.0]
t_ends = collect(0.0:0.1:12.0)


# ── Parameters (mirror run.jl) ────────────────────────────────────────────────

p = SDEParams(
    cx1=0.0, cy1=0.0, R1=4.0,   # inner cycle
    cx2=0.0, cy2=10.0, R2=4.0,   # outer cycle
    barrier=1.0,
    ω1=1.0,
    ω2=1.0,
    σ=1.5,
    tspan=(0.0, 200.0),
    u0=[1.0, 0.0],
    solver=:SOSRA,
    dt=1e-5,
)



# ── Grid ──────────────────────────────────────────────────────────────────────

Base.@kwdef struct FPGrid
    xmin::Float64 = -6.0
    xmax::Float64 = 6.0
    ymin::Float64 = -6.0
    ymax::Float64 = 16.0
    nx::Int = 120
    ny::Int = 220
    # Initial condition: Gaussian centred at (μx, μy)
    μx::Float64 = 1.0
    μy::Float64 = 0.0
    σ0x::Float64 = 0.5
    σ0y::Float64 = 0.5
end

g = FPGrid()

# ── Drift evaluated on a 2D grid ──────────────────────────────────────────────

function build_grid(g::FPGrid, p::SDEParams)
    xs = range(g.xmin, g.xmax, length=g.nx)
    ys = range(g.ymin, g.ymax, length=g.ny)
    dx = step(xs)
    dy = step(ys)

    Fx = zeros(g.nx, g.ny)
    Fy = zeros(g.nx, g.ny)

    Threads.@threads for j in 1:g.ny
        du = zeros(2)
        u = zeros(2)   # thread-local buffers
        for i in 1:g.nx
            u[1] = xs[i]
            u[2] = ys[j]
            drift!(du, u, p, 0.0)
            Fx[i, j] = du[1]
            Fy[i, j] = du[2]
        end
    end

    return xs, ys, dx, dy, Fx, Fy
end

# ── Initial PDF ───────────────────────────────────────────────────────────────

function initial_pdf(xs, ys, g::FPGrid)
    ρ = zeros(g.nx, g.ny)
    dx = step(xs)
    dy = step(ys)
    for j in 1:g.ny, i in 1:g.nx
        ρ[i, j] = exp(
            -0.5 * ((xs[i] - g.μx) / g.σ0x)^2
            -
            0.5 * ((ys[j] - g.μy) / g.σ0y)^2
        )
    end
    ρ ./= (sum(ρ) * dx * dy)
    return ρ
end

# ── Fokker-Planck RHS ─────────────────────────────────────────────────────────

function fp_rhs!(dρ, ρ, Fx, Fy, dx, dy, σ²_half, nx, ny)
    dρ .= 0.0
    inv_dx = 1.0 / dx
    inv_dy = 1.0 / dy
    inv_dx2 = inv_dx^2
    inv_dy2 = inv_dy^2
    # Conservative upwind: discretise -∂(Fx·ρ)/∂x - ∂(Fy·ρ)/∂y via donor-cell fluxes.
    # Internal fluxes cancel between neighbours → mass exactly conserved in the interior.
    Threads.@threads for j in 2:ny-1
        @inbounds for i in 2:nx-1
            flux_xR = Fx[i, j] >= 0.0 ? Fx[i, j] * ρ[i, j] : Fx[i+1, j] * ρ[i+1, j]
            flux_xL = Fx[i-1, j] >= 0.0 ? Fx[i-1, j] * ρ[i-1, j] : Fx[i, j] * ρ[i, j]
            flux_yT = Fy[i, j] >= 0.0 ? Fy[i, j] * ρ[i, j] : Fy[i, j+1] * ρ[i, j+1]
            flux_yB = Fy[i, j-1] >= 0.0 ? Fy[i, j-1] * ρ[i, j-1] : Fy[i, j] * ρ[i, j]

            diff = σ²_half * (
                (ρ[i+1, j] - 2ρ[i, j] + ρ[i-1, j]) * inv_dx2 +
                (ρ[i, j+1] - 2ρ[i, j] + ρ[i, j-1]) * inv_dy2
            )

            dρ[i, j] = -(flux_xR - flux_xL) * inv_dx - (flux_yT - flux_yB) * inv_dy + diff
        end
    end
end

function cfl_dt(Fx, Fy, dx, dy, σ; safety=0.4)
    max_fx = maximum(abs.(Fx))
    max_fy = maximum(abs.(Fy))
    dt_adv = safety / (max_fx / dx + max_fy / dy)
    dt_diff = safety * min(dx, dy)^2 / σ^2
    return min(dt_adv, dt_diff)
end

function rk4_step!(ρ, Fx, Fy, dx, dy, σ²_half, nx, ny, dt, k1, k2, k3, k4, ρ_tmp)
    fp_rhs!(k1, ρ, Fx, Fy, dx, dy, σ²_half, nx, ny)

    @. ρ_tmp = ρ + (dt / 2) * k1
    fp_rhs!(k2, ρ_tmp, Fx, Fy, dx, dy, σ²_half, nx, ny)

    @. ρ_tmp = ρ + (dt / 2) * k2
    fp_rhs!(k3, ρ_tmp, Fx, Fy, dx, dy, σ²_half, nx, ny)

    @. ρ_tmp = ρ + dt * k3
    fp_rhs!(k4, ρ_tmp, Fx, Fy, dx, dy, σ²_half, nx, ny)

    @. ρ += (dt / 6) * (k1 + 2k2 + 2k3 + k4)
end

# ── Solver ────────────────────────────────────────────────────────────────────

"""
    solve_fp(checkpoints, g, p; on_snapshot)

Integrate the FP equation from t=0 to `maximum(checkpoints)` in a single run.
`checkpoints` must be sorted ascending.
At each checkpoint, fires `on_snapshot(xs, ys, t, ρ_copy)` before continuing —
the copy is independent so on_snapshot can be dispatched asynchronously.
"""
function solve_fp(checkpoints::AbstractVector{Float64}, g::FPGrid, p::SDEParams;
    on_snapshot=nothing)
    xs, ys, dx, dy, Fx, Fy = build_grid(g, p)
    ρ = initial_pdf(xs, ys, g)
    σ²_half = p.σ^2 / 2
    dt = cfl_dt(Fx, Fy, dx, dy, p.σ)

    k1 = similar(ρ)
    k2 = similar(ρ)
    k3 = similar(ρ)
    k4 = similar(ρ)
    ρ_tmp = similar(ρ)

    t = 0.0
    steps = 0

    for target in checkpoints
        while t < target
            dt_actual = min(dt, target - t)
            rk4_step!(ρ, Fx, Fy, dx, dy, σ²_half, g.nx, g.ny, dt_actual, k1, k2, k3, k4, ρ_tmp)
            t += dt_actual
            steps += 1
        end
        mass = sum(ρ) * dx * dy
        println("t = $target  |  steps so far: $steps  |  dt ≈ $(round(dt, sigdigits=3))  |  mass = $(round(mass, sigdigits=5))")
        on_snapshot !== nothing && on_snapshot(xs, ys, target, copy(ρ))
    end
end

# ── Plot ──────────────────────────────────────────────────────────────────────

function plot_fp(xs, ys, ρ, g::FPGrid, p::SDEParams, t::Float64)
    fig = Figure(size=(650, 900))
    ax = Axis(fig[1, 1],
        xlabel="x",
        ylabel="y",
        title="Fokker-Planck density  t = $t",
        aspect=DataAspect(),
    )

    hm = heatmap!(ax, xs, ys, ρ,
        colormap=:inferno,
        interpolate=true,
    )
    Colorbar(fig[1, 2], hm, label="p(x, y, t)")

    θ = range(0, 2π, length=300)
    lines!(ax, p.cx1 .+ p.R1 .* cos.(θ), p.cy1 .+ p.R1 .* sin.(θ),
        color=:white, linewidth=1.5, linestyle=:dash, label="cycle 1")
    lines!(ax, p.cx2 .+ p.R2 .* cos.(θ), p.cy2 .+ p.R2 .* sin.(θ),
        color=:cyan, linewidth=1.5, linestyle=:dash, label="cycle 2")

    scatter!(ax, [p.cx1, p.cx2], [p.cy1, p.cy2],
        color=:white, marker=:cross, markersize=14, label="centres")

    axislegend(ax, position=:lt, framecolor=:grey)
    return fig
end

# ── I/O helpers ───────────────────────────────────────────────────────────────

"""
Format a Float64 as a compact filename stem: 1.0 → "1.0", 0.01 → "0.01", etc.
"""
snapshot_name(t::Float64) = string(t)

"""
Create the next fp_runs/run_NNN directory under base_dir and return its path.
Uses the same numbering scheme as run.jl (run_001, run_002, …).
"""
function fp_run_dir(base_dir::AbstractString)
    n = next_run_number(base_dir)
    dir = run_dir(base_dir, n)
    return dir
end

# ── Video ─────────────────────────────────────────────────────────────────────

function make_video(run_dir::AbstractString; fps::Int=10)
    pngs = filter(readdir(run_dir, join=true)) do f
        endswith(f, ".png") && tryparse(Float64, splitext(basename(f))[1]) !== nothing
    end
    isempty(pngs) && error("No PNGs found in $run_dir")
    sort!(pngs, by=f -> parse(Float64, splitext(basename(f))[1]))

    list_path = joinpath(run_dir, "frames.txt")
    open(list_path, "w") do io
        for p in pngs
            println(io, "file '$(abspath(p))'")
            println(io, "duration $(1/fps)")
        end
    end

    out_path = joinpath(run_dir, "density.mp4")
    run(`$(FFMPEG_jll.ffmpeg()) -y -f concat -safe 0 -i $list_path
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"
        -c:v libx264 -pix_fmt yuv420p $out_path`)
    rm(list_path)
    println("Video saved → ", out_path)
    return out_path
end

# ── Mass partition plot ────────────────────────────────────────────────────────

function plot_mass_partition(ts, m1, m2, run_dir::AbstractString, p::SDEParams)
    y_divide = (p.cy1 + p.cy2) / 2
    fig = Figure(size=(700, 400))
    ax = Axis(fig[1, 1],
        xlabel="t",
        ylabel="mass fraction",
        title="Basin mass partitioning  (divide: y = $y_divide)",
    )
    lines!(ax, ts, m1, color=:white, linewidth=2, label="cycle 1 (y < $y_divide)")
    lines!(ax, ts, m2, color=:cyan, linewidth=2, label="cycle 2 (y ≥ $y_divide)")
    axislegend(ax, position=:rt, framecolor=:grey)
    path = joinpath(run_dir, "mass_partition.png")
    save(path, fig)
    println("Mass partition plot saved → ", path)
end

# ── Run ───────────────────────────────────────────────────────────────────────

BASE_DIR = joinpath(@__DIR__, "fp_runs")
run_directory = fp_run_dir(BASE_DIR)
println("Saving snapshots to: ", run_directory)

plot_tasks = Task[]
mass_times = Float64[]
mass1_arr = Float64[]   # mass in cycle-1 basin (y < divide)
mass2_arr = Float64[]   # mass in cycle-2 basin (y ≥ divide)

solve_fp(sort(t_ends), g, p;
    on_snapshot=(xs, ys, t, ρ) -> begin
        # Mass partition (compute on calling thread before spawning plot)
        y_divide = (p.cy1 + p.cy2) / 2
        dx = step(xs)
        dy = step(ys)
        m1 = sum(ρ[i, j] for i in 1:length(xs), j in 1:length(ys) if ys[j] < y_divide) * dx * dy
        m2 = sum(ρ[i, j] for i in 1:length(xs), j in 1:length(ys) if ys[j] >= y_divide) * dx * dy
        push!(mass_times, t)
        push!(mass1_arr, m1)
        push!(mass2_arr, m2)

        push!(plot_tasks,
            Threads.@spawn begin
                fig = plot_fp(xs, ys, ρ, g, p, t)
                path = joinpath(run_directory, snapshot_name(t) * ".png")
                save(path, fig)
                println("  saved → ", path)
            end
        )
    end
)

# Block until all in-flight plot tasks have written their files
foreach(wait, plot_tasks)

# Mass partition plot
plot_mass_partition(mass_times, mass1_arr, mass2_arr, run_directory, p)

# Assemble video
make_video(run_directory; fps=10)
