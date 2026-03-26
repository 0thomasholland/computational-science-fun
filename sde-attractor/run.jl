# run.jl — Entry point. Execute blocks sequentially in VS Code (# %% markers)
# or run the whole file with: julia --project run.jl

# %% Imports
using Pkg; Pkg.activate(@__DIR__)
using DualLimitCycle
using CairoMakie
using StochasticDiffEq

# %% Parameters
# Each cycle is defined by its centre (cx, cy) and radius R.
# For concentric cycles set cx1=cy1=cx2=cy2=0.
# For offset cycles, move the centres independently.
p = SDEParams(
    cx1     = 0.0,   cy1 = 0.0,   R1 = 4.0,   # inner cycle
    cx2     = 0.0,   cy2 = 10.0,   R2 = 4.0,   # outer cycle
    barrier = 1.0,
    ω1      = 1.0,
    ω2      = 1.0,
    σ       = 1.5,
    tspan   = (0.0, 200.0),
    u0      = [1.0, 0.0],
    solver  = :SOSRA,
    dt      = 1e-5,
)

save_at_n = p.dt * 10

println("Kramers estimate of mean first passage time: ", round(kramers_estimate(p); sigdigits=3))

# %% Verify potential shape (scan along the inter-cycle path)
let
    ext = find_extrema(p)
    println("Potential extrema along inter-cycle path:")
    println("  minima : ", [(round(x; sigdigits=4), round(y; sigdigits=4)) for (x,y) in ext.minima])
    println("  maxima : ", [(round(x; sigdigits=4), round(y; sigdigits=4)) for (x,y) in ext.maxima])
    println("  barrier: ", round(barrier_height(p); sigdigits=4))

    # Plot V(x,0) along a 1D slice through the origin for a quick visual check
    xs = range(-p.R2 - 1, p.R2 + 1; length=400)
    vs = [V(x, 0.0, p) for x in xs]
    fig_v = Figure(size=(600, 350))
    ax = Axis(fig_v[1, 1]; xlabel="x  (y=0 slice)", ylabel="V(x,0)", title="Potential slice")
    lines!(ax, collect(xs), vs; color=:black)
    vlines!(ax, [-p.R1, p.R1]; color=:dodgerblue, linestyle=:dash, label="cycle 1")
    vlines!(ax, [-p.R2, p.R2]; color=:crimson,    linestyle=:dash, label="cycle 2")
    axislegend(ax; position=:rt, framevisible=false)
    save(joinpath(@__DIR__, "potential_check.png"), fig_v)
    println("Potential slice saved to potential_check.png")
end

# %% Solve
prob = make_problem(p)
cb, jump_log = saddle_callback(p)

solver_map = Dict(:EM => EM(), :SOSRA => SOSRA())
sol = solve(prob, solver_map[p.solver]; dt=p.dt, callback=cb, saveat=save_at_n, maxiters=Int(1e8))

# Extract to plain Float32 arrays immediately, then free the solution object.
# Float32 halves memory vs Float64 with no visible difference in plots.
t = Float32.(sol.t)
x = Float32.(sol[1, :])
y = Float32.(sol[2, :])
sol = nothing   # allow GC to reclaim the DiffEq solution object
GC.gc()

println("Integration complete. Steps saved: $(length(t)). Crossings: $(length(jump_log)).")

# %% Save and plot
BASE_DIR = joinpath(@__DIR__, "runs")

run_directory = save_run(t, x, y, jump_log, p, BASE_DIR)
fig = make_figure(t, x, y, p, jump_log, joinpath(run_directory, "figure.png"))

println("Figure saved to: ", joinpath(run_directory, "figure.png"))
