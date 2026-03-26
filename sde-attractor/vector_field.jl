# vector_field.jl — phase-plane vector field visualisation
# Run with: julia --project vector_field.jl

using Pkg; Pkg.activate(@__DIR__)
using DualLimitCycle
using CairoMakie

p = SDEParams(
    cx1=0.0, cy1=0.0, R1=4.0,
    cx2=0.0, cy2=10.0, R2=4.0,
    barrier=1.0, ω1=1.0, ω2=1.0, σ=0.0,
)

# ── Grid ───────────────────────────────────────────────────────────────────────
# Extent: largest circle centre + radius, with margin
lim_x = max(abs(p.cx1) + p.R1, abs(p.cx2) + p.R2) * 1.6
lim_y_lo = min(p.cy1 - p.R1, p.cy2 - p.R2) - 0.5
lim_y_hi = max(p.cy1 + p.R1, p.cy2 + p.R2) + 0.5

n  = 50
xs = range(-lim_x, lim_x;   length=n)
ys = range(lim_y_lo, lim_y_hi; length=n)

# ── Evaluate drift on grid ─────────────────────────────────────────────────────
du   = zeros(2)
segs = Point2f[]
mags = Float64[]
dirs = Tuple{Float64,Float64}[]

for x in xs, y in ys
    drift!(du, [x, y], p, 0.0)
    mag = sqrt(du[1]^2 + du[2]^2)
    mag < 1e-10 && continue
    push!(mags, mag)
    push!(dirs, (du[1]/mag, du[2]/mag))
    push!(segs, Point2f(x, y))
    push!(segs, Point2f(x, y))   # endpoint filled below
end

# Log-scale arrow lengths to compress the wide dynamic range
cell_x    = (last(xs) - first(xs)) / (n - 1)
cell_y    = (last(ys) - first(ys)) / (n - 1)
cellsize  = min(cell_x, cell_y)
log_mags  = log1p.(mags)
max_log   = maximum(log_mags)
arrowscale = cellsize * 0.75 / max_log

for k in eachindex(mags)
    s        = segs[2k - 1]
    len      = log_mags[k] * arrowscale
    d        = dirs[k]
    segs[2k] = Point2f(s[1] + d[1]*len, s[2] + d[2]*len)
end

seg_colors = vec([c for c in log_mags/max_log for _ in 1:2])

# ── Basin boundary contour (|d1| = |d2|, i.e. d1²-d2² = 0) ───────────────────
# Z[i,j] must correspond to (xf[i], yf[j]) for Makie's contour!
nfine = 300
xf = range(-lim_x, lim_x;     length=nfine)
yf = range(lim_y_lo, lim_y_hi; length=nfine)
Z  = [ let ρ1 = sqrt((x - p.cx1)^2 + (y - p.cy1)^2),
           ρ2 = sqrt((x - p.cx2)^2 + (y - p.cy2)^2)
           (ρ1 - p.R1)^2 - (ρ2 - p.R2)^2
       end
       for x in xf, y in yf ]   # (nfine × nfine): Z[i,j] = value at (xf[i], yf[j])

# ── Figure ─────────────────────────────────────────────────────────────────────
θ = range(0, 2π; length=400)

fig = Figure(size=(600, 800))
ax  = Axis(fig[1, 1];
    xlabel="x", ylabel="y",
    title="Phase-plane vector field",
    aspect=DataAspect())

linesegments!(ax, segs;
    color=seg_colors, colormap=:plasma, colorrange=(0,1), linewidth=1.5)

contour!(ax, xf, yf, Z;
    levels=[0.0], color=:orange, linestyle=:dash, linewidth=1.8,
    label="basin boundary")

lines!(ax, p.cx1 .+ p.R1 .* cos.(θ), p.cy1 .+ p.R1 .* sin.(θ);
    color=:dodgerblue, linewidth=2.5, label="cycle 1 (stable)")
lines!(ax, p.cx2 .+ p.R2 .* cos.(θ), p.cy2 .+ p.R2 .* sin.(θ);
    color=:crimson, linewidth=2.5, label="cycle 2 (stable)")

axislegend(ax; position=:rt, framevisible=false)

outpath = joinpath(@__DIR__, "vector_field.png")
save(outpath, fig)
println("Saved → ", outpath)
