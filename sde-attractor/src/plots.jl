# plots.jl — Phase portrait and time-series visualisation.
# Operates only on arrays; no imports from system.jl or potential.jl.

using CairoMakie

"""
    plot_phase_portrait(ax, x, y, params; jump_log=nothing)

Trajectory coloured by time. Overlays dashed circles for both limit cycles
at their specified centres and radii.
"""
function plot_phase_portrait(ax, x, y, params; jump_log=nothing)
    n = length(x)
    t_norm = collect(range(0f0, 1f0; length=n))
    lines!(ax, x, y;  linewidth=0.5)

    θ = range(0, 2π; length=300)
    cos_θ, sin_θ = cos.(θ), sin.(θ)

    for (cx, cy, R, label, col) in [
            (params.cx1, params.cy1, params.R1, "cycle 1", :dodgerblue),
            (params.cx2, params.cy2, params.R2, "cycle 2", :crimson),
        ]
        lines!(ax, cx .+ R .* cos_θ, cy .+ R .* sin_θ;
               color=col, linestyle=:dash, linewidth=1.2, label=label)
    end

    n_jumps = isnothing(jump_log) ? 0 : length(jump_log)
    title_str = n_jumps > 0 ? "Phase portrait  [$n_jumps crossings]" : "Phase portrait"
    ax.aspect = DataAspect()
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.title  = title_str
    axislegend(ax; position=:rt, framevisible=false)
end

"""
    plot_timeseries(ax, t, d1, d2, params; jump_log=nothing)

Draw the reaction coordinate q = d1 - d2 over time.
  q < 0  →  closer to cycle 1
  q > 0  →  closer to cycle 2
  q = 0  →  basin boundary
"""
function plot_timeseries(ax, t, d1, d2, params; jump_log=nothing)
    q = d1 .- d2
    lines!(ax, t, q; color=:black, linewidth=0.6)
    hlines!(ax, [0.0]; color=:grey50, linestyle=:dash, linewidth=1.0)

    if !isnothing(jump_log)
        for jmp in jump_log
            vlines!(ax, [jmp.time]; color=:grey70, linewidth=0.6, linestyle=:dot)
        end
    end

    ax.xlabel = "t"
    ax.ylabel = "d₁ − d₂"
    ax.title  = "Reaction coordinate  (−: cycle 1,  +: cycle 2)"
end

"""
    plot_xy_timeseries(ax_x, ax_y, t, x, y; jump_log=nothing)

Draw x(t) and y(t) on separate axes with optional crossing markers.
"""
function plot_xy_timeseries(ax_x, ax_y, t, x, y; jump_log=nothing)
    lines!(ax_x, t, x; color=:black, linewidth=0.5)
    lines!(ax_y, t, y; color=:black, linewidth=0.5)

    if !isnothing(jump_log)
        for jmp in jump_log
            vlines!(ax_x, [jmp.time]; color=:grey70, linewidth=0.6, linestyle=:dot)
            vlines!(ax_y, [jmp.time]; color=:grey70, linewidth=0.6, linestyle=:dot)
        end
    end

    ax_x.ylabel = "x(t)"
    ax_y.xlabel = "t"
    ax_y.ylabel = "y(t)"
    hidexdecorations!(ax_x; grid=false)
end

"""
    make_figure(t, x, y, params, jump_log, outpath)

Four-panel figure (2×2):
  top-left:  phase portrait
  top-right: reaction coordinate d1−d2
  bot-left:  x(t)
  bot-right: y(t)
Saves to `outpath` and returns the Figure.
"""
function make_figure(t, x, y, params, jump_log, outpath)
    ρ1 = sqrt.((x .- params.cx1).^2 .+ (y .- params.cy1).^2)
    ρ2 = sqrt.((x .- params.cx2).^2 .+ (y .- params.cy2).^2)
    d1 = ρ1 .- params.R1
    d2 = ρ2 .- params.R2

    fig = Figure(size=(1400, 900))

    # Phase portrait spans all three right-side rows on the left column
    ax_phase = Axis(fig[1:3, 1])
    ax_react = Axis(fig[1, 2])
    ax_x     = Axis(fig[2, 2])
    ax_y     = Axis(fig[3, 2])

    plot_phase_portrait(ax_phase, x, y, params; jump_log=jump_log)
    plot_timeseries(ax_react, t, d1, d2, params; jump_log=jump_log)
    plot_xy_timeseries(ax_x, ax_y, t, x, y; jump_log=jump_log)

    # Link the x-axes of all three right-column time series
    linkxaxes!(ax_react, ax_x, ax_y)
    hidexdecorations!(ax_react; grid=false)

    save(outpath, fig)
    return fig
end
