# replot.jl — Regenerate the figure for a saved run.
#
# Usage (from Julia REPL or terminal):
#   julia --project replot.jl 3        # replot run_003
#   julia --project replot.jl          # replot the most recent run

using Pkg; Pkg.activate(@__DIR__)
using DualLimitCycle
using CairoMakie

BASE_DIR = joinpath(@__DIR__, "runs")

# Determine run number: argument or latest
n = if !isempty(ARGS)
    parse(Int, ARGS[1])
else
    next_run_number(BASE_DIR) - 1
end

n < 1 && error("No runs found in $BASE_DIR")

println("Replotting run $(lpad(n, 3, '0')) ...")
data = load_run(n, BASE_DIR)

outpath = joinpath(BASE_DIR, "run_" * lpad(string(n), 3, '0'), "figure.png")
make_figure(data.t, data.x, data.y, data.params, data.jump_log, outpath)
println("Saved → ", outpath)
