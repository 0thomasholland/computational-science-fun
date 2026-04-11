using Pkg; Pkg.activate(@__DIR__)
using FFMPEG_jll

# ── Config ────────────────────────────────────────────────────────────────────

RUN = "run_002"          # ← change this
FPS = 10

# ── Run ───────────────────────────────────────────────────────────────────────

run_dir = joinpath(@__DIR__, "fp_runs", RUN)
isdir(run_dir) || error("Directory not found: $run_dir")

pngs = filter(readdir(run_dir, join=true)) do f
    endswith(f, ".png") && tryparse(Float64, splitext(basename(f))[1]) !== nothing
end
isempty(pngs) && error("No snapshot PNGs found in $run_dir")
sort!(pngs, by = f -> parse(Float64, splitext(basename(f))[1]))

println("Found $(length(pngs)) frames in $run_dir")

list_path = joinpath(run_dir, "frames.txt")
open(list_path, "w") do io
    for p in pngs
        println(io, "file '$(abspath(p))'")
        println(io, "duration $(1/FPS)")
    end
end

out_path = joinpath(run_dir, "density.mp4")
run(`$(FFMPEG_jll.ffmpeg()) -y -f concat -safe 0 -i $list_path
    -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"
    -c:v libx264 -pix_fmt yuv420p $out_path`)
rm(list_path)
println("Video saved → ", out_path)
