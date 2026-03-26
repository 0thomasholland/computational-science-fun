using JLD2
using Dates

"""
Scan base_dir for subdirectories matching run_NNN and return the next integer.
"""
function next_run_number(base_dir::AbstractString)
    if !isdir(base_dir)
        return 1
    end
    entries = readdir(base_dir)
    nums = Int[]
    for e in entries
        m = match(r"^run_(\d+)$", e)
        if m !== nothing
            push!(nums, parse(Int, m.captures[1]))
        end
    end
    return isempty(nums) ? 1 : maximum(nums) + 1
end

"""
Return path for run_NNN, creating the directory if needed.
"""
function run_dir(base_dir::AbstractString, n::Int)
    path = joinpath(base_dir, "run_" * lpad(string(n), 3, '0'))
    mkpath(path)
    return path
end

"""
Save a run to JLD2. Stores x, y, the signed distances d1/d2 to each cycle,
the jump log, parameters, and metadata.
"""
function save_run(t, x, y, jump_log, p, base_dir::AbstractString)
    n = next_run_number(base_dir)
    dir = run_dir(base_dir, n)
    data_path = joinpath(dir, "data.jld2")

    ρ1 = sqrt.((x .- p.cx1).^2 .+ (y .- p.cy1).^2)
    ρ2 = sqrt.((x .- p.cx2).^2 .+ (y .- p.cy2).^2)
    d1 = ρ1 .- p.R1
    d2 = ρ2 .- p.R2

    jldsave(data_path;
        t             = t,
        x             = x,
        y             = y,
        d1            = d1,
        d2            = d2,
        jump_log      = jump_log,
        params        = p,
        timestamp     = string(now()),
        julia_version = string(VERSION),
    )

    println("Saved run $n → $dir  ($(length(jump_log)) jumps detected)")
    return dir
end

"""
Load run_NNN/data.jld2 and return a named tuple matching what save_run wrote.
"""
function load_run(n::Int, base_dir::AbstractString)
    path = joinpath(base_dir, "run_" * lpad(string(n), 3, '0'), "data.jld2")
    d = load(path)
    return (
        t             = d["t"],
        x             = d["x"],
        y             = d["y"],
        d1            = d["d1"],
        d2            = d["d2"],
        jump_log      = d["jump_log"],
        params        = d["params"],
        timestamp     = d["timestamp"],
        julia_version = d["julia_version"],
    )
end
