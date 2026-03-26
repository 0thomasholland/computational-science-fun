# system.jl — SDE definition in Cartesian coordinates, parameter struct, and callbacks.

using DifferentialEquations

"""
    SDEParams

Parameter struct for the dual-limit-cycle system.
Each limit cycle is specified by its centre (cx, cy) and radius R.
The scale factor `_scale` is auto-computed from the geometry and barrier height.
"""
struct SDEParams
    cx1::Float64          # cycle 1 centre x
    cy1::Float64          # cycle 1 centre y
    R1::Float64           # cycle 1 radius
    cx2::Float64          # cycle 2 centre x
    cy2::Float64          # cycle 2 centre y
    R2::Float64           # cycle 2 radius
    barrier::Float64      # saddle height (max of V)
    ω1::Float64           # angular velocity on cycle 1
    ω2::Float64           # angular velocity on cycle 2
    σ::Float64            # noise amplitude
    tspan::Tuple{Float64,Float64}
    u0::Vector{Float64}
    solver::Symbol
    dt::Float64
    _scale::Float64       # auto-computed — do not set manually
end

"""
Keyword constructor. Computes `_scale` automatically from the geometry.
"""
function SDEParams(;
    cx1     = 0.0,
    cy1     = 0.0,
    R1      = 1.0,
    cx2     = 0.0,
    cy2     = 0.0,
    R2      = 3.0,
    barrier = 1.0,
    ω1      = 1.0,
    ω2      = 1.0,
    σ       = 0.3,
    tspan   = (0.0, 500.0),
    u0      = [1.0, 0.0],
    solver  = :SOSRA,
    dt      = 1e-3,
)
    scale = _compute_scale(cx1, cy1, R1, cx2, cy2, R2, barrier)
    return SDEParams(cx1, cy1, R1, cx2, cy2, R2, barrier, ω1, ω2, σ, tspan, u0, solver, dt, scale)
end

# ── Drift ──────────────────────────────────────────────────────────────────────

"""
    drift!(du, u, p, t)

In-place drift for the 2D SDE in Cartesian coordinates.

- Potential force: -∇V(x,y)
- Angular force: tangential velocity around each cycle centre, blended so
  that cycle i dominates when the system sits on that cycle (d_i ≈ 0).
  The blend weights are w_i ∝ d_j² (high near cycle i, low near cycle j).
"""
function drift!(du, u, p, t)
    x, y = u[1], u[2]

    dVdx, dVdy = ∇V(x, y, p)

    ρ1 = sqrt((x - p.cx1)^2 + (y - p.cy1)^2)
    ρ2 = sqrt((x - p.cx2)^2 + (y - p.cy2)^2)
    d1 = ρ1 - p.R1
    d2 = ρ2 - p.R2

    w1    = d2^2
    w2    = d1^2
    w_sum = w1 + w2

    if w_sum < 1e-16
        ωx, ωy = 0.0, 0.0
    else
        w1 /= w_sum;  w2 /= w_sum
        # Tangential velocity around centre i: ωi × (-(y-cyi), (x-cxi))
        ωx = w1 * p.ω1 * (-(y - p.cy1)) + w2 * p.ω2 * (-(y - p.cy2))
        ωy = w1 * p.ω1 * (  x - p.cx1)  + w2 * p.ω2 * (  x - p.cx2)
    end

    du[1] = -dVdx + ωx
    du[2] = -dVdy + ωy
end

# ── Diffusion ──────────────────────────────────────────────────────────────────

"""
    diffusion!(du, u, p, t)

Isotropic diagonal noise with amplitude `p.σ`.
"""
function diffusion!(du, u, p, t)
    du[1] = p.σ
    du[2] = p.σ
end

# ── Problem constructor ────────────────────────────────────────────────────────

"""
    make_problem(p)

Construct and return an `SDEProblem` ready to solve.
"""
function make_problem(p::SDEParams)
    return SDEProblem(drift!, diffusion!, p.u0, p.tspan, p)
end

# ── Saddle callback ────────────────────────────────────────────────────────────

"""
    saddle_callback(p)

Return a `ContinuousCallback` that fires when the system crosses the equidistant
surface between the two cycles: `d1² - d2² = 0` (i.e. |d1| = |d2|).

This surface is negative when the system is closer to cycle 1 and positive when
closer to cycle 2 — a natural generalisation of the saddle-radius condition for
offset cycles.  For concentric cycles it reduces to exactly the saddle radius.

Each crossing appends `(time, direction)` to the returned `jump_log`:
  +1 → cycle 1 → cycle 2
  -1 → cycle 2 → cycle 1
"""
function saddle_callback(p::SDEParams)
    jump_log = NamedTuple{(:time, :direction), Tuple{Float64,Int}}[]

    function condition(u, t, integrator)
        x, y = u[1], u[2]
        ρ1 = sqrt((x - p.cx1)^2 + (y - p.cy1)^2)
        ρ2 = sqrt((x - p.cx2)^2 + (y - p.cy2)^2)
        return (ρ1 - p.R1)^2 - (ρ2 - p.R2)^2
    end

    function affect!(integrator)
        x, y = integrator.u[1], integrator.u[2]
        du_tmp = zeros(2)
        drift!(du_tmp, integrator.u, integrator.p, integrator.t)

        ρ1 = sqrt((x - p.cx1)^2 + (y - p.cy1)^2)
        ρ2 = sqrt((x - p.cx2)^2 + (y - p.cy2)^2)
        d1 = ρ1 - p.R1;  d2 = ρ2 - p.R2

        inv_ρ1 = ρ1 < 1e-8 ? 0.0 : 1.0 / ρ1
        inv_ρ2 = ρ2 < 1e-8 ? 0.0 : 1.0 / ρ2

        dd1dt = ((x - p.cx1) * du_tmp[1] + (y - p.cy1) * du_tmp[2]) * inv_ρ1
        dd2dt = ((x - p.cx2) * du_tmp[1] + (y - p.cy2) * du_tmp[2]) * inv_ρ2

        # d/dt (d1² - d2²) = 2 d1 ḋ1 - 2 d2 ḋ2
        dcond_dt = 2 * d1 * dd1dt - 2 * d2 * dd2dt
        dir = dcond_dt ≥ 0 ? 1 : -1
        push!(jump_log, (time=integrator.t, direction=dir))
    end

    cb = ContinuousCallback(condition, affect!)
    return cb, jump_log
end
