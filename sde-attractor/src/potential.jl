# potential.jl — 2D Cartesian potential V(x,y) and related diagnostics.
# No dependencies on other project files.

"""
    _dists(x, y, p) → (d1, d2, ρ1, ρ2)

Signed distances from (x,y) to each limit cycle.
    d_i = ρ_i - R_i  (positive outside the circle, negative inside)
"""
@inline function _dists(x, y, p)
    ρ1 = sqrt((x - p.cx1)^2 + (y - p.cy1)^2)
    ρ2 = sqrt((x - p.cx2)^2 + (y - p.cy2)^2)
    return ρ1 - p.R1, ρ2 - p.R2, ρ1, ρ2
end

"""
    V(x, y, p)

2D potential with zero-minima along both limit cycles.
    V = _scale * d1² * d2²
where d_i is the signed distance from (x,y) to cycle i.
"""
function V(x, y, p)
    d1, d2, _, _ = _dists(x, y, p)
    return p._scale * d1^2 * d2^2
end

"""
    ∇V(x, y, p) → (dVdx, dVdy)

Analytical gradient of V(x,y).
    ∇(d1² d2²) = 2 d1 d2² ∇d1 + 2 d1² d2 ∇d2
"""
function ∇V(x, y, p)
    d1, d2, ρ1, ρ2 = _dists(x, y, p)

    inv_ρ1 = ρ1 < 1e-8 ? 0.0 : 1.0 / ρ1
    inv_ρ2 = ρ2 < 1e-8 ? 0.0 : 1.0 / ρ2

    dd1dx = (x - p.cx1) * inv_ρ1;  dd1dy = (y - p.cy1) * inv_ρ1
    dd2dx = (x - p.cx2) * inv_ρ2;  dd2dy = (y - p.cy2) * inv_ρ2

    c1 = 2.0 * p._scale * d1 * d2^2
    c2 = 2.0 * p._scale * d1^2 * d2
    return c1 * dd1dx + c2 * dd2dx, c1 * dd1dy + c2 * dd2dy
end

"""
    _compute_scale(cx1, cy1, R1, cx2, cy2, R2, barrier) → Float64

Numerically finds the peak of the unscaled potential d1²·d2² along the direct
path between the two cycles, then returns `barrier / peak` so that max(V) = barrier.

Called once during SDEParams construction — not on the hot path.
"""
function _compute_scale(cx1, cy1, R1, cx2, cy2, R2, barrier)
    D = sqrt((cx2 - cx1)^2 + (cy2 - cy1)^2)
    N = 500

    if D < 1e-8
        # Concentric: scan radially between the two circles
        r_a, r_b = min(R1, R2), max(R1, R2)
        max_unscaled = maximum(
            (r_a + (r_b - r_a) * i / (N - 1) - R1)^2 *
            (r_a + (r_b - r_a) * i / (N - 1) - R2)^2
            for i in 0:(N-1)
        )
    else
        nx, ny = (cx2 - cx1) / D, (cy2 - cy1) / D
        # Endpoints: nearest points on each circle facing the other
        x1e, y1e = cx1 + R1 * nx, cy1 + R1 * ny
        x2e, y2e = cx2 - R2 * nx, cy2 - R2 * ny
        max_unscaled = maximum(
            let s  = i / (N - 1),
                px = x1e + s * (x2e - x1e),
                py = y1e + s * (y2e - y1e),
                ρ1 = sqrt((px - cx1)^2 + (py - cy1)^2),
                ρ2 = sqrt((px - cx2)^2 + (py - cy2)^2)
                (ρ1 - R1)^2 * (ρ2 - R2)^2
            end
            for i in 0:(N-1)
        )
    end

    return barrier / max(max_unscaled, 1e-12)
end

"""
    find_extrema(p; n_samples=500)

Scan V along the direct inter-cycle path and return labelled extrema.
Returns `(minima=..., maxima=...)` as vectors of `(x, y)` tuples.
The endpoints (on the cycles themselves) have V=0 and are not included.
"""
function find_extrema(p; n_samples=500)
    D = sqrt((p.cx2 - p.cx1)^2 + (p.cy2 - p.cy1)^2)

    if D < 1e-8
        r_a, r_b = min(p.R1, p.R2), max(p.R1, p.R2)
        pts = [(r_a + (r_b - r_a) * i / (n_samples - 1), 0.0) for i in 0:(n_samples-1)]
    else
        nx, ny = (p.cx2 - p.cx1) / D, (p.cy2 - p.cy1) / D
        x1e, y1e = p.cx1 + p.R1 * nx, p.cy1 + p.R1 * ny
        x2e, y2e = p.cx2 - p.R2 * nx, p.cy2 - p.R2 * ny
        pts = [(x1e + (x2e - x1e) * i / (n_samples - 1),
                y1e + (y2e - y1e) * i / (n_samples - 1)) for i in 0:(n_samples-1)]
    end

    vs = [V(pt[1], pt[2], p) for pt in pts]

    # The cycle positions are the minima (V=0 at both endpoints by construction).
    # The saddle is the global maximum along the interior of the path.
    # Using argmax avoids the strict-inequality tie-breaking problem that arises
    # when the potential is symmetric and two adjacent samples have equal values.
    minima = [pts[1], pts[end]]
    maxima = Tuple{Float64,Float64}[]
    if length(pts) >= 3
        i_max = argmax(@view vs[2:end-1]) + 1   # +1: interior slice is offset by 1
        push!(maxima, pts[i_max])
    end
    return (minima=minima, maxima=maxima)
end

"""
    barrier_height(p)

Return V at the saddle between the two cycles (should equal p.barrier by construction).
"""
function barrier_height(p)
    ext = find_extrema(p)
    isempty(ext.maxima) && error("No saddle found — check geometry / barrier parameter")
    x_s, y_s = ext.maxima[1]
    return V(x_s, y_s, p)
end

"""
    kramers_estimate(p)

Rough Arrhenius estimate of mean first-passage time: τ ≈ exp(ΔV / σ²).
"""
function kramers_estimate(p)
    return exp(barrier_height(p) / p.σ^2)
end
