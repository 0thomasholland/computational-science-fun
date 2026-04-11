# Dual Limit Cycle SDE

A 2D stochastic dynamical system with two stable limit cycles at arbitrary positions
and radii. Noise drives rare (Kramers-type) transitions between them. Includes both
a stochastic trajectory integrator and a Fokker-Planck PDE solver that evolves the
probability density forward in time.

---

## Maths

### Potential

Each limit cycle is a circle with centre $(c_x, c_y)$ and radius $R$. The signed
distance from a point $(x, y)$ to cycle $i$ is

$$d_i = \sqrt{(x - c_{x,i})^2 + (y - c_{y,i})^2} - R_i$$

which is zero on the circle, negative inside, and positive outside.

The potential is a double-well in these distance coordinates:

$$V(x, y) = \lambda \, d_1^2 \, d_2^2$$

where $\lambda$ is chosen so that $\max V = \text{barrier}$. The minima (both at
$V = 0$) lie along the two circles; the saddle lies on the path between them.

**Why the gradient approach rather than coupled linear equations?**

- The cycles can be placed anywhere in the plane with independent centres and radii.
  A gradient potential $V(d_1, d_2)$ encodes the geometry exactly without any
  special-casing — the force always points toward the nearest cycle regardless of
  orientation or separation distance.
- The barrier height $\Delta V$ maps directly onto the Kramers escape rate
  $\tau \sim e^{\Delta V / \sigma^2}$, giving a single physically meaningful
  parameter to tune transition timescales. Coupled linear equations have no
  equivalent transparent link between coefficients and first-passage times.

The gradient is computed analytically:

$$\nabla V = 2\lambda \bigl( d_1 d_2^2 \,\nabla d_1 + d_1^2 d_2 \,\nabla d_2 \bigr), \qquad \nabla d_i = \frac{(x,y) - (c_{x,i}, c_{y,i})}{\rho_i}$$

### SDE

The system evolves as an Itô SDE in Cartesian coordinates:

$$d\mathbf{u} = f(\mathbf{u})\,dt + \sigma\,d\mathbf{W}$$

The drift has two components:

**Potential force** — pulls the trajectory onto the nearest cycle:

$$f_\text{pot} = -\nabla V(x, y)$$

**Angular force** — drives rotation around each cycle's centre. The contribution
from each cycle is weighted by proximity so that cycle $i$ dominates when
$d_i \approx 0$:

$$w_i = \frac{d_j^2}{d_1^2 + d_2^2}, \qquad f_\text{ang} = w_1 \,\omega_1 \begin{pmatrix} -(y - c_{y,1}) \\ x - c_{x,1} \end{pmatrix} + w_2 \,\omega_2 \begin{pmatrix} -(y - c_{y,2}) \\ x - c_{x,2} \end{pmatrix}$$

The diffusion is isotropic: $\sigma \mathbf{I}$.

### Transition detection

The equidistant surface $|d_1| = |d_2|$ separates the two basins of attraction.
A `ContinuousCallback` fires when the trajectory crosses $d_1^2 - d_2^2 = 0$,
logging the time and direction. For concentric cycles this reduces to the saddle
radius $r_s = (R_1 + R_2)/2$.

### Kramers estimate

The mean first-passage time between cycles scales as

$$\tau \approx \exp\!\left(\frac{\Delta V}{\sigma^2}\right)$$

where $\Delta V = \text{barrier}$. Use this to set `tspan` before running.

### Fokker-Planck equation

The evolution of the probability density $\rho(x, y, t)$ under the same drift and
diffusion satisfies

$$\frac{\partial \rho}{\partial t} = -\nabla \cdot (f \rho) + \frac{\sigma^2}{2} \nabla^2 \rho$$

This is solved numerically on a fixed Cartesian grid using conservative donor-cell
upwind advection and central-difference diffusion, time-stepped with explicit RK4.
Density snapshots at user-specified times are rendered asynchronously as the solver
advances.

---

## Implementation

```
src/
  potential.jl      — V(x,y), ∇V, _compute_scale, find_extrema, barrier_height, kramers_estimate
  system.jl         — SDEParams struct, drift!, diffusion!, make_problem, saddle_callback
  io.jl             — save_run / load_run (JLD2), run numbering
  plots.jl          — phase portrait, reaction coordinate, and x/y time series (CairoMakie)
  DualLimitCycle.jl — top-level module

run.jl              — SDE entry point; # %% blocks for interactive use in VS Code
sde_push.jl         — Fokker-Planck PDE solver; snapshots saved to fp_runs/
runs/               — SDE output: one numbered subdirectory per run
fp_runs/            — FP output:  one numbered subdirectory per run
```

### Parameters (`SDEParams`)

| Field | Meaning |
|-------|---------|
| `cx1, cy1, R1` | centre and radius of cycle 1 |
| `cx2, cy2, R2` | centre and radius of cycle 2 |
| `barrier` | height of the potential saddle (sets $\lambda$ automatically) |
| `ω1, ω2` | angular velocities on each cycle |
| `σ` | noise amplitude |
| `tspan, u0, solver, dt` | integration settings |

`_scale` ($\lambda$) is computed automatically from the geometry and `barrier`.

### Running the SDE solver

```julia
julia --project
] instantiate        # first time only

include("run.jl")    # or execute # %% blocks individually in VS Code
```

Output for each run is saved under `runs/run_NNN/`:
- `data.jld2` — trajectory arrays (`x`, `y`, `d1`, `d2`), `jump_log`, parameters
- `figure.png` — four-panel figure: phase portrait, reaction coordinate $d_1 - d_2$, $x(t)$, $y(t)$

### Running the Fokker-Planck solver

```julia
julia --project -t auto sde_push.jl
```

Edit `t_ends` in `sde_push.jl` to set snapshot times, e.g.:

```julia
t_ends = collect(range(0.01, 5.0, length=20))
```

Output for each run is saved under `fp_runs/run_NNN/` with one PNG per snapshot,
named by time (`0.01.png`, `0.5.png`, etc.). Plot tasks are dispatched
asynchronously so rendering overlaps with integration. The `-t auto` flag is
required for both the parallel RK4 inner loop and the async plot tasks.
