# Dual Limit Cycle SDE

A 2D stochastic dynamical system with two stable limit cycles. Noise drives rare
(Kramers-type) transitions between them.

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
A transition is detected when the trajectory crosses $d_1^2 - d_2^2 = 0$. For
concentric cycles this reduces to the saddle radius $r_s = (R_1 + R_2)/2$.

### Kramers estimate

The mean first-passage time between cycles scales as

$$\tau \approx \exp\!\left(\frac{\Delta V}{\sigma^2}\right)$$

where $\Delta V = \text{barrier}$. Use this to choose `tspan` before running.

---

## Implementation

```
src/
  potential.jl      — V(x,y), ∇V, _compute_scale, find_extrema, barrier_height, kramers_estimate
  system.jl         — SDEParams struct, drift!, diffusion!, make_problem, saddle_callback
  io.jl             — save_run / load_run (JLD2), run numbering
  plots.jl          — phase portrait + distance time series (CairoMakie)
  DualLimitCycle.jl — top-level module

run.jl              — entry point with # %% interactive blocks
runs/               — output: one numbered subdirectory per run
```

### Parameters (`SDEParams`)

| Field | Meaning |
|-------|---------|
| `cx1, cy1, R1` | centre and radius of cycle 1 |
| `cx2, cy2, R2` | centre and radius of cycle 2 |
| `barrier` | height of the potential saddle |
| `ω1, ω2` | angular velocities on each cycle |
| `σ` | noise amplitude |
| `tspan, u0, solver, dt` | integration settings |

`_scale` is computed automatically from the geometry and `barrier`.

### Running

```julia
julia --project
] instantiate        # first time only

include("run.jl")    # or execute # %% blocks individually in VS Code
```

Output for each run is saved under `runs/run_NNN/`:
- `data.jld2` — trajectory arrays (`x`, `y`, `d1`, `d2`), `jump_log`, parameters
- `figure.png` — phase portrait and signed-distance time series
