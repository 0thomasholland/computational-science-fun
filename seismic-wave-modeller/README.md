# Seismic Wave Forward Modeller

![Large simulation — 8 km × 8 km, 3 sources](output/20.gif)

![Minimal simulation — 3 km × 3 km, single source](output/10.gif)

A 2D elastic wave forward modeller written in Rust. Solves the velocity-stress formulation of the elastic wave equation on a staggered finite-difference grid, producing animated snapshots of seismic wavefield propagation.

## Features

- Staggered-grid finite differences (Virieux scheme) for P- and S-wave propagation
- Heterogeneous material properties: variable Vp, Vs, and density ρ
- Harmonic averaging of elastic moduli at interfaces; arithmetic averaging of density (Moczo et al. 2002)
- Ricker wavelet source with configurable frequency, position, and trigger delay
- Rayon-parallelised stress and velocity updates
- Energy monitoring (kinetic + potential) at each report step
- TOML config files; `dt` auto-computed from the 2D CFL condition if omitted
- Per-frame PNG output; stitch to mp4/gif with `make video`

## Quick Start

Requires Rust (stable), Cargo, and ffmpeg (for video).

```bash
make run                          # minimal config, 300×300 grid, 1.5 s
make run CONFIG=example_config.toml ITER=1   # 3-source, 500×500 grid, 2.0 s
make video ITER=0                 # stitch output/0/ frames → mp4 + gif
```

Or without Make:

```bash
cargo build --release
./target/release/seismic-wave-modeller --config minimal_config.toml --verbose
ffmpeg -framerate 30 -pattern_type glob -i "output/0/vmag_*.png" \
       -c:v libx264 -pix_fmt yuv420p output/0.mp4
```

## Makefile Targets

| Target | Description |
|--------|-------------|
| `make build` | Release build |
| `make run` | Run with `minimal_config.toml` (override: `CONFIG=`, `ITER=`) |
| `make run-example` | Run the 3-source example config |
| `make video` | Stitch frames → mp4 + gif (override: `ITER=`, `FIELD=`, `FPS=`) |
| `make clean` | Remove build artefacts |
| `make clean-output` | Delete all output frames and videos |

## Configuration

All parameters live in a TOML file. `minimal_config.toml` is a good starting point.

```toml
[grid]
nx = 300          # grid points in x
nz = 300          # grid points in z
dx = 10.0         # spacing (m)
dz = 10.0

[materials]
vp  = 3000.0      # P-wave velocity (m/s)
vs  = 1732.0      # S-wave velocity — must satisfy vs < vp
rho = 2700.0      # density (kg/m³)

[simulation]
total_time = 1.5  # seconds; dt auto-computed from CFL if omitted
cfl_safety = 0.5  # dt = cfl_safety * min(dx,dz) / (vp * √2)

[[sources]]
x         = 150   # grid index
z         = 150
amplitude = 25.0
frequency = 20.0  # Hz — Ricker wavelet dominant frequency
                  # wavelength = vp/f = 150 m → ~15 pts/wavelength at dx=10 m

[visualization]
field        = "vmag"   # vx | vz | vmag | divergence | curl | sigma_xx | sigma_zz | sigma_xz
video_length = 10.0     # animation duration (s)
fps          = 30.0
iteration    = 0        # controls output directory: output/<iteration>/
```

**Rule of thumb:** use at least 10 grid points per shortest wavelength. At `vp=3000 m/s`, `f=20 Hz`, wavelength = 150 m, so `dx ≤ 15 m`. `dt` is set automatically; the 2D CFL condition is `dt ≤ cfl_safety × min(dx,dz) / (vp × √2)`.

## Theory

Wavefields live on a staggered grid:

```
σxx,σzz ----vx---- σxx,σzz ----vx---- σxx,σzz
   |                 |                 |
   vz      σxz      vz      σxz      vz
   |                 |                 |
σxx,σzz ----vx---- σxx,σzz ----vx---- σxx,σzz
```

**Velocity updates** (momentum conservation):

$$\frac{\partial v_x}{\partial t} = \frac{1}{\rho}\!\left[\frac{\partial \sigma_{xx}}{\partial x} + \frac{\partial \sigma_{xz}}{\partial z}\right], \qquad \frac{\partial v_z}{\partial t} = \frac{1}{\rho}\!\left[\frac{\partial \sigma_{xz}}{\partial x} + \frac{\partial \sigma_{zz}}{\partial z}\right]$$

**Stress updates** (elastic constitutive relation):

$$\frac{\partial \sigma_{xx}}{\partial t} = (\lambda+2\mu)\frac{\partial v_x}{\partial x} + \lambda\frac{\partial v_z}{\partial z}, \qquad \frac{\partial \sigma_{zz}}{\partial t} = \lambda\frac{\partial v_x}{\partial x} + (\lambda+2\mu)\frac{\partial v_z}{\partial z}$$

$$\frac{\partial \sigma_{xz}}{\partial t} = \mu\!\left[\frac{\partial v_x}{\partial z} + \frac{\partial v_z}{\partial x}\right]$$

Time-stepping uses a leapfrog scheme: velocities at half-steps, stresses at full steps.

**Potential energy** uses the compliance-tensor form:

$$W = \frac{(\lambda+2\mu)(\sigma_{xx}^2+\sigma_{zz}^2) - 2\lambda\,\sigma_{xx}\sigma_{zz}}{8\mu(\lambda+\mu)} + \frac{\sigma_{xz}^2}{2\mu}$$

## Boundary Conditions

| Type | Status | Notes |
|------|--------|-------|
| Rigid (Dirichlet) | Implemented | `vx=vz=0` at all edges — reflects waves with wrong polarity relative to a free surface; no Rayleigh waves |
| Free surface | Planned | Zero traction (`σzz=σxz=0`) at top — required for realistic surface reflections and Rayleigh waves |
| Absorbing (PML/sponge) | Planned | Suppress artificial edge reflections |

## Todos

- [ ] Accept input velocity/density model files for heterogeneous media
- [ ] Free-surface boundary condition at top edge
- [ ] Absorbing boundary conditions (Clayton–Engquist or PML)
- [ ] Higher-order spatial stencils (4th-order FD)
