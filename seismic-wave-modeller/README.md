# Seismic Wave Forward Modeller

![An example output](output/3.gif)

A 2D elastic wave forward modeller written in Rust. Solves the velocity-stress formulation of the elastic wave equation on a staggered finite-difference grid, producing animated snapshots of seismic wavefield propagation through heterogeneous media.

## Features

- Staggered-grid finite differences (Virieux scheme) for P- and S-wave propagation
- Heterogeneous material properties: variable P-wave velocity (Vp), S-wave velocity (Vs), and density (ρ)
- Harmonic averaging at grid interfaces for stable material discontinuities
- TOML config files for easy parameter control
- Rayon-parallelised inner loops
- GIF output of wavefield snapshots

## Theory

Wavefields are stored on a staggered grid:

```
σxx,σzz ----vx---- σxx,σzz ----vx---- σxx,σzz
   |                 |                 |
   vz      σxz      vz      σxz      vz
   |                 |                 |
σxx,σzz ----vx---- σxx,σzz ----vx---- σxx,σzz
```

**Velocity updates** (from momentum conservation):

$$\frac{\partial v_x}{\partial t} = \frac{1}{\rho} \left[ \frac{\partial \sigma_{xx}}{\partial x} + \frac{\partial \sigma_{xz}}{\partial z} \right], \qquad \frac{\partial v_z}{\partial t} = \frac{1}{\rho} \left[ \frac{\partial \sigma_{xz}}{\partial x} + \frac{\partial \sigma_{zz}}{\partial z} \right]$$

**Stress updates** (from the elastic constitutive relation):

$$\frac{\partial \sigma_{xx}}{\partial t} = (\lambda + 2\mu) \frac{\partial v_x}{\partial x} + \lambda \frac{\partial v_z}{\partial z}, \qquad \frac{\partial \sigma_{zz}}{\partial t} = \lambda \frac{\partial v_x}{\partial x} + (\lambda + 2\mu) \frac{\partial v_z}{\partial z}$$

$$\frac{\partial \sigma_{xz}}{\partial t} = \mu \left[ \frac{\partial v_x}{\partial z} + \frac{\partial v_z}{\partial x} \right]$$

Time-stepping uses a leap-frog scheme: velocities at half-steps, stresses at full steps.

## Building and Running

Requires Rust (stable) and Cargo.

```bash
cargo build --release
cargo run --release -- --config example_config.toml
```

See `example_config.toml` and `minimal_config.toml` for parameter reference.

## Boundary Conditions

| Type | Status | Description |
|------|--------|-------------|
| Rigid | Implemented | Zero velocity at edges — causes reflections |
| Damping (sponge) | Planned | Absorbing boundary layer ~20–50 cells thick |
| Free surface | Planned | Zero normal stress at top surface |

## Todos

- [ ] Accept input velocity/density model files
- [ ] Implement absorbing (sponge) and free-surface boundary conditions
- [ ] Command-line interface for common parameters
- [ ] Loop over maximum Vp for stability CFL check
