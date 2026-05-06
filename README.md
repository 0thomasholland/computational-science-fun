# Computational Science Playground

Personal collection of computational physics and mathematics projects. Each is self-contained with its own build system and dependencies.

## Projects

### [Seismic Wave Modeller](seismic-wave-modeller/README.md) — Rust
2D elastic wave forward modeller using the Virieux staggered-grid finite-difference scheme. Simulates P- and S-wave propagation through heterogeneous media with configurable material properties. Outputs animated GIFs of the wavefield.

![Seismic wave model](seismic-wave-modeller/output/1.gif)

---

### [Critical Cellular Automaton](critical_cell/README.md) — Fortran
Sandpile model (Bak–Tang–Wiesenfeld) demonstrating self-organised criticality. Grains are added one at a time; cells exceeding a threshold topple and trigger avalanches. Three variants: 2D, 3D, and an OpenMP-parallelised 2D version using checkerboard domain decomposition for large (2000×2000) grids.

---

### [Dual Limit Cycle SDE](sde-attractor/README.md) — Julia
Stochastic dynamical system with two stable limit cycles connected by noise-driven (Kramers) transitions. Includes an Itô SDE trajectory integrator with automatic transition detection and a Fokker-Planck PDE solver that evolves the probability density forward in time. Barrier height controls transition rates via Kramers' formula.

---

### [N-Body Gravitational Solver](nBodySolver/README.md) — Julia
Early-stage N-body gravitational simulator. Core data structures (positions, velocities, masses in 3D) are defined; force calculation and integration are in progress.

---

## Languages & Tools

| Project | Language | Key dependencies |
|---------|----------|-----------------|
| seismic-wave-modeller | Rust | ndarray, rayon, plotters, clap |
| critical_cell | Fortran | gfortran, OpenMP |
| sde-attractor | Julia | DifferentialEquations.jl, CairoMakie, JLD2 |
| nBodySolver | Julia | (none yet) |
