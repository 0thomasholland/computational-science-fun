# N-Body Gravitational Solver

A Julia package for simulating N-body gravitational systems — early-stage work in progress.

## What it does (planned)

- Simulate N point masses interacting under Newtonian gravity in 3D.
- Track positions, velocities, and masses over time.
- Target use case: planetary systems, star clusters, and other gravitational N-body problems.

## Current state

The core data structure is defined (`State`: positions, velocities, masses, time). Integration and force-calculation routines are not yet implemented.

## Structure

```
src/nBodySolver.jl   — module with State struct (positions 3×N, velocities 3×N, masses N)
Project.toml         — Julia project manifest
```

## Running

```bash
julia --project
```

Requires Julia 1.x. No external dependencies beyond the standard library yet.
