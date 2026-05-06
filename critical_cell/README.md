# Critical Cellular Automaton (Sandpile Model)

A Fortran implementation of the [Bak–Tang–Wiesenfeld sandpile model](https://en.wikipedia.org/wiki/Abelian_sandpile_model) — a canonical example of self-organised criticality (SOC). Grains are added one at a time to a grid; cells that exceed a critical threshold topple, redistributing energy to their neighbours and potentially triggering avalanches of all sizes.

Three variants are provided:

| Directory | Description |
|-----------|-------------|
| `2d/` | 2D grid, criticality threshold 4, small grids (~20×20) |
| `3d/` | 3D grid, criticality threshold 6, redistributes to 6 face-adjacent neighbours |
| `multithreading/` | 2D OpenMP-parallelised version with checkerboard domain decomposition for large grids (2000×2000) |

## Physics

- Each cell holds a non-negative integer grain count.
- A cell is **critical** when its count exceeds the threshold (4 in 2D, 6 in 3D).
- Critical cells **topple**: they lose the threshold amount and each neighbour gains 1.
- After each full relaxation (no critical cells remain), one grain is added to a random cell.
- The system self-organises to a critical state where avalanche sizes follow a power law — a signature of SOC.

## Building and Running

Requires a Fortran compiler (`gfortran`) and `make`.

```bash
cd 2d/   # or 3d/ or multithreading/
make rebuild
make run
```

The multithreading variant requires OpenMP support (`gfortran -fopenmp`).

## Visualisation

Output is written to `output.csv` as a sequence of grid states. A Python script renders a video:

```bash
python3 visualise.py --output animation.mp4 --fps 15 --gridsize X Y [Z] --input output.csv
```

Requires `numpy`, `matplotlib`, and `ffmpeg`.

## Next Steps

- [ ] Record avalanche sizes and inter-event times to analyse the magnitude–frequency power law.
- [ ] Add configurable grid size and iteration count via command-line arguments.
