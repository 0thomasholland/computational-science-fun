#!/usr/bin/env python3
"""
2D Visualization of Critical Cellular Automaton
"""

import subprocess
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")


def parse_output_file(filename, grid_size=(20, 20)):
    """Parse output file efficiently using NumPy."""
    nx, ny = grid_size

    # First pass: count frames
    with open(filename, "r") as f:
        num_frames = sum(1 for line in f if line.startswith("#D"))

    if num_frames == 0:
        return [], []

    # Pre-allocate storage
    grid_states = np.zeros((num_frames + 1, nx, ny), dtype=np.int32)
    iterations = []
    current_frame = 0

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#INIT"):
                current_frame = 0
                grid_states[0, :, :] = 0
            elif line.startswith("#D"):
                current_frame += 1
                grid_states[current_frame] = grid_states[current_frame - 1]
                try:
                    iteration = int(line[2:])
                    iterations.append(iteration)
                    if iteration % 1000 == 0:
                        print(f"Found iteration: {iteration}")
                except ValueError:
                    iterations.append(current_frame)
            else:
                try:
                    i, j, value = map(int, line.split(",")[:3])
                    i, j = i - 1, j - 1

                    if 0 <= i < nx and 0 <= j < ny:
                        grid_states[current_frame, i, j] = value
                except (ValueError, IndexError):
                    continue

    return grid_states[: current_frame + 1], iterations


def create_visualization(
    grid_states, output_file="animation.mp4", fps=30, grid_size=(40, 40)
):
    """Render frames directly to ffmpeg pipe — bypasses FuncAnimation overhead."""
    dpi = 100
    fig_w, fig_h = 8, 8
    px_w, px_h = int(fig_w * dpi), int(fig_h * dpi)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    max_value = np.max(grid_states)

    im = ax.imshow(
        grid_states[0],
        cmap="Greys",
        vmin=0,
        vmax=max_value,
        origin="lower",
        interpolation="nearest",
    )
    plt.colorbar(im, ax=ax, label="Grain Count")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    title = ax.set_title("", fontsize=12, fontweight="bold")

    total_grains = grid_states.sum(axis=(1, 2))

    # Draw static elements once so they're in the buffer
    fig.canvas.draw()

    ffmpeg_cmd = [
        "ffmpeg", "-y",
        "-f", "rawvideo",
        "-vcodec", "rawvideo",
        "-s", f"{px_w}x{px_h}",
        "-pix_fmt", "rgba",
        "-r", str(fps),
        "-i", "pipe:0",
        "-vcodec", "libx264",
        "-preset", "ultrafast",
        "-pix_fmt", "yuv420p",
        "-b:v", "2000k",
        output_file,
    ]

    print(f"Creating animation with {len(grid_states)} frames...")
    print(f"Saving animation to '{output_file}'...")

    try:
        proc = subprocess.Popen(
            ffmpeg_cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        print("Error: ffmpeg not found. Install with: brew install ffmpeg")
        sys.exit(1)

    try:
        for frame in range(len(grid_states)):
            im.set_array(grid_states[frame])
            title.set_text(
                f"Iteration {frame + 1} | Total Grains: {total_grains[frame]}"
            )
            fig.canvas.draw()
            proc.stdin.write(fig.canvas.buffer_rgba().tobytes())

            if (frame + 1) % 500 == 0:
                print(f"  {frame + 1}/{len(grid_states)} frames written")

        proc.stdin.close()
        proc.wait()
    except BrokenPipeError:
        print("Error: ffmpeg pipe broke unexpectedly")
        proc.kill()
        sys.exit(1)

    plt.close(fig)
    print(f"✓ Animation saved successfully to '{output_file}'")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize 2D critical cellular automaton"
    )
    parser.add_argument("--input", "-i", default="output.csv")
    parser.add_argument("--output", "-o", default="animation.mp4")
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument(
        "--grid-size", type=int, nargs=2, default=[40, 40], metavar=("NX", "NY")
    )

    args = parser.parse_args()
    grid_size = tuple(args.grid_size)

    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    print(f"Reading simulation data from '{args.input}'...")
    grid_states, iterations = parse_output_file(args.input, grid_size=grid_size)

    if len(grid_states) == 0:
        print("Error: No valid data found in input file.")
        sys.exit(1)

    print(f"Parsed {len(grid_states)} simulation frames")
    print(f"Grid size: {grid_size[0]} x {grid_size[1]}")

    create_visualization(
        grid_states, output_file=args.output, fps=args.fps, grid_size=grid_size
    )


if __name__ == "__main__":
    main()
