import numpy as np
import matplotlib.pyplot as plt
import tqdm

# Geological constants
d = 1000.0  # Depth from surface to sill [meters]
rho = 2700.0  # Density of country rock [kg/m^3]
cp = 1000.0  # Constant specific heat capacity [J/(kg K)]
k = 2.5  # Thermal conductivity [W/(m K)]
kappa = k / (rho * cp)  # Thermal diffusivity [m^2/s]

# Boundary condition variables
Ts = 20.0  # Constant surface temperature [Celsius]
delta_T = 800.0  # Temperature scale of the magma [Celsius]
years = 365.25 * 24 * 3600  # Conversion factor: 1 year in seconds
tau = 5000.0 * years  # Timescale of magma recharge [seconds]


# Spatial grid (z = 0 is the sill, z = d is the Earth's surface)
N = 100  # Number of spatial nodes
dz = d / (N - 1)  # Distance between nodes [meters]
z = np.linspace(0, d, N)  # Array of depths

# Temporal grid (Ensuring numerical stability)
# For FTCS, dt must be <= dz^2 / (2 * kappa)
dt = (dz**2) / (2.5 * kappa)
total_time = 20000.0 * years  # Run simulation for 20,000 years
time_steps = int(total_time / dt)


def sill_temperature(t):
    # T_0(t) = T_s + delta_T * exp(-t/tau) * cos(t/4tau)
    return Ts + delta_T * np.exp(-t / tau) * np.cos(t / (4.0 * tau))


T = np.full((time_steps, N), Ts)


for n in tqdm.trange(0, time_steps - 1):
    current_time = n * dt

    # Step A: Update the interior rock nodes (Explicit FTCS scheme)
    # This vectorizes the loop for speed: T[n, 1 to N-1]
    T[n + 1, 1:-1] = T[n, 1:-1] + (kappa * dt / dz**2) * (
        T[n, 2:] - 2 * T[n, 1:-1] + T[n, :-2]
    )

    # Step B: Enforce Bottom Boundary (z = 0, the magmatic sill)
    T[n + 1, 0] = sill_temperature(current_time)

    # Step C: Enforce Top Boundary (z = d, Earth's surface)
    T[n + 1, -1] = Ts


print("Plotting")

fig, ax = plt.subplots(figsize=(10, 6))

# --- Hovmöller Diagram (Heatmap) ---
# Display the massive 2D array as a color map
# extent maps the pixel coordinates to our real physical units (meters and years)


extent = [0, total_time / years, 0, d]
im = ax.imshow(T.T, aspect="auto", cmap="inferno", extent=extent, origin="lower")

# Formatting
ax.set_title("Heatmap of Thermal Diffusion over Time")
ax.set_xlabel("Time (Years)")
ax.set_ylabel("Depth, z (meters)")

# Add the color scale to the side
fig.colorbar(im, ax=ax, label="Temperature (°Celsius)")


plt.tight_layout()
print("Saving figure")
plt.savefig("output_smol.png")
print("Saved")
