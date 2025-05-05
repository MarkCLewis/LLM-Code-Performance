import numpy as np
import time
import math
import numba # Import Numba

# Gravitational constant (can be set to 1 for simulation units)
G = 1.0

# --- Physics and Simulation Functions ---

# Use Numba for JIT compilation and parallel execution
# fastmath=True enables potentially unsafe floating point optimizations for speed
# parallel=True enables auto-parallelization (used with prange)
@numba.njit(parallel=True, fastmath=True)
def calculate_acceleration_parallel(masses, positions, G, epsilon):
    """
    Calculates the acceleration of each body due to gravitational interaction.
    Uses direct summation (O(N^2)) parallelized with Numba.

    Args:
        masses (np.ndarray): 1D array of masses (N).
        positions (np.ndarray): 2D array of positions (N, 3).
        G (float): Gravitational constant.
        epsilon (float): Softening length to avoid singularities.

    Returns:
        np.ndarray: 2D array of accelerations (N, 3).
    """
    n_bodies = masses.shape[0]
    # Initialize acceleration array (important for Numba)
    accelerations = np.zeros_like(positions)
    epsilon_sq = epsilon**2

    # Loop over each particle 'i' in parallel
    # prange signals Numba to parallelize this loop
    for i in numba.prange(n_bodies):
        acc_i = np.zeros(3) # Accumulator for particle i's acceleration
        pos_i = positions[i]
        # Inner loop: calculate force on particle 'i' from all other particles 'j'
        for j in range(n_bodies):
            if i == j:
                continue # Skip self-interaction

            # Vector difference
            r_ij = positions[j] - pos_i

            # Squared distance + softening
            dist_sq = np.sum(r_ij**2) + epsilon_sq

            # Inverse cube distance (avoiding sqrt for efficiency here)
            inv_dist_cubed = dist_sq**(-1.5)

            # Acceleration contribution from j onto i: G * m_j * r_ij / |r_ij|^3
            # Note: No need to divide by masses[i] here, as G * masses[j] is the relevant factor from F=ma
            acc_i += G * masses[j] * inv_dist_cubed * r_ij

        accelerations[i] = acc_i # Store final acceleration for particle i

    return accelerations

@numba.njit(parallel=True, fastmath=True)
def calculate_potential_energy_parallel(masses, positions, G):
    """
    Calculates total potential energy, parallelized with Numba.
    Sums over unique pairs (i < j).

    Args:
        masses (np.ndarray): 1D array of masses (N).
        positions (np.ndarray): 2D array of positions (N, 3).
        G (float): Gravitational constant.

    Returns:
        float: Total potential energy.
    """
    n_bodies = masses.shape[0]
    total_pe = 0.0 # Numba handles reductions for simple sums in parallel loops

    # Parallelize the outer loop over i
    for i in numba.prange(n_bodies):
        pe_i = 0.0 # Accumulator for pairs starting with i
        pos_i = positions[i]
        mass_i = masses[i]
        # Inner loop calculates interactions with j > i
        for j in range(i + 1, n_bodies):
            r_ij = positions[j] - pos_i
            dist = np.sqrt(np.sum(r_ij**2))
            # Avoid division by zero (softening isn't typically used in PE definition)
            if dist > 1e-12: # Use a very small threshold
                 pe_i -= G * mass_i * masses[j] / dist
        total_pe += pe_i # Safely add partial sum to total (Numba reduction)

    return total_pe

def calculate_energy_parallel(masses, positions, velocities, G):
    """
    Calculates the total kinetic, potential, and total energy of the system,
    using the parallelized potential energy calculation.

    Args:
        masses (np.ndarray): 1D array of masses (N).
        positions (np.ndarray): 2D array of positions (N, 3).
        velocities (np.ndarray): 2D array of velocities (N, 3).
        G (float): Gravitational constant.

    Returns:
        tuple: (total_kinetic_energy, total_potential_energy, total_energy)
    """
    # Kinetic Energy: 0.5 * m * v^2 (already fast with NumPy)
    vel_sq = np.sum(velocities**2, axis=1)
    kinetic_energy = 0.5 * np.sum(masses * vel_sq)

    # Use the Numba parallel potential energy calculation
    # Note: First call to a Numba function triggers compilation, which takes time.
    potential_energy = calculate_potential_energy_parallel(masses, positions, G)

    total_energy = kinetic_energy + potential_energy
    return kinetic_energy, potential_energy, total_energy


def initialize_circular_orbits(num_small_bodies, central_mass, min_radius, max_radius, G):
    """
    Initializes a system with a central massive body and N smaller bodies
    on roughly circular orbits with random orientations and radii.
    (No changes needed here for parallelization of the simulation itself)

    Args:
        num_small_bodies (int): Number of orbiting bodies.
        central_mass (float): Mass of the central body.
        min_radius (float): Minimum orbital radius for small bodies.
        max_radius (float): Maximum orbital radius for small bodies.
        G (float): Gravitational constant.

    Returns:
        tuple: (masses, positions, velocities) as np.ndarrays.
    """
    n_total = num_small_bodies + 1
    masses = np.zeros(n_total)
    positions = np.zeros((n_total, 3))
    velocities = np.zeros((n_total, 3))

    # Central Body (at origin, initially at rest)
    masses[0] = central_mass
    positions[0] = [0.0, 0.0, 0.0]
    velocities[0] = [0.0, 0.0, 0.0]

    # Assign very small mass to orbiting bodies relative to central mass
    small_mass = central_mass * 1e-9
    if small_mass == 0:
        small_mass = 1e-15 # Ensure non-zero mass

    masses[1:] = small_mass

    print(f"Initializing {num_small_bodies} orbiting bodies...")

    # Small Bodies
    for i in range(1, n_total):
        # Random radius
        r = np.random.uniform(min_radius, max_radius)

        # Random position on a sphere of radius r
        theta = np.random.uniform(0, 2 * np.pi) # Azimuthal angle
        phi = np.arccos(np.random.uniform(-1, 1)) # Inclination angle (cos(phi) uniform)
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        positions[i] = [x, y, z]

        # Velocity for circular orbit: v = sqrt(G * M_central / r)
        v_mag = np.sqrt(G * central_mass / r)

        # Velocity direction needs to be perpendicular to the position vector.
        while True:
            random_vec = np.random.rand(3) - 0.5
            plane_normal = np.cross(positions[i], random_vec)
            norm_mag = np.linalg.norm(plane_normal)
            if norm_mag > 1e-8:
                plane_normal /= norm_mag
                break

        vel_dir = np.cross(plane_normal, positions[i])
        vel_dir /= np.linalg.norm(vel_dir) # Normalize

        velocities[i] = v_mag * vel_dir

    print("Initialization complete.")
    return masses, positions, velocities


def run_simulation(masses, positions, velocities, G, dt, num_steps, epsilon, report_interval=10):
    """
    Runs the N-body simulation using the first-order kick-step method,
    calling the parallelized acceleration calculation.

    Args:
        masses (np.ndarray): Masses (N).
        positions (np.ndarray): Initial positions (N, 3).
        velocities (np.ndarray): Initial velocities (N, 3).
        G (float): Gravitational constant.
        dt (float): Time step.
        num_steps (int): Number of steps to simulate.
        epsilon (float): Softening length.
        report_interval (int): How many steps between progress reports.

    Returns:
       tuple: (final_positions, final_velocities) as np.ndarrays.
    """
    print(f"Starting PARALLEL simulation: {len(masses)} bodies, {num_steps} steps, dt={dt}, epsilon={epsilon}")
    # Trigger Numba compilation before timing the main loop (optional but good practice)
    print("Compiling Numba functions (if first run)...")
    _ = calculate_acceleration_parallel(masses[:2], positions[:2], G, epsilon)
    _ = calculate_energy_parallel(masses[:2], positions[:2], velocities[:2], G)
    print("Compilation complete.")

    start_sim_time = time.time()

    for step in range(num_steps):
        # 1. Calculate acceleration using the PARALLEL function
        accel = calculate_acceleration_parallel(masses, positions, G, epsilon)

        # 2. Update velocities (Kick)
        velocities += accel * dt

        # 3. Update positions (Step) using *new* velocities (Euler-Cromer)
        positions += velocities * dt

        # --- Reporting ---
        if (step + 1) % report_interval == 0 or step == num_steps - 1:
            elapsed_time = time.time() - start_sim_time
            # Estimate remaining time (very rough)
            time_per_step = elapsed_time / (step + 1)
            eta = time_per_step * (num_steps - (step + 1))
            print(f"Step {step + 1}/{num_steps} completed. Time elapsed: {elapsed_time:.2f} s. ETA: {eta:.2f} s")

    end_sim_time = time.time()
    print(f"Simulation finished. Total time: {end_sim_time - start_sim_time:.2f} s")
    return positions, velocities

# --- Main Execution ---
if __name__ == "__main__":

    # --- Simulation Parameters ---
    # N_SMALL_BODIES = 1_000_000 # Target number - Still VERY demanding, but feasible now
    N_SMALL_BODIES = 5000      # Use a smaller number for quicker testing (e.g., 5k-10k)
    CENTRAL_MASS = 1.0e6      # Mass of the central body (e.g., solar masses)
    MIN_RADIUS = 1.0          # Min orbital radius (e.g., AU)
    MAX_RADIUS = 100.0        # Max orbital radius (e.g., AU)
    DT = 0.001                # Time step
    N_STEPS = 1000            # Number of simulation steps
    SOFTENING = 0.05          # Softening length

    print("="*30)
    print(" PARALLEL N-Body Simulation Setup")
    print("="*30)
    print(f"Using Numba for parallel acceleration.")
    print(f"Number of CPU cores detected by Numba: {numba.config.NUMBA_DEFAULT_NUM_THREADS}")
    print(f"Number of orbiting bodies: {N_SMALL_BODIES}")
    print(f"Number of simulation steps: {N_STEPS}")
    print(f"Time step (dt): {DT}")
    print(f"Softening length (epsilon): {SOFTENING}")
    print(f"Gravitational Constant (G): {G}")
    print("\nNote: Even parallelized, N=1,000,000 is computationally")
    print("      intensive and may take significant time (hours/days).")
    print("      Consider available RAM and CPU cores.\n")


    # --- Initialization ---
    masses, initial_positions, initial_velocities = initialize_circular_orbits(
        num_small_bodies=N_SMALL_BODIES,
        central_mass=CENTRAL_MASS,
        min_radius=MIN_RADIUS,
        max_radius=MAX_RADIUS,
        G=G
    )

    # --- Initial Energy Calculation ---
    print("\nCalculating initial energy (using parallel PE)...")
    start_energy_time = time.time()
    # Call the parallel energy function
    KE_i, PE_i, E_i = calculate_energy_parallel(masses, initial_positions, initial_velocities, G)
    end_energy_time = time.time()
    print(f"Initial Kinetic Energy : {KE_i:.6e}")
    print(f"Initial Potential Energy: {PE_i:.6e}")
    print(f"Initial Total Energy    : {E_i:.6e}")
    print(f"(Energy calculation took {end_energy_time - start_energy_time:.2f} s)")


    # --- Run Simulation ---
    # Pass copies to the simulation function if you want to preserve originals
    final_positions, final_velocities = run_simulation(
        masses,
        initial_positions.copy(),
        initial_velocities.copy(),
        G=G,
        dt=DT,
        num_steps=N_STEPS,
        epsilon=SOFTENING,
        report_interval=max(1, N_STEPS // 20) # Report progress ~20 times
    )

    # --- Final Energy Calculation ---
    print("\nCalculating final energy (using parallel PE)...")
    start_energy_time = time.time()
    # Call the parallel energy function
    KE_f, PE_f, E_f = calculate_energy_parallel(masses, final_positions, final_velocities, G)
    end_energy_time = time.time()
    print(f"Final Kinetic Energy : {KE_f:.6e}")
    print(f"Final Potential Energy: {PE_f:.6e}")
    print(f"Final Total Energy    : {E_f:.6e}")
    print(f"(Energy calculation took {end_energy_time - start_energy_time:.2f} s)")


    # --- Energy Conservation Check ---
    print("\n--- Energy Conservation ---")
    energy_diff = E_f - E_i
    # Use initial energy for relative error, handle potential zero energy
    relative_error = np.abs(energy_diff / E_i) if abs(E_i) > 1e-12 else np.abs(energy_diff)
    print(f"Initial Total Energy: {E_i:.6e}")
    print(f"Final Total Energy  : {E_f:.6e}")
    print(f"Absolute Difference : {energy_diff:.6e}")
    print(f"Relative Error      : {relative_error:.6e} ({relative_error*100:.4f}%)")

    # Note: Parallel summation order might slightly change floating point results
    # compared to sequential version. fastmath=True can also impact precision.
    if relative_error > 0.01: # Example threshold (1%)
         print("\nWarning: Significant energy drift detected. Parallelism/fastmath might affect precision.")
         print("         Consider smaller dt or disabling fastmath if accuracy is paramount.")

    print("\nSimulation Complete.")