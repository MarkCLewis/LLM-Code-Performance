import numpy as np
import time
import math

# Gravitational constant (can be set to 1 for simulation units)
G = 1.0

# --- Physics and Simulation Functions ---

def calculate_acceleration(masses, positions, G, epsilon):
    """
    Calculates the acceleration of each body due to gravitational interaction.
    Uses direct summation (O(N^2)).

    Args:
        masses (np.ndarray): 1D array of masses (N).
        positions (np.ndarray): 2D array of positions (N, 3).
        G (float): Gravitational constant.
        epsilon (float): Softening length to avoid singularities.

    Returns:
        np.ndarray: 2D array of accelerations (N, 3).
    """
    n_bodies = masses.shape[0]
    accelerations = np.zeros_like(positions)
    epsilon_sq = epsilon**2

    # Loop over each unique pair of bodies (i, j) where i < j
    # This avoids double counting and self-interaction
    for i in range(n_bodies):
        for j in range(i + 1, n_bodies):
            # Vector difference
            r_ij = positions[j] - positions[i]

            # Squared distance + softening
            dist_sq = np.sum(r_ij**2) + epsilon_sq

            # Inverse cube distance (avoiding sqrt for efficiency)
            inv_dist_cubed = dist_sq**(-1.5)

            # Force magnitude scaled by G / dist^3
            force_scale = G * inv_dist_cubed

            # Calculate force vector components efficiently
            force_ij = force_scale * r_ij

            # Distribute acceleration according to F=ma (a = F/m)
            # Newton's 3rd Law: F_ji = -F_ij
            accelerations[i] += masses[j] * force_ij
            accelerations[j] -= masses[i] * force_ij

    return accelerations

def calculate_energy(masses, positions, velocities, G):
    """
    Calculates the total kinetic, potential, and total energy of the system.
    Potential energy calculation is O(N^2).

    Args:
        masses (np.ndarray): 1D array of masses (N).
        positions (np.ndarray): 2D array of positions (N, 3).
        velocities (np.ndarray): 2D array of velocities (N, 3).
        G (float): Gravitational constant.

    Returns:
        tuple: (total_kinetic_energy, total_potential_energy, total_energy)
    """
    # Kinetic Energy: 0.5 * m * v^2
    # sum(velocities**2, axis=1) calculates vx^2 + vy^2 + vz^2 for each body
    vel_sq = np.sum(velocities**2, axis=1)
    kinetic_energy = 0.5 * np.sum(masses * vel_sq)

    # Potential Energy: -G * m1 * m2 / r
    # Sum over all unique pairs (i < j)
    n_bodies = masses.shape[0]
    potential_energy = 0.0
    for i in range(n_bodies):
        for j in range(i + 1, n_bodies):
            r_ij = positions[j] - positions[i]
            dist = np.sqrt(np.sum(r_ij**2))
            # Avoid division by zero if particles somehow end up at the exact same spot
            # Though softening in dynamics should prevent this becoming extreme.
            if dist > 1e-10: # Add a small threshold
                 potential_energy -= G * masses[i] * masses[j] / dist

    total_energy = kinetic_energy + potential_energy
    return kinetic_energy, potential_energy, total_energy

def initialize_circular_orbits(num_small_bodies, central_mass, min_radius, max_radius, G):
    """
    Initializes a system with a central massive body and N smaller bodies
    on roughly circular orbits with random orientations and radii.

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
    # Avoid zero mass, but keep it small enough not to significantly affect
    # the central body's motion in this simplified setup.
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
        # Use uniform distribution on sphere surface
        theta = np.random.uniform(0, 2 * np.pi) # Azimuthal angle
        phi = np.arccos(np.random.uniform(-1, 1)) # Inclination angle (cos(phi) uniform)
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        positions[i] = [x, y, z]

        # Velocity for circular orbit: v = sqrt(G * M_central / r)
        v_mag = np.sqrt(G * central_mass / r)

        # Velocity direction needs to be perpendicular to the position vector.
        # Generate a random vector, ensure it's not parallel to position,
        # then use cross products to find a perpendicular direction in a
        # randomly oriented plane containing the position vector.
        while True:
            # Random vector for defining the orbital plane normal
            random_vec = np.random.rand(3) - 0.5
            # Calculate orbital plane normal (perpendicular to position and random vec)
            plane_normal = np.cross(positions[i], random_vec)
            norm_mag = np.linalg.norm(plane_normal)
            if norm_mag > 1e-8: # Ensure position and random vector weren't parallel
                plane_normal /= norm_mag
                break
            # If parallel, try a different random vector

        # Velocity direction (perpendicular to position and plane normal)
        vel_dir = np.cross(plane_normal, positions[i])
        vel_dir /= np.linalg.norm(vel_dir) # Normalize

        velocities[i] = v_mag * vel_dir

    print("Initialization complete.")
    return masses, positions, velocities


def run_simulation(masses, positions, velocities, G, dt, num_steps, epsilon, report_interval=10):
    """
    Runs the N-body simulation using the first-order kick-step method.

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
    print(f"Starting simulation: {len(masses)} bodies, {num_steps} steps, dt={dt}, epsilon={epsilon}")
    start_sim_time = time.time()

    for step in range(num_steps):
        # 1. Calculate acceleration based on current positions
        accel = calculate_acceleration(masses, positions, G, epsilon)

        # 2. Update velocities (Kick) using current acceleration
        # v(t + dt) = v(t) + a(t) * dt
        velocities += accel * dt

        # 3. Update positions (Step) using *new* velocities (Euler-Cromer)
        # r(t + dt) = r(t) + v(t + dt) * dt
        positions += velocities * dt

        # --- Reporting ---
        if (step + 1) % report_interval == 0 or step == num_steps - 1:
            elapsed_time = time.time() - start_sim_time
            print(f"Step {step + 1}/{num_steps} completed. Time elapsed: {elapsed_time:.2f} s")
            # Optional: Could calculate and report energy here periodically too

    end_sim_time = time.time()
    print(f"Simulation finished. Total time: {end_sim_time - start_sim_time:.2f} s")
    return positions, velocities

# --- Main Execution ---
if __name__ == "__main__":

    # --- Simulation Parameters ---
    N_SMALL_BODIES = 1_000_000 # Target number - BEWARE: EXTREMELY SLOW!
    # N_SMALL_BODIES = 100      # Use a much smaller number for testing (e.g., 100-1000)
    CENTRAL_MASS = 1.0e6      # Mass of the central body (e.g., solar masses)
    MIN_RADIUS = 1.0          # Min orbital radius (e.g., AU)
    MAX_RADIUS = 100.0        # Max orbital radius (e.g., AU)
    DT = 0.001                # Time step (smaller is generally more accurate but slower)
    N_STEPS = 1000            # Number of simulation steps
    SOFTENING = 0.05          # Softening length to prevent high forces

    print("="*30)
    print(" N-Body Simulation Setup")
    print("="*30)
    print(f"Number of orbiting bodies: {N_SMALL_BODIES}")
    print(f"Number of simulation steps: {N_STEPS}")
    print(f"Time step (dt): {DT}")
    print(f"Softening length (epsilon): {SOFTENING}")
    print(f"Gravitational Constant (G): {G}")
    print("\nWARNING: Running with N_SMALL_BODIES = 1,000,000 will be")
    print("         EXTREMELY slow due to O(N^2) complexity.")
    print("         Consider reducing N_SMALL_BODIES for testing.\n")


    # --- Initialization ---
    masses, initial_positions, initial_velocities = initialize_circular_orbits(
        num_small_bodies=N_SMALL_BODIES,
        central_mass=CENTRAL_MASS,
        min_radius=MIN_RADIUS,
        max_radius=MAX_RADIUS,
        G=G
    )

    # --- Initial Energy Calculation ---
    print("\nCalculating initial energy...")
    start_energy_time = time.time()
    KE_i, PE_i, E_i = calculate_energy(masses, initial_positions, initial_velocities, G)
    end_energy_time = time.time()
    print(f"Initial Kinetic Energy : {KE_i:.6e}")
    print(f"Initial Potential Energy: {PE_i:.6e}")
    print(f"Initial Total Energy    : {E_i:.6e}")
    print(f"(Energy calculation took {end_energy_time - start_energy_time:.2f} s)")


    # --- Run Simulation ---
    final_positions, final_velocities = run_simulation(
        masses,
        initial_positions.copy(), # Pass copies to avoid modifying originals if needed later
        initial_velocities.copy(),
        G=G,
        dt=DT,
        num_steps=N_STEPS,
        epsilon=SOFTENING,
        report_interval=max(1, N_STEPS // 20) # Report progress ~20 times
    )

    # --- Final Energy Calculation ---
    print("\nCalculating final energy...")
    start_energy_time = time.time()
    KE_f, PE_f, E_f = calculate_energy(masses, final_positions, final_velocities, G)
    end_energy_time = time.time()
    print(f"Final Kinetic Energy : {KE_f:.6e}")
    print(f"Final Potential Energy: {PE_f:.6e}")
    print(f"Final Total Energy    : {E_f:.6e}")
    print(f"(Energy calculation took {end_energy_time - start_energy_time:.2f} s)")


    # --- Energy Conservation Check ---
    print("\n--- Energy Conservation ---")
    energy_diff = E_f - E_i
    relative_error = np.abs(energy_diff / E_i) if abs(E_i) > 1e-12 else np.abs(energy_diff)
    print(f"Initial Total Energy: {E_i:.6e}")
    print(f"Final Total Energy  : {E_f:.6e}")
    print(f"Absolute Difference : {energy_diff:.6e}")
    print(f"Relative Error      : {relative_error:.6e} ({relative_error*100:.4f}%)")

    if relative_error > 0.01: # Example threshold (1%)
         print("\nWarning: Significant energy drift detected. Consider smaller dt or higher-order integrator.")

    print("\nSimulation Complete.")