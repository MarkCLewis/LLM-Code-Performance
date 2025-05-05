# -*- coding: utf-8 -*-
import numpy as np
import time
import math
import numba # Still useful for energy calculation potentially
from concurrent.futures import ProcessPoolExecutor # For parallelizing tree force calculation
import os # To get cpu count

# Gravitational constant
G = 1.0
# Barnes-Hut opening angle parameter
THETA_DEFAULT = 0.3

# --- k-D Tree Node Structure ---
# Using a simple Python class. Numba won't accelerate traversal directly
# unless we use jitclass (complex) or flatten the tree (complex).
class Node:
    def __init__(self, bounds, particles_indices, positions, masses):
        self.bounds = np.array(bounds, dtype=np.float64) # [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        self.particles_indices = particles_indices # Indices of particles in this node
        self.num_particles = len(particles_indices)
        self.children = [] # Child nodes
        self.is_leaf = self.num_particles <= 1 # Treat nodes with 0 or 1 particle as leaves

        if self.num_particles == 0:
            self.total_mass = 0.0
            self.center_of_mass = np.zeros(3, dtype=np.float64)
            # Size calculation needs care for empty nodes, maybe use parent's size?
            # Or set size to 0, it won't be used in force calcs anyway.
            self.size = 0.0

        elif self.num_particles == 1:
            self.is_leaf = True
            idx = self.particles_indices[0]
            self.total_mass = masses[idx]
            self.center_of_mass = positions[idx]
            # Size for a leaf node? Often treated differently or based on bounds.
            # Let's use the bounds size.
            self.size = np.max(self.bounds[1] - self.bounds[0])
        else:
            # Internal node or leaf with >1 particle (we'll subdivide if >1)
            self.is_leaf = False # Assume we will subdivide internal nodes > 1 particle
            node_masses = masses[self.particles_indices]
            node_positions = positions[self.particles_indices]
            self.total_mass = np.sum(node_masses)
            if self.total_mass > 1e-12: # Avoid division by zero
                 # Weighted average for COM
                self.center_of_mass = np.sum(node_masses[:, np.newaxis] * node_positions, axis=0) / self.total_mass
            else:
                # If total mass is effectively zero, place COM at geometric center
                self.center_of_mass = np.mean(self.bounds, axis=0)

            self.size = np.max(self.bounds[1] - self.bounds[0]) # Use max dimension of bounds

    def get_subdivision_axis_and_value(self, positions):
        """ Determine axis and value to split particles """
        # Choose the widest dimension of the current node's bounds
        dimensions = self.bounds[1] - self.bounds[0]
        axis = np.argmax(dimensions)

        # Choose the split value - median particle position along the axis
        node_positions = positions[self.particles_indices]
        median_value = np.median(node_positions[:, axis])

        # Alternative: split at the geometric center of the bounds
        # median_value = np.mean(self.bounds[:, axis])

        return axis, median_value

# --- k-D Tree Building Function ---
# This runs recursively in pure Python
def build_kdtree(particles_indices, positions, masses, bounds, max_leaf_size=1):
    """
    Recursively builds the k-D Tree.

    Args:
        particles_indices (list or np.array): Indices of particles for this node.
        positions (np.ndarray): Global positions array (N, 3).
        masses (np.ndarray): Global masses array (N).
        bounds (np.ndarray): [[xmin, ymin, zmin], [xmax, ymax, zmax]] for this node.
        max_leaf_size (int): Maximum particles in a leaf node (typically 1 for BH).

    Returns:
        Node: The root node of the (sub)tree.
    """
    node = Node(bounds, particles_indices, positions, masses)

    # Stop subdividing if node has few particles or is empty
    if node.num_particles <= max_leaf_size:
        node.is_leaf = True # Ensure it's marked as leaf
        return node

    # Subdivide if internal node
    node.is_leaf = False
    axis, split_value = node.get_subdivision_axis_and_value(positions)

    # Partition particles based on the split value
    left_indices = []
    right_indices = []
    node_positions = positions[particles_indices] # Get positions relevant to this node

    for i, idx in enumerate(particles_indices):
        if node_positions[i, axis] < split_value:
            left_indices.append(idx)
        else:
            right_indices.append(idx)

    # Handle cases where partition fails (all points on one side)
    # This can happen with duplicate coordinates or poor split choice.
    # A simple fix is to just put one particle in one child and the rest in the other,
    # or just stop subdividing. Let's stop subdividing for simplicity.
    if not left_indices or not right_indices:
        node.is_leaf = True
        # Recalculate properties as if it's a leaf with multiple particles
        node_masses = masses[particles_indices]
        node_positions = positions[particles_indices]
        node.total_mass = np.sum(node_masses)
        if node.total_mass > 1e-12:
             node.center_of_mass = np.sum(node_masses[:, np.newaxis] * node_positions, axis=0) / node.total_mass
        else:
             node.center_of_mass = np.mean(node.bounds, axis=0)
        # Size is already set in __init__
        return node


    # Create new bounds for children
    left_bounds = np.copy(node.bounds)
    left_bounds[1, axis] = split_value # Max value on split axis becomes split_value

    right_bounds = np.copy(node.bounds)
    right_bounds[0, axis] = split_value # Min value on split axis becomes split_value

    # Recursively build children
    left_child = build_kdtree(left_indices, positions, masses, left_bounds, max_leaf_size)
    right_child = build_kdtree(right_indices, positions, masses, right_bounds, max_leaf_size)
    node.children = [left_child, right_child]

    # Non-leaf nodes already had mass/COM calculated in __init__
    # based on all particles within them. No need to recalculate from children for k-D tree.

    return node


# --- Force Calculation using Tree ---
# This function runs for EACH particle, likely in a separate process.
# It recursively traverses the tree (in pure Python).
def calculate_force_tree_single_particle(p_idx, node, positions, masses, G, epsilon_sq, theta):
    """
    Calculates gravitational force on a single particle p_idx using the tree.
    This function is intended to be called in parallel for each particle.

    Args:
        p_idx (int): Index of the particle to calculate force for.
        node (Node): The current node in the tree traversal (start with root).
        positions (np.ndarray): Global positions array.
        masses (np.ndarray): Global masses array.
        G (float): Gravitational constant.
        epsilon_sq (float): Squared softening length.
        theta (float): Barnes-Hut opening angle parameter.

    Returns:
        np.ndarray: Force vector (3,) acting on particle p_idx.
    """
    force = np.zeros(3, dtype=np.float64)
    pos_p = positions[p_idx]

    # Stack for non-recursive traversal (can avoid Python recursion depth limits)
    stack = [node]
    while stack:
        current_node = stack.pop()

        if current_node is None or current_node.num_particles == 0:
            continue

        # Vector from node's COM to particle p
        dp = pos_p - current_node.center_of_mass
        dist_sq = np.sum(dp**2)

        # Check if p_idx is the only particle in a leaf node
        is_self_leaf = current_node.is_leaf and \
                       current_node.num_particles == 1 and \
                       current_node.particles_indices[0] == p_idx

        if is_self_leaf:
            continue # Don't interact particle with itself

        # Barnes-Hut criterion: s/d < theta  =>  (s/d)^2 < theta^2 => s^2 / d^2 < theta^2
        # Use squared values to avoid sqrt: size^2 / dist_sq < theta^2
        # Use dist_sq directly calculated. Add epsilon_sq to avoid division by zero if COM is at particle pos.
        if dist_sq < 1e-12: # Particle is effectively at the COM, handle carefully
            # This case is complex. If it's an internal node, we must open it.
            # If it's a leaf node (and not the particle itself), calculate direct force with softening.
             if current_node.is_leaf:
                 # Must be a different particle q_idx
                 q_idx = current_node.particles_indices[0]
                 # Direct force calculation (using epsilon_sq)
                 r_ij = positions[q_idx] - pos_p # Vector from p to q
                 # Distance uses softening
                 dist_sq_soft = np.sum(r_ij**2) + epsilon_sq
                 inv_dist_cubed = dist_sq_soft**(-1.5)
                 force += G * masses[p_idx] * masses[q_idx] * (-r_ij) * inv_dist_cubed # Force on p is toward q
                 continue # Finished with this leaf
             else:
                 # Internal node where particle is at COM - must open
                 for child in current_node.children:
                     if child: stack.append(child)
                 continue

        # Main BH criterion
        if (current_node.size**2 / dist_sq) < (theta**2):
            # Node is far enough away - treat as point mass
            # Force: -G * M_node * m_p * dp / |dp|^3
            # Add epsilon_sq to denominator distance calculation for force
            dist_sq_force = dist_sq + epsilon_sq # Use softening for force calc distance
            inv_dist_cubed = dist_sq_force**(-1.5)
            force -= G * current_node.total_mass * masses[p_idx] * dp * inv_dist_cubed

        else:
            # Node is too close or too large - open it
            if current_node.is_leaf:
                # Leaf node (must contain particles other than p_idx based on logic above)
                # Calculate direct force for all particles in the leaf
                for q_idx in current_node.particles_indices:
                    if p_idx == q_idx: continue # Should not happen based on is_self_leaf check, but safe check

                    r_pq = positions[q_idx] - pos_p # Vector from p to q
                    dist_sq_pq = np.sum(r_pq**2) + epsilon_sq # Softened distance squared
                    inv_dist_cubed_pq = dist_sq_pq**(-1.5)
                    # Force on p due to q
                    force += G * masses[p_idx] * masses[q_idx] * r_pq * inv_dist_cubed_pq

            else:
                # Internal node - add children to stack
                 for child in current_node.children:
                     if child: stack.append(child) # Add valid children to stack

    return force


# --- Wrapper Function for Parallel Acceleration Calculation ---
def calculate_acceleration_tree(positions, masses, G, epsilon, theta, num_workers):
    """
    Calculates accelerations for all particles using the k-D Tree and parallel processing.

    Args:
        positions (np.ndarray): Global positions array (N, 3).
        masses (np.ndarray): Global masses array (N).
        G (float): Gravitational constant.
        epsilon (float): Softening length.
        theta (float): Barnes-Hut opening angle.
        num_workers (int): Number of parallel processes to use.

    Returns:
        np.ndarray: Accelerations array (N, 3).
    """
    n_bodies = masses.shape[0]
    if n_bodies == 0:
        return np.zeros_like(positions)

    print("Building k-D Tree...")
    start_build = time.time()
    # Calculate initial bounds tightly
    min_coords = np.min(positions, axis=0) - 1e-3 # Add small buffer
    max_coords = np.max(positions, axis=0) + 1e-3
    initial_bounds = np.array([min_coords, max_coords])
    all_indices = list(range(n_bodies))
    root_node = build_kdtree(all_indices, positions, masses, initial_bounds, max_leaf_size=1)
    end_build = time.time()
    print(f"Tree build took {end_build - start_build:.3f} s")

    print(f"Calculating forces using tree on {num_workers} workers...")
    start_force = time.time()
    accelerations = np.zeros_like(positions)
    epsilon_sq = epsilon**2

    # Prepare arguments for parallel execution
    # Need to pass data that workers need. Tree structure (root_node), positions, masses are shared.
    # Using ProcessPoolExecutor - data is pickled and sent to workers. Large arrays can cause overhead.
    # Consider shared memory for optimization if this becomes a bottleneck.
    args = [(i, root_node, positions, masses, G, epsilon_sq, theta) for i in range(n_bodies)]

    forces = np.zeros_like(positions)
    # Use ProcessPoolExecutor for parallel map
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # map calculates results in order
        results = list(executor.map(calculate_force_tree_single_particle_unpack, args))

    # results should be a list of force vectors
    forces = np.array(results) # Convert list of forces back to Nx3 array

    # Acceleration a = F / m
    # Handle division by zero if any mass is zero
    non_zero_mass_indices = masses > 1e-15
    accelerations[non_zero_mass_indices] = forces[non_zero_mass_indices] / masses[non_zero_mass_indices, np.newaxis]

    end_force = time.time()
    print(f"Force calculation took {end_force - start_force:.3f} s")

    return accelerations

# Helper function to unpack arguments for executor.map, as it only takes one iterable
def calculate_force_tree_single_particle_unpack(args):
    return calculate_force_tree_single_particle(*args)


# --- Energy calculation using Numba (still O(N^2) but parallel) ---
# Keep this for accuracy check, although it will be slower than the force calc now
@numba.njit(parallel=True, fastmath=True)
def calculate_potential_energy_numba(masses, positions, G):
    n_bodies = masses.shape[0]
    total_pe = 0.0
    for i in numba.prange(n_bodies):
        pe_i = 0.0
        pos_i = positions[i]
        mass_i = masses[i]
        for j in range(i + 1, n_bodies):
            r_ij = positions[j] - pos_i
            dist = np.sqrt(np.sum(r_ij**2))
            if dist > 1e-12:
                 pe_i -= G * mass_i * masses[j] / dist
        total_pe += pe_i
    return total_pe

def calculate_energy_accurate(masses, positions, velocities, G):
    vel_sq = np.sum(velocities**2, axis=1)
    kinetic_energy = 0.5 * np.sum(masses * vel_sq)
    potential_energy = calculate_potential_energy_numba(masses, positions, G)
    total_energy = kinetic_energy + potential_energy
    return kinetic_energy, potential_energy, total_energy


# --- Initialization (Unchanged) ---
def initialize_circular_orbits(num_small_bodies, central_mass, min_radius, max_radius, G):
    n_total = num_small_bodies + 1
    masses = np.zeros(n_total)
    positions = np.zeros((n_total, 3))
    velocities = np.zeros((n_total, 3))
    masses[0] = central_mass
    positions[0] = [0.0, 0.0, 0.0]
    velocities[0] = [0.0, 0.0, 0.0]
    small_mass = central_mass * 1e-9
    if small_mass == 0: small_mass = 1e-15
    masses[1:] = small_mass
    print(f"Initializing {num_small_bodies} orbiting bodies...")
    for i in range(1, n_total):
        r = np.random.uniform(min_radius, max_radius)
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.arccos(np.random.uniform(-1, 1))
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        positions[i] = [x, y, z]
        v_mag = np.sqrt(G * central_mass / r)
        while True:
            random_vec = np.random.rand(3) - 0.5
            plane_normal = np.cross(positions[i], random_vec)
            norm_mag = np.linalg.norm(plane_normal)
            if norm_mag > 1e-8:
                plane_normal /= norm_mag
                break
        vel_dir = np.cross(plane_normal, positions[i])
        vel_dir /= np.linalg.norm(vel_dir)
        velocities[i] = v_mag * vel_dir
    print("Initialization complete.")
    return masses, positions, velocities

# --- Simulation Loop ---
def run_simulation(masses, positions, velocities, G, dt, num_steps, epsilon, theta, num_workers, report_interval=10):
    """
    Runs the N-body simulation using k-D Tree acceleration calculation.
    """
    print(f"Starting TREE simulation ({num_workers} workers): {len(masses)} bodies, {num_steps} steps, dt={dt}, epsilon={epsilon}, theta={theta}")

    # Pre-compile Numba energy function (optional)
    # print("Compiling Numba energy function (if first run)...")
    # _ = calculate_energy_accurate(masses[:2], positions[:2], velocities[:2], G)
    # print("Compilation complete.")

    start_sim_time = time.time()

    for step in range(num_steps):
        step_start_time = time.time()

        # 1. Calculate acceleration using the TREE method
        accel = calculate_acceleration_tree(positions, masses, G, epsilon, theta, num_workers)

        # 2. Update velocities (Kick)
        velocities += accel * dt

        # 3. Update positions (Step) using *new* velocities (Euler-Cromer)
        positions += velocities * dt

        # --- Reporting ---
        step_end_time = time.time()
        step_duration = step_end_time - step_start_time
        total_elapsed_time = step_end_time - start_sim_time

        if (step + 1) % report_interval == 0 or step == num_steps - 1:
            # Estimate remaining time (very rough)
            avg_step_time = total_elapsed_time / (step + 1)
            eta = avg_step_time * (num_steps - (step + 1))
            print(f"Step {step + 1}/{num_steps} completed. Step time: {step_duration:.3f} s. Total time: {total_elapsed_time:.2f} s. ETA: {eta:.2f} s")


    end_sim_time = time.time()
    print(f"Simulation finished. Total time: {end_sim_time - start_sim_time:.2f} s")
    return positions, velocities

# --- Main Execution ---
if __name__ == "__main__":

    # --- Simulation Parameters ---
    # N_SMALL_BODIES = 1_000_000 # Still very large, but O(N log N) makes it much more possible
    N_SMALL_BODIES = 50000      # A more reasonable number for testing (e.g., 50k)
    CENTRAL_MASS = 1.0e6
    MIN_RADIUS = 1.0
    MAX_RADIUS = 100.0
    DT = 0.001
    N_STEPS = 1000            # Reduce steps for faster testing initially?
    SOFTENING = 0.05
    THETA = 0.3               # Barnes-Hut opening angle
    NUM_WORKERS = os.cpu_count() # Use all available CPU cores

    print("="*40)
    print(" k-D Tree (Barnes-Hut) N-Body Simulation")
    print("="*40)
    print(f"Using k-D Tree (theta={THETA}) for O(N log N) acceleration.")
    print(f"Parallel force calculation using {NUM_WORKERS} workers.")
    print(f"Number of orbiting bodies: {N_SMALL_BODIES}")
    print(f"Number of simulation steps: {N_STEPS}")
    print(f"Time step (dt): {DT}")
    print(f"Softening length (epsilon): {SOFTENING}")
    print(f"Gravitational Constant (G): {G}")
    print("\nNote: Tree building and communication overhead exist.")
    print("      Energy calculation uses O(N^2) direct sum for accuracy check.\n")


    # --- Initialization ---
    masses, initial_positions, initial_velocities = initialize_circular_orbits(
        num_small_bodies=N_SMALL_BODIES,
        central_mass=CENTRAL_MASS,
        min_radius=MIN_RADIUS,
        max_radius=MAX_RADIUS,
        G=G
    )

    # --- Initial Energy Calculation (Accurate O(N^2)) ---
    print("\nCalculating initial energy (using O(N^2) parallel Numba sum)...")
    start_energy_time = time.time()
    # Make sure Numba PE function is compiled before timing if desired
    if N_SMALL_BODIES > 0:
       _ = calculate_potential_energy_numba(masses[:min(2, len(masses))], initial_positions[:min(2, len(masses))], G) # Pre-compile
    KE_i, PE_i, E_i = calculate_energy_accurate(masses, initial_positions, initial_velocities, G)
    end_energy_time = time.time()
    print(f"Initial Kinetic Energy : {KE_i:.6e}")
    print(f"Initial Potential Energy: {PE_i:.6e} (Numba O(N^2) calculation)")
    print(f"Initial Total Energy    : {E_i:.6e}")
    print(f"(Energy calculation took {end_energy_time - start_energy_time:.2f} s)")


    # --- Run Simulation ---
    final_positions, final_velocities = run_simulation(
        masses,
        initial_positions.copy(),
        initial_velocities.copy(),
        G=G,
        dt=DT,
        num_steps=N_STEPS,
        epsilon=SOFTENING,
        theta=THETA,
        num_workers=NUM_WORKERS,
        report_interval=max(1, N_STEPS // 20) # Report progress ~20 times
    )

    # --- Final Energy Calculation (Accurate O(N^2)) ---
    print("\nCalculating final energy (using O(N^2) parallel Numba sum)...")
    start_energy_time = time.time()
    KE_f, PE_f, E_f = calculate_energy_accurate(masses, final_positions, final_velocities, G)
    end_energy_time = time.time()
    print(f"Final Kinetic Energy : {KE_f:.6e}")
    print(f"Final Potential Energy: {PE_f:.6e} (Numba O(N^2) calculation)")
    print(f"Final Total Energy    : {E_f:.6e}")
    print(f"(Energy calculation took {end_energy_time - start_energy_time:.2f} s)")


    # --- Energy Conservation Check ---
    print("\n--- Energy Conservation (based on O(N^2) PE calculation) ---")
    energy_diff = E_f - E_i
    relative_error = np.abs(energy_diff / E_i) if abs(E_i) > 1e-12 else np.abs(energy_diff)
    print(f"Initial Total Energy: {E_i:.6e}")
    print(f"Final Total Energy  : {E_f:.6e}")
    print(f"Absolute Difference : {energy_diff:.6e}")
    print(f"Relative Error      : {relative_error:.6e} ({relative_error*100:.4f}%)")

    print(f"\nNote: Energy drift is expected with approximate methods (theta={THETA})")
    print( "      and the first-order integrator. Smaller theta or dt improves accuracy.")
    if relative_error > 0.05: # Example threshold (5%)
         print("\nWarning: Potentially significant energy drift detected.")

    print("\nSimulation Complete.")