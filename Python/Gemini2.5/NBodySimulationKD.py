import numpy as np
import time
import multiprocessing as mp
from functools import partial

# --- Global variables for worker processes ---
# Using globals is a common pattern with multiprocessing to avoid data serialization overhead.
g_pos = None
g_mass = None
g_kdtree_root = None
g_theta = 0.0

class KDNode:
    """A node in the k-D Tree for Barnes-Hut simulation."""
    def __init__(self, particle_indices, bbox):
        self.particle_indices = np.array(particle_indices, dtype=int)
        self.bbox = bbox
        self.is_leaf = len(self.particle_indices) <= 1
        
        self.left_child = None
        self.right_child = None
        
        self.total_mass = 0.0
        self.center_of_mass = np.zeros(3)
        self.size = np.max(self.bbox[1] - self.bbox[0])

def build_kdtree(particle_indices, bbox, depth=0):
    """Recursively builds the k-D tree."""
    node = KDNode(particle_indices, bbox)
    if node.is_leaf:
        return node

    # Determine split dimension and value
    split_dim = depth % 3
    particle_positions = g_pos[particle_indices]
    
    # Find the median to split the particles into two halves
    median_idx = np.argsort(particle_positions[:, split_dim])[len(particle_indices) // 2]
    median_particle_index = particle_indices[median_idx]
    split_val = g_pos[median_particle_index, split_dim]

    left_indices = particle_indices[particle_positions[:, split_dim] < split_val]
    right_indices = particle_indices[particle_positions[:, split_dim] >= split_val]

    # Create new bounding boxes for children
    left_bbox = np.copy(bbox)
    left_bbox[1, split_dim] = split_val
    node.left_child = build_kdtree(left_indices, left_bbox, depth + 1)

    right_bbox = np.copy(bbox)
    right_bbox[0, split_dim] = split_val
    node.right_child = build_kdtree(right_indices, right_bbox, depth + 1)
    
    return node

def annotate_tree(node):
    """Recursively calculates the center of mass and total mass for each node."""
    if node.is_leaf:
        if len(node.particle_indices) > 0:
            idx = node.particle_indices[0]
            node.total_mass = g_mass[idx]
            node.center_of_mass = g_pos[idx]
    else:
        annotate_tree(node.left_child)
        annotate_tree(node.right_child)
        
        m1 = node.left_child.total_mass
        m2 = node.right_child.total_mass
        node.total_mass = m1 + m2
        
        if node.total_mass > 0:
            node.center_of_mass = (m1 * node.left_child.center_of_mass + 
                                  m2 * node.right_child.center_of_mass) / node.total_mass
    return node

def calculate_force_on_particle(particle_idx, node):
    """Calculates the gravitational force on a single particle by traversing the tree."""
    accel = np.zeros(3)
    
    # Use a stack for non-recursive traversal
    stack = [node]
    while stack:
        current_node = stack.pop()
        
        # Do not interact with yourself
        if len(current_node.particle_indices) == 1 and current_node.particle_indices[0] == particle_idx:
            continue
            
        # Calculate distance and size
        r_vec = current_node.center_of_mass - g_pos[particle_idx]
        dist_sq = np.sum(r_vec**2)
        dist = np.sqrt(dist_sq)

        if dist > 0:
            # Barnes-Hut criterion: s/d < theta
            if current_node.size / dist < g_theta or current_node.is_leaf:
                # Treat as a macro-particle
                accel += current_node.total_mass * r_vec / (dist_sq * dist)
            else:
                # Node is too close, traverse its children
                if current_node.right_child: stack.append(current_node.right_child)
                if current_node.left_child: stack.append(current_node.left_child)
    return accel

def worker_calculate_force_chunk(particle_indices):
    """Worker function for multiprocessing. Calculates force for a chunk of particles."""
    chunk_accel = np.zeros((len(particle_indices), 3))
    for i, idx in enumerate(particle_indices):
        chunk_accel[i] = calculate_force_on_particle(idx, g_kdtree_root)
    return chunk_accel
    
def calculate_accelerations_parallel(pool, n_chunks):
    """Orchestrates parallel force calculation."""
    all_indices = np.arange(len(g_mass))
    index_chunks = np.array_split(all_indices, n_chunks)
    
    accel_results = pool.map(worker_calculate_force_chunk, index_chunks)
    
    # Stitch results back together
    return np.vstack(accel_results)

def init_worker(mass, pos, vel, root, theta):
    """Initializer for each worker process."""
    global g_mass, g_pos, g_vel, g_kdtree_root, g_theta
    g_mass, g_pos, g_vel = mass, pos, vel
    g_kdtree_root, g_theta = root, theta

# (The initialization and energy calculation functions from the previous example
#  can be reused here without modification. For brevity, they are omitted but
#  should be included in the final script.)
def initialize_circular_orbits(n_bodies, central_mass_val):
    total_bodies = n_bodies + 1; mass = np.ones(total_bodies); pos = np.zeros((total_bodies, 3)); vel = np.zeros((total_bodies, 3)); mass[0] = central_mass_val
    for i in range(1, total_bodies):
        radius = np.random.uniform(5.0, 150.0); theta_ang = np.arccos(2 * np.random.random() - 1); phi = 2 * np.pi * np.random.random()
        pos[i] = [radius * np.sin(theta_ang) * np.cos(phi), radius * np.sin(theta_ang) * np.sin(phi), radius * np.cos(theta_ang)]
        orbital_speed = np.sqrt(central_mass_val / radius); random_vec = np.random.randn(3)
        while np.allclose(np.cross(pos[i], random_vec), 0): random_vec = np.random.randn(3)
        vel[i] = np.cross(pos[i], random_vec); vel[i] = vel[i] / np.linalg.norm(vel[i]) * orbital_speed
    return mass, pos, vel
def calculate_energy(mass, pos, vel):
    ke = 0.5 * np.sum(mass[:, np.newaxis] * vel**2); pe = 0
    for i in range(pos.shape[0]):
        for j in range(i + 1, pos.shape[0]):
            r_vec = pos[j] - pos[i]; r = np.linalg.norm(r_vec)
            if r > 0: pe -= (mass[i] * mass[j]) / r
    return ke, pe
    
def main():
    # Simulation parameters
    # NOTE: 1M bodies is very memory intensive for k-D trees. 
    # A number like 100,000 is more feasible on consumer hardware.
    n_orbiting_bodies = 100_000 
    central_body_mass = 1_000_000
    n_steps = 10
    dt = 0.01
    theta = 0.3

    # Parallelization parameters
    num_processes = mp.cpu_count()
    n_chunks = num_processes * 4

    print("Initializing system...")
    mass, pos, vel = initialize_circular_orbits(n_orbiting_bodies, central_body_mass)
    
    # Calculate initial energy using direct summation for accuracy
    print("Calculating initial energy (direct summation)...")
    ke_initial, pe_initial = calculate_energy(mass, pos, vel)
    total_energy_initial = ke_initial + pe_initial
    print(f"Initial Total Energy: {total_energy_initial:.5e}")

    print(f"\nStarting simulation on {num_processes} cores with theta = {theta}...")
    sim_start_time = time.time()
    
    # Main simulation loop
    for step in range(n_steps):
        # 1. Build the k-D Tree (sequentially)
        min_coords = np.min(pos, axis=0)
        max_coords = np.max(pos, axis=0)
        bbox = np.array([min_coords, max_coords])
        root_node = build_kdtree(np.arange(len(mass)), bbox)

        # 2. Annotate the Tree with CoM info (sequentially)
        # We update the global position reference for the annotation function
        global g_pos, g_mass
        g_pos = pos
        g_mass = mass
        annotate_tree(root_node)
        
        # 3. Calculate forces in parallel using the tree
        with mp.Pool(processes=num_processes, initializer=init_worker, initargs=(mass, pos, vel, root_node, theta)) as pool:
            accel = calculate_accelerations_parallel(pool, n_chunks)

        # 4. Update velocities and positions (Kick-Step)
        vel += accel * dt
        pos += vel * dt

        if (step + 1) % 100 == 0:
            print(f"  - Step {step + 1}/{n_steps} completed.")
    
    sim_end_time = time.time()
    print(f"Simulation loop took: {sim_end_time - sim_start_time:.2f} seconds")

    print("\nCalculating final energy (direct summation)...")
    ke_final, pe_final = calculate_energy(mass, pos, vel)
    total_energy_final = ke_final + pe_final
    print(f"Final Total Energy:   {total_energy_final:.5e}")
    
    energy_change = abs((total_energy_final - total_energy_initial) / total_energy_initial)
    print(f"\nRelative change in total energy: {energy_change:.5%}")


if __name__ == "__main__":
    main()