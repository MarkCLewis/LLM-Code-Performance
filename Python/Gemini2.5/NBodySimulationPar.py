import numpy as np
import time
import multiprocessing as mp
from functools import partial

# --- Global variables for worker processes ---
# We define these as global so we don't need to pass them repeatedly,
# which can be inefficient with multiprocessing.
g_mass = None
g_pos = None

def initialize_circular_orbits(n_bodies, central_mass_val):
    """
    Initializes a system with a central body and smaller bodies in random
    circular orbits around it. (Identical to the sequential version)
    """
    total_bodies = n_bodies + 1
    mass = np.ones(total_bodies)
    pos = np.zeros((total_bodies, 3))
    vel = np.zeros((total_bodies, 3))

    mass[0] = central_mass_val
    pos[0] = [0, 0, 0]
    vel[0] = [0, 0, 0]

    for i in range(1, total_bodies):
        radius = np.random.uniform(5.0, 150.0)
        theta = np.arccos(2 * np.random.random() - 1)
        phi = 2 * np.pi * np.random.random()
        pos[i] = [
            radius * np.sin(theta) * np.cos(phi),
            radius * np.sin(theta) * np.sin(phi),
            radius * np.cos(theta)
        ]
        orbital_speed = np.sqrt(central_mass_val / radius)
        random_vec = np.random.randn(3)
        while np.allclose(np.cross(pos[i], random_vec), 0):
            random_vec = np.random.randn(3)
        vel[i] = np.cross(pos[i], random_vec)
        vel[i] = vel[i] / np.linalg.norm(vel[i]) * orbital_speed
    return mass, pos, vel

def calculate_energy_parallel(pool, n_chunks):
    """
    Calculates the total potential and kinetic energy in parallel.
    Uses the approximation that potential energy is only between the central
    body and the orbiting bodies.
    """
    # Kinetic energy is a fast numpy sum, no need to parallelize
    # g_mass and g_vel are accessible to the parent process
    ke = 0.5 * np.sum(g_mass[:, np.newaxis] * g_vel**2)

    # Parallelize the potential energy calculation
    # We create a partial function to pass the fixed central body position
    worker_func = partial(
        worker_calculate_pe_chunk,
        central_pos=g_pos[0],
        central_mass=g_mass[0]
    )
    
    # We only need to calculate PE for orbiting bodies (indices 1 to end)
    orbiting_indices = np.arange(1, len(g_mass))
    index_chunks = np.array_split(orbiting_indices, n_chunks)
    
    pe_results = pool.map(worker_func, index_chunks)
    
    return ke, np.sum(pe_results)

def worker_calculate_pe_chunk(indices, central_pos, central_mass):
    """Worker function to calculate potential energy for a chunk of bodies."""
    pe_chunk = 0.0
    # g_pos and g_mass are the global full arrays
    for i in indices:
        r_vec = g_pos[i] - central_pos
        r = np.linalg.norm(r_vec)
        if r > 0:
            pe_chunk -= (central_mass * g_mass[i]) / r
    return pe_chunk


def calculate_accelerations_parallel(pool, n_chunks):
    """
    Calculates the accelerations of all bodies in parallel using the
    central body approximation.
    """
    # Calculate acceleration for orbiting bodies in parallel
    worker_func = partial(worker_calculate_accel_chunk, central_pos=g_pos[0], central_mass=g_mass[0])
    
    orbiting_indices = np.arange(1, len(g_mass))
    index_chunks = np.array_split(orbiting_indices, n_chunks)

    # Run the parallel computation
    accel_results = pool.map(worker_func, index_chunks)

    # Combine results
    accel = np.zeros_like(g_pos)
    for i, chunk in enumerate(index_chunks):
        accel[chunk] = accel_results[i]

    # Calculate acceleration for the central body (sum of forces from others)
    # This is a reduction, so it's faster to do it sequentially after the parallel part.
    for i in range(1, len(g_mass)):
         r_vec = g_pos[i] - g_pos[0]
         r_mag_sq = np.sum(r_vec**2)
         r_mag = np.sqrt(r_mag_sq)
         # Force on central body is opposite to force on orbiting body
         force_on_central = -g_mass[0] * g_mass[i] * r_vec / (r_mag_sq * r_mag)
         accel[0] += force_on_central / g_mass[0]
         
    return accel

def worker_calculate_accel_chunk(indices, central_pos, central_mass):
    """Worker function that calculates acceleration for a chunk of bodies."""
    accel_chunk = np.zeros((len(indices), 3))
    for i, body_idx in enumerate(indices):
        r_vec = central_pos - g_pos[body_idx]
        r_mag_sq = np.sum(r_vec**2)
        r_mag = np.sqrt(r_mag_sq)
        accel_chunk[i] = central_mass * r_vec / (r_mag_sq * r_mag)
    return accel_chunk

def init_worker(mass, pos, vel):
    """Initializer for each worker process."""
    global g_mass, g_pos, g_vel
    g_mass = mass
    g_pos = pos
    g_vel = vel

def main():
    """
    Main function to run the parallel N-body simulation.
    """
    # Simulation parameters
    n_orbiting_bodies = 10_000
    central_body_mass = 1_000_000
    n_steps = 100
    dt = 0.01

    # Parallelization parameters
    # Use all available cores, but can be set to a specific number
    num_processes = mp.cpu_count()
    # Number of chunks to split the main array into for processing
    n_chunks = num_processes * 4 

    print("Initializing the system...")
    start_time = time.time()
    mass, pos, vel = initialize_circular_orbits(n_orbiting_bodies, central_body_mass)
    end_time = time.time()
    print(f"System initialization took: {end_time - start_time:.2f} seconds")

    # Use a multiprocessing Pool with an initializer
    # The initializer sets up global variables for each worker process to avoid
    # passing large arrays repeatedly.
    with mp.Pool(processes=num_processes, initializer=init_worker, initargs=(mass, pos, vel)) as pool:
        global g_vel # Make g_vel accessible to update it in the loop

        print(f"\nRunning simulation on {num_processes} cores.")
        print("Calculating initial energy...")
        start_time = time.time()
        ke_initial, pe_initial = calculate_energy_parallel(pool, n_chunks)
        total_energy_initial = ke_initial + pe_initial
        end_time = time.time()
        print(f"Energy calculation took: {end_time - start_time:.2f} seconds")
        print(f"Initial Total Energy: {total_energy_initial:.5e}")

        print(f"\nStarting simulation for {n_steps} steps...")
        sim_start_time = time.time()

        for step in range(n_steps):
            # The global g_pos is used by workers. Update g_vel which is also global.
            accel = calculate_accelerations_parallel(pool, n_chunks)
            vel += accel * dt
            pos += vel * dt
            g_vel = vel # Ensure the global velocity is updated for the next energy calc.
            
            if (step + 1) % 1 == 0:
                print(f"  - Step {step + 1}/{n_steps} completed.")
        
        sim_end_time = time.time()
        print(f"Simulation loop took: {sim_end_time - sim_start_time:.2f} seconds")

        print("\nCalculating final energy...")
        start_time = time.time()
        ke_final, pe_final = calculate_energy_parallel(pool, n_chunks)
        total_energy_final = ke_final + pe_final
        end_time = time.time()
        print(f"Energy calculation took: {end_time - start_time:.2f} seconds")
        print(f"Final Total Energy:   {total_energy_final:.5e}")

        energy_change = abs((total_energy_final - total_energy_initial) / total_energy_initial)
        print(f"\nRelative change in total energy: {energy_change:.5%}")

if __name__ == "__main__":
    # The __name__ == "__main__" guard is crucial for multiprocessing
    main()