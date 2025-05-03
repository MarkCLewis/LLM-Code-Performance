import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Array

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
mass_central = 1.989e30  # Mass of the central body (kg), e.g., mass of the Sun
time_step = 1e5  # Time step (s)
num_steps = 1000  # Number of steps to simulate
num_bodies = 1000000  # Number of small bodies
theta = 0.3  # Barnes-Hut approximation threshold (0.3)

# Initialize positions and velocities of bodies
def initialize_bodies(num_bodies, mass_central):
    radius = 1e11  # Distance from the central body (m), roughly 1 AU
    angles = np.linspace(0, 2 * np.pi, num_bodies, endpoint=False)
    x_positions = radius * np.cos(angles)
    y_positions = radius * np.sin(angles)
    z_positions = np.zeros(num_bodies)  # All bodies in the same plane
    
    velocities = np.sqrt(G * mass_central / radius)
    vx_positions = -velocities * np.sin(angles)
    vy_positions = velocities * np.cos(angles)
    vz_positions = np.zeros(num_bodies)  # No motion in the z-direction
    
    positions = np.array([x_positions, y_positions, z_positions]).T
    velocities = np.array([vx_positions, vy_positions, vz_positions]).T
    masses = np.full(num_bodies, 1e12)  # Small mass of each body
    
    return positions, velocities, masses

# Compute the total energy of the system (kinetic + potential)
def compute_energy(positions, velocities, masses, mass_central):
    kinetic_energy = 0.5 * np.sum(masses * np.linalg.norm(velocities, axis=1)**2)
    potential_energy = 0.0
    for i in range(len(masses)):
        r = np.linalg.norm(positions[i])
        potential_energy -= G * mass_central * masses[i] / r
        for j in range(i + 1, len(masses)):
            r_ij = np.linalg.norm(positions[i] - positions[j])
            potential_energy -= G * masses[i] * masses[j] / r_ij
    
    total_energy = kinetic_energy + potential_energy
    return total_energy

# Barnes-Hut Force Computation with KDTree
def compute_force(i, positions, masses, kdtree, theta, mass_central):
    xi, yi, zi = positions[i]
    mi = masses[i]
    
    # Search for nearby bodies using the k-d tree
    neighbors = kdtree.query_ball_point(positions[i], 2 * np.linalg.norm(positions[i]))
    
    # Calculate the total force on body i
    force = np.zeros(3)
    for j in neighbors:
        if i != j:
            r_ij = positions[i] - positions[j]
            r_ij_norm = np.linalg.norm(r_ij)
            force += G * masses[j] * r_ij / r_ij_norm**3
    
    # Barnes-Hut approximation: Apply when bodies are sufficiently far
    total_mass = np.sum(masses[neighbors])
    if total_mass > 0:
        center_of_mass = np.sum(positions[neighbors].T * masses[neighbors], axis=1) / total_mass
        distance_to_center = np.linalg.norm(center_of_mass - positions[i])
        if distance_to_center > theta * np.linalg.norm(positions[i]):
            r_center = center_of_mass - positions[i]
            force += G * total_mass * r_center / distance_to_center**3
    
    return i, force

# Perform the simulation using the kick-step method with Barnes-Hut
def run_simulation(positions, velocities, masses, mass_central, num_steps, time_step, theta):
    # Store energy before the simulation
    initial_energy = compute_energy(positions, velocities, masses, mass_central)
    
    # Create a KDTree for efficient neighbor search
    kdtree = KDTree(positions)
    
    # Use shared memory to store the force accumulations
    force_shared = Array('d', len(masses) * 3)  # 3 for x, y, z components per body
    
    with ProcessPoolExecutor() as executor:
        for step in range(num_steps):
            # Parallelize the force computation using Barnes-Hut with KDTree
            futures = []
            for i in range(len(masses)):
                futures.append(executor.submit(compute_force, i, positions, masses, kdtree, theta, mass_central))
            
            forces = np.zeros_like(positions)
            for future in as_completed(futures):
                i, force = future.result()
                forces[i] = force
            
            # First kick: Update velocities (half step)
            velocities += 0.5 * forces * time_step
            
            # Update positions: full step
            positions += velocities * time_step
            
            # Incrementally update the KDTree with new positions
            kdtree = KDTree(positions)
            
            # Recompute the forces for the next step
            futures = []
            for i in range(len(masses)):
                futures.append(executor.submit(compute_force, i, positions, masses, kdtree, theta, mass_central))
            
            forces = np.zeros_like(positions)
            for future in as_completed(futures):
                i, force = future.result()
                forces[i] = force
            
            # Second kick: Update velocities (half step)
            velocities += 0.5 * forces * time_step
    
    # Compute energy after the simulation
    final_energy = compute_energy(positions, velocities, masses, mass_central)
    
    return initial_energy, final_energy, positions, velocities

# Main function to initialize and run the simulation
def main():
    positions, velocities, masses = initialize_bodies(num_bodies, mass_central)
    
    # Perform the simulation
    initial_energy, final_energy, positions, velocities = run_simulation(positions, velocities, masses, mass_central, num_steps, time_step, theta)
    
    # Display energy before and after the simulation
    print(f"Initial Energy: {initial_energy:.2e} J")
    print(f"Final Energy: {final_energy:.2e} J")
    print(f"Energy difference: {final_energy - initial_energy:.2e} J")
    
    # Optional: Plot the positions of the bodies
    plt.figure(figsize=(8, 8))
    plt.scatter(positions[:, 0], positions[:, 1], s=0.1, color='blue', alpha=0.5)
    plt.title("Positions of Bodies after Simulation")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.axis('equal')
    plt.show()

if __name__ == "__main__":
    main()