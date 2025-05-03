import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
mass_central = 1.989e30  # Mass of the central body (kg), e.g., mass of the Sun
time_step = 1e5  # Time step (s)
num_steps = 1000  # Number of steps to simulate
num_bodies = 1000000  # Number of small bodies

# Initialize positions and velocities of bodies
def initialize_bodies(num_bodies, mass_central):
    # Place small bodies in circular orbits around the central body
    # Assume all bodies are at the same distance from the central body, in a circular orbit
    radius = 1e11  # Distance from the central body (m), roughly 1 AU
    angles = np.linspace(0, 2 * np.pi, num_bodies, endpoint=False)
    x_positions = radius * np.cos(angles)
    y_positions = radius * np.sin(angles)
    z_positions = np.zeros(num_bodies)  # All bodies in the same plane
    
    # Velocities based on circular orbits: v = sqrt(G * M / r)
    velocities = np.sqrt(G * mass_central / radius)
    vx_positions = -velocities * np.sin(angles)
    vy_positions = velocities * np.cos(angles)
    vz_positions = np.zeros(num_bodies)  # No motion in the z-direction
    
    # Combine positions and velocities
    positions = np.array([x_positions, y_positions, z_positions]).T
    velocities = np.array([vx_positions, vy_positions, vz_positions]).T
    
    # The mass of the small bodies (arbitrary small mass for demonstration purposes)
    masses = np.full(num_bodies, 1e12)  # Small mass of each body
    
    return positions, velocities, masses

# Compute the total energy of the system (kinetic + potential)
def compute_energy(positions, velocities, masses, mass_central):
    kinetic_energy = 0.5 * np.sum(masses * np.linalg.norm(velocities, axis=1)**2)
    potential_energy = 0.0
    # Gravitational potential energy between each pair of bodies (including central body)
    for i in range(len(masses)):
        # Potential energy with the central body
        r = np.linalg.norm(positions[i])
        potential_energy -= G * mass_central * masses[i] / r
        # Inter-body potential energy (optional, but usually not significant in large N-body simulations)
        for j in range(i + 1, len(masses)):
            r_ij = np.linalg.norm(positions[i] - positions[j])
            potential_energy -= G * masses[i] * masses[j] / r_ij
    
    total_energy = kinetic_energy + potential_energy
    return total_energy

# Compute the acceleration for a specific body due to the gravitational force
def compute_acceleration(i, positions, masses, mass_central):
    num_bodies = len(masses)
    acceleration = np.zeros(3)
    
    # Gravitational acceleration from the central body
    r = np.linalg.norm(positions[i])
    accel_central = -G * mass_central * positions[i] / r**3
    acceleration += accel_central
    
    # Gravitational accelerations due to other bodies
    for j in range(num_bodies):
        if i != j:
            r_ij = positions[i] - positions[j]
            r_ij_norm = np.linalg.norm(r_ij)
            accel_interbody = G * masses[j] * r_ij / r_ij_norm**3
            acceleration += accel_interbody
    
    return i, acceleration

# Perform the simulation using the kick-step method with parallelized force computation
def run_simulation(positions, velocities, masses, mass_central, num_steps, time_step):
    # Store energy before the simulation
    initial_energy = compute_energy(positions, velocities, masses, mass_central)
    
    # Perform the time integration
    with ProcessPoolExecutor() as executor:
        for step in range(num_steps):
            # Parallelize the acceleration computation
            futures = []
            for i in range(len(masses)):
                futures.append(executor.submit(compute_acceleration, i, positions, masses, mass_central))
            
            # Collect accelerations after all tasks are completed
            accelerations = np.zeros_like(positions)
            for future in as_completed(futures):
                i, acceleration = future.result()
                accelerations[i] = acceleration
            
            # First kick: Update velocities (half step)
            velocities += 0.5 * accelerations * time_step
            
            # Update positions: full step
            positions += velocities * time_step
            
            # Recompute accelerations at new positions
            futures = []
            for i in range(len(masses)):
                futures.append(executor.submit(compute_acceleration, i, positions, masses, mass_central))
            
            # Collect accelerations after all tasks are completed
            accelerations = np.zeros_like(positions)
            for future in as_completed(futures):
                i, acceleration = future.result()
                accelerations[i] = acceleration
            
            # Second kick: Update velocities (half step)
            velocities += 0.5 * accelerations * time_step
    
    # Compute energy after the simulation
    final_energy = compute_energy(positions, velocities, masses, mass_central)
    
    return initial_energy, final_energy, positions, velocities

# Main function to initialize and run the simulation
def main():
    positions, velocities, masses = initialize_bodies(num_bodies, mass_central)
    
    # Perform the simulation
    initial_energy, final_energy, positions, velocities = run_simulation(positions, velocities, masses, mass_central, num_steps, time_step)
    
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