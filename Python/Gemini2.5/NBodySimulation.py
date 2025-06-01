import numpy as np
import time

def calculate_energy(mass, pos, vel):
    """
    Calculates the total kinetic and potential energy of the N-body system.

    Args:
        mass (np.ndarray): Array of particle masses.
        pos (np.ndarray): Array of particle positions (shape: [N, 3]).
        vel (np.ndarray): Array of particle velocities (shape: [N, 3]).

    Returns:
        tuple: A tuple containing the total kinetic energy and total potential energy.
    """
    # Kinetic energy: 0.5 * m * v^2
    ke = 0.5 * np.sum(mass[:, np.newaxis] * vel**2)

    # Potential energy: -G * m1 * m2 / r
    pe = 0
    # Calculate potential energy between all unique pairs of particles
    for i in range(pos.shape[0]):
        for j in range(i + 1, pos.shape[0]):
            r_vec = pos[j] - pos[i]
            r = np.linalg.norm(r_vec)
            if r > 0:
                pe -= (mass[i] * mass[j]) / r
    return ke, pe

def initialize_circular_orbits(n_bodies, central_mass_val):
    """
    Initializes a system with a central body and smaller bodies in random
    circular orbits around it.

    Args:
        n_bodies (int): The number of smaller bodies to create.
        central_mass_val (float): The mass of the central body.

    Returns:
        tuple: A tuple containing arrays for mass, position, and velocity of all bodies.
    """
    # Total number of bodies is n_bodies + 1 (for the central body)
    total_bodies = n_bodies + 1
    mass = np.ones(total_bodies)
    pos = np.zeros((total_bodies, 3))
    vel = np.zeros((total_bodies, 3))

    # Central body properties
    mass[0] = central_mass_val
    pos[0] = [0, 0, 0]
    vel[0] = [0, 0, 0]

    # Create orbiting bodies
    for i in range(1, total_bodies):
        # Generate a random radius for the orbit
        radius = np.random.uniform(5.0, 150.0)
        
        # Generate a random 3D direction vector for the position
        theta = np.arccos(2 * np.random.random() - 1)
        phi = 2 * np.pi * np.random.random()
        pos[i] = [
            radius * np.sin(theta) * np.cos(phi),
            radius * np.sin(theta) * np.sin(phi),
            radius * np.cos(theta)
        ]

        # Calculate the velocity for a circular orbit: v = sqrt(G*M/r)
        # We assume G=1 for simplicity
        orbital_speed = np.sqrt(central_mass_val / radius)
        
        # To get a velocity vector perpendicular to the position vector,
        # we take the cross product of the position vector and a random vector.
        random_vec = np.random.randn(3)
        # Ensure the random vector is not parallel to the position vector
        while np.allclose(np.cross(pos[i], random_vec), 0):
            random_vec = np.random.randn(3)

        vel[i] = np.cross(pos[i], random_vec)
        vel[i] = vel[i] / np.linalg.norm(vel[i]) * orbital_speed

    return mass, pos, vel

def main():
    """
    Main function to run the N-body simulation.
    """
    # Simulation parameters
    n_orbiting_bodies = 1_000_000
    central_body_mass = 1_000_000
    n_steps = 1000
    dt = 0.01  # Time step

    print("Initializing the system...")
    start_time = time.time()
    mass, pos, vel = initialize_circular_orbits(n_orbiting_bodies, central_body_mass)
    end_time = time.time()
    print(f"System initialization took: {end_time - start_time:.2f} seconds")

    print("\nCalculating initial energy...")
    start_time = time.time()
    ke_initial, pe_initial = calculate_energy(mass, pos, vel)
    total_energy_initial = ke_initial + pe_initial
    end_time = time.time()
    print(f"Energy calculation took: {end_time - start_time:.2f} seconds")
    print(f"Initial Total Energy: {total_energy_initial:.5e}")
    print(f"  - Initial Kinetic Energy: {ke_initial:.5e}")
    print(f"  - Initial Potential Energy: {pe_initial:.5e}")

    print(f"\nStarting simulation for {n_steps} steps...")
    sim_start_time = time.time()

    # Main simulation loop
    for step in range(n_steps):
        # Calculate accelerations (a = F/m = -G * m_j * r_ij / |r_ij|^3)
        accel = np.zeros_like(pos)
        for i in range(pos.shape[0]):
            for j in range(pos.shape[0]):
                if i != j:
                    r_vec = pos[j] - pos[i]
                    r_mag_sq = np.sum(r_vec**2)
                    r_mag = np.sqrt(r_mag_sq)
                    accel[i] += mass[j] * r_vec / (r_mag_sq * r_mag)

        # First-order kick-step method
        # Kick: Update velocities
        vel += accel * dt
        # Drift: Update positions
        pos += vel * dt

        # Optional: Print progress
        if (step + 1) % 100 == 0:
            print(f"  - Step {step + 1}/{n_steps} completed.")
            
    sim_end_time = time.time()
    print(f"Simulation loop took: {sim_end_time - sim_start_time:.2f} seconds")

    print("\nCalculating final energy...")
    start_time = time.time()
    ke_final, pe_final = calculate_energy(mass, pos, vel)
    total_energy_final = ke_final + pe_final
    end_time = time.time()
    print(f"Energy calculation took: {end_time - start_time:.2f} seconds")
    print(f"Final Total Energy:   {total_energy_final:.5e}")
    print(f"  - Final Kinetic Energy:   {ke_final:.5e}")
    print(f"  - Final Potential Energy: {pe_final:.5e}")

    # Verify energy conservation
    energy_change = (total_energy_final - total_energy_initial) / total_energy_initial
    print(f"\nRelative change in total energy: {energy_change:.5%}")

if __name__ == "__main__":
    main()