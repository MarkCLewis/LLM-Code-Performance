import numpy as np
from concurrent.futures import ProcessPoolExecutor

G = 6.67430e-11  # Gravitational constant
dt = 1e-3        # Time step
num_bodies = 1000000  # Number of small bodies

class Body:
    def __init__(self, mass, pos, vel):
        self.mass = mass
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.force = np.zeros(3, dtype=float)

def initialize_bodies(n):
    bodies = []

    # Central body
    bodies.append(Body(1e30, [0, 0, 0], [0, 0, 0]))

    # Small bodies in circular orbits
    for _ in range(n):
        radius = 1e11 * np.random.rand()
        speed = np.sqrt(G * bodies[0].mass / radius)
        angle = 2 * np.pi * np.random.rand()

        pos = [radius * np.cos(angle), radius * np.sin(angle), 0]
        vel = [-speed * np.sin(angle), speed * np.cos(angle), 0]
        bodies.append(Body(1e24, pos, vel))

    return bodies

def calculate_forces(bodies, start, end):
    for i in range(start, end):
        bodies[i].force.fill(0)
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].pos - bodies[i].pos
                dist = np.linalg.norm(dx)
                force = G * bodies[i].mass * bodies[j].mass / dist**3 * dx
                bodies[i].force += force
    return bodies[start:end]

def simulate_step(bodies):
    num_workers = 4  # Adjust based on your system's capabilities
    chunk_size = len(bodies) // num_workers

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(calculate_forces, bodies, i, min(i + chunk_size, len(bodies)))
                   for i in range(0, len(bodies), chunk_size)]
        results = [future.result() for future in futures]

    for i, result in enumerate(results):
        bodies[i * chunk_size:(i + 1) * chunk_size] = result

    # Update velocities and positions
    for body in bodies:
        body.vel += body.force / body.mass * dt
        body.pos += body.vel * dt

def calculate_energy(bodies):
    kinetic_energy = 0
    potential_energy = 0

    for i, body_i in enumerate(bodies):
        kinetic_energy += 0.5 * body_i.mass * np.dot(body_i.vel, body_i.vel)

        for j in range(i + 1, len(bodies)):
            body_j = bodies[j]
            dx = body_j.pos - body_i.pos
            dist = np.linalg.norm(dx)
            potential_energy -= G * body_i.mass * body_j.mass / dist

    return kinetic_energy + potential_energy

def main():
    np.random.seed(1234)

    # Initialize bodies
    bodies = initialize_bodies(num_bodies)

    # Calculate initial energy
    initial_energy = calculate_energy(bodies)
    print(f"Initial Energy: {initial_energy:.6e}")

    # Run simulation for 1000 steps
    for _ in range(1000):
        simulate_step(bodies)

    # Calculate final energy
    final_energy = calculate_energy(bodies)
    print(f"Final Energy: {final_energy:.6e}")

if __name__ == "__main__":
    main()