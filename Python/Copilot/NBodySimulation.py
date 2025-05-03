import numpy as np

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

def simulate_step(bodies):
    # Reset forces
    for body in bodies:
        body.force.fill(0)

    # Calculate forces
    for i, body_i in enumerate(bodies):
        for j in range(i + 1, len(bodies)):
            body_j = bodies[j]
            dx = body_j.pos - body_i.pos
            dist = np.linalg.norm(dx)
            force = G * body_i.mass * body_j.mass / dist**3 * dx

            body_i.force += force
            body_j.force -= force

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