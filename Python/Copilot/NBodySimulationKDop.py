import numpy as np
from concurrent.futures import ProcessPoolExecutor

G = 6.67430e-11  # Gravitational constant
dt = 1e-3        # Time step
num_bodies = 1000000  # Number of small bodies
theta = 0.3      # Barnes-Hut approximation parameter

class Body:
    def __init__(self, mass, pos, vel):
        self.mass = mass
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.force = np.zeros(3, dtype=float)

class Node:
    def __init__(self, center, size):
        self.body = None
        self.center = np.array(center, dtype=float)
        self.size = size
        self.mass = 0.0
        self.com = np.zeros(3, dtype=float)
        self.children = [None] * 8

class NodePool:
    def __init__(self, size):
        self.pool = [Node([0, 0, 0], 0) for _ in range(size)]
        self.index = 0

    def allocate(self, center, size):
        if self.index >= len(self.pool):
            self.pool.append(Node(center, size))
        else:
            self.pool[self.index].center = np.array(center, dtype=float)
            self.pool[self.index].size = size
            self.pool[self.index].mass = 0.0
            self.pool[self.index].com.fill(0)
            self.pool[self.index].body = None
            self.pool[self.index].children = [None] * 8
        self.index += 1
        return self.pool[self.index - 1]

    def reset(self):
        self.index = 0

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

def insert_body(node, body, pool):
    if node.body is None and node.mass == 0.0:
        node.body = body
        node.mass = body.mass
        node.com = body.pos.copy()
        return

    if node.body is not None:
        old_body = node.body
        node.body = None
        insert_body(node, old_body, pool)

    node.mass += body.mass
    node.com = (node.com * (node.mass - body.mass) + body.pos * body.mass) / node.mass

    index = 0
    for i in range(3):
        if body.pos[i] > node.center[i]:
            index |= 1 << i

    if node.children[index] is None:
        new_center = node.center.copy()
        for i in range(3):
            new_center[i] += (1 if index & (1 << i) else -1) * node.size / 4
        node.children[index] = pool.allocate(new_center, node.size / 2)

    insert_body(node.children[index], body, pool)

def build_tree(bodies, pool):
    root = pool.allocate([0, 0, 0], 2 * 1e11)
    for body in bodies:
        insert_body(root, body, pool)
    return root

def calculate_force(body, node):
    if node.body is not None and node.body is not body:
        dx = node.body.pos - body.pos
        dist = np.linalg.norm(dx)
        force = G * body.mass * node.body.mass / dist**3 * dx
        body.force += force
        return

    dx = node.com - body.pos
    dist = np.linalg.norm(dx)

    if node.size / dist < theta:
        force = G * body.mass * node.mass / dist**3 * dx
        body.force += force
        return

    for child in node.children:
        if child is not None:
            calculate_force(body, child)

def calculate_forces(bodies, start, end, root):
    for i in range(start, end):
        bodies[i].force.fill(0)
        calculate_force(bodies[i], root)
    return bodies[start:end]

def simulate_step(bodies, pool):
    num_workers = 4  # Adjust based on your system's capabilities
    chunk_size = len(bodies) // num_workers

    pool.reset()
    root = build_tree(bodies, pool)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(calculate_forces, bodies, i, min(i + chunk_size, len(bodies)), root)
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

    # Create a node pool
    pool = NodePool(num_bodies)

    # Calculate initial energy
    initial_energy = calculate_energy(bodies)
    print(f"Initial Energy: {initial_energy:.6e}")

    # Run simulation for 1000 steps
    for _ in range(1000):
        simulate_step(bodies, pool)

    # Calculate final energy
    final_energy = calculate_energy(bodies)
    print(f"Final Energy: {final_energy:.6e}")

if __name__ == "__main__":
    main()