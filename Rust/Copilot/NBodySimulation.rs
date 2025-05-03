use rand::Rng;
use std::f64::consts::PI;

const G: f64 = 6.67430e-11; // Gravitational constant
const DT: f64 = 1e-3;       // Time step
const NUM_BODIES: usize = 1_000_000; // Number of small bodies

#[derive(Clone)]
struct Body {
    mass: f64,
    pos: [f64; 3],
    vel: [f64; 3],
    force: [f64; 3],
}

impl Body {
    fn new(mass: f64, pos: [f64; 3], vel: [f64; 3]) -> Self {
        Self {
            mass,
            pos,
            vel,
            force: [0.0, 0.0, 0.0],
        }
    }
}

fn initialize_bodies(n: usize) -> Vec<Body> {
    let mut bodies = Vec::with_capacity(n + 1);
    let mut rng = rand::thread_rng();

    // Central body
    bodies.push(Body::new(1e30, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]));

    // Small bodies in circular orbits
    for _ in 0..n {
        let radius = 1e11 * rng.gen::<f64>();
        let speed = (G * bodies[0].mass / radius).sqrt();
        let angle = 2.0 * PI * rng.gen::<f64>();

        bodies.push(Body::new(
            1e24,
            [radius * angle.cos(), radius * angle.sin(), 0.0],
            [-speed * angle.sin(), speed * angle.cos(), 0.0],
        ));
    }

    bodies
}

fn simulate_step(bodies: &mut [Body]) {
    // Reset forces
    for body in bodies.iter_mut() {
        body.force = [0.0, 0.0, 0.0];
    }

    // Calculate forces
    for i in 0..bodies.len() {
        for j in i + 1..bodies.len() {
            let dx = bodies[j].pos[0] - bodies[i].pos[0];
            let dy = bodies[j].pos[1] - bodies[i].pos[1];
            let dz = bodies[j].pos[2] - bodies[i].pos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let force = G * bodies[i].mass * bodies[j].mass / (dist * dist * dist);

            bodies[i].force[0] += force * dx;
            bodies[i].force[1] += force * dy;
            bodies[i].force[2] += force * dz;

            bodies[j].force[0] -= force * dx;
            bodies[j].force[1] -= force * dy;
            bodies[j].force[2] -= force * dz;
        }
    }

    // Update velocities and positions
    for body in bodies.iter_mut() {
        body.vel[0] += body.force[0] / body.mass * DT;
        body.vel[1] += body.force[1] / body.mass * DT;
        body.vel[2] += body.force[2] / body.mass * DT;

        body.pos[0] += body.vel[0] * DT;
        body.pos[1] += body.vel[1] * DT;
        body.pos[2] += body.vel[2] * DT;
    }
}

fn calculate_energy(bodies: &[Body]) -> f64 {
    let mut kinetic_energy = 0.0;
    let mut potential_energy = 0.0;

    for i in 0..bodies.len() {
        kinetic_energy += 0.5 * bodies[i].mass * (bodies[i].vel[0] * bodies[i].vel[0] + bodies[i].vel[1] * bodies[i].vel[1] + bodies[i].vel[2] * bodies[i].vel[2]);

        for j in i + 1..bodies.len() {
            let dx = bodies[j].pos[0] - bodies[i].pos[0];
            let dy = bodies[j].pos[1] - bodies[i].pos[1];
            let dz = bodies[j].pos[2] - bodies[i].pos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            potential_energy -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }

    kinetic_energy + potential_energy
}

fn main() {
    // Initialize bodies
    let mut bodies = initialize_bodies(NUM_BODIES);

    // Calculate initial energy
    let initial_energy = calculate_energy(&bodies);
    println!("Initial Energy: {}", initial_energy);

    // Run simulation for 1000 steps
    for _ in 0..1000 {
        simulate_step(&mut bodies);
    }

    // Calculate final energy
    let final_energy = calculate_energy(&bodies);
    println!("Final Energy: {}", final_energy);
}