use rand::Rng;
use rayon::prelude::*;
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
    bodies.par_iter_mut().for_each(|body| {
        body.force = [0.0, 0.0, 0.0];
    });

    // Calculate forces
    bodies.par_iter_mut().enumerate().for_each(|(i, body_i)| {
        for j in i + 1..bodies.len() {
            let body_j = &mut bodies[j];
            let dx = body_j.pos[0] - body_i.pos[0];
            let dy = body_j.pos[1] - body_i.pos[1];
            let dz = body_j.pos[2] - body_i.pos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let force = G * body_i.mass * body_j.mass / (dist * dist * dist);

            body_i.force[0] += force * dx;
            body_i.force[1] += force * dy;
            body_i.force[2] += force * dz;

            body_j.force[0] -= force * dx;
            body_j.force[1] -= force * dy;
            body_j.force[2] -= force * dz;
        }
    });

    // Update velocities and positions
    bodies.par_iter_mut().for_each(|body| {
        body.vel[0] += body.force[0] / body.mass * DT;
        body.vel[1] += body.force[1] / body.mass * DT;
        body.vel[2] += body.force[2] / body.mass * DT;

        body.pos[0] += body.vel[0] * DT;
        body.pos[1] += body.vel[1] * DT;
        body.pos[2] += body.vel[2] * DT;
    });
}

fn calculate_energy(bodies: &[Body]) -> f64 {
    let kinetic_energy: f64 = bodies.par_iter().map(|body| {
        0.5 * body.mass * (body.vel[0] * body.vel[0] + body.vel[1] * body.vel[1] + body.vel[2] * body.vel[2])
    }).sum();

    let potential_energy: f64 = bodies.par_iter().enumerate().map(|(i, body_i)| {
        (i + 1..bodies.len()).map(|j| {
            let body_j = &bodies[j];
            let dx = body_j.pos[0] - body_i.pos[0];
            let dy = body_j.pos[1] - body_i.pos[1];
            let dz = body_j.pos[2] - body_i.pos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            -G * body_i.mass * body_j.mass / dist
        }).sum::<f64>()
    }).sum();

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