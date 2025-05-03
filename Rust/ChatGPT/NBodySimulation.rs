use rand::Rng;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
struct Body {
    position: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

const G: f64 = 6.67430e-11;
const DT: f64 = 1.0;

impl Body {
    fn new(position: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Self { position, velocity, mass }
    }
}

fn compute_energy(bodies: &[Body]) -> f64 {
    let mut energy = 0.0;
    let n = bodies.len();

    for i in 0..n {
        let bi = &bodies[i];
        let kinetic = 0.5 * bi.mass * (bi.velocity.iter().map(|v| v * v).sum::<f64>());
        energy += kinetic;
        
        for j in (i + 1)..n {
            let bj = &bodies[j];
            let r = ((0..3).map(|k| (bi.position[k] - bj.position[k]).powi(2)).sum::<f64>()).sqrt();
            if r > 0.0 {
                energy -= G * bi.mass * bj.mass / r;
            }
        }
    }
    energy
}

fn update_positions(bodies: &mut Vec<Body>) {
    for body in bodies.iter_mut() {
        for i in 0..3 {
            body.position[i] += body.velocity[i] * DT;
        }
    }
}

fn update_velocities(bodies: &mut Vec<Body>) {
    let n = bodies.len();
    let mut accelerations = vec![[0.0; 3]; n];
    
    for i in 0..n {
        let mut acc = [0.0; 3];
        for j in 0..n {
            if i != j {
                let r_vec: [f64; 3] = [bodies[i].position[0] - bodies[j].position[0], bodies[i].position[1] - bodies[j].position[1], bodies[i].position[2] - bodies[j].position[2]];
                let r_mag = (r_vec.iter().map(|v| v * v).sum::<f64>()).sqrt();
                if r_mag > 0.0 {
                    let factor = -G * bodies[j].mass / (r_mag * r_mag * r_mag);
                    for k in 0..3 {
                        acc[k] += factor * r_vec[k];
                    }
                }
            }
        }
        accelerations[i] = acc;
    }

    for i in 0..n {
        for k in 0..3 {
            bodies[i].velocity[k] += accelerations[i][k] * DT;
        }
    }
}

fn initialize_orbiting_bodies(n: usize, central_mass: f64) -> Vec<Body> {
    let mut bodies = Vec::with_capacity(n + 1);
    bodies.push(Body::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], central_mass));
    let mut rng = rand::thread_rng();
    let radius = 1.0e9;

    for i in 0..n {
        let angle = 2.0 * PI * (i as f64 / n as f64);
        let pos = [radius * angle.cos(), radius * angle.sin(), 0.0];
        let speed = (G * central_mass / radius).sqrt();
        let vel = [-speed * angle.sin(), speed * angle.cos(), 0.0];
        bodies.push(Body::new(pos, vel, 1.0));
    }
    bodies
}

fn main() {
    let n_bodies = 1_000_000;
    let central_mass = 1.989e30;
    let mut bodies = initialize_orbiting_bodies(n_bodies, central_mass);
    
    let initial_energy = compute_energy(&bodies);
    println!("Initial Energy: {:.5e}", initial_energy);

    for _ in 0..1000 {
        update_velocities(&mut bodies);
        update_positions(&mut bodies);
    }

    let final_energy = compute_energy(&bodies);
    println!("Final Energy: {:.5e}", final_energy);
}
