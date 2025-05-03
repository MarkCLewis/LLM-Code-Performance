use rayon::prelude::*;
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
    bodies.par_iter().enumerate().map(|(i, bi)| {
        let kinetic = 0.5 * bi.mass * (bi.velocity.iter().map(|v| v * v).sum::<f64>());
        let potential: f64 = bodies.iter().enumerate()
            .filter(|(j, _)| *j > i)
            .map(|(_, bj)| {
                let r = ((0..3).map(|k| (bi.position[k] - bj.position[k]).powi(2)).sum::<f64>()).sqrt();
                if r > 0.0 { -G * bi.mass * bj.mass / r } else { 0.0 }
            }).sum();
        kinetic + potential
    }).sum()
}

fn update_positions(bodies: &mut Vec<Body>) {
    bodies.par_iter_mut().for_each(|body| {
        for i in 0..3 {
            body.position[i] += body.velocity[i] * DT;
        }
    });
}

fn update_velocities(bodies: &mut Vec<Body>) {
    let n = bodies.len();
    let mut accelerations = vec![[0.0; 3]; n];
    
    bodies.par_iter().enumerate().for_each(|(i, bi)| {
        let mut acc = [0.0; 3];
        for (j, bj) in bodies.iter().enumerate() {
            if i != j {
                let r_vec: [f64; 3] = [bi.position[0] - bj.position[0], bi.position[1] - bj.position[1], bi.position[2] - bj.position[2]];
                let r_mag = (r_vec.iter().map(|v| v * v).sum::<f64>()).sqrt();
                if r_mag > 0.0 {
                    let factor = -G * bj.mass / (r_mag * r_mag * r_mag);
                    acc.iter_mut().zip(r_vec.iter()).for_each(|(a, r)| *a += factor * r);
                }
            }
        }
        accelerations[i] = acc;
    });

    bodies.par_iter_mut().enumerate().for_each(|(i, body)| {
        for k in 0..3 {
            body.velocity[k] += accelerations[i][k] * DT;
        }
    });
}

fn initialize_orbiting_bodies(n: usize, central_mass: f64) -> Vec<Body> {
    let mut bodies = Vec::with_capacity(n + 1);
    bodies.push(Body::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], central_mass));
    let radius = 1.0e9;

    bodies.par_extend((0..n).into_par_iter().map(|i| {
        let angle = 2.0 * PI * (i as f64 / n as f64);
        let pos = [radius * angle.cos(), radius * angle.sin(), 0.0];
        let speed = (G * central_mass / radius).sqrt();
        let vel = [-speed * angle.sin(), speed * angle.cos(), 0.0];
        Body::new(pos, vel, 1.0)
    }));
    
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
