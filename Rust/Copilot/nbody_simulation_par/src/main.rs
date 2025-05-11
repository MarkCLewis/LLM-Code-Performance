use rayon::prelude::*;
use std::f64::consts::PI;

const G: f64 = 6.67430e-11; // Gravitational constant
const DT: f64 = 1e-3 * 360.0 * 24.0 * 365.0; // Time step
const NUM_BODIES: usize = 10000; // Number of bodies

#[derive(Clone, Copy)]
struct Body {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    mass: f64,
}

fn initialize_bodies() -> Vec<Body> {
    let mut bodies = Vec::with_capacity(NUM_BODIES);
    bodies.push(Body {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: 0.0,
        vz: 0.0,
        mass: 1e30, // Central body mass
    });

    for i in 1..NUM_BODIES {
        let angle = 2.0 * PI * i as f64 / (NUM_BODIES - 1) as f64;
        bodies.push(Body {
            x: angle.cos() * 1e11,
            y: angle.sin() * 1e11,
            z: 0.0,
            vx: -angle.sin() * (G * bodies[0].mass / 1e11).sqrt(),
            vy: angle.cos() * (G * bodies[0].mass / 1e11).sqrt(),
            vz: 0.0,
            mass: 1e24, // Small body mass
        });
    }
    bodies
}

fn calculate_energy(bodies: &[Body]) -> f64 {
    bodies.par_iter().map(|body| {
        let kinetic = 0.5 * body.mass * (body.vx.powi(2) + body.vy.powi(2) + body.vz.powi(2));
        let potential: f64 = bodies.par_iter().filter(|&&other| body != &other).map(|other| {
            let dx = body.x - other.x;
            let dy = body.y - other.y;
            let dz = body.z - other.z;
            let distance = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
            -G * body.mass * other.mass / distance
        }).sum();
        kinetic + 0.5 * potential
    }).sum()
}

fn kick_step(bodies: &mut [Body]) {
    let accelerations: Vec<(f64, f64, f64)> = bodies.par_iter().map(|body| {
        let mut ax = 0.0;
        let mut ay = 0.0;
        let mut az = 0.0;
        for other in bodies.iter() {
            if body != other {
                let dx = other.x - body.x;
                let dy = other.y - body.y;
                let dz = other.z - body.z;
                let distance = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
                let force = G * other.mass / distance.powi(3);
                ax += force * dx;
                ay += force * dy;
                az += force * dz;
            }
        }
        (ax, ay, az)
    }).collect();

    bodies.par_iter_mut().zip(accelerations.par_iter()).for_each(|(body, &(ax, ay, az))| {
        body.vx += ax * DT;
        body.vy += ay * DT;
        body.vz += az * DT;
    });

    bodies.par_iter_mut().for_each(|body| {
        body.x += body.vx * DT;
        body.y += body.vy * DT;
        body.z += body.vz * DT;
    });
}

fn main() {
    let mut bodies = initialize_bodies();

    let initial_energy = calculate_energy(&bodies);
    println!("Initial energy: {}", initial_energy);

    for _ in 0..100 {
        kick_step(&mut bodies);
    }

    let final_energy = calculate_energy(&bodies);
    println!("Final energy: {}", final_energy);
}
