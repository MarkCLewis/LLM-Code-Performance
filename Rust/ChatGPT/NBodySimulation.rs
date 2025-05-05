use rand::prelude::*;
use std::f64::consts::PI;

const G: f64 = 6.67430e-11;
const DT: f64 = 1.0;
const EPS2: f64 = 1e-10;
const N_BODIES: usize = 1_000_000;
const STEPS: usize = 1000;

#[derive(Clone, Copy, Debug)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    fn zero() -> Self {
        Vec3 { x: 0.0, y: 0.0, z: 0.0 }
    }

    fn norm_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    fn add(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    fn sub(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    fn scale(&self, s: f64) -> Vec3 {
        Vec3 {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

#[derive(Clone, Debug)]
struct Body {
    pos: Vec3,
    vel: Vec3,
    acc: Vec3,
    mass: f64,
}

fn initialize_system(n: usize) -> Vec<Body> {
    let mut bodies = Vec::with_capacity(n + 1);
    bodies.push(Body {
        pos: Vec3::zero(),
        vel: Vec3::zero(),
        acc: Vec3::zero(),
        mass: 1e20,
    });

    let mut rng = StdRng::seed_from_u64(42);
    let radius = 1e7;

    for i in 1..=n {
        let angle = 2.0 * PI * i as f64 / n as f64;
        let r = radius * (1.0 + 0.1 * rng.gen::<f64>());
        let x = r * angle.cos();
        let y = r * angle.sin();
        let z = 0.0;

        let v = (G * bodies[0].mass / r).sqrt();
        let vx = -v * angle.sin();
        let vy = v * angle.cos();
        let vz = 0.0;

        bodies.push(Body {
            pos: Vec3 { x, y, z },
            vel: Vec3 { x: vx, y: vy, z: vz },
            acc: Vec3::zero(),
            mass: 1.0,
        });
    }

    bodies
}

fn compute_forces(bodies: &mut Vec<Body>) {
    let n = bodies.len();
    for i in 0..n {
        bodies[i].acc = Vec3::zero();
        for j in 0..n {
            if i != j {
                let dx = bodies[j].pos.sub(&bodies[i].pos);
                let dist2 = dx.norm_squared() + EPS2;
                let dist = dist2.sqrt();
                let force = G * bodies[j].mass / (dist2 * dist);
                bodies[i].acc = bodies[i].acc.add(&dx.scale(force));
            }
        }
    }
}

fn update_bodies(bodies: &mut Vec<Body>) {
    for body in bodies.iter_mut() {
        body.vel = body.vel.add(&body.acc.scale(DT));
        body.pos = body.pos.add(&body.vel.scale(DT));
    }
}

fn compute_energy(bodies: &Vec<Body>) -> f64 {
    let mut kinetic = 0.0;
    let mut potential = 0.0;
    let n = bodies.len();
    for i in 0..n {
        let v2 = bodies[i].vel.norm_squared();
        kinetic += 0.5 * bodies[i].mass * v2;
        for j in (i + 1)..n {
            let dx = bodies[i].pos.sub(&bodies[j].pos);
            let dist = (dx.norm_squared() + EPS2).sqrt();
            potential -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }
    kinetic + potential
}

fn main() {
    println!("Initializing system...");
    let mut bodies = initialize_system(N_BODIES);

    println!("Calculating initial energy...");
    let initial_energy = compute_energy(&bodies);
    println!("Initial energy: {:.6e}", initial_energy);

    for step in 0..STEPS {
        compute_forces(&mut bodies);
        update_bodies(&mut bodies);
        if step % 100 == 0 {
            println!("Step {}", step);
        }
    }

    println!("Calculating final energy...");
    let final_energy = compute_energy(&bodies);
    println!("Final energy: {:.6e}", final_energy);
    println!("Energy difference: {:.6e}", (final_energy - initial_energy).abs());
}
