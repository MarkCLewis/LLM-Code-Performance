use rand::prelude::*;
use rayon::prelude::*;
use std::f64::consts::PI;

const G: f64 = 6.67430e-11;
const DT: f64 = 1.0;
const EPS2: f64 = 1e-10;
const THETA: f64 = 0.3;
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

#[derive(Clone)]
struct KDNode<'a> {
    min: Vec3,
    max: Vec3,
    cm: Vec3,
    mass: f64,
    left: Option<Box<KDNode<'a>>>,
    right: Option<Box<KDNode<'a>>>,
    indices: &'a [usize],
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

fn build_kdtree<'a>(bodies: &'a [Body], indices: &'a [usize], depth: usize) -> Option<Box<KDNode<'a>>> {
    if indices.is_empty() {
        return None;
    }

    let axis = depth % 3;
    let mut sorted_indices = indices.to_vec();
    sorted_indices.sort_by(|&i, &j| bodies[i].pos[axis].partial_cmp(&bodies[j].pos[axis]).unwrap());
    let mid = sorted_indices.len() / 2;

    let bounds = |f: fn(&Vec3) -> f64, cmp: fn(f64, f64) -> f64| -> f64 {
        sorted_indices.iter().map(|&i| f(&bodies[i].pos)).fold(f64::INFINITY, cmp)
    };

    let min = Vec3 {
        x: sorted_indices.iter().map(|&i| bodies[i].pos.x).fold(f64::INFINITY, f64::min),
        y: sorted_indices.iter().map(|&i| bodies[i].pos.y).fold(f64::INFINITY, f64::min),
        z: sorted_indices.iter().map(|&i| bodies[i].pos.z).fold(f64::INFINITY, f64::min),
    };
    let max = Vec3 {
        x: sorted_indices.iter().map(|&i| bodies[i].pos.x).fold(f64::NEG_INFINITY, f64::max),
        y: sorted_indices.iter().map(|&i| bodies[i].pos.y).fold(f64::NEG_INFINITY, f64::max),
        z: sorted_indices.iter().map(|&i| bodies[i].pos.z).fold(f64::NEG_INFINITY, f64::max),
    };
    let total_mass: f64 = sorted_indices.iter().map(|&i| bodies[i].mass).sum();
    let cm = sorted_indices
        .iter()
        .map(|&i| bodies[i].pos.scale(bodies[i].mass))
        .fold(Vec3::zero(), |a, b| a.add(&b))
        .scale(1.0 / total_mass);

    let left = build_kdtree(bodies, &sorted_indices[..mid], depth + 1);
    let right = build_kdtree(bodies, &sorted_indices[mid..], depth + 1);

    Some(Box::new(KDNode {
        min,
        max,
        cm,
        mass: total_mass,
        left,
        right,
        indices: sorted_indices.as_slice(),
    }))
}

fn compute_force<'a>(b: &Body, node: &KDNode<'a>, bodies: &'a [Body]) -> Vec3 {
    if node.mass == 0.0 || node.indices.iter().any(|&i| std::ptr::eq(&bodies[i], b)) {
        return Vec3::zero();
    }

    let dx = node.cm.sub(&b.pos);
    let dist2 = dx.norm_squared() + EPS2;
    let size = (node.max.x - node.min.x).max((node.max.y - node.min.y).max(node.max.z - node.min.z));

    if node.indices.len() <= 1 || (size * size / dist2 < THETA * THETA) {
        let dist = dist2.sqrt();
        let force = G * node.mass / (dist2 * dist);
        dx.scale(force)
    } else {
        let left = node.left.as_ref().map(|l| compute_force(b, l, bodies)).unwrap_or(Vec3::zero());
        let right = node.right.as_ref().map(|r| compute_force(b, r, bodies)).unwrap_or(Vec3::zero());
        left.add(&right)
    }
}

fn compute_forces(bodies: &mut Vec<Body>) {
    let indices: Vec<usize> = (0..bodies.len()).collect();
    let tree = build_kdtree(bodies, &indices, 0).unwrap();
    bodies.par_iter_mut().for_each(|b| {
        b.acc = compute_force(b, &tree, bodies);
    });
}

fn update_bodies(bodies: &mut Vec<Body>) {
    bodies.par_iter_mut().for_each(|body| {
        body.vel = body.vel.add(&body.acc.scale(DT));
        body.pos = body.pos.add(&body.vel.scale(DT));
    });
}

fn compute_energy(bodies: &Vec<Body>) -> f64 {
    let kinetic: f64 = bodies.par_iter().map(|b| 0.5 * b.mass * b.vel.norm_squared()).sum();
    let potential: f64 = (0..bodies.len()).into_par_iter().map(|i| {
        let mut pot = 0.0;
        let bi = &bodies[i];
        for j in (i + 1)..bodies.len() {
            let bj = &bodies[j];
            let dx = bi.pos.sub(&bj.pos);
            let dist = (dx.norm_squared() + EPS2).sqrt();
            pot -= G * bi.mass * bj.mass / dist;
        }
        pot
    }).sum();
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

use std::ops::Index;
impl Index<usize> for Vec3 {
    type Output = f64;
    fn index(&self, i: usize) -> &Self::Output {
        match i {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds"),
        }
    }
}