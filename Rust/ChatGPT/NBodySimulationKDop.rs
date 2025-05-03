use rand::Rng;
use std::f64::consts::PI;
use rayon::prelude::*;
use std::cmp::Ordering;

#[derive(Clone, Copy, Debug)]
struct Body {
    position: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

const G: f64 = 6.67430e-11;
const DT: f64 = 1.0;
const THETA: f64 = 0.3;

impl Body {
    fn new(position: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Self { position, velocity, mass }
    }
}

#[derive(Debug)]
struct Node {
    center_of_mass: [f64; 3],
    total_mass: f64,
    boundary: [[f64; 3]; 2],
    body: Option<Body>,
    children: Vec<Node>,
}

impl Node {
    fn new(boundary: [[f64; 3]; 2]) -> Self {
        Self {
            center_of_mass: [0.0; 3],
            total_mass: 0.0,
            boundary,
            body: None,
            children: Vec::new(),
        }
    }

    fn insert(&mut self, body: Body) {
        if self.body.is_none() && self.children.is_empty() {
            self.body = Some(body);
            self.center_of_mass = body.position;
            self.total_mass = body.mass;
        } else {
            if self.children.is_empty() {
                self.subdivide();
            }
            let idx = self.find_quadrant(body.position);
            self.children[idx].insert(body);
            self.update_center_of_mass();
        }
    }

    fn subdivide(&mut self) {
        let min = self.boundary[0];
        let max = self.boundary[1];
        let mid = [(min[0] + max[0]) / 2.0, (min[1] + max[1]) / 2.0, (min[2] + max[2]) / 2.0];

        for i in 0..8 {
            let mut new_min = min;
            let mut new_max = mid;
            if i & 1 != 0 { new_min[0] = mid[0]; new_max[0] = max[0]; }
            if i & 2 != 0 { new_min[1] = mid[1]; new_max[1] = max[1]; }
            if i & 4 != 0 { new_min[2] = mid[2]; new_max[2] = max[2]; }
            self.children.push(Node::new([new_min, new_max]));
        }
    }

    fn find_quadrant(&self, position: [f64; 3]) -> usize {
        let mid = [(self.boundary[0][0] + self.boundary[1][0]) / 2.0,
                   (self.boundary[0][1] + self.boundary[1][1]) / 2.0,
                   (self.boundary[0][2] + self.boundary[1][2]) / 2.0];
        ((position[0] > mid[0]) as usize) +
        (((position[1] > mid[1]) as usize) << 1) +
        (((position[2] > mid[2]) as usize) << 2)
    }

    fn update_center_of_mass(&mut self) {
        let mut total_mass = 0.0;
        let mut weighted_position = [0.0; 3];
        for child in &self.children {
            total_mass += child.total_mass;
            for i in 0..3 {
                weighted_position[i] += child.center_of_mass[i] * child.total_mass;
            }
        }
        if total_mass > 0.0 {
            self.total_mass = total_mass;
            for i in 0..3 {
                self.center_of_mass[i] = weighted_position[i] / total_mass;
            }
        }
    }

    // Optimized force computation using parallelism for subtree traversal
    fn compute_force(&self, body: &Body) -> [f64; 3] {
        let mut force = [0.0; 3];

        // Early return for the single body case
        if let Some(b) = self.body {
            if b.position != body.position {
                let r_vec: [f64; 3] = [b.position[0] - body.position[0], b.position[1] - body.position[1], b.position[2] - body.position[2]];
                let r_mag = (r_vec.iter().map(|v| v * v).sum::<f64>()).sqrt();
                if r_mag > 0.0 {
                    let factor = G * b.mass / (r_mag * r_mag * r_mag);
                    for i in 0..3 {
                        force[i] += factor * r_vec[i];
                    }
                }
            }
        } else {
            let size = self.boundary[1][0] - self.boundary[0][0];
            let r_mag = ((0..3).map(|i| (self.center_of_mass[i] - body.position[i]).powi(2)).sum::<f64>()).sqrt();
            if size / r_mag < THETA {
                let factor = G * self.total_mass / (r_mag * r_mag * r_mag);
                for i in 0..3 {
                    force[i] += factor * (self.center_of_mass[i] - body.position[i]);
                }
            } else {
                // Use Rayon parallelism to compute forces on children nodes
                force = self.children.par_iter()
                    .map(|child| child.compute_force(body))
                    .fold(|| [0.0; 3], |mut acc, child_force| {
                        for i in 0..3 {
                            acc[i] += child_force[i];
                        }
                        acc
                    })
                    .reduce(|| [0.0; 3], |mut acc, child_force| {
                        for i in 0..3 {
                            acc[i] += child_force[i];
                        }
                        acc
                    });
            }
        }
        force
    }
}

fn main() {
    println!("KD-tree based N-body simulation with Barnes-Hut method");

    // Example initialization and usage
    let num_bodies = 1000;
    let mut bodies: Vec<Body> = Vec::with_capacity(num_bodies);
    let mut rng = rand::thread_rng();
    for _ in 0..num_bodies {
        let position = [rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)];
        let velocity = [rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)];
        let mass = rng.gen_range(1e10..1e15);
        bodies.push(Body::new(position, velocity, mass));
    }

    let boundary = [[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]];
    let mut root_node = Node::new(boundary);
    
    // Insert bodies into the KD-tree (root node)
    for body in &bodies {
        root_node.insert(*body);
    }

    // Example force computation (this is just for testing and should be done in a simulation loop)
    let body = bodies[0];
    let force = root_node.compute_force(&body);
    println!("{:?}", force);
}