use rayon::prelude::*;
use std::f64::consts::PI;

const G: f64 = 6.67430e-11; // Gravitational constant
const DT: f64 = 1e-3; // Time step
const NUM_BODIES: usize = 1000000; // Number of bodies
const THETA: f64 = 0.3; // Theta value for approximation

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

struct KDNode {
    body: Body,
    left: Option<Box<KDNode>>,
    right: Option<Box<KDNode>>,
    min: [f64; 3],
    max: [f64; 3],
}

struct KDTreePool {
    nodes: Vec<KDNode>,
    index: usize,
}

impl KDTreePool {
    fn new(capacity: usize) -> Self {
        KDTreePool {
            nodes: Vec::with_capacity(capacity),
            index: 0,
        }
    }

    fn allocate_node(&mut self, body: Body) -> &mut KDNode {
        if self.index >= self.nodes.capacity() {
            panic!("Pool capacity exceeded");
        }
        self.nodes.push(KDNode {
            body,
            left: None,
            right: None,
            min: [body.x, body.y, body.z],
            max: [body.x, body.y, body.z],
        });
        let node = &mut self.nodes[self.index];
        self.index += 1;
        node
    }
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

fn build_kdtree(bodies: &mut [Body], depth: usize, pool: &mut KDTreePool) -> Option<&mut KDNode> {
    if bodies.is_empty() {
        return None;
    }

    let axis = depth % 3;
    bodies.sort_by(|a, b| a[axis].partial_cmp(&b[axis]).unwrap());

    let median = bodies.len() / 2;
    let body = bodies[median];
    let node = pool.allocate_node(body);
    node.left = build_kdtree(&mut bodies[..median], depth + 1, pool).map(Box::new);
    node.right = build_kdtree(&mut bodies[median + 1..], depth + 1, pool).map(Box::new);

    for i in 0..3 {
        if let Some(ref left_node) = node.left {
            node.min[i] = node.min[i].min(left_node.min[i]);
            node.max[i] = node.max[i].max(left_node.max[i]);
        }
        if let Some(ref right_node) = node.right {
            node.min[i] = node.min[i].min(right_node.min[i]);
            node.max[i] = node.max[i].max(right_node.max[i]);
        }
    }

    Some(node)
}

fn calculate_force(node: &KDNode, body: &Body, ax: &mut f64, ay: &mut f64, az: &mut f64) {
    let dx = node.body.x - body.x;
    let dy = node.body.y - body.y;
    let dz = node.body.z - body.z;
    let distance = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();

    let size = node.max.iter().zip(node.min.iter()).map(|(max, min)| max - min).fold(0.0, f64::max);

    if size / distance < THETA || (node.left.is_none() && node.right.is_none()) {
        let force = G * node.body.mass / distance.powi(3);
        *ax += force * dx;
        *ay += force * dy;
        *az += force * dz;
    } else {
        if let Some(ref left) = node.left {
            calculate_force(left, body, ax, ay, az);
        }
        if let Some(ref right) = node.right {
            calculate_force(right, body, ax, ay, az);
        }
    }
}

fn kick_step(bodies: &mut [Body], root: &KDNode) {
    let accelerations: Vec<(f64, f64, f64)> = bodies.par_iter().map(|body| {
        let mut ax = 0.0;
        let mut ay = 0.0;
        let mut az = 0.0;
        calculate_force(root, body, &mut ax, &mut ay, &mut az);
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
    let mut pool = KDTreePool::new(NUM_BODIES * 2);

    let initial_energy = calculate_energy(&bodies);
    println!("Initial energy: {}", initial_energy);

    for _ in 0..1000 {
        pool.index = 0; // Reset pool index for reuse
        let root = build_kdtree(&mut bodies, 0, &mut pool).unwrap();
        kick_step(&mut bodies, &root);
    }

    let final_energy = calculate_energy(&bodies);
    println!("Final energy: {}", final_energy);
}
