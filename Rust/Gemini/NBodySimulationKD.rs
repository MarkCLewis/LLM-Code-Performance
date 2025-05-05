use rayon::prelude::*;
use std::f64::consts::PI;
use std::time::Instant;

const GRAVITATIONAL_CONSTANT: f64 = 6.67430e-11;
const THETA: f64 = 0.3;

#[derive(Clone, Debug)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    fn zero() -> Self {
        Vec3::new(0.0, 0.0, 0.0)
    }

    fn distance_squared(&self, other: &Self) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)
    }

    fn distance(&self, other: &Self) -> f64 {
        self.distance_squared(other).sqrt()
    }

    fn sub(&self, other: &Self) -> Self {
        Vec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    fn add(&self, other: &Self) -> Self {
        Vec3::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    fn mul(&self, scalar: f64) -> Self {
        Vec3::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }

    fn div(&self, scalar: f64) -> Self {
        Vec3::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }

    fn norm(&self) -> f64 {
        self.distance(&Vec3::zero())
    }
}

#[derive(Clone, Debug)]
struct Body {
    mass: f64,
    position: Vec3,
    velocity: Vec3,
    force: Vec3,
}

impl Body {
    fn new(mass: f64, position: Vec3, velocity: Vec3) -> Self {
        Body {
            mass,
            position,
            velocity,
            force: Vec3::zero(),
        }
    }
}

struct System {
    bodies: Vec<Body>,
}

impl System {
    fn new() -> Self {
        System { bodies: Vec::new() }
    }

    fn add_body(&mut self, body: Body) {
        self.bodies.push(body);
    }

    fn calculate_force_element(body_i: &mut Body, mass_element_cm: &Vec3, mass_element_mass: f64) {
        let r_vec = mass_element_cm.sub(&body_i.position);
        let r_sq = body_i.position.distance_squared(mass_element_cm);
        if r_sq > 1e-9 {
            let r = r_sq.sqrt();
            let magnitude = (GRAVITATIONAL_CONSTANT * body_i.mass * mass_element_mass) / r_sq;
            body_i.force.x += r_vec.x * magnitude / r;
            body_i.force.y += r_vec.y * magnitude / r;
            body_i.force.z += r_vec.z * magnitude / r;
        }
    }

    fn calculate_total_energy(&self) -> f64 {
        let kinetic_energy: f64 = self
            .bodies
            .par_iter()
            .map(|body| 0.5 * body.mass * body.velocity.distance_squared(&Vec3::zero()))
            .sum();

        let potential_energy = self
            .bodies
            .par_iter()
            .enumerate()
            .map(|(i, body_i)| {
                self.bodies
                    .iter()
                    .skip(i + 1)
                    .map(|body_j| {
                        let r = body_i.position.distance(&body_j.position);
                        -(GRAVITATIONAL_CONSTANT * body_i.mass * body_j.mass) / r
                    })
                    .sum::<f64>()
            })
            .sum();

        kinetic_energy + potential_energy
    }

    fn initialize_circular_orbits(
        num_orbiting: usize,
        central_mass: f64,
        orbit_radius: f64,
        orbiting_mass: f64,
    ) -> Self {
        let mut system = System::new();

        // Initialize the central body
        system.add_body(Body::new(central_mass, Vec3::zero(), Vec3::zero()));

        // Initialize the orbiting bodies
        for i in 0..num_orbiting {
            let angle = 2.0 * PI * (i as f64) / (num_orbiting as f64);
            let x = orbit_radius * angle.cos();
            let y = orbit_radius * angle.sin();
            let z = 0.0;
            let position = Vec3::new(x, y, z);

            // Calculate the orbital velocity for a circular orbit
            let orbital_speed = (GRAVITATIONAL_CONSTANT * central_mass / orbit_radius).sqrt();
            let vx = -orbital_speed * angle.sin();
            let vy = orbital_speed * angle.cos();
            let vz = 0.0;
            let velocity = Vec3::new(vx, vy, vz);

            system.add_body(Body::new(orbiting_mass, position, velocity));
        }
        system
    }

    // kD-tree Node
    #[derive(Clone, Debug)]
    struct KDNode {
        body_index: Option<usize>, // None for internal nodes
        center_of_mass: Vec3,
        total_mass: f64,
        min_bound: Vec3,
        max_bound: Vec3,
        left: Option<Box<KDNode>>,
        right: Option<Box<KDNode>>,
    }

    fn build_kd_tree(system: &System) -> Option<Box<KDNode>> {
        let num_bodies = system.bodies.len();
        if num_bodies == 0 {
            return None;
        }

        let mut indices: Vec<usize> = (0..num_bodies).collect();
        let mut min_bound = Vec3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut max_bound = Vec3::new(f64::MIN, f64::MIN, f64::MIN);

        for body in &system.bodies {
            min_bound.x = min_bound.x.min(body.position.x);
            min_bound.y = min_bound.y.min(body.position.y);
            min_bound.z = min_bound.z.min(body.position.z);
            max_bound.x = max_bound.x.max(body.position.x);
            max_bound.y = max_bound.y.max(body.position.y);
            max_bound.z = max_bound.z.max(body.position.z);
        }

        System::build_kd_tree_recursive(system, &mut indices, min_bound, max_bound, 0)
    }

    fn build_kd_tree_recursive(
        system: &System,
        indices: &mut [usize],
        min_bound: Vec3,
        max_bound: Vec3,
        depth: usize,
    ) -> Option<Box<KDNode>> {
        let num_indices = indices.len();

        if num_indices == 0 {
            return None;
        }

        if num_indices == 1 {
            let index = indices[0];
            let body = &system.bodies[index];
            return Some(Box::new(KDNode {
                body_index: Some(index),
                center_of_mass: body.position.clone(),
                total_mass: body.mass,
                min_bound,
                max_bound,
                left: None,
                right: None,
            }));
        }

        let mut node = KDNode {
            body_index: None,
            center_of_mass: Vec3::zero(),
            total_mass: 0.0,
            min_bound: min_bound.clone(),
            max_bound: max_bound.clone(),
            left: None,
            right: None,
        };

        for &index in indices.iter() {
            let body = &system.bodies[index];
            node.total_mass += body.mass;
            node.center_of_mass.x += body.mass * body.position.x;
            node.center_of_mass.y += body.mass * body.position.y;
            node.center_of_mass.z += body.mass * body.position.z;
        }

        node.center_of_mass.x /= node.total_mass;
        node.center_of_mass.y /= node.total_mass;
        node.center_of_mass.z /= node.total_mass;

        let split_dim = depth % 3;
        indices.sort_by(|&a, &b| {
            let pos_a = &system.bodies[a].position;
            let pos_b = &system.bodies[b].position;
            match split_dim {
                0 => pos_a.x.partial_cmp(&pos_b.x).unwrap(),
                1 => pos_a.y.partial_cmp(&pos_b.y).unwrap(),
                2 => pos_a.z.partial_cmp(&pos_b.z).unwrap(),
                _ => unreachable!(),
            }
        });

        let median_index = num_indices / 2;
        let left_indices = &mut indices[0..median_index];
        let right_indices = &mut indices[median_index..num_indices];

        let mut left_max_bound = max_bound.clone();
        let mut right_min_bound = min_bound.clone();

        match split_dim {
            0 => {
                left_max_bound.x = system.bodies[left_indices.last().unwrap()].position.x;
                right_min_bound.x = system.bodies[right_indices.first().unwrap()].position.x;
            }
            1 => {
                left_max_bound.y = system.bodies[left_indices.last().unwrap()].position.y;
                right_min_bound.y = system.bodies[right_indices.first().unwrap()].position.y;
            }
            2 => {
                left_max_bound.z = system.bodies[left_indices.last().unwrap()].position.z;
                right_min_bound.z = system.bodies[right_indices.first().unwrap()].position.z;
            }
            _ => unreachable!(),
        }

        node.left = System::build_kd_tree_recursive(system, left_indices, min_bound.clone(), left_max_bound, depth + 1);
        node.right = System::build_kd_tree_recursive(system, right_indices, right_min_bound, max_bound.clone(), depth + 1);

        Some(Box::new(node))
    }

    fn calculate_force_kd_tree(body_i: &mut Body, node: &KDNode, system: &System) {  // Pass the system
        if let Some(index) = node.body_index {
            if index != system.bodies.as_ptr_range().start as usize { // Correctly compare body pointers
                System::calculate_force_element(body_i, &node.center_of_mass, node.total_mass);
            }
            return;
        }

        let s = node.max_bound.sub(&node.min_bound).norm();
        let d = body_i.position.distance(&node.center_of_mass);

        if d == 0.0 || s / d < THETA {
            System::calculate_force_element(body_i, &node.center_of_mass, node.total_mass);
        } else {
            if let Some(left) = &node.left {
                System::calculate_force_kd_tree(body_i, left, system); // Pass the system
            }
            if let Some(right) = &node.right {
                System::calculate_force_kd_tree(body_i, right, system); // Pass the system
            }
        }
    }

    fn kick_step(&mut self, dt: f64) {
        let kd_tree_option = System::build_kd_tree(self);

        // Reset forces in parallel
        self.bodies.par_iter_mut().for_each(|body| {
            body.force.x = 0.0;
            body.force.y = 0.0;
            body.force.z = 0.0;
        });

        // Calculate forces using the kD-tree in parallel
        if let Some(kd_tree) = &kd_tree_option {
            self.bodies.par_iter_mut().for_each(|body| {
                System::calculate_force_kd_tree(body, kd_tree, self); // Pass self (system)
            });
        }

        // Update velocities (kick) in parallel
        self.bodies.par_iter_mut().for_each(|body| {
            body.velocity.x += (body.force.x / body.mass) * (dt / 2.0);
            body.velocity.y += (body.force.y / body.mass) * (dt / 2.0);
            body.velocity.z += (body.force.z / body.mass) * (dt / 2.0);
        });
    }

    fn drift_step(&mut self, dt: f64) {
        self.bodies.par_iter_mut().for_each(|body| {
            body.position.x += body.velocity.x * dt;
            body.position.y += body.velocity.y * dt;
            body.position.z += body.velocity.z * dt;
        });
    }

    fn first_order_kick_step(&mut self, dt: f64) {
        self.kick_step(dt);
        self.drift_step(dt);
        self.kick_step(dt);
    }
}

fn main() {
    let num_orbiting_bodies = 1_000_000;
    let central_mass = 1.989e30; // Mass of the Sun (kg)
    let orbit_radius = 1.496e11; // 1 AU (m)
    let orbiting_mass = 5.972e24; // Mass of the Earth (kg)
    let num_steps = 1000;
    let time_step = 3600.0 * 24.0 * 7.0; // 1 week in seconds

    // Initialize the system
    let mut system =
        System::initialize_circular_orbits(num_orbiting_bodies, central_mass, orbit_radius, orbiting_mass);
    let initial_system = system.clone();

    println!("Initial number of bodies: {}", system.bodies.len());

    // Calculate initial energy
    let initial_energy = initial_system.calculate_total_energy();
    println!("Initial total energy: {:.e} J", initial_energy);

    println!("Running simulation for {} steps...", num_steps);
    let start_time = Instant::now();
    for step in 0..num_steps {
        system.first_order_kick_step(time_step);
        if (step + 1) % 100 == 0 {
            println!("Step {} completed.", step + 1);
        }
    }
    let elapsed_time = start_time.elapsed();
    println!("Simulation finished in {:.2} seconds.", elapsed_time.as_secs_f64());

    // Calculate final energy
    let final_energy = system.calculate_total_energy();
    println!("Final total energy: {:.e} J", final_energy);

    // Calculate the energy difference
    let energy_difference = (final_energy - initial_energy).abs();
    let relative_energy_difference = energy_difference / initial_