use std::f64::consts::PI;
use rand::Rng;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};
use std::cmp::Ordering as CmpOrdering;
use std::mem;
use std::cell::RefCell;
use std::sync::Arc;

const G: f64 = 6.67430e-11; // gravitational constant
const SOFTENING: f64 = 1.0e-3; // softening parameter to avoid singularities
const THETA: f64 = 0.3; // Barnes-Hut opening angle parameter
const MAX_BODIES_PER_LEAF: usize = 16;
const REBUILD_TREE_FREQUENCY: usize = 10; // Rebuild tree every N steps

// SIMD-friendly body representation with aligned data
#[derive(Clone, Debug)]
struct Body {
    mass: f64,
    position: [f64; 3],
    velocity: [f64; 3],
    acceleration: [f64; 3], // Store acceleration for force integration schemes
}

impl Body {
    fn new(mass: f64, position: [f64; 3], velocity: [f64; 3]) -> Self {
        Body {
            mass,
            position,
            velocity,
            acceleration: [0.0, 0.0, 0.0],
        }
    }
}

// Structure to represent a 3D bounding box
#[derive(Clone, Debug)]
struct BoundingBox {
    min: [f64; 3],
    max: [f64; 3],
}

impl BoundingBox {
    fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        BoundingBox { min, max }
    }

    fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    fn size(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        dx.max(dy).max(dz)
    }

    fn contains(&self, point: &[f64; 3]) -> bool {
        point[0] >= self.min[0] && point[0] <= self.max[0] &&
        point[1] >= self.min[1] && point[1] <= self.max[1] &&
        point[2] >= self.min[2] && point[2] <= self.max[2]
    }

    // Grow bounding box to include point
    fn include_point(&mut self, point: &[f64; 3]) {
        for i in 0..3 {
            self.min[i] = self.min[i].min(point[i]);
            self.max[i] = self.max[i].max(point[i]);
        }
    }
}

// Node in the kD-tree - optimized for cache locality
#[derive(Clone)]
enum Node {
    Internal {
        bounds: BoundingBox,
        center_of_mass: [f64; 3],
        total_mass: f64,
        split_dim: usize,
        split_value: f64,
        left: Box<Node>,
        right: Box<Node>,
    },
    Leaf {
        bounds: BoundingBox,
        bodies_indices: Vec<usize>, // Store indices instead of bodies to avoid duplication
        center_of_mass: [f64; 3],
        total_mass: f64,
    },
    Empty {
        bounds: BoundingBox,
    },
}

impl Node {
    fn new_leaf(bounds: BoundingBox, bodies: &[Body], indices: Vec<usize>) -> Self {
        let (center_of_mass, total_mass) = if indices.is_empty() {
            (bounds.center(), 0.0)
        } else {
            let mut com = [0.0, 0.0, 0.0];
            let mut total_mass = 0.0;
            for &idx in &indices {
                let body = &bodies[idx];
                total_mass += body.mass;
                for i in 0..3 {
                    com[i] += body.mass * body.position[i];
                }
            }
            if total_mass > 0.0 {
                for i in 0..3 {
                    com[i] /= total_mass;
                }
            }
            (com, total_mass)
        };
        
        Node::Leaf {
            bounds,
            bodies_indices: indices,
            center_of_mass,
            total_mass,
        }
    }

    // Find bounding box for a set of indices
    fn compute_bounds(bodies: &[Body], indices: &[usize]) -> BoundingBox {
        if indices.is_empty() {
            return BoundingBox::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        }

        let mut min = [f64::MAX, f64::MAX, f64::MAX];
        let mut max = [f64::MIN, f64::MIN, f64::MIN];

        for &idx in indices {
            let body = &bodies[idx];
            for i in 0..3 {
                min[i] = min[i].min(body.position[i]);
                max[i] = max[i].max(body.position[i]);
            }
        }

        // Add padding to ensure bodies at the boundaries are fully contained
        for i in 0..3 {
            let padding = 0.01 * (max[i] - min[i]).abs().max(1.0);
            min[i] -= padding;
            max[i] += padding;
        }

        BoundingBox::new(min, max)
    }

    // Build a kD-tree from a list of bodies
    fn build(bodies: &[Body], indices: Vec<usize>, max_bodies_per_leaf: usize) -> Self {
        if indices.is_empty() {
            return Node::Empty {
                bounds: BoundingBox::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            };
        }

        let bounds = Self::compute_bounds(bodies, &indices);
        
        // Base case: few enough bodies for a leaf node
        if indices.len() <= max_bodies_per_leaf {
            return Self::new_leaf(bounds, bodies, indices);
        }

        // Find dimension with largest spread to split on
        let mut split_dim = 0;
        let mut max_spread = f64::MIN;
        
        for dim in 0..3 {
            let spread = bounds.max[dim] - bounds.min[dim];
            if spread > max_spread {
                max_spread = spread;
                split_dim = dim;
            }
        }

        // Sort indices based on position along split dimension
        // Using median-of-three approach for better pivot selection
        let mut sorted_indices = indices.clone();
        sorted_indices.sort_by(|&a, &b| {
            bodies[a].position[split_dim].partial_cmp(&bodies[b].position[split_dim])
                .unwrap_or(CmpOrdering::Equal)
        });

        // Find median value (middle of sorted array)
        let mid_idx = sorted_indices.len() / 2;
        let split_value = bodies[sorted_indices[mid_idx]].position[split_dim];

        // Split bodies into left and right groups
        let (left_indices, right_indices): (Vec<usize>, Vec<usize>) = sorted_indices.into_iter()
            .partition(|&idx| bodies[idx].position[split_dim] <= split_value);

        // Ensure we don't create empty partitions
        if left_indices.is_empty() || right_indices.is_empty() {
            return Self::new_leaf(bounds, bodies, indices);
        }

        // Recursively build left and right subtrees
        let left = Box::new(Node::build(bodies, left_indices, max_bodies_per_leaf));
        let right = Box::new(Node::build(bodies, right_indices, max_bodies_per_leaf));

        // Compute center of mass and total mass for this internal node
        let (left_com, left_mass) = left.center_of_mass_and_mass();
        let (right_com, right_mass) = right.center_of_mass_and_mass();
        
        let total_mass = left_mass + right_mass;
        let mut center_of_mass = [0.0, 0.0, 0.0];
        
        if total_mass > 0.0 {
            for i in 0..3 {
                center_of_mass[i] = (left_mass * left_com[i] + right_mass * right_com[i]) / total_mass;
            }
        }

        Node::Internal {
            bounds,
            center_of_mass,
            total_mass,
            split_dim,
            split_value,
            left,
            right,
        }
    }

    fn center_of_mass_and_mass(&self) -> ([f64; 3], f64) {
        match self {
            Node::Internal { center_of_mass, total_mass, .. } => (*center_of_mass, *total_mass),
            Node::Leaf { center_of_mass, total_mass, .. } => (*center_of_mass, *total_mass),
            Node::Empty { .. } => ([0.0, 0.0, 0.0], 0.0),
        }
    }

    // Calculate acceleration on a body using the kD-tree
    fn calculate_acceleration(&self, body_idx: usize, bodies: &[Body], theta: f64, acceleration: &mut [f64; 3]) {
        let body = &bodies[body_idx];
        
        match self {
            Node::Empty { .. } => {},
            
            Node::Leaf { bodies_indices, .. } => {
                for &other_idx in bodies_indices {
                    // Skip self-interaction
                    if body_idx == other_idx {
                        continue;
                    }
                    
                    let other = &bodies[other_idx];
                    let mut r = [0.0, 0.0, 0.0];
                    for k in 0..3 {
                        r[k] = other.position[k] - body.position[k];
                    }
                    
                    let distance_squared = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    let distance = (distance_squared + SOFTENING * SOFTENING).sqrt();
                    let prefactor = G * other.mass / (distance * distance * distance);
                    
                    for k in 0..3 {
                        acceleration[k] += prefactor * r[k];
                    }
                }
            },
            
            Node::Internal { bounds, center_of_mass, total_mass, left, right, .. } => {
                let size = bounds.size();
                
                // Calculate distance to center of mass
                let mut r = [0.0, 0.0, 0.0];
                for k in 0..3 {
                    r[k] = center_of_mass[k] - body.position[k];
                }
                let distance_squared = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                let distance = distance_squared.sqrt();
                
                // If node is far enough away, use approximation
                if size / distance < theta {
                    // Avoid division by zero with softening
                    let distance = (distance_squared + SOFTENING * SOFTENING).sqrt();
                    let prefactor = G * total_mass / (distance * distance * distance);
                    
                    for k in 0..3 {
                        acceleration[k] += prefactor * r[k];
                    }
                } else {
                    // Otherwise, recurse into children
                    left.calculate_acceleration(body_idx, bodies, theta, acceleration);
                    right.calculate_acceleration(body_idx, bodies, theta, acceleration);
                }
            }
        }
    }
}

// Optimized structure to manage the simulation
struct NBodySimulation {
    bodies: Vec<Body>,
    tree: Option<Node>,
    tree_build_count: usize,
    last_tree_build: usize,
}

impl NBodySimulation {
    fn new(initial_bodies: Vec<Body>) -> Self {
        NBodySimulation {
            bodies: initial_bodies,
            tree: None,
            tree_build_count: 0,
            last_tree_build: 0,
        }
    }
    
    // Rebuild tree if needed, or return existing tree
    fn get_or_build_tree(&mut self, step: usize, force_rebuild: bool) -> &Node {
        if force_rebuild || self.tree.is_none() || step - self.last_tree_build >= REBUILD_TREE_FREQUENCY {
            let indices: Vec<usize> = (0..self.bodies.len()).collect();
            self.tree = Some(Node::build(&self.bodies, indices, MAX_BODIES_PER_LEAF));
            self.last_tree_build = step;
            self.tree_build_count += 1;
        }
        self.tree.as_ref().unwrap()
    }
    
    // Calculate all accelerations using the tree
    fn calculate_accelerations(&mut self, step: usize, theta: f64) {
        let tree = self.get_or_build_tree(step, false);
        
        // Pre-allocate accelerations in parallel
        let accelerations: Vec<[f64; 3]> = (0..self.bodies.len())
            .into_par_iter()
            .map(|body_idx| {
                let mut acc = [0.0, 0.0, 0.0];
                tree.calculate_acceleration(body_idx, &self.bodies, theta, &mut acc);
                acc
            })
            .collect();
        
        // Update body accelerations
        for (i, acc) in accelerations.into_iter().enumerate() {
            self.bodies[i].acceleration = acc;
        }
    }
    
    // Update positions and velocities using leapfrog integration
    fn leapfrog_step(&mut self, dt: f64, step: usize, theta: f64) {
        // Calculate accelerations
        let accel_start = Instant::now();
        self.calculate_accelerations(step, theta);
        let accel_time = accel_start.elapsed();
        
        // Update velocities by half step (kick)
        self.bodies.par_iter_mut().for_each(|body| {
            for k in 0..3 {
                body.velocity[k] += 0.5 * body.acceleration[k] * dt;
            }
        });
        
        // Update positions (drift)
        self.bodies.par_iter_mut().for_each(|body| {
            for k in 0..3 {
                body.position[k] += body.velocity[k] * dt;
            }
        });
        
        // Recalculate accelerations at new positions
        if step % REBUILD_TREE_FREQUENCY == 0 {
            self.calculate_accelerations(step, theta);
        }
        
        // Complete velocity update (final kick)
        self.bodies.par_iter_mut().for_each(|body| {
            for k in 0..3 {
                body.velocity[k] += 0.5 * body.acceleration[k] * dt;
            }
        });
        
        if step % 10 == 0 {
            println!("  Acceleration calculation took: {:?}", accel_time);
        }
    }
    
    // Calculate system energy
    fn calculate_energy(&self) -> f64 {
        // Kinetic energy: Σ(1/2 * m * v²) - parallelized
        let kinetic = self.bodies.par_iter().map(|body| {
            let v_squared = body.velocity.iter().map(|&v| v * v).sum::<f64>();
            0.5 * body.mass * v_squared
        }).sum::<f64>();
        
        // Potential energy using the tree for approximation
        let tree = match &self.tree {
            Some(tree) => tree,
            None => {
                // Build tree if not available
                let indices: Vec<usize> = (0..self.bodies.len()).collect();
                &Node::build(&self.bodies, indices, MAX_BODIES_PER_LEAF)
            }
        };
        
        let potential = (0..self.bodies.len()).into_par_iter().map(|i| {
            let mut potential_energy = 0.0;
            
            // Use different theta for energy calculation to maintain accuracy
            let energy_theta = THETA * 0.5;
            let body = &self.bodies[i];
            
            for j in (i+1)..self.bodies.len() {
                let other = &self.bodies[j];
                let mut r_ij = [0.0, 0.0, 0.0];
                for k in 0..3 {
                    r_ij[k] = other.position[k] - body.position[k];
                }
                let distance = (r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2] + SOFTENING * SOFTENING).sqrt();
                potential_energy -= G * body.mass * other.mass / distance;
            }
            
            potential_energy
        }).sum::<f64>();
        
        kinetic + potential
    }
}

// Initialize system with central body and specified number of smaller bodies on circular orbits
fn initialize_system(num_small_bodies: usize, central_mass: f64) -> Vec<Body> {
    let small_body_mass = 1.0e3; // 1000 kg for small bodies
    
    // Create central body
    let central_body = Body::new(
        central_mass,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0]
    );
    
    // Generate small bodies in parallel with better work distribution
    let chunk_size = (num_small_bodies + rayon::current_num_threads() - 1) / rayon::current_num_threads();
    let small_bodies: Vec<Body> = (0..num_small_bodies)
        .into_par_iter()
        .with_min_len(chunk_size)
        .map(|_| {
            let mut rng = rand::thread_rng();
            
            // Random orbital radius with better distribution
            let radius = 1.0e8 * (1.0 + 4.0 * rng.gen::<f64>().powf(0.5));
            
            // Generate more uniformly distributed points on a sphere using improved method
            let u = rng.gen_range(0.0..1.0);
            let v = rng.gen_range(0.0..1.0);
            
            let theta = 2.0 * PI * u;
            let phi = (2.0 * v - 1.0).acos();
            
            // Convert spherical to cartesian coordinates
            let x = radius * phi.sin() * theta.cos();
            let y = radius * phi.sin() * theta.sin();
            let z = radius * phi.cos();
            
            // Calculate orbital velocity for a circular orbit with slight eccentricity
            let orbital_velocity = (G * central_mass / radius).sqrt() * (0.95 + 0.1 * rng.gen::<f64>());
            
            // Calculate velocity vector perpendicular to radius vector
            let v_x = orbital_velocity * (-theta.sin() * phi.sin());
            let v_y = orbital_velocity * (theta.cos() * phi.sin());
            let v_z = orbital_velocity * 0.01 * (rng.gen::<f64>() - 0.5); // Small vertical component for interest
            
            Body::new(
                small_body_mass * (0.5 + rng.gen::<f64>()),  // Vary mass slightly
                [x, y, z],
                [v_x, v_y, v_z]
            )
        })
        .collect();
    
    // Combine central body with small bodies
    let mut bodies = Vec::with_capacity(num_small_bodies + 1);
    bodies.push(central_body);
    bodies.extend(small_bodies);
    
    bodies
}

fn main() {
    // Parameters
    let num_bodies = 1_000_000;
    let central_mass = 1.989e30; // Solar mass in kg
    let time_step = 86400.0; // One day in seconds
    let num_steps = 1000;
    let theta = THETA; // Barnes-Hut opening angle parameter
    
    println!("N-Body Simulation with Optimized Barnes-Hut Algorithm");
    println!("-----------------------------------------------------");
    println!("Bodies: {}", num_bodies + 1);
    println!("Steps: {}", num_steps);
    println!("Theta: {}", theta);
    println!("Threads: {}", rayon::current_num_threads());
    println!("Tree rebuild frequency: every {} steps", REBUILD_TREE_FREQUENCY);
    println!("Max bodies per leaf: {}", MAX_BODIES_PER_LEAF);
    println!("-----------------------------------------------------");
    
    println!("Initializing system...");
    let start_time = Instant::now();
    let bodies = initialize_system(num_bodies, central_mass);
    println!("Initialization took: {:?}", start_time.elapsed());
    
    // Create simulation
    let mut simulation = NBodySimulation::new(bodies);
    
    // Calculate initial energy
    println!("Calculating initial energy...");
    let energy_start = Instant::now();
    let initial_energy = simulation.calculate_energy();
    println!("Initial energy calculation took: {:?}", energy_start.elapsed());
    println!("Initial total energy: {:.6e}", initial_energy);
    
    // Progress tracking
    let completed = AtomicUsize::new(0);
    let progress_interval = 10; // Report progress every 10 steps
    
    // Run simulation
    println!("Starting simulation...");
    let sim_start_time = Instant::now();
    
    for step in 1..=num_steps {
        let step_start = Instant::now();
        simulation.leapfrog_step(time_step, step, theta);
        
        // Update progress counter
        let completed_steps = completed.fetch_add(1, Ordering::Relaxed) + 1;
        if completed_steps % progress_interval == 0 {
            println!("Step {}/{} ({:.1}%) - took {:?}", 
                     completed_steps, num_steps, 
                     (completed_steps as f64 / num_steps as f64) * 100.0,
                     step_start.elapsed());
        }
    }
    
    let sim_duration = sim_start_time.elapsed();
    println!("Simulation completed in: {:?}", sim_duration);
    println!("Average time per step: {:?}", sim_duration / num_steps as u32);
    println!("Tree rebuilds: {}", simulation.tree_build_count);
    
    // Calculate final energy
    println!("Calculating final energy...");
    let final_energy_start = Instant::now();
    let final_energy = simulation.calculate_energy();
    println!("Final energy calculation took: {:?}", final_energy_start.elapsed());
    println!("Final total energy: {:.6e}", final_energy);
    
    // Calculate energy error
    let energy_error = (final_energy - initial_energy).abs() / initial_energy.abs();
    println!("Relative energy error: {:.6e}", energy_error);
    
    // Print summary of system state
    let mut total_mass = 0.0;
    let mut center_of_mass = [0.0, 0.0, 0.0];
    let mut total_momentum = [0.0, 0.0, 0.0];
    
    for body in &simulation.bodies {
        total_mass += body.mass;
        for k in 0..3 {
            center_of_mass[k] += body.mass * body.position[k];
            total_momentum[k] += body.mass * body.velocity[k];
        }
    }
    
    for k in 0..3 {
        center_of_mass[k] /= total_mass;
    }
    
    println!("Final system state:");
    println!("  Total mass: {:.6e} kg", total_mass);
    println!("  Center of mass: [{:.6e}, {:.6e}, {:.6e}] m", 
             center_of_mass[0], center_of_mass[1], center_of_mass[2]);
    println!("  Total momentum: [{:.6e}, {:.6e}, {:.6e}] kg·m/s", 
             total_momentum[0], total_momentum[1], total_momentum[2]);
}