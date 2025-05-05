use std::f64::consts::PI;
use rand::Rng;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};
use std::cmp::Ordering as CmpOrdering;

const G: f64 = 6.67430e-11; // gravitational constant
const SOFTENING: f64 = 1.0e-3; // softening parameter to avoid singularities
const THETA: f64 = 0.3; // Barnes-Hut opening angle parameter

#[derive(Clone, Debug)]
struct Body {
    mass: f64,
    position: [f64; 3],
    velocity: [f64; 3],
}

impl Body {
    fn new(mass: f64, position: [f64; 3], velocity: [f64; 3]) -> Self {
        Body {
            mass,
            position,
            velocity,
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
}

// Node in the kD-tree
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
        bodies: Vec<Body>,
        center_of_mass: [f64; 3],
        total_mass: f64,
    },
    Empty {
        bounds: BoundingBox,
    },
}

impl Node {
    fn new_leaf(bounds: BoundingBox, bodies: Vec<Body>) -> Self {
        let (center_of_mass, total_mass) = if bodies.is_empty() {
            (bounds.center(), 0.0)
        } else {
            let mut com = [0.0, 0.0, 0.0];
            let mut total_mass = 0.0;
            for body in &bodies {
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
            bodies,
            center_of_mass,
            total_mass,
        }
    }

    // Build a kD-tree from a list of bodies
    fn build(bodies: &[Body], max_bodies_per_leaf: usize) -> Self {
        if bodies.is_empty() {
            return Node::Empty {
                bounds: BoundingBox::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            };
        }

        // Find bounding box for all bodies
        let mut min = [f64::MAX, f64::MAX, f64::MAX];
        let mut max = [f64::MIN, f64::MIN, f64::MIN];

        for body in bodies {
            for i in 0..3 {
                min[i] = min[i].min(body.position[i]);
                max[i] = max[i].max(body.position[i]);
            }
        }

        // Add some padding to ensure bodies at the boundaries are fully contained
        for i in 0..3 {
            let padding = 0.01 * (max[i] - min[i]).abs().max(1.0);
            min[i] -= padding;
            max[i] += padding;
        }

        let bounds = BoundingBox::new(min, max);
        
        // Base case: few enough bodies for a leaf node
        if bodies.len() <= max_bodies_per_leaf {
            return Node::new_leaf(bounds, bodies.to_vec());
        }

        // Find dimension with largest spread to split on
        let mut split_dim = 0;
        let mut max_spread = f64::MIN;
        
        for dim in 0..3 {
            let spread = max[dim] - min[dim];
            if spread > max_spread {
                max_spread = spread;
                split_dim = dim;
            }
        }

        // Find median value along that dimension
        let mut values: Vec<f64> = bodies.iter().map(|b| b.position[split_dim]).collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(CmpOrdering::Equal));
        let split_value = if values.len() % 2 == 0 {
            (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
        } else {
            values[values.len() / 2]
        };

        // Split bodies into left and right groups
        let (left_bodies, right_bodies): (Vec<Body>, Vec<Body>) = bodies.iter()
            .cloned()
            .partition(|b| b.position[split_dim] <= split_value);

        // Ensure we don't create empty partitions
        if left_bodies.is_empty() || right_bodies.is_empty() {
            return Node::new_leaf(bounds, bodies.to_vec());
        }

        // Recursively build left and right subtrees
        let left = Box::new(Node::build(&left_bodies, max_bodies_per_leaf));
        let right = Box::new(Node::build(&right_bodies, max_bodies_per_leaf));

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
    fn calculate_acceleration(&self, body: &Body, theta: f64) -> [f64; 3] {
        match self {
            Node::Empty { .. } => [0.0, 0.0, 0.0],
            
            Node::Leaf { bodies, .. } => {
                let mut acc = [0.0, 0.0, 0.0];
                for other in bodies {
                    // Skip self-interaction
                    if body.position == other.position {
                        continue;
                    }
                    
                    let mut r = [0.0, 0.0, 0.0];
                    for k in 0..3 {
                        r[k] = other.position[k] - body.position[k];
                    }
                    
                    let distance_squared = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    let distance = (distance_squared + SOFTENING * SOFTENING).sqrt();
                    let prefactor = G * other.mass / (distance * distance * distance);
                    
                    for k in 0..3 {
                        acc[k] += prefactor * r[k];
                    }
                }
                acc
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
                    let distance = (distance_squared + SOFTENING * SOFTENING).sqrt();
                    let prefactor = G * total_mass / (distance * distance * distance);
                    
                    [prefactor * r[0], prefactor * r[1], prefactor * r[2]]
                } else {
                    // Otherwise, recurse into children
                    let left_acc = left.calculate_acceleration(body, theta);
                    let right_acc = right.calculate_acceleration(body, theta);
                    
                    [
                        left_acc[0] + right_acc[0],
                        left_acc[1] + right_acc[1],
                        left_acc[2] + right_acc[2],
                    ]
                }
            }
        }
    }
}

// Calculate total system energy (kinetic + potential)
fn calculate_energy(bodies: &[Body]) -> f64 {
    // Kinetic energy: Σ(1/2 * m * v²) - parallelized
    let kinetic = bodies.par_iter().map(|body| {
        let v_squared = body.velocity.iter().map(|&v| v * v).sum::<f64>();
        0.5 * body.mass * v_squared
    }).sum::<f64>();
    
    // Potential energy: -Σ(G * m1 * m2 / r) - parallelized
    let potential = (0..bodies.len()).into_par_iter().map(|i| {
        ((i+1)..bodies.len()).map(|j| {
            let mut r_ij = [0.0; 3];
            for k in 0..3 {
                r_ij[k] = bodies[j].position[k] - bodies[i].position[k];
            }
            let distance = (r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2] + SOFTENING * SOFTENING).sqrt();
            -G * bodies[i].mass * bodies[j].mass / distance
        }).sum::<f64>()
    }).sum::<f64>();
    
    kinetic + potential
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
    
    // Generate small bodies in parallel
    let small_bodies: Vec<Body> = (0..num_small_bodies)
        .into_par_iter()
        .map(|_| {
            let mut rng = rand::thread_rng();
            
            // Random orbital radius between 1e8 and 5e8 meters
            let radius = rng.gen_range(1.0e8..5.0e8);
            
            // Random angle for position on orbital plane
            let theta = rng.gen_range(0.0..(2.0 * PI));
            
            // Random inclination angle
            let phi = rng.gen_range(0.0..PI);
            
            // Convert spherical to cartesian coordinates
            let x = radius * phi.sin() * theta.cos();
            let y = radius * phi.sin() * theta.sin();
            let z = radius * phi.cos();
            
            // Calculate orbital velocity for a circular orbit
            let orbital_velocity = (G * central_mass / radius).sqrt();
            
            // Calculate velocity vector perpendicular to radius vector
            let v_x = orbital_velocity * (-theta.sin() * phi.sin());
            let v_y = orbital_velocity * (theta.cos() * phi.sin());
            let v_z = 0.0; // No vertical velocity component for circular orbit
            
            Body::new(
                small_body_mass,
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

// Kick-step method (first-order) time integration using kD-tree for force calculation
fn kick_step(bodies: &mut [Body], dt: f64, theta: f64) {
    let n = bodies.len();
    let old_positions: Vec<[f64; 3]> = bodies.iter().map(|b| b.position).collect();
    
    // Update positions based on current velocities (kick) - parallelized
    bodies.par_iter_mut().for_each(|body| {
        for k in 0..3 {
            body.position[k] += body.velocity[k] * dt;
        }
    });
    
    // Build kD-tree from current positions
    let tree_start = Instant::now();
    let old_bodies: Vec<Body> = bodies.iter().enumerate().map(|(i, body)| {
        Body {
            mass: body.mass,
            position: old_positions[i],
            velocity: body.velocity.clone(),
        }
    }).collect();
    
    let tree = Node::build(&old_bodies, 16); // Max 16 bodies per leaf
    println!("Tree building took: {:?}", tree_start.elapsed());
    
    // Calculate accelerations and update velocities (step) - parallelized
    let accel_start = Instant::now();
    bodies.par_iter_mut().for_each(|body| {
        // Calculate acceleration using the tree
        let acceleration = tree.calculate_acceleration(body, theta);
        
        // Update velocity
        for k in 0..3 {
            body.velocity[k] += acceleration[k] * dt;
        }
    });
    println!("Acceleration calculation took: {:?}", accel_start.elapsed());
}

fn main() {
    // Parameters
    let num_bodies = 1_000_000;
    let central_mass = 1.989e30; // Solar mass in kg
    let time_step = 86400.0; // One day in seconds
    let num_steps = 1000;
    let theta = THETA; // Barnes-Hut opening angle parameter
    
    println!("Initializing system with {} bodies...", num_bodies + 1);
    let start_time = Instant::now();
    let mut bodies = initialize_system(num_bodies, central_mass);
    println!("Initialization took: {:?}", start_time.elapsed());
    
    // Calculate initial energy
    println!("Calculating initial energy...");
    let energy_start = Instant::now();
    let initial_energy = calculate_energy(&bodies);
    println!("Initial energy calculation took: {:?}", energy_start.elapsed());
    println!("Initial total energy: {:.6e}", initial_energy);
    
    // Progress tracking
    let completed = AtomicUsize::new(0);
    let progress_interval = 10; // Report progress every 10 steps
    
    // Run simulation
    println!("Starting simulation for {} steps with theta={} and {} threads...", 
             num_steps, theta, rayon::current_num_threads());
    let sim_start_time = Instant::now();
    
    for step in 1..=num_steps {
        let step_start = Instant::now();
        kick_step(&mut bodies, time_step, theta);
        
        // Update progress counter
        let completed_steps = completed.fetch_add(1, Ordering::Relaxed) + 1;
        if completed_steps % progress_interval == 0 {
            println!("Completed step {}/{} ({:.1}%) - took {:?}", 
                     completed_steps, num_steps, 
                     (completed_steps as f64 / num_steps as f64) * 100.0,
                     step_start.elapsed());
        }
    }
    
    let sim_duration = sim_start_time.elapsed();
    println!("Simulation completed in: {:?}", sim_duration);
    println!("Average time per step: {:?}", sim_duration / num_steps as u32);
    
    // Calculate final energy
    println!("Calculating final energy...");
    let final_energy_start = Instant::now();
    let final_energy = calculate_energy(&bodies);
    println!("Final energy calculation took: {:?}", final_energy_start.elapsed());
    println!("Final total energy: {:.6e}", final_energy);
    
    // Calculate energy error
    let energy_error = (final_energy - initial_energy).abs() / initial_energy.abs();
    println!("Relative energy error: {:.6e}", energy_error);
    
    // Print summary of system state
    let mut total_mass = 0.0;
    let mut center_of_mass = [0.0, 0.0, 0.0];
    let mut total_momentum = [0.0, 0.0, 0.0];
    
    for body in &bodies {
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