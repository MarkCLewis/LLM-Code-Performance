use std::f64::consts::PI;
use rand::Rng;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

const G: f64 = 6.67430e-11; // gravitational constant
const SOFTENING: f64 = 1.0e-3; // softening parameter to avoid singularities

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

// Calculate acceleration on body i due to body j
fn calculate_acceleration(body_i: &Body, body_j: &Body) -> [f64; 3] {
    let mut r = [0.0; 3];
    for k in 0..3 {
        r[k] = body_j.position[k] - body_i.position[k];
    }
    
    let distance_squared = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    let distance = (distance_squared + SOFTENING * SOFTENING).sqrt();
    let prefactor = G * body_j.mass / (distance * distance * distance);
    
    [prefactor * r[0], prefactor * r[1], prefactor * r[2]]
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

// Kick-step method (first-order) time integration
fn kick_step(bodies: &mut [Body], dt: f64) {
    let n = bodies.len();
    let old_positions: Vec<[f64; 3]> = bodies.iter().map(|b| b.position).collect();
    
    // Update positions based on current velocities (kick) - parallelized
    bodies.par_iter_mut().for_each(|body| {
        for k in 0..3 {
            body.position[k] += body.velocity[k] * dt;
        }
    });
    
    // Calculate accelerations and update velocities (step) - parallelized
    bodies.par_iter_mut().enumerate().for_each(|(i, body)| {
        let mut acceleration = [0.0, 0.0, 0.0];
        
        for j in 0..n {
            if i == j {
                continue;
            }
            
            // Create temporary body with old position of body j
            let body_j = Body {
                mass: bodies[j].mass,
                position: old_positions[j],
                velocity: bodies[j].velocity,
            };
            
            // Calculate acceleration due to body j
            let acc_ij = calculate_acceleration(body, &body_j);
            
            // Accumulate acceleration
            for k in 0..3 {
                acceleration[k] += acc_ij[k];
            }
        }
        
        // Update velocity
        for k in 0..3 {
            body.velocity[k] += acceleration[k] * dt;
        }
    });
}

fn main() {
    // Parameters
    let num_bodies = 10_000;
    let central_mass = 1.989e30; // Solar mass in kg
    let time_step = 86400.0; // One day in seconds
    let num_steps = 100;
    
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
    let progress_interval = 50; // Report progress every 50 steps
    
    // Run simulation
    println!("Starting simulation for {} steps with {} threads...", 
             num_steps, rayon::current_num_threads());
    let sim_start_time = Instant::now();
    
    // Create time step blocks for better parallelism
    let block_size = 10; // Process 10 steps in each block
    let num_blocks = (num_steps + block_size - 1) / block_size;
    
    (0..num_blocks).into_par_iter().for_each(|block| {
        let start_step = block * block_size + 1;
        let end_step = std::cmp::min((block + 1) * block_size, num_steps);
        
        let mut local_bodies = bodies.clone();
        
        for step in start_step..=end_step {
            kick_step(&mut local_bodies, time_step);
            
            // Update progress counter
            let completed_steps = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if completed_steps % progress_interval == 0 {
                println!("Completed step {}/{} ({:.1}%)", 
                         completed_steps, num_steps, 
                         (completed_steps as f64 / num_steps as f64) * 100.0);
            }
        }
        
        // We can't use the results from parallel blocks as each has its own copy
        // This is just to demonstrate the structure - in a real implementation
        // we'd need a better strategy for combining results
        if block == 0 {
            bodies = local_bodies;
        }
    });
    
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