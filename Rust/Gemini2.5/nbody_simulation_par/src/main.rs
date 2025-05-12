// Add Rayon prelude
use rayon::prelude::*;

use rand::Rng;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign};
use std::time::Instant;

// --- Constants --- (Same as before)
const G: f64 = 1.0;
const DT: f64 = 0.001;
const SOFTENING: f64 = 1e-5;

// --- Vec3 Struct --- (Same as before)
#[derive(Debug, Clone, Copy, PartialEq, Default)] // Added Default
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    fn new(x: f64, y: f64, z: f64) -> Self { Vec3 { x, y, z } }
    fn zero() -> Self { Vec3 { x: 0.0, y: 0.0, z: 0.0 } }
    fn magnitude_squared(&self) -> f64 { self.x * self.x + self.y * self.y + self.z * self.z }
    fn magnitude(&self) -> f64 { self.magnitude_squared().sqrt() }
}
impl Add for Vec3 { type Output = Self; fn add(self, other: Self) -> Self { Vec3 { x: self.x + other.x, y: self.y + other.y, z: self.z + other.z } } }
impl AddAssign for Vec3 { fn add_assign(&mut self, other: Self) { self.x += other.x; self.y += other.y; self.z += other.z; } }
impl Sub for Vec3 { type Output = Self; fn sub(self, other: Self) -> Self { Vec3 { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z } } }
impl SubAssign for Vec3 { fn sub_assign(&mut self, other: Self) { self.x -= other.x; self.y -= other.y; self.z -= other.z; } }
impl Mul<f64> for Vec3 { type Output = Self; fn mul(self, scalar: f64) -> Self { Vec3 { x: self.x * scalar, y: self.y * scalar, z: self.z * scalar } } }
impl Div<f64> for Vec3 { type Output = Self; fn div(self, scalar: f64) -> Self { Vec3 { x: self.x / scalar, y: self.y / scalar, z: self.z / scalar } } }


// --- Body Struct --- (Same as before)
// Make Send + Sync explicit for Rayon, though derive(Clone, Copy) should suffice
#[derive(Debug, Clone, Copy)]
struct Body {
    mass: f64,
    pos: Vec3,
    vel: Vec3,
}
unsafe impl Send for Body {}
unsafe impl Sync for Body {}


impl Body {
    fn new(mass: f64, pos: Vec3, vel: Vec3) -> Self {
        Body { mass, pos, vel }
    }
}

// --- Simulation Functions ---

/// Calculates the total energy of the system (Kinetic + Potential) in parallel.
// fn calculate_total_energy(bodies: &[Body]) -> f64 {
//     let n = bodies.len();

//     // Parallel calculation of Kinetic Energy
//     let kinetic_energy: f64 = bodies
//         .par_iter() // Use parallel iterator
//         .map(|b| 0.5 * b.mass * b.vel.magnitude_squared())
//         .sum(); // Parallel sum

//     // Parallel calculation of Potential Energy
//     // Iterate over indices `i`, and for each `i`, sum contributions from `j > i`
//     let potential_energy: f64 = (0..n)
//         .into_par_iter() // Parallelize the outer loop index `i`
//         .map(|i| {
//             let mut local_potential = 0.0;
//             // Inner loop remains sequential for each `i`
//             for j in (i + 1)..n {
//                 let dr = bodies[j].pos - bodies[i].pos;
//                 let dist_sq = dr.magnitude_squared();
//                 let dist = (dist_sq + SOFTENING * SOFTENING).sqrt();
//                 local_potential -= G * bodies[i].mass * bodies[j].mass / dist;
//             }
//             local_potential
//         })
//         .sum(); // Parallel sum of results from each `i`

//     kinetic_energy + potential_energy
// }


/// Performs one step of the simulation using the kick-step method in parallel.
fn simulate_step(bodies: &mut [Body], dt: f64) {
    let n = bodies.len();
    // We need accelerations based on positions *before* updates.
    // Calculate accelerations in parallel.
    // Each thread calculates acceleration for a subset of bodies.
    let accelerations: Vec<Vec3> = (0..n)
        .into_par_iter() // Parallel iteration over body indices
        .map(|i| {
            let mut acc_i = Vec3::zero();
            // Inner loop is sequential for each body 'i', but needs access to all 'j'
            // bodies slice is captured immutably here.
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dr = bodies[j].pos - bodies[i].pos; // Read positions
                let dist_sq = dr.magnitude_squared();
                let dist_cubed = (dist_sq + SOFTENING * SOFTENING).powf(1.5);
                let force_factor = G * bodies[j].mass / dist_cubed;
                acc_i += dr * force_factor;
            }
            acc_i // Return calculated acceleration for body i
        })
        .collect(); // Collect results into a Vec

    // 2. Kick: Update Velocities in parallel using calculated accelerations
    bodies
        .par_iter_mut() // Parallel mutable iterator over bodies
        .zip(accelerations.par_iter()) // Zip with parallel iterator over accelerations
        .for_each(|(body, acc)| {
            body.vel += *acc * dt; // Update velocity
        });

    // 3. Step: Update Positions in parallel using the *new* velocities
    bodies
        .par_iter_mut() // Parallel mutable iterator
        .for_each(|body| {
            body.pos += body.vel * dt; // Update position
        });
}


/// Initializes a system with one central body and N smaller bodies in parallel.
fn initialize_circular_orbits(
    num_small_bodies: usize,
    central_mass: f64,
    small_mass: f64,
    avg_radius: f64,
    radius_spread: f64,
) -> Vec<Body> {
    let mut bodies = Vec::with_capacity(num_small_bodies + 1);

    // Add central body (sequentially, it's just one)
    let central_body = Body::new(central_mass, Vec3::zero(), Vec3::zero());
    bodies.push(central_body);

    // Generate orbiting bodies in parallel
    let small_bodies: Vec<Body> = (0..num_small_bodies)
        .into_par_iter() // Parallel iterator over the count
        .map(|_| {
            // Each thread gets its own RNG instance
            let mut rng = rand::thread_rng();

            let r_mag = avg_radius + radius_spread * (rng.gen::<f64>() - 0.5) * 2.0;
            let theta = rng.gen::<f64>() * 2.0 * PI;
            let phi = (rng.gen::<f64>() * 2.0 - 1.0).acos();

            let x = r_mag * phi.sin() * theta.cos();
            let y = r_mag * phi.sin() * theta.sin();
            let z = r_mag * phi.cos();
            let pos = Vec3::new(x, y, z);

            let speed = (G * central_mass / r_mag).sqrt();
            let mut vel_dir = Vec3::new(-y, x, 0.0);
            if vel_dir.magnitude_squared() < 1e-12 {
                 vel_dir = Vec3::new(1.0, 0.0, 0.0);
            }
            let vel = vel_dir * (speed / vel_dir.magnitude());

            Body::new(small_mass, pos, vel) // Return the generated body
        })
        .collect(); // Collect the parallel results into a Vec

    // Extend the main vector with the generated small bodies
    bodies.extend(small_bodies);

    bodies
}


// --- Main Execution ---
fn main() {
    let n_small_bodies = 10_000;
    // let n_small_bodies = 5_000; // Use a smaller number for faster testing/debugging
    let n_steps = 100;
    let central_mass = 1_000_000.0;
    let small_mass = 1.0;
    let orbital_radius = 100.0;
    let radius_spread = 20.0;

    println!("N-Body Simulation (Multithreaded with Rayon)"); // Updated title
    println!("-----------------");
    println!("Number of small bodies: {}", n_small_bodies);
    println!("Number of steps: {}", n_steps);
    println!("Time step (dt): {}", DT);
    println!("Gravitational Constant (G): {}", G);
    println!("Softening Factor: {}", SOFTENING);
    println!("Central Mass: {}", central_mass);
    println!("Small Body Mass: {}", small_mass);
    println!("Target Orbital Radius: {}", orbital_radius);
    println!("-----------------");
    println!("Using {} CPU threads (default for Rayon)", rayon::current_num_threads());


    println!("Initializing system (parallel)..."); // Updated message
    let start_init = Instant::now();
    let mut bodies = initialize_circular_orbits(
        n_small_bodies,
        central_mass,
        small_mass,
        orbital_radius,
        radius_spread,
    );
    let init_duration = start_init.elapsed();
    println!("Initialization complete ({:.2?} total bodies) in {:.2?}", bodies.len(), init_duration);


    // println!("Calculating initial energy (parallel)..."); // Updated message
    // let start_energy = Instant::now();
    // let initial_energy = calculate_total_energy(&bodies);
    // let energy_duration = start_energy.elapsed();
    // println!("Initial Total Energy: {:.6e} (calculated in {:.2?})", initial_energy, energy_duration);


    println!("Starting simulation...");
    let start_sim = Instant::now();
    for step in 0..n_steps {
        let step_start = Instant::now(); // Time each step
        simulate_step(&mut bodies, DT);
        let step_duration = step_start.elapsed();

        // Optional: Print progress more frequently for long sims
        if (step + 1) % 10 == 0 || step == 0 { // Print every 10 steps + first step
             let elapsed_total = start_sim.elapsed();
             let avg_step_time = elapsed_total.as_secs_f64() / (step + 1) as f64;
             println!(
                 "Step {:>4}/{}: Step time: {:>8.3?} | Avg step time: {:>6.3} s | Total elapsed: {:>8.2?}",
                 step + 1,
                 n_steps,
                 step_duration,
                 avg_step_time,
                 elapsed_total
             );
             // Calculating energy every step is VERY expensive, do it less often if needed
             // if (step + 1) % 100 == 0 {
             //     let current_energy = calculate_total_energy(&bodies);
             //     println!("  Current Energy (@step {}): {:.6e}", step + 1, current_energy);
             // }
        }
    }
    let sim_duration = start_sim.elapsed();
    println!("Simulation finished in {:.2?}", sim_duration);


    // println!("Calculating final energy (parallel)..."); // Updated message
    // let start_energy = Instant::now();
    // let final_energy = calculate_total_energy(&bodies);
    // let energy_duration = start_energy.elapsed();
    // println!("Final Total Energy:   {:.6e} (calculated in {:.2?})", final_energy, energy_duration);


    // let energy_diff = final_energy - initial_energy;
    // let relative_energy_diff = if initial_energy.abs() > 1e-12 {
    //     (energy_diff / initial_energy).abs()
    // } else {
    //     energy_diff.abs()
    // };
    // println!("Absolute Energy Change: {:.6e}", energy_diff);
    // println!("Relative Energy Change: {:.6e} ({:.4}%)", relative_energy_diff, relative_energy_diff * 100.0);
    // Optional: Print final position of a few bodies
    println!("\nFinal state of first few bodies:");
    for i in 0..std::cmp::min(5, bodies.len()) {
        println!("Body {}: Pos=({:.2}, {:.2}, {:.2}), Vel=({:.2}, {:.2}, {:.2})",
            i,
            bodies[i].pos.x, bodies[i].pos.y, bodies[i].pos.z,
            bodies[i].vel.x, bodies[i].vel.y, bodies[i].vel.z);
    }
}