use rand::Rng;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign};
use std::time::Instant;

// --- Constants ---
const G: f64 = 1.0; // Gravitational constant (set to 1 for simplicity)
const DT: f64 = 0.001; // Time step
const SOFTENING: f64 = 1e-5; // Softening factor to prevent division by zero/huge forces

// --- Vector3D Struct ---
#[derive(Debug, Clone, Copy, PartialEq)]
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
        Vec3 { x: 0.0, y: 0.0, z: 0.0 }
    }

    fn magnitude_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn magnitude(&self) -> f64 {
        self.magnitude_squared().sqrt()
    }

    // // Optional: Normalize method if needed elsewhere
    // fn normalize(&self) -> Self {
    //     let mag = self.magnitude();
    //     if mag == 0.0 {
    //         Vec3::zero()
    //     } else {
    //         *self / mag
    //     }
    // }

    // // Optional: Cross product if needed for more complex initializations
    // fn cross(&self, other: &Self) -> Self {
    //     Vec3 {
    //         x: self.y * other.z - self.z * other.y,
    //         y: self.z * other.x - self.x * other.z,
    //         z: self.x * other.y - self.y * other.x,
    //     }
    // }
}

// Operator Overloading for Vec3
impl Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Vec3 { x: self.x + other.x, y: self.y + other.y, z: self.z + other.z }
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Vec3 { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, other: Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}


impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Vec3 { x: self.x * scalar, y: self.y * scalar, z: self.z * scalar }
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        Vec3 { x: self.x / scalar, y: self.y / scalar, z: self.z / scalar }
    }
}

// --- Body Struct ---
#[derive(Debug, Clone, Copy)]
struct Body {
    mass: f64,
    pos: Vec3,
    vel: Vec3,
}

impl Body {
    fn new(mass: f64, pos: Vec3, vel: Vec3) -> Self {
        Body { mass, pos, vel }
    }
}

// --- Simulation Functions ---

/// Calculates the total energy of the system (Kinetic + Potential).
// fn calculate_total_energy(bodies: &[Body]) -> f64 {
//     let mut kinetic_energy = 0.0;
//     let mut potential_energy = 0.0;
//     let n = bodies.len();

//     for i in 0..n {
//         // Kinetic Energy: 0.5 * m * v^2
//         kinetic_energy += 0.5 * bodies[i].mass * bodies[i].vel.magnitude_squared();

//         // Potential Energy: -G * m_i * m_j / |r_i - r_j|
//         // Sum over pairs (i < j) to avoid double counting and self-interaction
//         for j in (i + 1)..n {
//             let dr = bodies[j].pos - bodies[i].pos;
//             let dist_sq = dr.magnitude_squared();
//             let dist = (dist_sq + SOFTENING * SOFTENING).sqrt(); // Add softening
//             potential_energy -= G * bodies[i].mass * bodies[j].mass / dist;
//         }
//     }

//     kinetic_energy + potential_energy
// }

/// Performs one step of the simulation using the kick-step (Euler-Cromer) method.
fn simulate_step(bodies: &mut [Body], dt: f64) {
    let n = bodies.len();
    let mut accelerations = vec![Vec3::zero(); n]; // Store accelerations for this step

    // 1. Calculate Accelerations (Force Calculation) based on current positions
    // F_i = sum_{j!=i} G * m_i * m_j * (r_j - r_i) / |r_j - r_i|^3
    // a_i = F_i / m_i = sum_{j!=i} G * m_j * (r_j - r_i) / |r_j - r_i|^3
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dr = bodies[j].pos - bodies[i].pos;
            let dist_sq = dr.magnitude_squared();
            // Add softening to prevent division by zero and extreme forces at close range
            let dist_cubed = (dist_sq + SOFTENING * SOFTENING).powf(1.5);
            let force_factor = G * bodies[j].mass / dist_cubed;
            accelerations[i] += dr * force_factor;
        }
    }

    // 2. Kick: Update Velocities using calculated accelerations
    // v_new = v_old + a * dt
    for i in 0..n {
        bodies[i].vel += accelerations[i] * dt;
    }

    // 3. Step: Update Positions using the *new* velocities
    // r_new = r_old + v_new * dt
    for i in 0..n {
        bodies[i].pos += bodies[i].vel * dt;
    }
}

/// Initializes a system with one central body and N smaller bodies in circular orbits.
fn initialize_circular_orbits(
    num_small_bodies: usize,
    central_mass: f64,
    small_mass: f64,
    avg_radius: f64,
    radius_spread: f64, // Allow some randomness in orbital radius
) -> Vec<Body> {
    let mut bodies = Vec::with_capacity(num_small_bodies + 1);
    let mut rng = rand::thread_rng();

    // Add central body
    let central_body = Body::new(central_mass, Vec3::zero(), Vec3::zero());
    bodies.push(central_body);

    // Add orbiting bodies
    for _ in 0..num_small_bodies {
        // Random radius around the average
        let r_mag = avg_radius + radius_spread * (rng.r#gen::<f64>() - 0.5) * 2.0; // Uniform spread
        
        // Random position on a sphere of radius r_mag
        // Use spherical coordinates to distribute somewhat evenly
        let theta = rng.r#gen::<f64>() * 2.0 * PI; // Azimuthal angle (0 to 2pi)
        let phi = (rng.r#gen::<f64>() * 2.0 - 1.0).acos(); // Polar angle (0 to pi) - ensures uniform spherical distribution

        let x = r_mag * phi.sin() * theta.cos();
        let y = r_mag * phi.sin() * theta.sin();
        let z = r_mag * phi.cos();
        let pos = Vec3::new(x, y, z);

        // Calculate circular orbit speed: v = sqrt(G * M_central / r)
        let speed = (G * central_mass / r_mag).sqrt();

        // Calculate velocity vector perpendicular to position vector (for circular orbit)
        // A simple way: create a vector roughly perpendicular in the xy-plane projection,
        // then ensure it's truly perpendicular using cross product logic if needed,
        // but for random orientations, a simpler tangential velocity often suffices for starting.
        // Let's choose a simple perpendicular vector in the xy-plane projection first.
        // If pos = (x, y, z), a perpendicular vector in xy plane is (-y, x, 0). Normalize and scale.
        let mut vel_dir = Vec3::new(-y, x, 0.0); 
        if vel_dir.magnitude_squared() < 1e-12 { // Avoid division by zero if body is near z-axis
             vel_dir = Vec3::new(1.0, 0.0, 0.0); // Assign arbitrary perpendicular if needed
        }
       
        // Normalize the direction and scale by speed
        // Note: This simplified velocity might not be perfectly in the plane defined by pos and Z-axis,
        // leading to initially slightly non-circular/non-planar orbits, which is often acceptable.
        // For perfectly planar circular orbits, one would define an orbital plane first.
        let vel = vel_dir * (speed / vel_dir.magnitude());


        bodies.push(Body::new(small_mass, pos, vel));
    }

    bodies
}

// --- Main Execution ---
fn main() {
    let n_small_bodies = 10_000;
    // let n_small_bodies = 1_000; // Use a smaller number for faster testing
    let n_steps = 100;
    let central_mass = 1_000_000.0;
    let small_mass = 1.0;
    let orbital_radius = 100.0;
    let radius_spread = 20.0; // +/- 10 units from orbital_radius

    println!("N-Body Simulation");
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


    println!("Initializing system...");
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


    // println!("Calculating initial energy...");
    // let initial_energy = calculate_total_energy(&bodies);
    // println!("Initial Total Energy: {:.6e}", initial_energy);


    println!("Starting simulation...");
    let start_sim = Instant::now();
    for step in 0..n_steps {
        simulate_step(&mut bodies, DT);

        // Optional: Print progress
        if (step + 1) % 100 == 0 {
            let elapsed = start_sim.elapsed();
            let avg_step_time = elapsed.as_secs_f64() / (step + 1) as f64;
            println!(
                "Step {}/{} completed. Elapsed: {:.2?}. Avg time/step: {:.3} s",
                step + 1,
                n_steps,
                elapsed,
                avg_step_time
            );
             // Optional: Calculate energy mid-simulation (can be slow)
             // let current_energy = calculate_total_energy(&bodies);
             // println!("  Current Energy: {:.6e}", current_energy);
        }
    }
    let sim_duration = start_sim.elapsed();
    println!("Simulation finished in {:.2?}", sim_duration);


    // println!("Calculating final energy...");
    // let final_energy = calculate_total_energy(&bodies);
    // println!("Final Total Energy:   {:.6e}", final_energy);

    // let energy_diff = final_energy - initial_energy;
    // let relative_energy_diff = if initial_energy.abs() > 1e-12 {
    //     (energy_diff / initial_energy).abs()
    // } else {
    //     energy_diff.abs() // Avoid division by zero if initial energy is near zero
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