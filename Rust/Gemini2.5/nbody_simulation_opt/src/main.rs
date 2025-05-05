// (Keep previous use statements)
use rand::Rng;
use rayon::prelude::*;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div};
use std::time::Instant;
use std::fmt;
use typed_arena::Arena; // Import Arena

// --- Constants ---
const G: f64 = 1.0;
const EPSILON_SQ: f64 = 1e-6;
const THETA: f64 = 0.3;
const THETA_SQ: f64 = THETA * THETA; // Precompute theta squared

// --- Data Structures ---

#[derive(Debug, Clone, Copy, PartialEq)]
struct Vec3d { x: f64, y: f64, z: f64 }

// --- Vec3d Implementation (Add #[inline]) ---
impl Vec3d {
    #[inline] fn new(x: f64, y: f64, z: f64) -> Self { Vec3d { x, y, z } }
    #[inline] fn zero() -> Self { Vec3d { x: 0.0, y: 0.0, z: 0.0 } }
    #[inline] fn norm_sq(&self) -> f64 { self.x.mul_add(self.x, self.y.mul_add(self.y, self.z * self.z)) } // Use FMA
    #[inline] fn norm(&self) -> f64 { self.norm_sq().sqrt() }
    #[inline] fn dist_sq(&self, other: &Vec3d) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        dx.mul_add(dx, dy.mul_add(dy, dz * dz)) // Use FMA
    }
}
// --- Operator Overloads (Add #[inline]) ---
impl Add for Vec3d { type Output = Self; #[inline] fn add(self, other: Self) -> Self { Vec3d::new(self.x + other.x, self.y + other.y, self.z + other.z) } }
impl Sub for Vec3d { type Output = Self; #[inline] fn sub(self, other: Self) -> Self { Vec3d::new(self.x - other.x, self.y - other.y, self.z - other.z) } }
impl Mul<f64> for Vec3d { type Output = Self; #[inline] fn mul(self, scalar: f64) -> Self { Vec3d::new(self.x * scalar, self.y * scalar, self.z * scalar) } }
impl Mul<Vec3d> for f64 { type Output = Vec3d; #[inline] fn mul(self, vec: Vec3d) -> Vec3d { vec * self } }
impl Div<f64> for Vec3d { type Output = Self; #[inline] fn div(self, scalar: f64) -> Self { let inv_scalar = 1.0 / scalar; Vec3d::new(self.x * inv_scalar, self.y * inv_scalar, self.z * inv_scalar) } }


#[derive(Debug, Clone, Copy)]
struct Body {
    id: usize,
    position: Vec3d,
    velocity: Vec3d,
    mass: f64,
}

// --- Bounding Box Helper (Add #[inline]) ---
#[derive(Debug, Clone, Copy)]
struct BoundingBox { min: Vec3d, max: Vec3d }
impl BoundingBox {
    fn compute_global(bodies: &[Body]) -> Self {
        // (Same implementation as before)
        if bodies.is_empty() { return BoundingBox { min: Vec3d::zero(), max: Vec3d::zero() }; }
        let mut min_p = bodies[0].position; let mut max_p = bodies[0].position;
        for body in bodies.iter().skip(1) {
            min_p.x = min_p.x.min(body.position.x); min_p.y = min_p.y.min(body.position.y); min_p.z = min_p.z.min(body.position.z);
            max_p.x = max_p.x.max(body.position.x); max_p.y = max_p.y.max(body.position.y); max_p.z = max_p.z.max(body.position.z);
        }
        let size = max_p - min_p; let buffer = size * 0.01 + Vec3d::new(1e-5, 1e-5, 1e-5);
        BoundingBox { min: min_p - buffer, max: max_p + buffer }
    }
    #[inline] fn center(&self) -> Vec3d { (self.min + self.max) * 0.5 }
    #[inline] fn size(&self) -> Vec3d { self.max - self.min }
    #[inline] fn max_side_length(&self) -> f64 { let size = self.size(); size.x.max(size.y).max(size.z) }
    #[inline] fn max_side_length_sq(&self) -> f64 { let s = self.max_side_length(); s * s } // Calculate squared size directly
    fn split(&self, dim: usize, val: f64) -> (BoundingBox, BoundingBox) {
        // (Same implementation as before)
        let mut left_max = self.max; let mut right_min = self.min;
        match dim { 0 => { left_max.x = val; right_min.x = val; } 1 => { left_max.y = val; right_min.y = val; } _ => { left_max.z = val; right_min.z = val; } }
        (BoundingBox { min: self.min, max: left_max }, BoundingBox { min: right_min, max: self.max })
    }
}


// --- k-d Tree Node Structure (Using Arena References) ---
// Note: Lifetimes 'a indicate references tied to the Arena's lifetime
enum KdNode<'a> {
    Empty, // Still useful as a base case / placeholder
    Leaf {
        body_id: usize,
        position: Vec3d,
        mass: f64,
        // No bounds needed in leaf if using parent bounds during traversal
    },
    Internal {
        split_dim: usize,
        split_val: f64,
        bounds: BoundingBox, // Store bounds for 's' calculation
        total_mass: f64,
        center_of_mass: Vec3d,
        left: &'a KdNode<'a>,  // Use reference tied to arena lifetime
        right: &'a KdNode<'a>, // Use reference tied to arena lifetime
    },
}

// Manual Debug impl needed as derive won't work easily with lifetimes/recursion here
impl<'a> fmt::Debug for KdNode<'a> {
     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KdNode::Empty => write!(f, "Empty"),
            KdNode::Leaf { body_id, position, mass } => f.debug_struct("Leaf")
                .field("body_id", body_id)
                .field("position", position)
                .field("mass", mass)
                .finish(),
            KdNode::Internal { split_dim, split_val, bounds, total_mass, center_of_mass, .. } => f.debug_struct("Internal")
                .field("split_dim", split_dim)
                .field("split_val", split_val)
                .field("bounds", bounds)
                .field("total_mass", total_mass)
                .field("center_of_mass", center_of_mass)
                .field("left", &"...") // Avoid deep print
                .field("right", &"...")
                .finish(),
        }
    }
}


impl<'a> KdNode<'a> {
    #[inline]
    fn get_coord(body: &Body, dim: usize) -> f64 {
        match dim { 0 => body.position.x, 1 => body.position.y, _ => body.position.z, }
    }

    // Build function now takes and returns references tied to the arena
    fn build(
        arena: &'a Arena<KdNode<'a>>, // Pass arena reference
        mut body_indices: &mut [usize],
        all_bodies: &[Body],
        bounds: BoundingBox,
        depth: usize,
    ) -> &'a KdNode<'a> { // Return arena reference
        let n = body_indices.len();

        if n == 0 {
            // Allocate an Empty node in the arena
            return arena.alloc(KdNode::Empty);
        }

        if n == 1 {
            let body_idx = body_indices[0];
            let body = &all_bodies[body_idx];
            // Allocate a Leaf node in the arena
            return arena.alloc(KdNode::Leaf {
                body_id: body.id,
                position: body.position,
                mass: body.mass,
                // bounds are implicit via parent in traversal now
            });
        }

        let split_dim = depth % 3;
        let median_index = n / 2;

        // Use partitioning (O(N)) - requires ordered_float::NotNan
         body_indices.select_nth_unstable_by_key(median_index, |&idx| {
             ordered_float::NotNan::new(Self::get_coord(&all_bodies[idx], split_dim)).unwrap()
         });

        let median_body_idx = body_indices[median_index];
        let split_val = Self::get_coord(&all_bodies[median_body_idx], split_dim);

        // Split indices (median goes left)
        let (left_indices, right_indices) = body_indices.split_at_mut(median_index + 1);

        // Split bounds
        let (left_bounds, right_bounds) = bounds.split(split_dim, split_val);

        // Build children recursively, getting references back
        let left_child_ref = Self::build(arena, left_indices, all_bodies, left_bounds, depth + 1);
        let right_child_ref = Self::build(arena, right_indices, all_bodies, right_bounds, depth + 1);

        // Calculate CoM and total mass for this internal node
        let (total_mass, center_of_mass) = Self::combine_node_properties(left_child_ref, right_child_ref);

        // Allocate an Internal node in the arena
        arena.alloc(KdNode::Internal {
            split_dim,
            split_val,
            bounds, // Store bounds for this node
            total_mass,
            center_of_mass,
            left: left_child_ref, // Store references
            right: right_child_ref, // Store references
        })
    }

    // Helper to combine mass/CoM (using FMA)
    #[inline]
    fn combine_node_properties(left: &KdNode, right: &KdNode) -> (f64, Vec3d) {
        let (m_l, com_l) = left.get_mass_com();
        let (m_r, com_r) = right.get_mass_com();
        let total_mass = m_l + m_r;

        if total_mass < 1e-99 { // Check against small number instead of exact zero
            (0.0, Vec3d::zero())
        } else {
             // Use FMA: com_l * m_l + com_r * m_r
             let weighted_sum_pos = Vec3d::new(
                 com_l.x.mul_add(m_l, com_r.x * m_r),
                 com_l.y.mul_add(m_l, com_r.y * m_r),
                 com_l.z.mul_add(m_l, com_r.z * m_r),
             );
             (total_mass, weighted_sum_pos / total_mass)
        }
    }

    // Helper to get mass and CoM
    #[inline]
    fn get_mass_com(&self) -> (f64, Vec3d) {
        match self {
            KdNode::Empty => (0.0, Vec3d::zero()),
            KdNode::Leaf { mass, position, .. } => (*mass, *position),
            KdNode::Internal { total_mass, center_of_mass, .. } => (*total_mass, *center_of_mass),
        }
    }


    // Calculate acceleration using optimized BH criterion and math
    fn calculate_acceleration_on(
        &self,
        target_body: &Body,
        // theta_sq is now used instead of theta
    ) -> Vec3d {
        match self {
            KdNode::Empty => Vec3d::zero(),
            KdNode::Leaf { body_id, position, mass, .. } => {
                if *body_id == target_body.id { return Vec3d::zero(); } // Self-interaction check

                let rij = *position - target_body.position;
                let dist_sq = rij.norm_sq();
                // Add small value to dist_sq before sqrt/division to prevent issues near zero,
                // in addition to EPSILON_SQ in denominator power calc.
                if dist_sq < 1e-99 { return Vec3d::zero(); }

                // Optimized inverse cube calculation
                let r_soft_sq = dist_sq + EPSILON_SQ;
                let r_soft = r_soft_sq.sqrt();
                let inv_r_soft_cubed = 1.0 / (r_soft_sq * r_soft); // r^3 = r^2 * r

                // Acceleration = G * mass * rij / r_soft^3 = rij * (G * mass * inv_r_soft_cubed)
                rij * (G * (*mass) * inv_r_soft_cubed)
            }
            KdNode::Internal {
                bounds,
                total_mass,
                center_of_mass,
                left,
                right,
                ..
            } => {
                 let d_sq = target_body.position.dist_sq(center_of_mass);
                 // Add small value to d_sq before sqrt/division checks
                 if d_sq < 1e-99 {
                     // If exactly at CoM, must recurse (should be rare)
                      return left.calculate_acceleration_on(target_body)
                           + right.calculate_acceleration_on(target_body);
                 }

                let s_sq = bounds.max_side_length_sq(); // Use squared size

                // Optimized Barnes-Hut criterion: s*s < theta^2 * d*d
                if s_sq < THETA_SQ * d_sq {
                    // Node is sufficiently far away, approximate using CoM
                    let r_soft_sq = d_sq + EPSILON_SQ;
                    let r_soft = r_soft_sq.sqrt();
                    let inv_r_soft_cubed = 1.0 / (r_soft_sq * r_soft);
                    let force_magnitude_part = G * (*total_mass) * inv_r_soft_cubed;

                    (*center_of_mass - target_body.position) * force_magnitude_part
                } else {
                    // Node is too close, recurse into children
                    left.calculate_acceleration_on(target_body)
                    + right.calculate_acceleration_on(target_body)
                }
            }
        }
    }
}


// --- Initialization Function --- (Assign unique IDs - Same as before)
fn initialize_solar_system(/* ... */) -> Vec<Body> {
    // (No changes needed here)
    let num_small_bodies: usize = 100_000; // Example value
    let central_mass: f64 = 1_000_000.0;
    let small_mass: f64 = 1.0;
    let max_radius: f64 = 100.0;

    let mut bodies = Vec::with_capacity(num_small_bodies + 1);
    let mut rng = rand::thread_rng();
    bodies.push(Body { id: 0, position: Vec3d::zero(), velocity: Vec3d::zero(), mass: central_mass });
    for i in 0..num_small_bodies {
        let r = rng.gen::<f64>().sqrt() * max_radius;
        let theta_angle = rng.gen::<f64>() * 2.0 * PI;
        let phi = (rng.gen::<f64>() * 2.0 - 1.0).acos();
        let x = r * phi.sin() * theta_angle.cos(); let y = r * phi.sin() * theta_angle.sin(); let z = r * phi.cos();
        let pos = Vec3d::new(x, y, z);
        let v_mag = if r > 1e-6 { (G * central_mass / r).sqrt() } else { 0.0 };
        let mut vel = Vec3d::new(-pos.y, pos.x, 0.0);
        let vel_norm = vel.norm();
         if vel_norm > 1e-9 { vel = vel * (v_mag / vel_norm); }
         else { let random_angle = rng.gen::<f64>() * 2.0 * PI; vel = Vec3d::new(v_mag * random_angle.cos(), v_mag * random_angle.sin(), 0.0); }
        bodies.push(Body { id: i + 1, position: pos, velocity: vel, mass: small_mass });
    }
    bodies
}


// --- Simulation Step (Using Arena k-d Tree) ---
fn simulation_step(bodies: &mut [Body], dt: f64) {
    let n = bodies.len();
    if n == 0 { return; }

    // --- 1. Setup Arena and Build Tree ---
    let tree_node_arena = Arena::new(); // Create arena for this step
    let mut indices: Vec<usize> = (0..n).collect();
    let global_bounds = BoundingBox::compute_global(bodies);
    // Build tree using the arena, get reference to root
    let kdtree_root = KdNode::build(&tree_node_arena, &mut indices, bodies, global_bounds, 0);


    // --- 2. Calculate Accelerations using Tree (Parallel) ---
    let accelerations: Vec<Vec3d> = bodies
        .par_iter()
        .map(|target_body| {
            // Pass theta_sq implicitly via constant
            kdtree_root.calculate_acceleration_on(target_body)
        })
        .collect();

    // --- 3. Update Velocities and Positions (Parallel) ---
    bodies
        .par_iter_mut()
        .zip(accelerations.par_iter())
        .for_each(|(body, &acc)| {
            // Use FMA for velocity update? v = v + a*dt
            // body.velocity = body.velocity + acc * dt;
             body.velocity.x = acc.x.mul_add(dt, body.velocity.x);
             body.velocity.y = acc.y.mul_add(dt, body.velocity.y);
             body.velocity.z = acc.z.mul_add(dt, body.velocity.z);

            // Use FMA for position update? p = p + v*dt
            // body.position = body.position + body.velocity * dt;
             body.position.x = body.velocity.x.mul_add(dt, body.position.x);
             body.position.y = body.velocity.y.mul_add(dt, body.position.y);
             body.position.z = body.velocity.z.mul_add(dt, body.position.z);
        });

    // Arena `tree_node_arena` goes out of scope here, memory is freed.
}

// --- Energy Calculation (Accurate O(N^2) Parallel - Same as before) ---
fn calculate_energy(bodies: &[Body]) -> f64 {
    // (No changes needed here - keep for verification)
    let kinetic_energy = bodies.par_iter().map(|b| 0.5 * b.mass * b.velocity.norm_sq()).sum::<f64>();
    let potential_energy = bodies.par_iter().enumerate().map(|(i, body_i)| {
        let mut pe_i = 0.0;
        for j in (i + 1)..bodies.len() {
             let body_j = &bodies[j]; let dist_sq = body_i.position.dist_sq(&body_j.position);
             let dist = (dist_sq + EPSILON_SQ).sqrt();
             if dist > 1e-9 { pe_i -= G * body_i.mass * body_j.mass / dist; }
        } pe_i
    }).sum::<f64>();
    kinetic_energy + potential_energy
}

// --- Main Simulation ---
fn main() {
    // --- Parameters ---
    const NUM_SMALL_BODIES: usize = 100_000; // Keep N manageable for testing
    // ... (other parameters same as before)
    const NUM_STEPS: usize = 1000;
    const DT: f64 = 0.001;

    println!("--- N-Body Simulation (Optimized Parallel Barnes-Hut k-d Tree) ---");
    println!("Theta (θ): {}, Theta^2: {:.3}", THETA, THETA_SQ);
    println!("Number of bodies: {}", NUM_SMALL_BODIES + 1);
    // ... (rest of print statements)
    let num_threads = rayon::current_num_threads();
    println!("Using {} CPU threads", num_threads);


    // --- Initialization ---
    println!("Initializing system...");
    let start_init = Instant::now();
    let mut bodies = initialize_solar_system(NUM_SMALL_BODIES, 1_000_000.0, 1.0, 100.0);
    println!("Initialization took: {:?}", start_init.elapsed());


    // --- Initial Energy (Accurate) ---
    println!("Calculating initial energy (accurate O(N^2) parallel)...");
    let start_e_init = Instant::now();
    let initial_energy = calculate_energy(&bodies); // Accurate O(N^2) calc
    println!("Initial Energy: {:.6e} (calculated in {:?})", initial_energy, start_e_init.elapsed());


    // --- Simulation Loop ---
    println!("Starting simulation...");
    let start_sim = Instant::now();
    // Remove detailed timing block for cleaner main loop, rely on avg step time
    for step in 0..NUM_STEPS {
        simulation_step(&mut bodies, DT);

        if (step + 1) % 50 == 0 || step == 0 {
            let elapsed_total = start_sim.elapsed();
            let avg_step_time = elapsed_total / (step + 1) as u32;
             println!(
                 "Step [{}/{}] | Sim Time: {:.3} | Avg Step: {:?} | Total Elapsed: {:?}",
                 step + 1, NUM_STEPS, (step + 1) as f64 * DT, avg_step_time, elapsed_total
             );
        }
    }
    let sim_duration = start_sim.elapsed();
    println!("Simulation finished.");
    println!("Total simulation loop took: {:?}", sim_duration);
    println!("Average time per step: {:?}", sim_duration / NUM_STEPS as u32);


    // --- Final Energy (Accurate) ---
    println!("Calculating final energy (accurate O(N^2) parallel)...");
    let start_e_final = Instant::now();
    let final_energy = calculate_energy(&bodies); // Accurate O(N^2) calc
    println!("Final Energy:   {:.6e} (calculated in {:?})", final_energy, start_e_final.elapsed());


    // --- Energy Conservation Check ---
    let energy_change = (final_energy - initial_energy).abs();
    let relative_energy_change = if initial_energy.abs() > 1e-9 { energy_change / initial_energy.abs() } else { energy_change };
    println!("Absolute Energy Change: {:.6e}", energy_change);
    println!("Relative Energy Change: {:.6e} ({:.4}%)", relative_energy_change, relative_energy_change * 100.0);
    println!("Note: Energy conservation reflects accuracy of Barnes-Hut (theta={}).", THETA);

}

// --- ordered-float module (keep as before or use crate) ---
// (Make sure this module is included if not using the external crate)
mod ordered_float { /* ... same as before ... */
    use std::cmp::Ordering;
    #[derive(Debug, Clone, Copy, PartialEq)] pub struct NotNan<T>(T);
    impl<T: PartialOrd> NotNan<T> { pub fn new(val: T) -> Result<Self, &'static str> { Ok(NotNan(val)) } }
    impl<T: PartialOrd> PartialOrd for NotNan<T> { fn partial_cmp(&self, other: &Self) -> Option<Ordering> { self.0.partial_cmp(&other.0) } }
    impl<T: PartialOrd> Eq for NotNan<T> {} impl<T: PartialOrd> Ord for NotNan<T> { fn cmp(&self, other: &Self) -> Ordering { self.partial_cmp(other).unwrap_or(Ordering::Equal) } }
}