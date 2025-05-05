use rayon::prelude::*;
use rand::Rng;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, Index, IndexMut}; // Added Index, IndexMut
use std::time::Instant;
use std::cmp::Ordering; // For sorting

// --- Constants ---
const G: f64 = 1.0;
const DT: f64 = 0.001;
const SOFTENING: f64 = 1e-5; // Softening for both direct and approximate forces
const THETA: f64 = 0.3; // Opening angle criterion for Barnes-Hut
const THETA2: f64 = THETA * THETA; // Pre-calculate squared theta

// --- Vec3 Struct ---
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

// Implement Indexing for Vec3 to access x, y, z by 0, 1, 2
impl Index<usize> for Vec3 {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Invalid Vec3 index"),
        }
    }
}

impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Invalid Vec3 index"),
        }
    }
}

// Vec3 methods and operator overloading (same as before)
impl Vec3 { /* ... methods new, zero, magnitude_squared, magnitude ... */
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

// --- Body Struct --- (same as before)
#[derive(Debug, Clone, Copy)]
struct Body {
    mass: f64,
    pos: Vec3,
    vel: Vec3,
}
unsafe impl Send for Body {}
unsafe impl Sync for Body {}
impl Body { fn new(mass: f64, pos: Vec3, vel: Vec3) -> Self { Body { mass, pos, vel } } }


// --- kD-Tree Structures ---

#[derive(Debug, Clone, Copy)]
struct BoundingBox {
    min: Vec3,
    max: Vec3,
}

impl BoundingBox {
    // Create an empty bounding box (min=infinity, max=-infinity)
    fn empty() -> Self {
        BoundingBox {
            min: Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            max: Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }

    // Extend the box to include a point
    fn extend(&mut self, point: Vec3) {
        self.min.x = self.min.x.min(point.x);
        self.min.y = self.min.y.min(point.y);
        self.min.z = self.min.z.min(point.z);
        self.max.x = self.max.x.max(point.x);
        self.max.y = self.max.y.max(point.y);
        self.max.z = self.max.z.max(point.z);
    }

    // Extend the box to include another box
    fn extend_box(&mut self, other: &BoundingBox) {
        self.extend(other.min);
        self.extend(other.max);
    }

    // Calculate the maximum dimension (width, height, or depth)
    fn max_dim(&self) -> f64 {
        (self.max.x - self.min.x)
            .max(self.max.y - self.min.y)
            .max(self.max.z - self.min.z)
    }
}


#[derive(Debug)]
enum KdNode {
    Internal {
        bounds: BoundingBox,
        axis: usize, // 0=x, 1=y, 2=z
        // split_val: f64, // Not strictly needed if we store children
        total_mass: f64,
        center_of_mass: Vec3,
        left: Box<KdNode>,
        right: Box<KdNode>,
    },
    Leaf {
        body_index: usize, // Index into the main bodies array
        pos: Vec3, // Store position for convenience
        mass: f64, // Store mass for convenience
    },
    Empty,
}

impl KdNode {
    fn bounds(&self) -> Option<BoundingBox> {
        match self {
            KdNode::Internal { bounds, .. } => Some(*bounds),
            KdNode::Leaf { pos, .. } => Some(BoundingBox { min: *pos, max: *pos }),
            KdNode::Empty => None,
        }
    }

     fn total_mass(&self) -> f64 {
        match self {
            KdNode::Internal { total_mass, .. } => *total_mass,
            KdNode::Leaf { mass, .. } => *mass,
            KdNode::Empty => 0.0,
        }
    }

    fn center_of_mass(&self) -> Option<Vec3> {
         match self {
            KdNode::Internal { center_of_mass, .. } => Some(*center_of_mass),
            KdNode::Leaf { pos, .. } => Some(*pos), // CoM of a leaf is its position
            KdNode::Empty => None,
        }
    }
}


// --- kD-Tree Construction ---

// Helper function to build the tree recursively
// Takes mutable slice of *indices* into the bodies array
fn build_kdtree_recursive(
    indices: &mut [usize],
    bodies: &[Body],
    depth: usize,
) -> KdNode {
    let n = indices.len();

    if n == 0 {
        return KdNode::Empty;
    }

    if n == 1 {
        let idx = indices[0];
        return KdNode::Leaf { body_index: idx, pos: bodies[idx].pos, mass: bodies[idx].mass };
    }

    // Select axis based on depth
    let axis = depth % 3;

    // Sort indices based on the position along the current axis
    // This is the most expensive part of the build
    indices.sort_unstable_by(|&a, &b| {
        bodies[a].pos[axis]
            .partial_cmp(&bodies[b].pos[axis])
            .unwrap_or(Ordering::Equal)
    });

    // Choose median index for partitioning
    let median_idx = n / 2;
    // let median_body_idx = indices[median_idx]; // Body index at the median split

    // Recursively build left and right subtrees
    // Note: We split the *indices* slice
    let (left_indices, right_indices_with_median) = indices.split_at_mut(median_idx);
    let (median_slice, right_indices) = right_indices_with_median.split_at_mut(1); // Take median out for now
    let median_body_idx = median_slice[0];


    // Build subtrees (could potentially be parallelized with Rayon's join, but adds complexity)
    let left_child = build_kdtree_recursive(left_indices, bodies, depth + 1);
    let right_child = build_kdtree_recursive(right_indices, bodies, depth + 1);

    // Create a leaf node for the median body we split on
    let median_node = KdNode::Leaf {
        body_index: median_body_idx,
        pos: bodies[median_body_idx].pos,
        mass: bodies[median_body_idx].mass
    };

    // Combine left + median, then combine (left+median) + right to get CoM and bounds
    // This structure is slightly different from pure kD-split but common in BH kD-trees

    // Combine left and median first
    let (combined_mass_lm, combined_com_lm, combined_bounds_lm) =
        combine_nodes(&left_child, &median_node);

    // Now combine the (left+median) result with the right child
    let total_mass: f64;
    let center_of_mass: Vec3;
    let mut bounds: BoundingBox;

    match combine_nodes_precalculated(combined_mass_lm, combined_com_lm, combined_bounds_lm, &right_child) {
         (m, com, b) => {
            total_mass = m;
            center_of_mass = com;
            bounds = b;
         }
    }

    // If bounds are still invalid (e.g., only empty children), compute from median
     if !bounds.min.x.is_finite() {
        bounds = BoundingBox { min: bodies[median_body_idx].pos, max: bodies[median_body_idx].pos };
        if let Some(lb) = left_child.bounds() { bounds.extend_box(&lb); }
        if let Some(rb) = right_child.bounds() { bounds.extend_box(&rb); }
    }


    // Need to handle the median body. A common approach is to rebuild the tree
    // including the median body in one of the children during the recursive call,
    // or to handle 3-way split (left, median, right).
    // Let's reconstruct the node using the children generated *including* the median in the right branch
    // Simpler approach: median splits the space, bodies <= go left, > go right.
    // Let's revert to the simpler split approach:
    let median_spatial_idx = n / 2; // Index in the *sorted* indices slice
    let (left_indices_split, right_indices_split) = indices.split_at_mut(median_spatial_idx);

    // Recursively build left and right subtrees using the simple split
    let left_child_final = build_kdtree_recursive(left_indices_split, bodies, depth + 1);
    let right_child_final = build_kdtree_recursive(right_indices_split, bodies, depth + 1);


    // Calculate combined properties
    let (total_mass_final, center_of_mass_final, bounds_final) = combine_nodes(&left_child_final, &right_child_final);

     // Ensure bounds are valid if children were empty
    let final_bounds = if total_mass_final > 0.0 && bounds_final.min.x.is_finite() {
        bounds_final
    } else if let Some(lb) = left_child_final.bounds() { // If only left exists
         lb
    } else if let Some(rb) = right_child_final.bounds() { // If only right exists
        rb
    } else {
        BoundingBox::empty() // Should not happen if n > 0
    };


    KdNode::Internal {
        bounds: final_bounds,
        axis,
        // split_val: bodies[indices[median_spatial_idx]].pos[axis], // Value used for split
        total_mass: total_mass_final,
        center_of_mass: center_of_mass_final,
        left: Box::new(left_child_final),
        right: Box::new(right_child_final),
    }
}

// Helper to combine mass, CoM, and bounds of two nodes
fn combine_nodes(node1: &KdNode, node2: &KdNode) -> (f64, Vec3, BoundingBox) {
    let m1 = node1.total_mass();
    let m2 = node2.total_mass();
    let com1_opt = node1.center_of_mass();
    let com2_opt = node2.center_of_mass();
    let bounds1_opt = node1.bounds();
    let bounds2_opt = node2.bounds();

    let total_mass = m1 + m2;
    let center_of_mass = if total_mass == 0.0 {
        Vec3::zero()
    } else {
        let com1 = com1_opt.unwrap_or(Vec3::zero());
        let com2 = com2_opt.unwrap_or(Vec3::zero());
        (com1 * m1 + com2 * m2) / total_mass
    };

    let mut bounds = bounds1_opt.unwrap_or_else(BoundingBox::empty);
    if let Some(b2) = bounds2_opt {
        bounds.extend_box(&b2);
    }

    (total_mass, center_of_mass, bounds)
}

// Helper when one side is pre-calculated
fn combine_nodes_precalculated(m1: f64, com1: Vec3, bounds1: BoundingBox, node2: &KdNode) -> (f64, Vec3, BoundingBox) {
     let m2 = node2.total_mass();
     let com2_opt = node2.center_of_mass();
     let bounds2_opt = node2.bounds();

     let total_mass = m1 + m2;
     let center_of_mass = if total_mass == 0.0 {
         Vec3::zero()
     } else {
         let com2 = com2_opt.unwrap_or(Vec3::zero());
         (com1 * m1 + com2 * m2) / total_mass
     };

     let mut bounds = bounds1;
     if let Some(b2) = bounds2_opt {
         bounds.extend_box(&b2);
     }
     (total_mass, center_of_mass, bounds)
}


// --- Force Calculation using kD-Tree ---

// Recursive helper function to calculate force on a target body by traversing the tree
fn calculate_force_recursive(
    node: &KdNode,
    target_body_idx: usize,
    target_pos: Vec3,
    target_mass: f64, // Mass of the target body (for softening maybe?)
    bodies: &[Body],
    theta2: f64,
) -> Vec3 {
    match node {
        KdNode::Empty => Vec3::zero(),
        KdNode::Leaf { body_index, pos, mass } => {
            if *body_index == target_body_idx {
                Vec3::zero() // Don't interact with self
            } else {
                // Direct force calculation
                let dr = *pos - target_pos;
                let dist_sq = dr.magnitude_squared() + SOFTENING * SOFTENING; // Softening
                let dist_pow_1_5 = dist_sq.powf(1.5);
                if dist_pow_1_5 == 0.0 { return Vec3::zero(); } // Avoid division by zero

                // let force_factor = G * (*mass) / dist_pow_1_5; // Force = G * m_leaf * dr / dist^3
                // dr * force_factor * target_mass // We need ACCELERATION so multiply by G*m_leaf / dist^3
                // F_target = G * m_target * m_leaf * dr / dist^3
                // a_target = F_target / m_target = G * m_leaf * dr / dist^3
                // Corrected: We want acceleration on target, so don't multiply by target_mass here
                let accel_factor = G * (*mass) / dist_pow_1_5;
                dr * accel_factor
            }
        }
        KdNode::Internal { bounds, total_mass, center_of_mass, left, right, .. } => {
            let dr = *center_of_mass - target_pos;
            let dist_sq = dr.magnitude_squared();

            // If target is exactly at CoM (unlikely but possible), handle carefully.
            // Maybe open the node? Or apply softening... softened distance:
            let dist_sq_softened = dist_sq + SOFTENING * SOFTENING;
             if dist_sq_softened == 0.0 { return Vec3::zero(); } // Avoid division by zero


            let s = bounds.max_dim();
            let s_squared = s * s;

            // Barnes-Hut criterion: s^2 / d^2 < theta^2
            if s_squared / dist_sq_softened < theta2 {
                // Use approximation: treat node as single particle
                let dist_pow_1_5 = dist_sq_softened.powf(1.5);
                 if dist_pow_1_5 == 0.0 { return Vec3::zero(); } // Avoid division by zero
                // accel = G * M_node * dr / dist^3
                let accel_factor = G * (*total_mass) / dist_pow_1_5;
                dr * accel_factor
            } else {
                // Node is too close/large, traverse children
                let force_left = calculate_force_recursive(&*left, target_body_idx, target_pos, target_mass, bodies, theta2);
                let force_right = calculate_force_recursive(&*right, target_body_idx, target_pos, target_mass, bodies, theta2);
                force_left + force_right
            }
        }
    }
}


// --- Simulation Step using kD-Tree ---
fn simulate_step_kdtree(bodies: &mut [Body], dt: f64) {
    let n = bodies.len();
    if n == 0 { return; }

    // 0. Build kD-Tree (needs to be rebuilt every step)
    let build_start = Instant::now();
    let mut indices: Vec<usize> = (0..n).collect();
    let tree = build_kdtree_recursive(&mut indices, bodies, 0);
    let build_duration = build_start.elapsed();
    // println!("Tree build time: {:?}", build_duration); // Optional: timing


    // 1. Calculate Accelerations using the kD-Tree (Parallelized)
    let force_start = Instant::now();
    let accelerations: Vec<Vec3> = (0..n)
        .into_par_iter()
        .map(|i| {
            calculate_force_recursive(&tree, i, bodies[i].pos, bodies[i].mass, bodies, THETA2)
        })
        .collect();
     let force_duration = force_start.elapsed();
     // println!("Force calc time: {:?}", force_duration); // Optional: timing


    // 2. Kick: Update Velocities (Parallelized)
    bodies
        .par_iter_mut()
        .zip(accelerations.par_iter())
        .for_each(|(body, acc)| {
            body.vel += *acc * dt;
        });

    // 3. Step: Update Positions (Parallelized)
    bodies
        .par_iter_mut()
        .for_each(|body| {
            body.pos += body.vel * dt;
        });
}


// --- Exact Energy Calculation (O(N^2), unchanged, for verification) ---
fn calculate_total_energy(bodies: &[Body]) -> f64 {
    let n = bodies.len();
    let kinetic_energy: f64 = bodies
        .par_iter()
        .map(|b| 0.5 * b.mass * b.vel.magnitude_squared())
        .sum();

    let potential_energy: f64 = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut local_potential = 0.0;
            for j in (i + 1)..n {
                let dr = bodies[j].pos - bodies[i].pos;
                let dist_sq = dr.magnitude_squared();
                let dist = (dist_sq + SOFTENING * SOFTENING).sqrt();
                 if dist == 0.0 { continue } // Avoid division by zero if softening is zero
                local_potential -= G * bodies[i].mass * bodies[j].mass / dist;
            }
            local_potential
        })
        .sum();

    kinetic_energy + potential_energy
}


// --- Initialization (Unchanged) ---
fn initialize_circular_orbits( /* ... */ ) -> Vec<Body> { /* ... Same as before ... */
    let num_small_bodies: usize = 1_000_000; // Example parameter
    let central_mass: f64 = 1_000_000.0;
    let small_mass: f64 = 1.0;
    let avg_radius: f64 = 100.0;
    let radius_spread: f64 = 20.0;


    let mut bodies = Vec::with_capacity(num_small_bodies + 1);
    let central_body = Body::new(central_mass, Vec3::zero(), Vec3::zero());
    bodies.push(central_body);

    let small_bodies: Vec<Body> = (0..num_small_bodies)
        .into_par_iter() // Parallel iterator over the count
        .map(|_| {
            let mut rng = rand::thread_rng();
            let r_mag = avg_radius + radius_spread * (rng.gen::<f64>() - 0.5) * 2.0;
            let theta = rng.gen::<f64>() * 2.0 * PI;
            let phi = (rng.gen::<f64>() * 2.0 - 1.0).acos();
            let x = r_mag * phi.sin() * theta.cos();
            let y = r_mag * phi.sin() * theta.sin();
            let z = r_mag * phi.cos();
            let pos = Vec3::new(x, y, z);
            let speed = (G * central_mass / r_mag.max(1e-6)).sqrt(); // Avoid division by zero if r_mag is tiny
            let mut vel_dir = Vec3::new(-y, x, 0.0);
            if vel_dir.magnitude_squared() < 1e-12 { vel_dir = Vec3::new(1.0, 0.0, 0.0); }
            let vel = vel_dir * (speed / vel_dir.magnitude().max(1e-6)); // Avoid division by zero
            Body::new(small_mass, pos, vel)
        })
        .collect();
    bodies.extend(small_bodies);
    bodies
}


// --- Main Execution ---
fn main() {
    let n_small_bodies = 1_000_000;
    // let n_small_bodies = 10_000; // Use a smaller number for faster testing
    let n_steps = 1000;
    let central_mass = 1_000_000.0;
    let small_mass = 1.0;
    let orbital_radius = 100.0;
    let radius_spread = 20.0;

    println!("N-Body Simulation (Multithreaded with kD-Tree / Barnes-Hut)"); // Updated title
    println!("-----------------");
    println!("Theta: {}", THETA);
    // ... (rest of parameter printing) ...
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


    println!("Initializing system (parallel)...");
    let start_init = Instant::now();
    // Pass parameters directly to initializer
    let mut bodies = initialize_circular_orbits(
        // n_small_bodies, central_mass, small_mass, orbital_radius, radius_spread
    );
    let init_duration = start_init.elapsed();
    println!("Initialization complete ({:.2?} total bodies) in {:.2?}", bodies.len(), init_duration);


    println!("Calculating initial energy (parallel O(N^2))..."); // Still exact calculation
    let start_energy = Instant::now();
    let initial_energy = calculate_total_energy(&bodies);
    let energy_duration = start_energy.elapsed();
    println!("Initial Total Energy: {:.6e} (calculated in {:.2?})", initial_energy, energy_duration);


    println!("Starting simulation using kD-Tree...");
    let start_sim = Instant::now();
    for step in 0..n_steps {
        let step_start = Instant::now();
        simulate_step_kdtree(&mut bodies, DT); // Use the kD-Tree step function
        let step_duration = step_start.elapsed();

        if (step + 1) % 10 == 0 || step == 0 {
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
             // Energy check less frequently, as it's O(N^2) and simulation is approximate
             if (step + 1) % 100 == 0 {
                  let energy_check_start = Instant::now();
                  let current_energy = calculate_total_energy(&bodies);
                  let energy_check_duration = energy_check_start.elapsed();
                  println!(
                      "  Energy (@step {}): {:.6e} (Check took: {:?})",
                      step + 1, current_energy, energy_check_duration
                 );
             }
        }
    }
    let sim_duration = start_sim.elapsed();
    println!("Simulation finished in {:.2?}", sim_duration);


    println!("Calculating final energy (parallel O(N^2))..."); // Still exact calculation
    let start_energy = Instant::now();
    let final_energy = calculate_total_energy(&bodies);
    let energy_duration = start_energy.elapsed();
    println!("Final Total Energy:   {:.6e} (calculated in {:.2?})", final_energy, energy_duration);

    // Energy comparison will now show drift due to approximation!
    let energy_diff = final_energy - initial_energy;
    let relative_energy_diff = if initial_energy.abs() > 1e-12 {
        (energy_diff / initial_energy).abs()
    } else {
        energy_diff.abs()
    };
    println!("Absolute Energy Change: {:.6e}", energy_diff);
    println!("Relative Energy Change: {:.6e} ({:.4}%)", relative_energy_diff, relative_energy_diff * 100.0);
    println!("NOTE: Energy drift is expected with the approximate Barnes-Hut method.");
}