/* Changes in this version:

build_kd_tree_parallel(): This is the entry point for the parallel kD-tree construction. It initializes the indices and bounding box and then calls the recursive parallel function.
build_kd_tree_recursive_parallel(): This is the recursive function that builds the kD-tree in parallel using rayon::join(). It splits the indices at the median and recursively calls itself for the left and right subtrees in parallel.
To use this parallel kD-tree construction:

Make sure you have the rayon dependency in your Cargo.toml.
Replace the line in kick_step() that calls the sequential build_kd_tree() with Self::build_kd_tree_parallel(self).
Now, the kD-tree will be constructed in parallel, which should significantly reduce the time spent building the tree, especially for a large number of bodies. The rest of the parallelization using Rayon for force calculation and the time-stepping remains the same.

Remember to compile and run in release mode (cargo build --release and then execute the binary) to see the performance benefits.
*/

impl System {
    // ... (rest of the System implementation from the previous response)

    fn build_kd_tree_parallel(&self) -> Option<Box<Self::KDNode>> {
        let num_bodies = self.bodies.len();
        if num_bodies == 0 {
            return None;
        }

        let mut indices: Vec<usize> = (0..num_bodies).collect();
        let min_bound = self.bodies.iter().fold(Vec3::new(f64::MAX, f64::MAX, f64::MAX), |acc, b| Vec3::new(acc.x.min(b.position.x), acc.y.min(b.position.y), acc.z.min(b.position.z)));
        let max_bound = self.bodies.iter().fold(Vec3::new(f64::MIN, f64::MIN, f64::MIN), |acc, b| Vec3::new(acc.x.max(b.position.x), acc.y.max(b.position.y), acc.z.max(b.position.z)));

        Self::build_kd_tree_recursive_parallel(self, &mut indices, min_bound, max_bound, 0)
    }

    fn build_kd_tree_recursive_parallel(
        system: &System,
        indices: &mut [usize],
        min_bound: Vec3,
        max_bound: Vec3,
        depth: usize,
    ) -> Option<Box<Self::KDNode>> {
        let num_indices = indices.len();

        if num_indices == 0 {
            return None;
        }

        let mut node = Self::KDNode {
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

        if node.total_mass > 0.0 {
            node.center_of_mass.x /= node.total_mass;
            node.center_of_mass.y /= node.total_mass;
            node.center_of_mass.z /= node.total_mass;
        }

        if num_indices == 1 {
            node.body_index = Some(indices[0]);
            return Some(Box::new(node));
        }

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
        let (mut left_indices, mut right_indices) = indices.split_at_mut(median_index);

        let mut left_max_bound = max_bound.clone();
        let mut right_min_bound = min_bound.clone();

        match split_dim {
            0 => {
                if let Some(&idx) = left_indices.last() {
                    left_max_bound.x = system.bodies[idx].position.x;
                }
                if let Some(&idx) = right_indices.first() {
                    right_min_bound.x = system.bodies[idx].position.x;
                }
            }
            1 => {
                if let Some(&idx) = left_indices.last() {
                    left_max_bound.y = system.bodies[idx].position.y;
                }
                if let Some(&idx) = right_indices.first() {
                    right_min_bound.y = system.bodies[idx].position.y;
                }
            }
            2 => {
                if let Some(&idx) = left_indices.last() {
                    left_max_bound.z = system.bodies[idx].position.z;
                }
                if let Some(&idx) = right_indices.first() {
                    right_min_bound.z = system.bodies[idx].position.z;
                }
            }
            _ => unreachable!(),
        }

        let (left_result, right_result) = rayon::join(
            || Self::build_kd_tree_recursive_parallel(system, left_indices, min_bound.clone(), left_max_bound, depth + 1),
            || Self::build_kd_tree_recursive_parallel(system, right_indices, right_min_bound, max_bound.clone(), depth + 1),
        );

        node.left = left_result;
        node.right = right_result;

        Some(Box::new(node))
    }

    fn kick_step(&mut self, dt: f64) {
        let kd_tree_option = Self::build_kd_tree_parallel(self);

        // Reset forces in parallel
        self.bodies.par_iter_mut().for_each(|body| {
            body.force.x = 0.0;
            body.force.y = 0.0;
            body.force.z = 0.0;
        });

        // Calculate forces using the kD-tree in parallel
        if let Some(kd_tree) = &kd_tree_option {
            self.bodies.par_iter_mut().for_each(|body| {
                Self::calculate_force_kd_tree(body, kd_tree, self);
            });
        }

        // Update velocities (kick) in parallel
        self.bodies.par_iter_mut().for_each(|body| {
            body.velocity.x += (body.force.x / body.mass) * (dt / 2.0);
            body.velocity.y += (body.force.y / body.mass) * (dt / 2.0);
            body.velocity.z += (body.force.z / body.mass) * (dt / 2.0);
        });
    }

    // ... (rest of the System implementation)
}

fn main() {
    // ... (rest of the main function)
}