use core::num;
use std::char::MAX;
use std::{fs::File, io::Write};

use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use rayon::{join, prelude::*};

use crate::array_particle::*;

use crate::quickstat::*;

//use super::array_particle::Particle;

// Add portable_simd feature when available
#[cfg(feature = "portable_simd")]
use std::simd::*;

const MAX_PARTS: usize = 8;
const THETA: f64 = 0.3;
const NEGS: [usize; MAX_PARTS] = [usize::MAX; MAX_PARTS];

// Adjust chunk size for better cache performance and parallel efficiency
const OPTIMAL_CHUNK_SIZE: usize = 8192;

#[derive(Clone, Copy)]
pub enum KDTree {
    Leaf {
        num_parts: usize,
        leaf_parts: [usize; MAX_PARTS]
    },

    Internal {
        split_dim: usize,
        split_val: f64,
        m: f64,
        cm: [f64; 3],
        size: f64,
        left: usize,
        right: usize
    }
}

impl KDTree {
    pub fn leaf<'a>(num_parts: usize, particles: [usize; MAX_PARTS]) -> KDTree {
        KDTree::Leaf {
            num_parts,
            leaf_parts: particles,
        }
    }
}

// Optimized SOA-based KD-tree
#[derive(Clone)]
pub struct KDTreeSoA {
    // Node storage using SOA pattern
    node_types: Vec<u8>, // 0 for leaf, 1 for internal
    split_dims: Vec<usize>,
    split_vals: Vec<f64>,
    masses: Vec<f64>,
    centers_of_mass: Vec<[f64; 3]>,
    sizes: Vec<f64>,
    lefts: Vec<usize>,
    rights: Vec<usize>,
    
    // Leaf storage
    leaf_counts: Vec<usize>,
    leaf_particles: Vec<Vec<usize>>,
    
    // Total number of nodes
    node_count: usize,
}

impl KDTreeSoA {
    pub fn new(capacity: usize) -> Self {
        let node_capacity = nodes_needed_for_particles(capacity);
        
        KDTreeSoA {
            node_types: vec![0; node_capacity],
            split_dims: vec![0; node_capacity],
            split_vals: vec![0.0; node_capacity],
            masses: vec![0.0; node_capacity],
            centers_of_mass: vec![[0.0, 0.0, 0.0]; node_capacity],
            sizes: vec![0.0; node_capacity],
            lefts: vec![0; node_capacity],
            rights: vec![0; node_capacity],
            
            leaf_counts: vec![0; node_capacity],
            leaf_particles: vec![Vec::new(); node_capacity],
            
            node_count: 0,
        }
    }
    
    // Build tree from system
    pub fn build(&mut self, system: &ParticleSystem) -> usize {
        let mut indices: Vec<usize> = (0..system.count).collect();
        self.node_count = 0;
        self.build_node(&mut indices, system)
    }
    
    fn build_node(&mut self, indices: &mut [usize], system: &ParticleSystem) -> usize {
        let node_idx = self.node_count;
        self.node_count += 1;
        
        if indices.len() <= MAX_PARTS {
            // Create leaf node
            self.node_types[node_idx] = 0;
            self.leaf_counts[node_idx] = indices.len();
            self.leaf_particles[node_idx] = indices.to_vec();
            return node_idx;
        }
        
        // Create internal node
        self.node_types[node_idx] = 1;
        
        // Calculate bounds and center of mass
        let mut min = [f64::MAX, f64::MAX, f64::MAX];
        let mut max = [f64::MIN, f64::MIN, f64::MIN];
        let mut m = 0.0;
        let mut cm = [0.0, 0.0, 0.0];
        
        for &i in indices.iter() {
            m += system.masses[i];
            cm[0] += system.masses[i] * system.positions[i][0];
            cm[1] += system.masses[i] * system.positions[i][1];
            cm[2] += system.masses[i] * system.positions[i][2];
            
            for d in 0..3 {
                min[d] = f64::min(min[d], system.positions[i][d]);
                max[d] = f64::max(max[d], system.positions[i][d]);
            }
        }
        
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        
        // Determine split dimension
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim;
            }
        }
        
        // Store node data
        self.masses[node_idx] = m;
        self.centers_of_mass[node_idx] = cm;
        self.split_dims[node_idx] = split_dim;
        self.sizes[node_idx] = max[split_dim] - min[split_dim];
        
        // Partition indices based on positions
        let mid = indices.len() / 2;
        quickstat_index(indices, mid, |i1, i2| {
            system.positions[i1][split_dim] < system.positions[i2][split_dim]
        });
        
        self.split_vals[node_idx] = system.positions[indices[mid]][split_dim];
        
        // Recursively build child nodes
        let (left_indices, right_indices) = indices.split_at_mut(mid);
        
        self.lefts[node_idx] = self.build_node(left_indices, system);
        self.rights[node_idx] = self.build_node(right_indices, system);
        
        node_idx
    }
    
    // Compute acceleration for a particle
    pub fn calc_accel(&self, p_idx: usize, system: &ParticleSystem) -> [f64; 3] {
        self.accel_recur(0, p_idx, system)
    }
    
    fn accel_recur(&self, node_idx: usize, p_idx: usize, system: &ParticleSystem) -> [f64; 3] {
        if node_idx >= self.node_count {
            return [0.0, 0.0, 0.0];
        }
        
        if self.node_types[node_idx] == 0 {
            // Leaf node - compute direct interactions
            let mut acc = [0.0, 0.0, 0.0];
            
            for &i in &self.leaf_particles[node_idx][0..self.leaf_counts[node_idx]] {
                if i != p_idx {
                    let pp_acc = calc_pp_accel_soa(p_idx, i, system);
                    acc[0] += pp_acc[0];
                    acc[1] += pp_acc[1];
                    acc[2] += pp_acc[2];
                }
            }
            
            acc
        } else {
            // Internal node - use Barnes-Hut approximation when possible
            let p = &system.positions[p_idx];
            let cm = &self.centers_of_mass[node_idx];
            
            let dx = p[0] - cm[0];
            let dy = p[1] - cm[1];
            let dz = p[2] - cm[2];
            let dist_sqr = dx*dx + dy*dy + dz*dz;
            
            if self.sizes[node_idx] * self.sizes[node_idx] < THETA * THETA * dist_sqr {
                // Far enough - use approximation
                let dist = f64::sqrt(dist_sqr);
                let magi = -self.masses[node_idx] / (dist_sqr * dist);
                [dx * magi, dy * magi, dz * magi]
            } else {
                // Too close - traverse children
                let left_acc = self.accel_recur(self.lefts[node_idx], p_idx, system);
                let right_acc = self.accel_recur(self.rights[node_idx], p_idx, system);
                [
                    left_acc[0] + right_acc[0],
                    left_acc[1] + right_acc[1],
                    left_acc[2] + right_acc[2]
                ]
            }
        }
    }
}

// Original functionality
fn nodes_needed_for_particles(num_parts: usize) -> usize {
    if num_parts <= MAX_PARTS {
        1
    } else {
        let min_num_leaves = num_parts / (MAX_PARTS / 2);
        let num_leaves = usize::pow(2, f64::log2(min_num_leaves as f64).ceil() as u32);
        2 * num_leaves - 1
    }
}

pub fn allocate_node_vec(num_parts: usize) -> Vec<KDTree> {
    let num_nodes = nodes_needed_for_particles(num_parts);
    let mut ret = Vec::new();
    ret.resize(num_nodes, KDTree::leaf(0, NEGS));
    ret
}

// Returns the index of the last Node used in the construction.
pub fn build_tree<'a>(
    indices: &mut Vec<usize>,
    start: usize,
    end: usize,
    particles: &Vec<Particle>,
    cur_node: usize,
    nodes: &mut Vec<KDTree>,
) -> usize {
    // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
    let np = end - start;
    // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
    if np <= MAX_PARTS {
        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[start + i]
        }
        nodes[cur_node] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
        cur_node
    } else {
        // Pick split dim and value
        let mut min = [1e100, 1e100, 1e100];
        let mut max = [-1e100, -1e100, -1e100];
        let mut m = 0.0;
        let mut cm = [0.0, 0.0, 0.0];
        for i in start..end {
            m += particles[indices[i]].m;
            cm[0] += particles[indices[i]].m * particles[indices[i]].p[0];
            cm[1] += particles[indices[i]].m * particles[indices[i]].p[1];
            cm[2] += particles[indices[i]].m * particles[indices[i]].p[2];
            min[0] = f64::min(min[0], particles[indices[i]].p[0]);
            min[1] = f64::min(min[1], particles[indices[i]].p[1]);
            min[2] = f64::min(min[2], particles[indices[i]].p[2]);
            max[0] = f64::max(max[0], particles[indices[i]].p[0]);
            max[1] = f64::max(max[1], particles[indices[i]].p[1]);
            max[2] = f64::max(max[2], particles[indices[i]].p[2]);
        }
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = (start + end) / 2;
        quickstat_index(&mut indices[start..end], mid - start, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let left = build_tree(indices, start, mid, particles, cur_node + 1, nodes);
        let right = build_tree(indices, mid, end, particles, left + 1, nodes);

        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        nodes[cur_node] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: left+1 };

        right
    }
}

// Returns the index of the last Node used in the construction.
pub fn build_tree_par1<'a>(
    indices: &mut Vec<usize>,
    start: usize,
    end: usize,
    particles: &Vec<Particle>,
    cur_node: usize,
    nodes: &mut Vec<KDTree>,
) -> usize {
    // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
    let np = end - start;
    // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
    if np <= MAX_PARTS {
        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[start + i]
        }
        nodes[cur_node] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
        cur_node
    } else {
        // Pick split dim and value
        let (m, cm_tmp, min, max) = (start..end).into_par_iter().fold(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
            |(m, cm, min, max), i| {
                (m + particles[indices[i]].m,
                [cm[0] + particles[indices[i]].m * particles[indices[i]].p[0],
                cm[1] + particles[indices[i]].m * particles[indices[i]].p[1],
                cm[2] + particles[indices[i]].m * particles[indices[i]].p[2]],
                [f64::min(min[0], particles[indices[i]].p[0]),
                f64::min(min[1], particles[indices[i]].p[1]),
                f64::min(min[2], particles[indices[i]].p[2])],
                [f64::max(max[0], particles[indices[i]].p[0]),
                f64::max(max[1], particles[indices[i]].p[1]),
                f64::max(max[2], particles[indices[i]].p[2])])
        }).reduce(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
            |(m1, cm1, min1, max1),(m2, cm2, min2, max2)| {
                (m1 + m2, 
                [cm1[0] + cm2[0], cm1[1] + cm2[1], cm1[2] + cm2[2]], 
                [f64::min(min1[0], min2[0]), f64::min(min1[1], min2[1]), f64::min(min1[2], min2[2])], 
                [f64::max(max1[0], max2[0]), f64::min(max1[1], max2[1]), f64::max(max1[2], max2[2])])
        });
        let mut cm = [cm_tmp[0], cm_tmp[1], cm_tmp[2]];
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = (start + end) / 2;
        quickstat_index(&mut indices[start..end], mid - start, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let left = build_tree_par1(indices, start, mid, particles, cur_node + 1, nodes);
        let right = build_tree_par1(indices, mid, end, particles, left + 1, nodes);

        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        nodes[cur_node] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: left+1 };

        right
    }
}

// Returns the index of the last Node used in the construction.
pub fn build_tree_par1_chunk<'a>(
    indices: &mut Vec<usize>,
    start: usize,
    end: usize,
    particles: &Vec<Particle>,
    cur_node: usize,
    nodes: &mut Vec<KDTree>,
) -> usize {
    // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
    let np = end - start;
    // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
    if np <= MAX_PARTS {
        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[start + i]
        }
        nodes[cur_node] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
        cur_node
    } else {
        // Pick split dim and value
        let (m, cm_tmp, min, max) = indices[start..end].par_chunks(500000).fold(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
            |(m, cm, min, max), chunk| {
                let mut min = [1e100, 1e100, 1e100];
                let mut max = [-1e100, -1e100, -1e100];
                let mut m = 0.0;
                let mut cm = [0.0, 0.0, 0.0];
                for i in chunk {
                    m += particles[*i].m;
                    cm[0] += particles[*i].m * particles[*i].p[0];
                    cm[1] += particles[*i].m * particles[*i].p[1];
                    cm[2] += particles[*i].m * particles[*i].p[2];
                    min[0] = f64::min(min[0], particles[*i].p[0]);
                    min[1] = f64::min(min[1], particles[*i].p[1]);
                    min[2] = f64::min(min[2], particles[*i].p[2]);
                    max[0] = f64::max(max[0], particles[*i].p[0]);
                    max[1] = f64::max(max[1], particles[*i].p[1]);
                    max[2] = f64::max(max[2], particles[*i].p[2]);
                }
                (m, cm, min, max)
        }).reduce(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
            |(m1, cm1, min1, max1),(m2, cm2, min2, max2)| {
                (m1 + m2, 
                [cm1[0] + cm2[0], cm1[1] + cm2[1], cm1[2] + cm2[2]], 
                [f64::min(min1[0], min2[0]), f64::min(min1[1], min2[1]), f64::min(min1[2], min2[2])], 
                [f64::max(max1[0], max2[0]), f64::min(max1[1], max2[1]), f64::max(max1[2], max2[2])])
        });
        let mut cm = [cm_tmp[0], cm_tmp[1], cm_tmp[2]];
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = (start + end) / 2;
        quickstat_index(&mut indices[start..end], mid - start, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let left = build_tree_par1_chunk(indices, start, mid, particles, cur_node + 1, nodes);
        let right = build_tree_par1_chunk(indices, mid, end, particles, left + 1, nodes);

        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        nodes[cur_node] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: left+1 };

        right
    }
}

// Returns the index of the last Node used in the construction.
pub fn build_tree_par2<'a>(
    indices: &mut Vec<usize>,
    buffer: &mut Vec<usize>,
    start: usize,
    end: usize,
    particles: &Vec<Particle>,
    cur_node: usize,
    nodes: &mut Vec<KDTree>,
) -> usize {
    // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
    let np = end - start;
    // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
    if np <= MAX_PARTS {
        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[start + i]
        }
        nodes[cur_node] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
        cur_node
    } else {
        // Pick split dim and value
        let mut min = [1e100, 1e100, 1e100];
        let mut max = [-1e100, -1e100, -1e100];
        let mut m = 0.0;
        let mut cm = [0.0, 0.0, 0.0];
        for i in start..end {
            m += particles[indices[i]].m;
            cm[0] += particles[indices[i]].m * particles[indices[i]].p[0];
            cm[1] += particles[indices[i]].m * particles[indices[i]].p[1];
            cm[2] += particles[indices[i]].m * particles[indices[i]].p[2];
            min[0] = f64::min(min[0], particles[indices[i]].p[0]);
            min[1] = f64::min(min[1], particles[indices[i]].p[1]);
            min[2] = f64::min(min[2], particles[indices[i]].p[2]);
            max[0] = f64::max(max[0], particles[indices[i]].p[0]);
            max[1] = f64::max(max[1], particles[indices[i]].p[1]);
            max[2] = f64::max(max[2], particles[indices[i]].p[2]);
        }
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = (start + end) / 2;
        quickstat_index_par(&mut indices[start..end], &mut buffer[start..end], mid - start, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let left = build_tree_par2(indices, buffer, start, mid, particles, cur_node + 1, nodes);
        let right = build_tree_par2(indices, buffer, mid, end, particles, left + 1, nodes);

        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        nodes[cur_node] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: left+1 };

        right
    }
}

pub fn build_tree_par3<'a>(
    indices: &mut [usize],
    cur_node: usize,
    particles: &Vec<Particle>,
    nodes: &mut [KDTree],
    thread_cnt: usize,
) {
    // println!("indices: {} cur: {} nodes: {}", indices.len(), cur_node, nodes.len());
    let np = indices.len();
    if np <= MAX_PARTS {
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[i]
        }
        nodes[0] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
    } else {
        // TODO: Test if this is faster chuncked.
        // Pick split dim and value
        let (m, cm_tmp, min, max) = (0..indices.len()).into_par_iter().fold(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
            |(m, cm, min, max), i| {
                (m + particles[indices[i]].m,
                [cm[0] + particles[indices[i]].m * particles[indices[i]].p[0],
                cm[1] + particles[indices[i]].m * particles[indices[i]].p[1],
                cm[2] + particles[indices[i]].m * particles[indices[i]].p[2]],
                [f64::min(min[0], particles[indices[i]].p[0]),
                f64::min(min[1], particles[indices[i]].p[1]),
                f64::min(min[2], particles[indices[i]].p[2])],
                [f64::max(max[0], particles[indices[i]].p[0]),
                f64::max(max[1], particles[indices[i]].p[1]),
                f64::max(max[2], particles[indices[i]].p[2])])
        }).reduce(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
            |(m1, cm1, min1, max1),(m2, cm2, min2, max2)| {
                (m1 + m2, 
                [cm1[0] + cm2[0], cm1[1] + cm2[1], cm1[2] + cm2[2]], 
                [f64::min(min1[0], min2[0]), f64::min(min1[1], min2[1]), f64::min(min1[2], min2[2])], 
                [f64::max(max1[0], max2[0]), f64::min(max1[1], max2[1]), f64::max(max1[2], max2[2])])
        });
        let mut cm = [cm_tmp[0], cm_tmp[1], cm_tmp[2]];
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = indices.len() / 2;
        quickstat_index(indices, mid, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let (left_indices, right_indices) = indices.split_at_mut(mid);
        let num_nodes = nodes_needed_for_particles(left_indices.len());
        let (_this_node, other_nodes) = nodes.split_at_mut(1);
        let (left_nodes, right_nodes) = other_nodes.split_at_mut(num_nodes);
        // println!("num_nodes = {}", num_nodes);
        let (_, _) = 
        if thread_cnt < num_cpus::get() {
            join(|| build_tree_par3(left_indices, cur_node + 1, particles, left_nodes, thread_cnt * 2),
            || build_tree_par3(right_indices, cur_node + 1 + num_nodes, particles, right_nodes, thread_cnt * 2))
        } else {
            (build_tree_par3(left_indices, cur_node + 1, particles, left_nodes, thread_cnt * 2),
            build_tree_par3(right_indices, cur_node + 1 + num_nodes, particles, right_nodes, thread_cnt * 2))
        }
        ;

        nodes[0] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: cur_node + 1 + num_nodes };
    }
}

pub fn build_tree_par3_chunk<'a>(
    indices: &mut [usize],
    cur_node: usize,
    particles: &Vec<Particle>,
    nodes: &mut [KDTree],
    thread_cnt: usize,
) {
    // println!("indices: {} cur: {} nodes: {}", indices.len(), cur_node, nodes.len());
    let np = indices.len();
    if np <= MAX_PARTS {
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[i]
        }
        nodes[0] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
    } else {
        // TODO: Test if this is faster chuncked.
        // Pick split dim and value
        let (m, cm_tmp, min, max) = indices.par_chunks(500000).fold(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
        |(m, cm, min, max), chunk| {
            let mut min = [1e100, 1e100, 1e100];
            let mut max = [-1e100, -1e100, -1e100];
            let mut m = 0.0;
            let mut cm = [0.0, 0.0, 0.0];
            for i in chunk {
                m += particles[*i].m;
                cm[0] += particles[*i].m * particles[*i].p[0];
                cm[1] += particles[*i].m * particles[*i].p[1];
                cm[2] += particles[*i].m * particles[*i].p[2];
                min[0] = f64::min(min[0], particles[*i].p[0]);
                min[1] = f64::min(min[1], particles[*i].p[1]);
                min[2] = f64::min(min[2], particles[*i].p[2]);
                max[0] = f64::max(max[0], particles[*i].p[0]);
                max[1] = f64::max(max[1], particles[*i].p[1]);
                max[2] = f64::max(max[2], particles[*i].p[2]);
            }
            (m, cm, min, max)
    }).reduce(|| (0.0, [0.0,0.0,0.0],[1e100, 1e100, 1e100],[-1e100, -1e100, -1e100]), 
        |(m1, cm1, min1, max1),(m2, cm2, min2, max2)| {
            (m1 + m2, 
            [cm1[0] + cm2[0], cm1[1] + cm2[1], cm1[2] + cm2[2]], 
            [f64::min(min1[0], min2[0]), f64::min(min1[1], min2[1]), f64::min(min1[2], min2[2])], 
            [f64::max(max1[0], max2[0]), f64::min(max1[1], max2[1]), f64::max(max1[2], max2[2])])
    });
    let mut cm = [cm_tmp[0], cm_tmp[1], cm_tmp[2]];
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = indices.len() / 2;
        quickstat_index(indices, mid, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let (left_indices, right_indices) = indices.split_at_mut(mid);
        let num_nodes = nodes_needed_for_particles(left_indices.len());
        let (_this_node, other_nodes) = nodes.split_at_mut(1);
        let (left_nodes, right_nodes) = other_nodes.split_at_mut(num_nodes);
        // println!("num_nodes = {}", num_nodes);
        let (_, _) = 
        if thread_cnt < num_cpus::get() {
            join(|| build_tree_par3_chunk(left_indices, cur_node + 1, particles, left_nodes, thread_cnt * 2),
            || build_tree_par3_chunk(right_indices, cur_node + 1 + num_nodes, particles, right_nodes, thread_cnt * 2))
        } else {
            (build_tree_par3_chunk(left_indices, cur_node + 1, particles, left_nodes, thread_cnt * 2),
            build_tree_par3_chunk(right_indices, cur_node + 1 + num_nodes, particles, right_nodes, thread_cnt * 2))
        }
        ;

        nodes[0] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: cur_node + 1 + num_nodes };
    }
}

pub fn build_tree_par4<'a>(
    indices: &mut [usize],
    cur_node: usize,
    particles: &Vec<Particle>,
    nodes: &mut [KDTree],
    thread_cnt: usize,
) {
    // println!("indices: {} cur: {} nodes: {}", indices.len(), cur_node, nodes.len());
    let np = indices.len();
    if np <= MAX_PARTS {
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[i]
        }
        nodes[0] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
    } else {
        // Pick split dim and value
        let mut min = [1e100, 1e100, 1e100];
        let mut max = [-1e100, -1e100, -1e100];
        let mut m = 0.0;
        let mut cm = [0.0, 0.0, 0.0];
        for i in 0..indices.len() {
            m += particles[indices[i]].m;
            cm[0] += particles[indices[i]].m * particles[indices[i]].p[0];
            cm[1] += particles[indices[i]].m * particles[indices[i]].p[1];
            cm[2] += particles[indices[i]].m * particles[indices[i]].p[2];
            min[0] = f64::min(min[0], particles[indices[i]].p[0]);
            min[1] = f64::min(min[1], particles[indices[i]].p[1]);
            min[2] = f64::min(min[2], particles[indices[i]].p[2]);
            max[0] = f64::max(max[0], particles[indices[i]].p[0]);
            max[1] = f64::max(max[1], particles[indices[i]].p[1]);
            max[2] = f64::max(max[2], particles[indices[i]].p[2]);
        }
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = indices.len() / 2;
        quickstat_index(indices, mid, 
            |i1, i2| particles[i1].p[split_dim] < particles[i2].p[split_dim]);
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let (left_indices, right_indices) = indices.split_at_mut(mid);
        let num_nodes = nodes_needed_for_particles(left_indices.len());
        let (_this_node, other_nodes) = nodes.split_at_mut(1);
        let (left_nodes, right_nodes) = other_nodes.split_at_mut(num_nodes);
        // println!("num_nodes = {}", num_nodes);
        let (_, _) = 
        if thread_cnt < num_cpus::get() {
            join(|| build_tree_par4(left_indices, cur_node + 1, particles, left_nodes, thread_cnt * 2),
            || build_tree_par4(right_indices, cur_node + 1 + num_nodes, particles, right_nodes, thread_cnt * 2))
        } else {
            (build_tree_par4(left_indices, cur_node + 1, particles, left_nodes, thread_cnt * 2),
            build_tree_par4(right_indices, cur_node + 1 + num_nodes, particles, right_nodes, thread_cnt * 2))
        }
        ;

        nodes[0] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: cur_node + 1 + num_nodes };
    }
}

fn accel_recur(cur_node: usize, p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> [f64; 3] {
    // println!("accel {}", cur_node);
    match nodes[cur_node] {
        KDTree::Leaf { num_parts, leaf_parts} => {
            let mut acc = [0.0, 0.0, 0.0];
            for i in 0..(num_parts) {
                if leaf_parts[i] != p {
                    let pp_acc = calc_pp_accel(&particles[p], &particles[leaf_parts[i]]);
                    acc[0] += pp_acc[0];
                    acc[1] += pp_acc[1];
                    acc[2] += pp_acc[2];
                }
            }
            acc
        }
        KDTree::Internal { m, cm, size, left, right, .. } => {
            let dx = particles[p].p[0] - cm[0];
            let dy = particles[p].p[1] - cm[1];
            let dz = particles[p].p[2] - cm[2];
            let dist_sqr = dx * dx + dy * dy + dz * dz;
            // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
            if size * size < THETA * THETA * dist_sqr {
                let dist = f64::sqrt(dist_sqr);
                let magi = -m / (dist_sqr * dist);
                [dx * magi, dy * magi, dz * magi]
            } else {
                let left_acc = accel_recur(left, p, particles, nodes);
                let right_acc = accel_recur(right, p, particles, nodes);
                [left_acc[0] + right_acc[0], left_acc[1] + right_acc[1], left_acc[2] + right_acc[2]]
            }
        }
    }
}

pub fn calc_accel(p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> [f64; 3] {
    accel_recur(0, p, particles, nodes)
}

pub fn simple_sim(bodies: &mut Vec<Particle>, dt: f64, steps: i64) {
    let initial_energy = calc_total_energy(bodies);
    println!("Initial total energy: {:.6e}", initial_energy);

    for step in 0..steps {
        let mut indices = Vec::new();
        for i in 0..bodies.len() {
            indices.push(i);
        }
        let mut nodes = allocate_node_vec(bodies.len());
        build_tree(&mut indices, 0, bodies.len(), bodies, 0, &mut nodes);

        // Calculate accelerations
        let mut accels = vec![[0.0, 0.0, 0.0]; bodies.len()];
        for i in 0..bodies.len() {
            accels[i] = calc_accel(i, bodies, &nodes);
        }

        // Update positions and velocities
        for i in 0..bodies.len() {
            bodies[i].v[0] += accels[i][0] * dt;
            bodies[i].v[1] += accels[i][1] * dt;
            bodies[i].v[2] += accels[i][2] * dt;
            bodies[i].p[0] += bodies[i].v[0] * dt;
            bodies[i].p[1] += bodies[i].v[1] * dt;
            bodies[i].p[2] += bodies[i].v[2] * dt;
        }
    }

    let final_energy = calc_total_energy(bodies);
    println!("Final total energy: {:.6e}", final_energy);
    println!("Energy change: {:.6e} ({:.2}%)", 
        final_energy - initial_energy,
        (final_energy - initial_energy).abs() / initial_energy.abs() * 100.0);
}

fn print_tree(step: i64, tree: &Vec<KDTree>, particles: &Vec<Particle>) -> std::io::Result<()> {
    let mut file = File::create(format!("tree{}.txt", step))?;

    file.write_fmt(format_args!("{}\n", tree.len()))?;
    for n in tree {
        match n {
            KDTree::Leaf { num_parts, leaf_parts } => {
                file.write_fmt(format_args!("L {}\n", num_parts))?;
                for i in 0..*num_parts {
                    let p = leaf_parts[i];
                    file.write_fmt(format_args!(
                        "{} {} {}\n",
                        particles[p].p[0], particles[p].p[1], particles[p].p[2]
                    ))?;
                }
            }
            KDTree::Internal { split_dim, split_val, left, right, .. } => {
                file.write_fmt(format_args!(
                    "I {} {} {} {}\n",
                    split_dim, split_val, left, right
                ))?;
            }
        }
    }

    Ok(())
}

// New optimized simulation function using SOA approach
pub fn simple_sim_soa(system: &mut ParticleSystem, dt: f64, steps: i64) {
    let initial_energy = calc_total_energy_soa(system);
    println!("Initial total energy: {:.6e}", initial_energy);
    
    // Create KD tree - only once, reuse for each step
    let mut kdtree = KDTreeSoA::new(system.count);
    
    // Temporary acceleration storage
    let mut accels = vec![[0.0, 0.0, 0.0]; system.count];
    
    for step in 0..steps {
        // Build tree
        kdtree.build(system);
        
        // Calculate accelerations in parallel
        accels.par_iter_mut().enumerate().for_each(|(i, acc)| {
            *acc = kdtree.calc_accel(i, system);
        });
        
        // Update velocities and positions
        for i in 0..system.count {
            // Update velocity first (using old position)
            system.velocities[i][0] += accels[i][0] * dt;
            system.velocities[i][1] += accels[i][1] * dt;
            system.velocities[i][2] += accels[i][2] * dt;
            
            // Then update position (using new velocity)
            system.positions[i][0] += system.velocities[i][0] * dt;
            system.positions[i][1] += system.velocities[i][1] * dt;
            system.positions[i][2] += system.velocities[i][2] * dt;
        }
    }
    
    let final_energy = calc_total_energy_soa(system);
    println!("Final total energy: {:.6e}", final_energy);
    println!("Energy change: {:.6e} ({:.2}%)", 
        final_energy - initial_energy,
        (final_energy - initial_energy).abs() / initial_energy.abs() * 100.0);
}

#[cfg(test)]
mod tests {
    use crate::{array_kd_tree, array_particle};

    #[test]
    fn single_node() {
        let parts = array_particle::two_bodies();
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        assert_eq!(node_vec.len(), 2);
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        match node_vec[0] {
            array_kd_tree::KDTree::Leaf { num_parts, .. } => assert_eq!(num_parts, parts.len()),
            _ => assert!(false, "Root isn't leaf of right size when small.")
        };
    }

    #[test]
    fn two_leaves() {
        let parts = array_particle::circular_orbits(11);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
        assert!(std::matches!(node_vec[0], array_kd_tree::KDTree::Internal { .. }));
        match (node_vec[1], node_vec[2]) {
            (array_kd_tree::KDTree::Leaf { num_parts: n1, ..}, array_kd_tree::KDTree::Leaf {num_parts: n2, ..}) => {
                assert_eq!(n1 + n2, 12);
            }
            _ => assert!(false, "Node vectors weren't leaves.")
        }
    }

    #[test]
    fn two_leaves_par3() {
        let parts = array_particle::circular_orbits(11);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree_par3(&mut indices, 0, &parts, &mut node_vec, 1);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
        assert!(std::matches!(node_vec[0], array_kd_tree::KDTree::Internal { .. }));
        match (node_vec[1], node_vec[2]) {
            (array_kd_tree::KDTree::Leaf { num_parts: n1, ..}, array_kd_tree::KDTree::Leaf {num_parts: n2, ..}) => {
                assert_eq!(n1 + n2, 12);
            }
            _ => assert!(false, "Node vectors weren't leaves.")
        }
    }

    #[test]
    fn big_solar() {
        let parts = array_particle::circular_orbits(5000);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }

    #[test]
    fn big_solar_par2() {
        let parts = array_particle::circular_orbits(5000);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        let mut buffer: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree_par2(&mut indices, &mut buffer, 0, parts.len(), &parts, 0, &mut node_vec);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }

    #[test]
    fn big_solar_par3() {
        let parts = array_particle::circular_orbits(5000);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree_par3(&mut indices,  0, &parts, &mut node_vec, 1);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }

    #[test]
    fn big_solar_par4() {
        let parts = array_particle::circular_orbits(5000);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree_par4(&mut indices,  0, &parts, &mut node_vec, 1);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }

    #[test]
    fn big_solar_with_steps() {
        let mut parts = array_particle::circular_orbits(5000);
        array_kd_tree::simple_sim(&mut parts, 1e-3, 10);

        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        let mut buffer: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree_par2(&mut indices, &mut buffer, 0, parts.len(), &parts, 0, &mut node_vec);
        recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }

    fn recur_test_tree_struct(
        node: usize,
        nodes: &Vec<array_kd_tree::KDTree>,
        particles: &Vec<array_particle::Particle>,
        mut min: [f64; 3],
        mut max: [f64; 3],
    ) {
        match nodes[node] {
            array_kd_tree::KDTree::Leaf { num_parts, leaf_parts } => {
                for index in 0..num_parts {
                    let i = leaf_parts[index];
                    for dim in 0..2 {
                        assert!(
                            particles[i].p[dim] >= min[dim],
                            "Particle dim {} is below min. i={} p={} min={}",
                            dim,
                            i,
                            particles[i].p[dim],
                            min[dim]
                        );
                        assert!(
                            particles[i].p[dim] < max[dim],
                            "Particle dim {} is above max. i={} p={} max={}",
                            dim,
                            i,
                            particles[i].p[dim],
                            max[dim]
                        );
                    }
                }
            }
            array_kd_tree::KDTree::Internal { split_dim, split_val, left, right, .. } => {
                let split_dim = split_dim;
                let tmin = min[split_dim];
                let tmax = max[split_dim];
                max[split_dim] = split_val;
                recur_test_tree_struct(left, nodes, particles, min, max);
                max[split_dim] = tmax;
                min[split_dim] = split_val;
                recur_test_tree_struct(right, nodes, particles, min, max);
                min[split_dim] = tmin;
            }
        }
    }
}
