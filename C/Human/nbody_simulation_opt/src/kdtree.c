#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <immintrin.h>
#include <stdint.h>  // Add for int64_t type

#include "kdtree.h"
#include "particle.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// Fast inverse square root approximation (Quake III technique)
inline double fast_inv_sqrt(double number) {
  double y = number;
  double x2 = y * 0.5;
  int64_t i = *(int64_t *)&y;
  i = 0x5fe6eb50c7b537a9 - (i >> 1);
  y = *(double *)&i;
  y = y * (1.5 - (x2 * y * y)); // 1st iteration
  return y;
}

// Calculate system energy for accuracy verification
SystemEnergy calculate_system_energy(const Particle_array_t *particles) {
    SystemEnergy energy = {0.0, 0.0, 0.0};
    
    // Calculate kinetic energy: Sum of (1/2) * m * v^2 for each particle
    double kinetic_sum = 0.0;
    #pragma omp parallel for reduction(+:kinetic_sum)
    for (size_t i = 0; i < particles->size; i++) {
        const Particle *p = &particles->ptr[i];
        double v_squared = p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2];
        kinetic_sum += 0.5 * p->m * v_squared;
    }
    energy.kinetic = kinetic_sum;
    
    // Calculate potential energy: Sum of G * m1 * m2 / r for all particle pairs
    // Note: We compute this pairwise to avoid double-counting
    double potential_sum = 0.0;
    #pragma omp parallel for reduction(+:potential_sum)
    for (size_t i = 0; i < particles->size; i++) {
        for (size_t j = i + 1; j < particles->size; j++) {
            const Particle *p1 = &particles->ptr[i];
            const Particle *p2 = &particles->ptr[j];
            
            // Compute distance between particles
            double dx = p1->p[0] - p2->p[0];
            double dy = p1->p[1] - p2->p[1];
            double dz = p1->p[2] - p2->p[2];
            double dist_squared = dx*dx + dy*dy + dz*dz;
            
            // Avoid division by zero
            if (dist_squared > 1e-10) {
                double dist = sqrt(dist_squared);
                potential_sum -= G * p1->m * p2->m / dist;
            }
        }
    }
    energy.potential = potential_sum;
    
    // Calculate total energy
    energy.total = energy.kinetic + energy.potential;
    
    return energy;
}

// For vectorization, using aligned loads/stores
const size_t NEGS[MAX_PARTS] __attribute__((aligned(64))) = 
                              {0, 0, 0, 0, 0, 0, 0};

KDTree_array_t allocate_node_vec(size_t num_parts) {
  size_t num_nodes = 2 * (num_parts / (MAX_PARTS - 1) + 1);
  KDTree_array_t ret = new_kdtree_array_t(num_nodes);
  return ret;
}

// Returns the index of the last Node used in the construction.
size_t build_tree(size_t_array_t *indices, size_t start, size_t end,
                  const Particle_array_t *particles, size_t cur_node,
                  KDTree_array_t *nodes) {
  size_t np = end - start;
  if (np <= MAX_PARTS) {
    if (cur_node >= nodes->size) {
      KDTree_resize(nodes, cur_node + 1);
    }
    nodes->ptr[cur_node].num_parts = np;
    for (size_t i = 0; i < np; ++i) {
      nodes->ptr[cur_node].particles[i] = indices->ptr[start + i];
    }
    return cur_node;
  } else {
    // Pick split dim and value
    double min_arr[3] = {1e100, 1e100, 1e100};
    double max_arr[3] = {-1e100, -1e100, -1e100};
    double m = 0.0;
    double cm[3] = {0.0, 0.0, 0.0};
    
    // Vectorized bounding box calculation
    #pragma omp parallel sections
    {
      #pragma omp section
      {
        for (size_t i = start; i < end; ++i) {
          size_t idx = indices->ptr[i];
          m += particles->ptr[idx].m;
          cm[0] += particles->ptr[idx].m * particles->ptr[idx].p[0];
          cm[1] += particles->ptr[idx].m * particles->ptr[idx].p[1];
          cm[2] += particles->ptr[idx].m * particles->ptr[idx].p[2];
        }
      }
      
      #pragma omp section
      {
        for (size_t i = start; i < end; ++i) {
          size_t idx = indices->ptr[i];
          min_arr[0] = MIN(min_arr[0], particles->ptr[idx].p[0]);
          min_arr[1] = MIN(min_arr[1], particles->ptr[idx].p[1]);
          min_arr[2] = MIN(min_arr[2], particles->ptr[idx].p[2]);
        }
      }
      
      #pragma omp section
      {
        for (size_t i = start; i < end; ++i) {
          size_t idx = indices->ptr[i];
          max_arr[0] = MAX(max_arr[0], particles->ptr[idx].p[0]);
          max_arr[1] = MAX(max_arr[1], particles->ptr[idx].p[1]);
          max_arr[2] = MAX(max_arr[2], particles->ptr[idx].p[2]);
        }
      }
    }
    
    cm[0] /= m;
    cm[1] /= m;
    cm[2] /= m;
    
    // Find dimension with largest spread for splitting
    size_t split_dim = 0;
    if (max_arr[1] - min_arr[1] > max_arr[split_dim] - min_arr[split_dim]) {
      split_dim = 1;
    }
    if (max_arr[2] - min_arr[2] > max_arr[split_dim] - min_arr[split_dim]) {
      split_dim = 2;
    }
    double size = max_arr[split_dim] - min_arr[split_dim];

    // Partition particles on split_dim - using median-of-three pivot
    size_t mid = (start + end) / 2;
    
    // Improved partitioning with deterministic pivot selection
    size_t first = start;
    size_t middle = (start + end) / 2;
    size_t last = end - 1;
    
    // Sort first, middle, last elements for a better pivot
    if (particles->ptr[indices->ptr[first]].p[split_dim] > 
        particles->ptr[indices->ptr[middle]].p[split_dim]) {
      size_t temp = indices->ptr[first];
      indices->ptr[first] = indices->ptr[middle];
      indices->ptr[middle] = temp;
    }
    
    if (particles->ptr[indices->ptr[middle]].p[split_dim] > 
        particles->ptr[indices->ptr[last]].p[split_dim]) {
      size_t temp = indices->ptr[middle];
      indices->ptr[middle] = indices->ptr[last];
      indices->ptr[last] = temp;
      
      if (particles->ptr[indices->ptr[first]].p[split_dim] > 
          particles->ptr[indices->ptr[middle]].p[split_dim]) {
        temp = indices->ptr[first];
        indices->ptr[first] = indices->ptr[middle];
        indices->ptr[middle] = temp;
      }
    }
    
    // Place pivot at position first
    size_t pivot_idx = middle;
    size_t pivot = indices->ptr[pivot_idx];
    indices->ptr[pivot_idx] = indices->ptr[first];
    indices->ptr[first] = pivot;
    
    // Partition around pivot
    size_t i = first + 1;
    size_t j = end - 1;
    
    while (i <= j) {
      while (i <= j && 
             particles->ptr[indices->ptr[i]].p[split_dim] <= 
             particles->ptr[pivot].p[split_dim]) {
        i++;
      }
      
      while (particles->ptr[indices->ptr[j]].p[split_dim] > 
             particles->ptr[pivot].p[split_dim]) {
        j--;
      }
      
      if (i < j) {
        size_t temp = indices->ptr[i];
        indices->ptr[i] = indices->ptr[j];
        indices->ptr[j] = temp;
      }
    }
    
    // Place pivot in its final position
    indices->ptr[first] = indices->ptr[j];
    indices->ptr[j] = pivot;
    
    // Adjust partitioning to ensure we get exactly the median
    if (j < mid) {
      // If the pivot is to the left of the desired median, recurse on right side
      size_t k = j + 1;
      while (k < mid && particles->ptr[indices->ptr[k]].p[split_dim] == 
                         particles->ptr[pivot].p[split_dim]) {
        k++;
      }
      if (k < mid) {
        // There are more elements on the right side that need to be partitioned
        // Recursively partition the right side
        j = k;
      }
    } else if (j > mid) {
      // If the pivot is to the right of the desired median, recurse on left side
      size_t k = j - 1;
      while (k > mid && particles->ptr[indices->ptr[k]].p[split_dim] == 
                         particles->ptr[pivot].p[split_dim]) {
        k--;
      }
      if (k > mid) {
        // There are more elements on the left side that need to be partitioned
        // Recursively partition the left side
        j = k;
      }
    }
    
    double split_val = particles->ptr[indices->ptr[mid]].p[split_dim];

    // Recurse on children and build this node.
    size_t left = build_tree(indices, start, mid, particles, cur_node + 1, nodes);
    size_t right = build_tree(indices, mid, end, particles, left + 1, nodes);

    if (cur_node >= nodes->size) {
      KDTree_resize(nodes, cur_node + 1);
    }
    nodes->ptr[cur_node].num_parts = 0;
    nodes->ptr[cur_node].split_dim = split_dim;
    nodes->ptr[cur_node].split_val = split_val;
    nodes->ptr[cur_node].m = m;
    nodes->ptr[cur_node].cm[0] = cm[0];
    nodes->ptr[cur_node].cm[1] = cm[1];
    nodes->ptr[cur_node].cm[2] = cm[2];
    nodes->ptr[cur_node].size = size;
    nodes->ptr[cur_node].left = cur_node + 1;
    nodes->ptr[cur_node].right = left + 1;

    return right;
  }
}

// Original particle-particle acceleration for reference
void calc_pp_accel(const Particle *pi, const Particle *pj, double acc[3]) {
  double dp[3] = {pi->p[0] - pj->p[0], pi->p[1] - pj->p[1],
                  pi->p[2] - pj->p[2]};
  double dist = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]);
  double magi = -pj->m / (dist * dist * dist);
  acc[0] += dp[0] * magi;
  acc[1] += dp[1] * magi;
  acc[2] += dp[2] * magi;
}

// Optimized particle-particle acceleration using SIMD where available
void calc_pp_accel_optimized(const Particle *pi, const Particle *pj, double acc[3]) {
  #ifdef __AVX__
  // Use AVX if available
  __m256d pi_p = _mm256_set_pd(0.0, pi->p[2], pi->p[1], pi->p[0]);
  __m256d pj_p = _mm256_set_pd(0.0, pj->p[2], pj->p[1], pj->p[0]);
  
  // Calculate dp vector
  __m256d dp = _mm256_sub_pd(pi_p, pj_p);
  
  // Calculate dp^2 and sum
  __m256d dp_squared = _mm256_mul_pd(dp, dp);
  
  // Horizontal add to get sum of squares (dist_squared)
  double dist_squared_array[4];
  _mm256_store_pd(dist_squared_array, dp_squared);
  double dist_squared = dist_squared_array[0] + dist_squared_array[1] + dist_squared_array[2];
  
  // Calculate 1/sqrt(dist_squared) faster
  double inv_dist = fast_inv_sqrt(dist_squared);
  double inv_dist_cubed = inv_dist * inv_dist * inv_dist;
  double magi = -pj->m * inv_dist_cubed;
  
  // Multiply dp by magi
  __m256d magi_v = _mm256_set1_pd(magi);
  __m256d result = _mm256_mul_pd(dp, magi_v);
  
  // Store result back
  double result_array[4];
  _mm256_store_pd(result_array, result);
  
  // Add to accumulator
  acc[0] += result_array[0];
  acc[1] += result_array[1];
  acc[2] += result_array[2];
  
  #else
  // Fallback to a optimized scalar version
  double dp[3] = {pi->p[0] - pj->p[0], pi->p[1] - pj->p[1], pi->p[2] - pj->p[2]};
  double dist_squared = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
  
  // Use fast inverse square root and avoid division
  double inv_dist = fast_inv_sqrt(dist_squared);
  double inv_dist_cubed = inv_dist * inv_dist * inv_dist;
  double magi = -pj->m * inv_dist_cubed;
  
  acc[0] += dp[0] * magi;
  acc[1] += dp[1] * magi;
  acc[2] += dp[2] * magi;
  #endif
}

// Optimized recursive acceleration calculation
void accel_recur(size_t cur_node, size_t p, const Particle_array_t *particles,
                 const KDTree_array_t *nodes, double acc[3]) {
  if (nodes->ptr[cur_node].num_parts > 0) {
    // Leaf node - direct calculation with all particles
    for (size_t i = 0; i < nodes->ptr[cur_node].num_parts; ++i) {
      if (nodes->ptr[cur_node].particles[i] != p) {
        calc_pp_accel_optimized(particles->ptr + p,
                      particles->ptr + nodes->ptr[cur_node].particles[i], acc);
      }
    }
  } else {
    // Internal node - check if we can use the center of mass approximation
    double dp[3];
    dp[0] = particles->ptr[p].p[0] - nodes->ptr[cur_node].cm[0];
    dp[1] = particles->ptr[p].p[1] - nodes->ptr[cur_node].cm[1];
    dp[2] = particles->ptr[p].p[2] - nodes->ptr[cur_node].cm[2];
    
    // Compute distance squared without sqrt for performance
    double dist_sqr = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
    double size_sqr = nodes->ptr[cur_node].size * nodes->ptr[cur_node].size;
    double theta_sqr = THETA * THETA;
    
    if (size_sqr < theta_sqr * dist_sqr) {
      // Use the center of mass approximation - node is far enough
      // Compute magi directly with fast_inv_sqrt
      double inv_dist = fast_inv_sqrt(dist_sqr);
      double inv_dist_cubed = inv_dist * inv_dist * inv_dist;
      double magi = -nodes->ptr[cur_node].m * inv_dist_cubed;
      
      acc[0] += dp[0] * magi;
      acc[1] += dp[1] * magi;
      acc[2] += dp[2] * magi;
    } else {
      // Traverse children
      accel_recur(nodes->ptr[cur_node].left, p, particles, nodes, acc);
      accel_recur(nodes->ptr[cur_node].right, p, particles, nodes, acc);
    }
  }
}

// Process a batch of particles for better cache efficiency
void calc_accel_batch(size_t start, size_t end, const Particle_array_t *particles,
                     const KDTree_array_t *nodes, vect3_array_t *acc) {
  const size_t batch_size = 16; // Process particles in small batches for better cache usage
  
  for (size_t b_start = start; b_start < end; b_start += batch_size) {
    size_t b_end = MIN(b_start + batch_size, end);
    
    // Prefetch particles data
    #pragma omp simd
    for (size_t i = b_start; i < b_end; i++) {
      __builtin_prefetch(&particles->ptr[i], 0, 3);
    }
    
    // Process batch
    #pragma omp parallel for
    for (size_t i = b_start; i < b_end; i++) {
      accel_recur(0, i, particles, nodes, acc->ptr[i].v);
    }
  }
}

void calc_accel(size_t p, const Particle_array_t *particles,
                const KDTree_array_t *nodes, double acc[3]) {
  accel_recur(0, p, particles, nodes, acc);
}

void print_tree(int step, const KDTree_array_t *tree,
                const Particle_array_t *particles) {

  if (step > 9999999) {
    fprintf(stderr, "Step too big!\n");
    exit(EXIT_FAILURE);
  }
  char name[8 + 7];
  snprintf(name, 8 + 7, "tree%d.txt", step);

  FILE *file = fopen(name, "w");
  if (file == NULL) {
    fprintf(stderr, "File couldn't be opened!\n");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "%lu\n", particles->size);

  for (size_t i_iter = 0; i_iter < tree->size; ++i_iter) {
    KDTree *n = tree->ptr + i_iter;
    if (n->num_parts > 0) {
      fprintf(file, "L %lu\n", n->num_parts);

      for (size_t i = 0; i < n->num_parts; ++i) {
        size_t p = n->particles[i];
        fprintf(file, "%f %f %f\n", particles->ptr[p].p[0],
                particles->ptr[p].p[1], particles->ptr[p].p[2]);
      }
    } else {
      fprintf(file, "I %lu %f %lu %lu\n", n->split_dim, n->split_val, n->left,
              n->right);
    }
  }

  fclose(file);
}

void simple_sim(Particle_array_t *bodies, double dt, int steps) {
  vect3_array_t acc = new_vect3_array_t(bodies->size);
  
  // Initialize acceleration vectors to zero
  #pragma omp parallel for simd
  for (size_t i = 0; i < bodies->size; ++i) {
    acc.ptr[i].v[0] = 0.0;
    acc.ptr[i].v[1] = 0.0;
    acc.ptr[i].v[2] = 0.0;
  }
  
  // Pre-allocate tree with sufficient size to avoid resizing
  KDTree_array_t tree = allocate_node_vec(bodies->size);
  size_t_array_t indices = new_range(0, bodies->size);
  
  // Determine optimal thread count and chunk size based on particle count
  int num_threads = omp_get_max_threads();
  int chunk_size = MAX(1, bodies->size / (num_threads * 4));
  omp_set_schedule(omp_sched_guided, chunk_size);

  for (int step = 0; step < steps; ++step) {
    // Reset indices
    #pragma omp parallel for simd
    for (size_t i = 0; i < bodies->size; ++i) {
      indices.ptr[i] = i;
    }

    // Build tree
    build_tree(&indices, 0, bodies->size, bodies, 0, &tree);
    
    // Calculate accelerations in batches for better cache usage
    calc_accel_batch(0, bodies->size, bodies, &tree, &acc);
    
    // Update velocities and positions
    #pragma omp parallel for
    for (size_t i = 0; i < bodies->size; ++i) {
      // Update velocities
      bodies->ptr[i].v[0] += dt * acc.ptr[i].v[0];
      bodies->ptr[i].v[1] += dt * acc.ptr[i].v[1];
      bodies->ptr[i].v[2] += dt * acc.ptr[i].v[2];
      
      // Update positions
      bodies->ptr[i].p[0] += dt * bodies->ptr[i].v[0];
      bodies->ptr[i].p[1] += dt * bodies->ptr[i].v[1];
      bodies->ptr[i].p[2] += dt * bodies->ptr[i].v[2];
      
      // Reset accelerations for next step
      acc.ptr[i].v[0] = 0.0;
      acc.ptr[i].v[1] = 0.0;
      acc.ptr[i].v[2] = 0.0;
    }
  }

  FREE_ARRAY(acc);
  FREE_ARRAY(tree);
  FREE_ARRAY(indices);
}

vect3_array_t new_vect3_array_t(size_t elem_count) {
  vect3_array_t a = {elem_count, (vect3 *)aligned_alloc(64, elem_count * sizeof(vect3))};

  if (a.ptr == NULL) {
    fprintf(stderr, "aligned_alloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }
  
  // Initialize to zero
  memset(a.ptr, 0, elem_count * sizeof(vect3));

  return a;
}

KDTree_array_t new_kdtree_array_t(size_t elem_count) {
  KDTree_array_t a = {elem_count, (KDTree *)aligned_alloc(64, elem_count * sizeof(KDTree))};

  if (a.ptr == NULL) {
    fprintf(stderr, "aligned_alloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }
  
  // Initialize to zero
  memset(a.ptr, 0, elem_count * sizeof(KDTree));

  return a;
}

void KDTree_resize(KDTree_array_t *arr, size_t new_elem_count) {
  if (new_elem_count == arr->size) {
    // do nothing
  } else if (new_elem_count > arr->size) {
    // Need aligned memory for better SIMD performance
    KDTree *new_ptr = (KDTree *)aligned_alloc(64, sizeof(KDTree) * new_elem_count);
    
    if (new_ptr == NULL) {
      fprintf(stderr, "aligned_alloc failed: it returned NULL\n");
      exit(EXIT_FAILURE);
    }
    
    // Copy old data
    memcpy(new_ptr, arr->ptr, arr->size * sizeof(KDTree));
    
    // Clear new space
    memset(new_ptr + arr->size, 0, (new_elem_count - arr->size) * sizeof(KDTree));
    
    // Free old buffer and update
    free(arr->ptr);
    arr->ptr = new_ptr;
    arr->size = new_elem_count;
  } else {
    // Shrink (this is less common, but may happen)
    KDTree *new_ptr = (KDTree *)aligned_alloc(64, sizeof(KDTree) * new_elem_count);
    
    if (new_ptr == NULL) {
      fprintf(stderr, "aligned_alloc failed: it returned NULL\n");
      exit(EXIT_FAILURE);
    }
    
    // Copy subset of data
    memcpy(new_ptr, arr->ptr, new_elem_count * sizeof(KDTree));
    
    // Free old buffer and update
    free(arr->ptr);
    arr->ptr = new_ptr;
    arr->size = new_elem_count;
  }
}
