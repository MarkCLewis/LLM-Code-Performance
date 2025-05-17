#ifndef KDTREE
#define KDTREE

#include "particle.h"
#include <immintrin.h> // For SIMD intrinsics

#define MAX_PARTS ((size_t)7)
#define THETA ((double)0.3)
#define ALIGN_TO_CACHE __attribute__((aligned(64))) // 64-byte alignment for cache lines

// Gravitational constant (in simulation units)
#define G ((double)1.0)

// Allocate memory aligned to 64-byte boundary for better SIMD performance
#define ALIGNED_ALLOC(size) aligned_alloc(64, (size))

typedef struct {
  // For leaves
  size_t num_parts;
  size_t particles[MAX_PARTS];

  // For internal nodes
  size_t split_dim;
  double split_val;
  double m;
  double cm[3];
  double size;
  size_t left;
  size_t right;
} ALIGN_TO_CACHE KDTree;

// Vector structure with SIMD-friendly alignment
typedef struct {
  double v[3] ALIGN_TO_CACHE;
} vect3;

void simple_sim(Particle_array_t *bodies, double dt, int steps);

typedef struct {
  size_t size;
  KDTree *ptr;
} KDTree_array_t;

KDTree_array_t new_kdtree_array_t(size_t elem_count);

typedef struct {
  size_t size;
  vect3 *ptr;
} vect3_array_t;

vect3_array_t new_vect3_array_t(size_t elem_count);

void KDTree_resize(KDTree_array_t *arr, size_t new_elem_count);
KDTree_array_t allocate_node_vec(size_t num_parts);

size_t build_tree(size_t_array_t *indices, size_t start, size_t end,
                  const Particle_array_t *particles, size_t cur_node,
                  KDTree_array_t *nodes);

// Optimized versions of key functions
void calc_pp_accel_optimized(const Particle *pi, const Particle *pj, double acc[3]);
void calc_accel_batch(size_t start, size_t end, const Particle_array_t *particles,
                     const KDTree_array_t *nodes, vect3_array_t *acc);

// Energy calculation for accuracy verification
typedef struct {
  double kinetic;    // Total kinetic energy
  double potential;  // Total potential energy
  double total;      // Total energy (kinetic + potential)
} SystemEnergy;

// Calculate total system energy (kinetic + potential)
SystemEnergy calculate_system_energy(const Particle_array_t *particles);

#endif
