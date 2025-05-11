#ifndef PARTICLE
#define PARTICLE

#include <stddef.h>

#define FREE_ARRAY(arr)                                                        \
  if ((arr).ptr != NULL) {                                                     \
    free((arr).ptr);                                                           \
    (arr).size = 0;                                                            \
    (arr).ptr = NULL;                                                          \
  }

// Define cache line size for alignment
#define CACHE_LINE_SIZE 64

// Particle structure aligned to cache line size
typedef struct {
  // Position and velocity vectors kept together for vectorization
  double p[3];
  double v[3];
  double r;
  double m;
  double pad[1]; // Padding to ensure 64-byte alignment
} __attribute__((aligned(CACHE_LINE_SIZE))) Particle;

typedef struct {
  size_t size;
  Particle *ptr;
} Particle_array_t;

Particle_array_t new_array_t(size_t elem_count);

typedef struct {
  size_t size;
  size_t *ptr;
} size_t_array_t;

size_t_array_t new_size_t_array_t(size_t elem_count);
size_t_array_t new_range(size_t start, size_t end);

Particle_array_t two_bodies();
Particle_array_t circular_orbits(size_t n);

// Optimized memory allocation for particles
Particle_array_t new_array_aligned_t(size_t elem_count);

#endif
