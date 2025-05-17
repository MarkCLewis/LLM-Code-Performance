#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle.h"

Particle_array_t two_bodies() {
  Particle_array_t bodies = new_array_aligned_t(2);
  Particle p0 = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1.0, {0.0}};
  Particle p1 = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1e-4, 1e-20, {0.0}};
  bodies.ptr[0] = p0;
  bodies.ptr[1] = p1;
  return bodies;
}

Particle_array_t new_array_t(size_t elem_count) {
  Particle_array_t a = {elem_count,
                        (Particle *)calloc(elem_count, sizeof(Particle))};

  if (a.ptr == NULL) {
    fprintf(stderr, "calloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }

  return a;
}

// Aligned allocation for better SIMD performance
Particle_array_t new_array_aligned_t(size_t elem_count) {
  Particle_array_t a;
  a.size = elem_count;
  
  // Allocate aligned memory
  a.ptr = (Particle *)aligned_alloc(CACHE_LINE_SIZE, 
                                   elem_count * sizeof(Particle));
  
  if (a.ptr == NULL) {
    fprintf(stderr, "aligned_alloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }
  
  // Initialize memory to zero
  memset(a.ptr, 0, elem_count * sizeof(Particle));
  
  return a;
}

size_t_array_t new_size_t_array_t(size_t elem_count) {
  size_t_array_t a = {elem_count, (size_t *)calloc(elem_count, sizeof(size_t))};

  if (a.ptr == NULL) {
    fprintf(stderr, "calloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }

  return a;
}

size_t_array_t new_range(size_t start, size_t end) {
  if (end <= start) {
    size_t_array_t a = {0, (size_t *)NULL};
    return a;
  }

  size_t_array_t a = new_size_t_array_t(end - start);

  for (size_t i = 0; i < end - start; ++i) {
    a.ptr[i] = i + start;
  }

  return a;
}

double random_d() { return (double)rand() / RAND_MAX; }

Particle_array_t circular_orbits(size_t n) {
  Particle_array_t particle_buf = new_array_aligned_t(n + 1);

  Particle star = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.00465047, 1.0, {0.0}};
  particle_buf.ptr[0] = star;

  for (size_t i = 1; i <= n; ++i) {
    double d = 0.1 + (i * 5.0 / n);
    double v = sqrt(1.0 / d);
    double theta = random_d() * 6.28;
    double x = d * cos(theta);
    double y = d * sin(theta);
    double vx = -v * sin(theta);
    double vy = v * cos(theta);
    Particle planet = {{x, y, 0.0}, {vx, vy, 0.0}, 1e-14, 1e-7, {0.0}};
    particle_buf.ptr[i] = planet;
  }
  return particle_buf;
}
