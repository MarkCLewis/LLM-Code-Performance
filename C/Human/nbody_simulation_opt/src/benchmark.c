#define _POSIX_C_SOURCE 199309L  // Required for clock_gettime and CLOCK_MONOTONIC

#include "kdtree.h"
#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <string.h>

// Measure elapsed time using high-resolution clock
double get_time() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Simple run function to benchmark different particle counts
void run_benchmark(int n_particles, int steps) {
  printf("Benchmark: %d particles, %d steps\n", n_particles, steps);
  
  double dt = 1e-3;
  double start_time;
  double elapsed;
  
  // Generate particles
  Particle_array_t particles = circular_orbits(n_particles);
  
  // Run simulation and measure time
  start_time = get_time();
  simple_sim(&particles, dt, steps);
  elapsed = get_time() - start_time;
  
  printf("Completed in %.6f seconds\n", elapsed);
  printf("Performance: %.2f particle-steps per second\n", 
         (double)n_particles * steps / elapsed);
  printf("-------------------------\n");
  
  // Clean up
  FREE_ARRAY(particles);
}

int main(int argc, char *argv[]) {
  // Initialize random seed
  srand(42);
  
  // Run benchmarks with increasing particle counts
  printf("Running benchmarks with %d threads\n", omp_get_max_threads());
  printf("=========================\n");
  
  // Small test
  run_benchmark(1000, 10);
  
  // Medium test
  run_benchmark(5000, 10);
  
  // Larger test if command line argument is provided
  if (argc > 1 && strcmp(argv[1], "large") == 0) {
    run_benchmark(20000, 5);
  }
  
  return 0;
} 