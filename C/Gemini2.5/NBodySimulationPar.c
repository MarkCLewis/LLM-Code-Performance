#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // For srand and timing
#include <omp.h>  // <<< Include OpenMP header

// --- Constants ---
// (Keep constants the same as the sequential version)
#define G              1.0
#define DT             0.001
#define NUM_STEPS      1000
#define N_ORBITING     1000000 // <<< WARNING: Still computationally intensive, but faster with threads
#define N_TOTAL        (N_ORBITING + 1)
#define EPSILON_SQ     1e-6
#define M_PI           3.1415927 // Gemini 2.5 didn't create this.

#define CENTRAL_MASS   1000.0
#define ORBITING_MASS  0.01
#define MIN_RADIUS     5.0
#define MAX_RADIUS     50.0

// --- Data Structure ---
typedef struct {
    double *posX, *posY, *posZ;
    double *velX, *velY, *velZ;
    double *mass;
    int N;
} System;

// --- Function Prototypes ---
void allocate_system(System *sys, int N);
void free_system(System *sys);
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting);
double calculate_energy(const System *sys);
void calculate_forces(const System *sys, double *forceX, double *forceY, double *forceZ);
void simulate_step(System *sys);

// --- Main Function ---
int main() {
    printf("N-Body Simulation (OpenMP Multithreaded)\n");
    printf("----------------------------------------\n");
    printf("Bodies: %d (%d orbiting + 1 central)\n", N_TOTAL, N_ORBITING);
    printf("Steps:  %d\n", NUM_STEPS);
    printf("DT:     %.4f\n", DT);
    printf("G:      %.2f\n", G);
    // Print OpenMP info
    printf("OpenMP Max Threads: %d\n", omp_get_max_threads());
    printf("WARNING: N=%d is still large. Execution time depends heavily on core count.\n\n", N_TOTAL);


    System sys;
    allocate_system(&sys, N_TOTAL);

    printf("Initializing system...\n");
    initialize_system(&sys, CENTRAL_MASS, ORBITING_MASS, MIN_RADIUS, MAX_RADIUS, N_ORBITING);
    printf("Initialization complete.\n");

    printf("Calculating initial energy...\n");
    double initial_energy = calculate_energy(&sys);
    printf("Initial System Energy: %.6e\n", initial_energy);

    printf("Starting simulation...\n");
    // Use OpenMP timer for more precision if needed: double start_omp = omp_get_wtime();
    time_t start_time = time(NULL);

    for (int step = 0; step < NUM_STEPS; ++step) {
        simulate_step(&sys);

        if ((step + 1) % (NUM_STEPS / 10) == 0 || step == NUM_STEPS - 1) {
             time_t current_time = time(NULL);
             // double elapsed_omp = omp_get_wtime() - start_omp;
             double elapsed_sec = difftime(current_time, start_time);
             printf("Step %d / %d completed. Elapsed time: %.2f s\n",
                    step + 1, NUM_STEPS, elapsed_sec);
        }
    }

    time_t end_time = time(NULL);
    // double end_omp = omp_get_wtime();
    double total_time = difftime(end_time, start_time); // end_omp - start_omp;
    printf("Simulation finished.\n");
    printf("Total simulation time: %.2f seconds\n", total_time);


    printf("Calculating final energy...\n");
    double final_energy = calculate_energy(&sys);
    printf("Final System Energy:   %.6e\n", final_energy);

    double energy_diff = final_energy - initial_energy;
    double energy_rel_diff = (initial_energy != 0.0) ? (energy_diff / fabs(initial_energy)) : 0.0; // Use fabs for denominator
    printf("Energy Difference (Final - Initial): %.6e\n", energy_diff);
    printf("Relative Energy Difference: %.6e (%.4f%%)\n", energy_rel_diff, energy_rel_diff * 100.0);

    free_system(&sys);

    return 0;
}

// --- Function Implementations ---

// Allocate memory (no changes needed for OpenMP)
void allocate_system(System *sys, int N) {
    sys->N = N;
    sys->posX = (double*)malloc(N * sizeof(double));
    sys->posY = (double*)malloc(N * sizeof(double));
    sys->posZ = (double*)malloc(N * sizeof(double));
    sys->velX = (double*)malloc(N * sizeof(double));
    sys->velY = (double*)malloc(N * sizeof(double));
    sys->velZ = (double*)malloc(N * sizeof(double));
    sys->mass = (double*)malloc(N * sizeof(double));

    if (!sys->posX || !sys->posY || !sys->posZ || !sys->velX || !sys->velY || !sys->velZ || !sys->mass) {
        fprintf(stderr, "Error: Memory allocation failed!\n");
        free(sys->posX); free(sys->posY); free(sys->posZ);
        free(sys->velX); free(sys->velY); free(sys->velZ);
        free(sys->mass);
        exit(EXIT_FAILURE);
    }
}

// Free memory (no changes needed for OpenMP)
void free_system(System *sys) {
    free(sys->posX);
    free(sys->posY);
    free(sys->posZ);
    free(sys->velX);
    free(sys->velY);
    free(sys->velZ);
    free(sys->mass);
    sys->N = 0;
}

// Initialize the system (Parallelized loop)
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting) {
    if (sys->N != n_orbiting + 1) {
        fprintf(stderr, "Error: System size mismatch during initialization!\n");
        exit(EXIT_FAILURE);
    }

    // Seed random number generator (once, before parallel region)
    // Note: Using rand() inside parallel region is generally not thread-safe.
    // Consider thread-safe alternatives for production code.
    srand(time(NULL));

    // --- Central Body (index 0) - done sequentially ---
    sys->mass[0] = central_mass;
    sys->posX[0] = 0.0;
    sys->posY[0] = 0.0;
    sys->posZ[0] = 0.0;
    sys->velX[0] = 0.0;
    sys->velY[0] = 0.0;
    sys->velZ[0] = 0.0;

    // --- Orbiting Bodies (indices 1 to N-1) - Parallelized ---
    #pragma omp parallel for // <<< Parallelize this loop
    for (int i = 1; i < sys->N; ++i) {
        // WARNING: rand() called within parallel region - potential issues.
        // For demo purposes only. Use thread-safe RNG in real applications.
        // Local seeding (less ideal statistically): unsigned int seed = time(NULL) ^ omp_get_thread_num();
        // double rand_val = rand_r(&seed) / (double)RAND_MAX;

        sys->mass[i] = orbiting_mass;

        double r = min_radius + (max_radius - min_radius) * ((double)rand() / RAND_MAX);
        double angle = 2.0 * M_PI * ((double)rand() / RAND_MAX);

        sys->posX[i] = r * cos(angle);
        sys->posY[i] = r * sin(angle);
        sys->posZ[i] = 0.0;

        if (r > 0) {
            double orbital_speed = sqrt(G * sys->mass[0] / r); // Use mass[0] directly
            sys->velX[i] = -orbital_speed * sin(angle);
            sys->velY[i] =  orbital_speed * cos(angle);
            sys->velZ[i] =  0.0;
        } else {
            sys->velX[i] = 0.0;
            sys->velY[i] = 0.0;
            sys->velZ[i] = 0.0;
        }
    }
}

// Calculate total energy (Parallelized loops with reduction)
double calculate_energy(const System *sys) {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;

    // Kinetic Energy (Parallelized)
    #pragma omp parallel for reduction(+:kinetic_energy) // <<< Parallelize KE loop
    for (int i = 0; i < sys->N; ++i) {
        double v_sq = sys->velX[i] * sys->velX[i] +
                      sys->velY[i] * sys->velY[i] +
                      sys->velZ[i] * sys->velZ[i];
        kinetic_energy += 0.5 * sys->mass[i] * v_sq;
    }

    // Potential Energy (Parallelized outer loop)
    #pragma omp parallel for reduction(+:potential_energy) // <<< Parallelize PE loop (outer)
    for (int i = 0; i < sys->N; ++i) {
        for (int j = i + 1; j < sys->N; ++j) { // Inner loop remains sequential for each 'i'
            double dx = sys->posX[j] - sys->posX[i];
            double dy = sys->posY[j] - sys->posY[i];
            double dz = sys->posZ[j] - sys->posZ[i];

            double r_sq = dx * dx + dy * dy + dz * dz;
            double r = sqrt(r_sq + EPSILON_SQ);

            if (r > sqrt(EPSILON_SQ)) {
                 potential_energy -= G * sys->mass[i] * sys->mass[j] / r;
            }
        }
    }

    return kinetic_energy + potential_energy;
}

// Calculate forces (Parallelized O(N^2) calculation)
void calculate_forces(const System *sys, double *forceX, double *forceY, double *forceZ) {
    // Initialize forces to zero (can be done sequentially or in parallel)
    // Doing it inside the parallel loop implicitly initializes the sums per thread.
    // Let's explicitly zero the arrays first (parallelized).
    #pragma omp parallel for
    for(int i=0; i<sys->N; ++i) {
        forceX[i] = 0.0;
        forceY[i] = 0.0;
        forceZ[i] = 0.0;
    }


    // Parallelize the loop calculating force ON particle 'i'
    #pragma omp parallel for schedule(dynamic) // <<< Parallelize outer loop, dynamic schedule might be better for uneven load
    for (int i = 0; i < sys->N; ++i) {
        double fx_i = 0.0; // Use local accumulators for particle i's force components
        double fy_i = 0.0;
        double fz_i = 0.0;

        // Inner loop: iterate through all OTHER particles 'j' to calculate force ON 'i'
        for (int j = 0; j < sys->N; ++j) {
            if (i == j) continue; // Skip self-interaction

            // Calculate vector from i to j: r_ij = pos_j - pos_i
            double dx = sys->posX[j] - sys->posX[i];
            double dy = sys->posY[j] - sys->posY[i];
            double dz = sys->posZ[j] - sys->posZ[i];

            double r_sq = dx * dx + dy * dy + dz * dz;
            double dist = sqrt(r_sq + EPSILON_SQ); // Apply softening

            if (dist > 1e-10) { // Avoid division by ~zero
                // Force magnitude F = G * m_i * m_j / (dist^2 + epsilon^2)
                // Force vector direction is along r_ij (from i to j)
                // F_vector = F * (r_ij / dist)
                double force_mag_over_dist = G * sys->mass[i] * sys->mass[j] / (r_sq + EPSILON_SQ) / dist;

                // Add force component exerted BY j ON i to the accumulators for i
                fx_i += force_mag_over_dist * dx;
                fy_i += force_mag_over_dist * dy;
                fz_i += force_mag_over_dist * dz;
            }
        }
        // Write the total calculated force for particle 'i' to the shared array.
        // This is safe because each thread writes to a unique force[i].
        forceX[i] = fx_i;
        forceY[i] = fy_i;
        forceZ[i] = fz_i;
    }
    // Implicit barrier at the end of the parallel for loop ensures all forces are computed
}

// Perform one simulation step (Parallelized O(N) updates)
void simulate_step(System *sys) {
    // Allocate force arrays (outside parallel regions if possible, but here is fine too)
    // These arrays will be written to by calculate_forces
    double *forceX = (double*)malloc(sys->N * sizeof(double));
    double *forceY = (double*)malloc(sys->N * sizeof(double));
    double *forceZ = (double*)malloc(sys->N * sizeof(double));
    if (!forceX || !forceY || !forceZ) {
        fprintf(stderr, "Error: Force array allocation failed in simulate_step!\n");
        exit(EXIT_FAILURE);
    }

    // 1. Calculate forces (Parallelized inside the function)
    calculate_forces(sys, forceX, forceY, forceZ);

    // 2. Update velocities ("Kick") - Parallelized
    #pragma omp parallel for // <<< Parallelize velocity update
    for (int i = 0; i < sys->N; ++i) {
        if (sys->mass[i] > 0) {
            double inv_mass = 1.0 / sys->mass[i]; // Compute inverse mass once
            sys->velX[i] += forceX[i] * inv_mass * DT;
            sys->velY[i] += forceY[i] * inv_mass * DT;
            sys->velZ[i] += forceZ[i] * inv_mass * DT;
        }
    }

    // 3. Update positions ("Drift") - Parallelized
    #pragma omp parallel for // <<< Parallelize position update
    for (int i = 0; i < sys->N; ++i) {
        sys->posX[i] += sys->velX[i] * DT;
        sys->posY[i] += sys->velY[i] * DT;
        sys->posZ[i] += sys->velZ[i] * DT;
    }

    // Free temporary force arrays
    free(forceX);
    free(forceY);
    free(forceZ);
}