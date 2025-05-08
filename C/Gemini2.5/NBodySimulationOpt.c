#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string.h> // For memcpy, memset

// --- Constants ---
#define G              1.0
#define DT             0.001
#define NUM_STEPS      10
#define N_ORBITING     1000000 // Increase N slightly to see optimization benefits
#define N_TOTAL        (N_ORBITING + 1)
#define EPSILON_SQ     1e-8    // Softening factor squared (can be tuned)
#define THETA          0.3     // Barnes-Hut opening angle
const double THETA_SQ = THETA * THETA; // Precompute theta squared
#define M_PI           3.1415927 // Gemini 2.5 didn't create this.


#define CENTRAL_MASS   1000.0
#define ORBITING_MASS  0.01
#define MIN_RADIUS     5.0
#define MAX_RADIUS     50.0

// --- Data Structures ---
typedef struct {
    // Use doubles for potentially large coordinate ranges
    double *posX, *posY, *posZ;
    double *velX, *velY, *velZ;
    double *mass;
    int N;
} System;

typedef struct {
    double minX, maxX, minY, maxY, minZ, maxZ;
} BoundingBox;

typedef struct KDTreeNode {
    BoundingBox bbox;
    double centerOfMassX, centerOfMassY, centerOfMassZ;
    double totalMass;
    double nodeSize;   // Characteristic size (max dimension width)
    double nodeSizeSq; // Size squared for optimized check

    int particleIndex; // Index in the original System arrays, -1 if internal node
    int numParticles;  // Number of particles within or below this node

    struct KDTreeNode *left;
    struct KDTreeNode *right;
} KDTreeNode;

// --- Node Pool (Arena Allocator) ---
KDTreeNode *g_node_pool = NULL;
size_t g_node_pool_size = 0;
size_t g_node_pool_next_idx = 0;

void init_node_pool(size_t estimated_nodes) {
    // Estimate nodes needed: ~2*N for a balanced tree, add buffer
    g_node_pool_size = estimated_nodes; // Usually around 2*N - 1 nodes
    g_node_pool = (KDTreeNode*)malloc(g_node_pool_size * sizeof(KDTreeNode));
    if (!g_node_pool) {
        fprintf(stderr, "Fatal: Failed to allocate node pool!\n");
        exit(EXIT_FAILURE);
    }
    g_node_pool_next_idx = 0;
    printf("Node pool initialized with capacity for %zu nodes.\n", g_node_pool_size);
}

void free_node_pool() {
    free(g_node_pool);
    g_node_pool = NULL;
    g_node_pool_size = 0;
    g_node_pool_next_idx = 0;
}

// Reset pool usage for the next tree build
void reset_node_pool() {
    g_node_pool_next_idx = 0;
}

// Allocate a node from the pool
KDTreeNode* alloc_node_from_pool() {
    if (g_node_pool_next_idx >= g_node_pool_size) {
        // Option 1: Error out (simplest)
        fprintf(stderr, "Fatal: Node pool exhausted! (Requested %zu, Size %zu)\n",
                g_node_pool_next_idx + 1, g_node_pool_size);
        // Option 2: Try to resize (more complex)
        // size_t new_size = g_node_pool_size * 2;
        // KDTreeNode* new_pool = (KDTreeNode*)realloc(g_node_pool, new_size * sizeof(KDTreeNode));
        // if (!new_pool) { ... error ... }
        // g_node_pool = new_pool;
        // g_node_pool_size = new_size;
        // printf("Warning: Resized node pool to %zu\n", new_size);
        exit(EXIT_FAILURE); // Stick with Option 1 for now
    }
    // No need to initialize memory if relying on struct assignments later
    // KDTreeNode* node = &g_node_pool[g_node_pool_next_idx++];
    // memset(node, 0, sizeof(KDTreeNode)); // Optional: zero out memory
    // return node;
     return &g_node_pool[g_node_pool_next_idx++];
}

// Structure to hold particle index and coordinate for sorting/partitioning
typedef struct {
    int index;
    double coord; // Coordinate value along the current splitting dimension
} SortParticle;


// --- Function Prototypes ---
void allocate_system(System *sys, int N);
void free_system(System *sys);
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting);
double calculate_energy(const System *sys);

// Tree-related functions
BoundingBox compute_bounding_box(const System *sys, const int *indices, int count);
KDTreeNode* build_kdtree(const System *sys, int *indices, SortParticle *aux_sort_array, int count, int depth); // Added aux array
// removed: void free_kdtree(KDTreeNode *node); // Replaced by reset_node_pool
void compute_node_properties(KDTreeNode *node, const System *sys);
void calculate_force_from_node(int targetParticleIdx, KDTreeNode *node, const System *sys, double *fx, double *fy, double *fz); // Removed theta (uses THETA_SQ)
void calculate_forces_tree(const System *sys, KDTreeNode *root, double *forceX, double *forceY, double *forceZ); // Removed theta

// Simulation step using the tree
void simulate_step_tree(System *sys, double *forceX, double *forceY, double *forceZ, int *particle_indices, SortParticle *sort_particles); // Pass reusable arrays

// Partitioning/Median finding functions
void swap_sort_particles(SortParticle *a, SortParticle *b);
int partition_sort_particles(SortParticle *arr, int low, int high, int pivot_idx);
int find_kth_smallest(SortParticle *arr, int low, int high, int k); // k is 0-based index


// --- Main Function ---
int main() {
    printf("N-Body Simulation (OpenMP + kD-Tree Barnes-Hut - Optimized)\n");
    printf("----------------------------------------------------------\n");
    printf("Bodies: %d (%d orbiting + 1 central)\n", N_TOTAL, N_ORBITING);
    printf("Steps:  %d\n", NUM_STEPS);
    printf("DT:     %.4f\n", DT);
    printf("G:      %.2f\n", G);
    printf("Theta:  %.2f (Theta^2 = %.4f)\n", THETA, THETA_SQ);
    printf("OpenMP Max Threads: %d\n", omp_get_max_threads());
    printf("Optimizations: Node Pool, Quickselect Partitioning, Squared BH Criterion.\n");
    printf("NOTE: Tree building remains sequential. Energy check uses O(N^2) potential.\n\n");

    System sys;
    allocate_system(&sys, N_TOTAL);

    // --- Pre-allocate reusable arrays ---
    double *forceX = (double*)malloc(sys.N * sizeof(double));
    double *forceY = (double*)malloc(sys.N * sizeof(double));
    double *forceZ = (double*)malloc(sys.N * sizeof(double));
    int *particle_indices = (int*)malloc(sys.N * sizeof(int)); // For building tree
    SortParticle *sort_particles = (SortParticle*)malloc(sys.N * sizeof(SortParticle)); // Aux array for partitioning
    if (!forceX || !forceY || !forceZ || !particle_indices || !sort_particles) {
        fprintf(stderr, "Fatal: Failed to allocate reusable arrays!\n");
        exit(EXIT_FAILURE);
    }

    // --- Initialize Node Pool ---
    // Estimate: Max depth log2(N), total nodes < 2*N. Add buffer.
    size_t estimated_nodes = (size_t)(2.5 * sys.N);
    init_node_pool(estimated_nodes);

    printf("Initializing system...\n");
    initialize_system(&sys, CENTRAL_MASS, ORBITING_MASS, MIN_RADIUS, MAX_RADIUS, N_ORBITING);
    printf("Initialization complete.\n");

    printf("Calculating initial energy (using O(N^2) potential)...\n");
    double initial_energy = calculate_energy(&sys);
    printf("Initial System Energy: %.6e\n", initial_energy);

    printf("Starting simulation with kD-Tree...\n");
    double start_omp = omp_get_wtime();
    double tree_build_time = 0.0;
    double force_calc_time = 0.0;

    for (int step = 0; step < NUM_STEPS; ++step) {
        // Simulate step, passing reusable arrays
        simulate_step_tree(&sys, forceX, forceY, forceZ, particle_indices, sort_particles);

        // Optional: More detailed timing (adds overhead)
        // double step_start = omp_get_wtime();
        // reset_node_pool();
        // // Build Tree...
        // double build_end = omp_get_wtime();
        // tree_build_time += (build_end - step_start);
        // // Calc Forces...
        // double force_end = omp_get_wtime();
        // force_calc_time += (force_end - build_end);
        // // Update vel/pos...

        if ((step + 1) % (NUM_STEPS / 10) == 0 || step == NUM_STEPS - 1) {
             double current_omp = omp_get_wtime();
             printf("Step %d / %d completed. Cumulative time: %.2f s\n",
                    step + 1, NUM_STEPS, current_omp - start_omp);
        }
    }

    double end_omp = omp_get_wtime();
    double total_time = end_omp - start_omp;
    printf("Simulation finished.\n");
    // printf("Approx Tree Build Time: %.2f s\n", tree_build_time); // If timing enabled
    // printf("Approx Force Calc Time: %.2f s\n", force_calc_time); // If timing enabled
    printf("Total simulation time: %.2f seconds\n", total_time);

    printf("Calculating final energy (using O(N^2) potential)...\n");
    double final_energy = calculate_energy(&sys);
    printf("Final System Energy:   %.6e\n", final_energy);

    double energy_diff = final_energy - initial_energy;
    double energy_rel_diff = (initial_energy != 0.0) ? (energy_diff / fabs(initial_energy)) : 0.0;
    printf("Energy Difference (Final - Initial): %.6e\n", energy_diff);
    printf("Relative Energy Difference: %.6e (%.4f%%)\n", energy_rel_diff, energy_rel_diff * 100.0);
    printf("NOTE: Energy difference reflects approximation error from tree method.\n");

    // --- Free resources ---
    free_node_pool();
    free(forceX);
    free(forceY);
    free(forceZ);
    free(particle_indices);
    free(sort_particles);
    free_system(&sys);

    return 0;
}

// --- System Init/Energy/Free (Mostly Unchanged) ---
void allocate_system(System *sys, int N) { /* ... unchanged ... */
    sys->N = N;
    sys->posX = (double*)malloc(N * sizeof(double));
    sys->posY = (double*)malloc(N * sizeof(double));
    sys->posZ = (double*)malloc(N * sizeof(double));
    sys->velX = (double*)malloc(N * sizeof(double));
    sys->velY = (double*)malloc(N * sizeof(double));
    sys->velZ = (double*)malloc(N * sizeof(double));
    sys->mass = (double*)malloc(N * sizeof(double));
    // Error checking...
    if (!sys->posX || !sys->posY || !sys->posZ || !sys->velX || !sys->velY || !sys->velZ || !sys->mass) {
        fprintf(stderr, "Error: Memory allocation failed!\n"); exit(EXIT_FAILURE);
    }
}
void free_system(System *sys) { /* ... unchanged ... */
    free(sys->posX); free(sys->posY); free(sys->posZ);
    free(sys->velX); free(sys->velY); free(sys->velZ);
    free(sys->mass); sys->N = 0;
}
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting) { /* ... unchanged ... */
     if (sys->N != n_orbiting + 1) { fprintf(stderr, "Error: System size mismatch!\n"); exit(EXIT_FAILURE); }
    srand(time(NULL));
    sys->mass[0] = central_mass;
    sys->posX[0] = 0.0; sys->posY[0] = 0.0; sys->posZ[0] = 0.0;
    sys->velX[0] = 0.0; sys->velY[0] = 0.0; sys->velZ[0] = 0.0;
    #pragma omp parallel for
    for (int i = 1; i < sys->N; ++i) {
        sys->mass[i] = orbiting_mass;
        double r = min_radius + (max_radius - min_radius) * ((double)rand() / RAND_MAX);
        double angle = 2.0 * M_PI * ((double)rand() / RAND_MAX);
        sys->posX[i] = r * cos(angle); sys->posY[i] = r * sin(angle); sys->posZ[i] = 0.0;
        if (r > 1e-9) { // Avoid division by zero for speed calculation
            double orbital_speed = sqrt(G * sys->mass[0] / r);
            sys->velX[i] = -orbital_speed * sin(angle); sys->velY[i] =  orbital_speed * cos(angle); sys->velZ[i] =  0.0;
        } else { sys->velX[i] = 0.0; sys->velY[i] = 0.0; sys->velZ[i] = 0.0; }
    }
}
double calculate_energy(const System *sys) { /* ... unchanged O(N^2) version ... */
    double kinetic_energy = 0.0; double potential_energy = 0.0;
    #pragma omp parallel for reduction(+:kinetic_energy)
    for (int i = 0; i < sys->N; ++i) {
        double v_sq = sys->velX[i] * sys->velX[i] + sys->velY[i] * sys->velY[i] + sys->velZ[i] * sys->velZ[i];
        kinetic_energy += 0.5 * sys->mass[i] * v_sq;
    }
    #pragma omp parallel for reduction(+:potential_energy)
    for (int i = 0; i < sys->N; ++i) {
        for (int j = i + 1; j < sys->N; ++j) {
            double dx = sys->posX[j] - sys->posX[i]; double dy = sys->posY[j] - sys->posY[i]; double dz = sys->posZ[j] - sys->posZ[i];
            double r_sq = dx * dx + dy * dy + dz * dz; double r = sqrt(r_sq + EPSILON_SQ);
            if (r > sqrt(EPSILON_SQ)) { potential_energy -= G * sys->mass[i] * sys->mass[j] / r; }
        }
    } return kinetic_energy + potential_energy;
}

// --- Optimized kD-Tree Functions ---

// Helper: Swap two SortParticle structs
void swap_sort_particles(SortParticle *a, SortParticle *b) {
    SortParticle temp = *a;
    *a = *b;
    *b = temp;
}

// Helper: Partition subarray around a pivot element (Lomuto scheme variation)
// Returns the final index of the pivot element after partitioning.
int partition_sort_particles(SortParticle *arr, int low, int high, int pivot_idx) {
    double pivot_value = arr[pivot_idx].coord;
    // Move pivot to end
    swap_sort_particles(&arr[pivot_idx], &arr[high]);
    int store_idx = low;
    for (int i = low; i < high; ++i) {
        if (arr[i].coord < pivot_value) {
            swap_sort_particles(&arr[store_idx], &arr[i]);
            store_idx++;
        }
    }
    // Move pivot to its final place
    swap_sort_particles(&arr[store_idx], &arr[high]);
    return store_idx;
}

// Helper: Find the k-th smallest element (0-based index k) using Quickselect
// Modifies the array `arr` in place between low and high.
int find_kth_smallest(SortParticle *arr, int low, int high, int k) {
    if (low == high) {
        return low;
    }

    // Select a pivot index (e.g., median of 3 or random) - simple middle for now
    int pivot_idx = low + (high - low) / 2;

    pivot_idx = partition_sort_particles(arr, low, high, pivot_idx);

    // Check if pivot is the k-th element
    if (k == pivot_idx) {
        return k;
    } else if (k < pivot_idx) {
        return find_kth_smallest(arr, low, pivot_idx - 1, k);
    } else {
        return find_kth_smallest(arr, pivot_idx + 1, high, k);
    }
}

// Compute bounding box (unchanged logic, maybe optimize later if needed)
BoundingBox compute_bounding_box(const System *sys, const int *indices, int count) { /* ... unchanged ... */
    BoundingBox bbox = { INFINITY, -INFINITY, INFINITY, -INFINITY, INFINITY, -INFINITY };
    if (count == 0) return (BoundingBox){ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    for (int i = 0; i < count; ++i) {
        int p_idx = indices[i];
        const double px = sys->posX[p_idx], py = sys->posY[p_idx], pz = sys->posZ[p_idx];
        if (px < bbox.minX) bbox.minX = px; if (px > bbox.maxX) bbox.maxX = px;
        if (py < bbox.minY) bbox.minY = py; if (py > bbox.maxY) bbox.maxY = py;
        if (pz < bbox.minZ) bbox.minZ = pz; if (pz > bbox.maxZ) bbox.maxZ = pz;
    } return bbox;
}


// Compute node properties (added nodeSizeSq)
void compute_node_properties(KDTreeNode *node, const System *sys) {
    if (node == NULL || node->numParticles == 0) {
        // Initialize fields for empty node
        node->totalMass = 0.0; node->centerOfMassX = 0.0; node->centerOfMassY = 0.0; node->centerOfMassZ = 0.0;
        node->nodeSize = 0.0; node->nodeSizeSq = 0.0;
        return;
    }

    if (node->particleIndex != -1) { // Leaf node
        int p_idx = node->particleIndex;
        node->totalMass = sys->mass[p_idx];
        node->centerOfMassX = sys->posX[p_idx];
        node->centerOfMassY = sys->posY[p_idx];
        node->centerOfMassZ = sys->posZ[p_idx];
    } else { // Internal node
        double m_left = (node->left) ? node->left->totalMass : 0.0;
        double m_right = (node->right) ? node->right->totalMass : 0.0;
        node->totalMass = m_left + m_right;

        if (node->totalMass > 1e-20) { // Use a small tolerance for mass
            double comX_left = (node->left) ? node->left->centerOfMassX : 0.0;
            double comY_left = (node->left) ? node->left->centerOfMassY : 0.0;
            double comZ_left = (node->left) ? node->left->centerOfMassZ : 0.0;
            double comX_right = (node->right) ? node->right->centerOfMassX : 0.0;
            double comY_right = (node->right) ? node->right->centerOfMassY : 0.0;
            double comZ_right = (node->right) ? node->right->centerOfMassZ : 0.0;

            node->centerOfMassX = (m_left * comX_left + m_right * comX_right) / node->totalMass;
            node->centerOfMassY = (m_left * comY_left + m_right * comY_right) / node->totalMass;
            node->centerOfMassZ = (m_left * comZ_left + m_right * comZ_right) / node->totalMass;
        } else {
             // Center of mass is ill-defined for zero mass; use geometric center
             node->centerOfMassX = (node->bbox.minX + node->bbox.maxX) * 0.5;
             node->centerOfMassY = (node->bbox.minY + node->bbox.maxY) * 0.5;
             node->centerOfMassZ = (node->bbox.minZ + node->bbox.maxZ) * 0.5;
             node->totalMass = 0.0; // Ensure mass is exactly zero
        }
    }

    // Calculate node size and size squared
    double sizeX = node->bbox.maxX - node->bbox.minX;
    double sizeY = node->bbox.maxY - node->bbox.minY;
    double sizeZ = node->bbox.maxZ - node->bbox.minZ;
    node->nodeSize = fmax(sizeX, fmax(sizeY, sizeZ));
    node->nodeSizeSq = node->nodeSize * node->nodeSize;
}

// Build kD-Tree using Quickselect partitioning and node pool
KDTreeNode* build_kdtree(const System *sys, int *indices, SortParticle *aux_sort_array, int count, int depth) {
    KDTreeNode *node = alloc_node_from_pool(); // Get node from pool

    node->left = node->right = NULL;
    node->particleIndex = -1;
    node->numParticles = count;
    node->bbox = compute_bounding_box(sys, indices, count);

    if (count == 0) { compute_node_properties(node, sys); return node; }
    if (count == 1) { node->particleIndex = indices[0]; compute_node_properties(node, sys); return node; }

    int split_dim = depth % 3;

    // Populate auxiliary array for partitioning this subset
    for(int i = 0; i < count; ++i) {
        aux_sort_array[i].index = indices[i];
        switch(split_dim) {
            case 0: aux_sort_array[i].coord = sys->posX[indices[i]]; break;
            case 1: aux_sort_array[i].coord = sys->posY[indices[i]]; break;
            case 2: aux_sort_array[i].coord = sys->posZ[indices[i]]; break;
        }
    }

    // Find the median element's final index using Quickselect (k is 0-based median index)
    int k_median = count / 2; // Target rank for median
    int median_final_idx = find_kth_smallest(aux_sort_array, 0, count - 1, k_median);
    // The element at aux_sort_array[median_final_idx] is the true median.
    // The array aux_sort_array is now partitioned around this median value,
    // BUT the element itself might not be exactly at index k_median if duplicates exist.
    // For tree building, we need to split near the middle index 'k_median'.
    // A simple approach is to partition around the element found at index k_median after find_kth_smallest.
    median_final_idx = partition_sort_particles(aux_sort_array, 0, count - 1, k_median);


    // Update the original indices array based on the partitioned auxiliary array
    for(int i = 0; i < count; ++i) {
        indices[i] = aux_sort_array[i].index;
    }

    // Determine split counts (handle edge case where partition might be skewed)
    int left_count = median_final_idx;
    int right_count = count - (median_final_idx + 1);

     // Recurse, passing down the corresponding slice of the indices array AND the aux array
    node->left = build_kdtree(sys, indices, aux_sort_array, left_count, depth + 1);
    node->right = build_kdtree(sys, indices + median_final_idx + 1, aux_sort_array + median_final_idx + 1, right_count, depth + 1);


    compute_node_properties(node, sys); // Compute properties after children are built
    return node;
}


// Calculate force from node (Optimized criterion)
void calculate_force_from_node(int targetParticleIdx, KDTreeNode *node, const System *sys, double *fx, double *fy, double *fz) {
    if (node == NULL || node->numParticles == 0 || node->particleIndex == targetParticleIdx) {
        return; // Skip empty nodes, self-interaction
    }

    // Vector from target particle to node's COM (note direction change for force)
    double dx = node->centerOfMassX - sys->posX[targetParticleIdx];
    double dy = node->centerOfMassY - sys->posY[targetParticleIdx];
    double dz = node->centerOfMassZ - sys->posZ[targetParticleIdx];

    double distSq = dx * dx + dy * dy + dz * dz;

    // Check Barnes-Hut criterion using squared values: nodeSize^2 < dist^2 * theta^2
    // Also check if the node is a leaf. Need dist > 0 check.
    int is_leaf = (node->particleIndex != -1);
    // Consider node "far" if distSq is reasonably large AND criterion met
    int treat_as_single = is_leaf || (distSq * THETA_SQ > node->nodeSizeSq);

    if (treat_as_single) {
         // Calculate force using softened distance squared
        double distSq_soft = distSq + EPSILON_SQ;
        // Check if distSq_soft is significantly greater than zero before division
         if (distSq_soft > 1e-20 && node->totalMass > 1e-20) {
            double dist_soft = sqrt(distSq_soft);
            double force_mag_over_dist = G * node->totalMass * sys->mass[targetParticleIdx] / distSq_soft / dist_soft;

            // Force pulls target particle towards node COM
            *fx += force_mag_over_dist * dx;
            *fy += force_mag_over_dist * dy;
            *fz += force_mag_over_dist * dz;
         }
    } else { // Internal node too close: recurse
        calculate_force_from_node(targetParticleIdx, node->left, sys, fx, fy, fz);
        calculate_force_from_node(targetParticleIdx, node->right, sys, fx, fy, fz);
    }
}

// Calculate forces for all particles using the kD-Tree (Parallelized)
void calculate_forces_tree(const System *sys, KDTreeNode *root, double *forceX, double *forceY, double *forceZ) {
    #pragma omp parallel for schedule(dynamic) // Dynamic schedule often good for tree traversals
    for (int i = 0; i < sys->N; ++i) {
        double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;
        // Skip force calculation if mass is zero or negligible (optional optimization)
        if (sys->mass[i] > 1e-20) {
             calculate_force_from_node(i, root, sys, &fx_i, &fy_i, &fz_i);
        }
        forceX[i] = fx_i;
        forceY[i] = fy_i;
        forceZ[i] = fz_i;
    }
}

// --- Simulation Step using Optimized kD-Tree ---
void simulate_step_tree(System *sys, double *forceX, double *forceY, double *forceZ, int *particle_indices, SortParticle *sort_particles) {

    // --- 1. Build kD-Tree (Sequential, using node pool & Quickselect) ---
    reset_node_pool(); // Reset arena allocator for this step
    for (int i = 0; i < sys->N; ++i) { particle_indices[i] = i; } // Initialize indices

    // Pass sort_particles as auxiliary space to build_kdtree
    KDTreeNode *tree_root = build_kdtree(sys, particle_indices, sort_particles, sys->N, 0);

    // --- 2. Calculate Forces using Tree (Parallel) ---
    calculate_forces_tree(sys, tree_root, forceX, forceY, forceZ);

    // Tree structure (in node pool) is implicitly discarded on next reset_node_pool()

    // --- 3. Update Velocities (Parallel) ---
    #pragma omp parallel for
    for (int i = 0; i < sys->N; ++i) {
        // Precompute inverse mass if mass is non-zero
        if (sys->mass[i] > 1e-20) {
            double inv_mass = 1.0 / sys->mass[i];
            sys->velX[i] += forceX[i] * inv_mass * DT;
            sys->velY[i] += forceY[i] * inv_mass * DT;
            sys->velZ[i] += forceZ[i] * inv_mass * DT;
        }
         // Optional: Zero out velocity if mass is zero? Or handle earlier.
    }

    // --- 4. Update Positions (Parallel) ---
    #pragma omp parallel for
    for (int i = 0; i < sys->N; ++i) {
        sys->posX[i] += sys->velX[i] * DT;
        sys->posY[i] += sys->velY[i] * DT;
        sys->posZ[i] += sys->velZ[i] * DT;
    }
}