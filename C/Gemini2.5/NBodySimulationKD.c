#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string.h> // For memcpy

// --- Constants ---
#define G              1.0
#define DT             0.001
#define NUM_STEPS      10
#define N_ORBITING     1000000 // Reduced N for faster testing with tree overhead
#define N_TOTAL        (N_ORBITING + 1)
#define EPSILON_SQ     1e-6
#define THETA          0.3     // Barnes-Hut opening angle <<< NEW
#define M_PI           3.1415927 // Gemini 2.5 didn't create this.


#define CENTRAL_MASS   1000.0
#define ORBITING_MASS  0.01
#define MIN_RADIUS     5.0
#define MAX_RADIUS     50.0

// --- Data Structures ---
typedef struct {
    double *posX, *posY, *posZ;
    double *velX, *velY, *velZ;
    double *mass;
    int N;
} System;

// Structure for bounding box
typedef struct {
    double minX, maxX, minY, maxY, minZ, maxZ;
} BoundingBox;

// Structure for kD-Tree Node
typedef struct KDTreeNode {
    BoundingBox bbox;
    double centerOfMassX, centerOfMassY, centerOfMassZ;
    double totalMass;
    double nodeSize; // Characteristic size (e.g., max dimension width) for Barnes-Hut

    int particleIndex; // Index in the original System arrays, -1 if internal node
    int numParticles;  // Number of particles within or below this node

    struct KDTreeNode *left;
    struct KDTreeNode *right;
} KDTreeNode;

// Structure to hold particle index and coordinate for sorting
typedef struct {
    int index;
    double coord;
} SortParticle;


// --- Function Prototypes ---
void allocate_system(System *sys, int N);
void free_system(System *sys);
void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting);
double calculate_energy(const System *sys); // Still O(N^2) potential for check

// Tree-related functions
BoundingBox compute_bounding_box(const System *sys, const int *indices, int count);
KDTreeNode* build_kdtree(const System *sys, int *indices, int count, int depth);
void free_kdtree(KDTreeNode *node);
void compute_node_properties(KDTreeNode *node, const System *sys);
void calculate_force_from_node(int targetParticleIdx, KDTreeNode *node, const System *sys, double theta, double *fx, double *fy, double *fz);
void calculate_forces_tree(const System *sys, KDTreeNode *root, double theta, double *forceX, double *forceY, double *forceZ);

// Simulation step using the tree
void simulate_step_tree(System *sys);

// qsort comparison functions
int compare_sort_particles(const void *a, const void *b);

// --- Global Comparison Dimension for qsort ---
// Hacky way to pass dimension to comparison function for qsort
// A better way involves passing context, but this is simpler for demo
const System *g_sys_ptr = NULL;
int g_sort_dim = 0;


// --- Main Function ---
int main() {
    printf("N-Body Simulation (OpenMP + kD-Tree Barnes-Hut)\n");
    printf("------------------------------------------------\n");
    printf("Bodies: %d (%d orbiting + 1 central)\n", N_TOTAL, N_ORBITING);
    printf("Steps:  %d\n", NUM_STEPS);
    printf("DT:     %.4f\n", DT);
    printf("G:      %.2f\n", G);
    printf("Theta:  %.2f\n", THETA); // <<< Print theta
    printf("OpenMP Max Threads: %d\n", omp_get_max_threads());
    printf("NOTE: Tree building is sequential. Energy check uses O(N^2) potential.\n\n");


    System sys;
    allocate_system(&sys, N_TOTAL);
    g_sys_ptr = &sys; // Set global pointer for qsort comparison

    printf("Initializing system...\n");
    initialize_system(&sys, CENTRAL_MASS, ORBITING_MASS, MIN_RADIUS, MAX_RADIUS, N_ORBITING);
    printf("Initialization complete.\n");

    printf("Calculating initial energy (using O(N^2) potential)...\n");
    double initial_energy = calculate_energy(&sys);
    printf("Initial System Energy: %.6e\n", initial_energy);

    printf("Starting simulation with kD-Tree...\n");
    time_t start_time = time(NULL);
    double start_omp = omp_get_wtime(); // Use OpenMP timer

    for (int step = 0; step < NUM_STEPS; ++step) {
        simulate_step_tree(&sys); // <<< Use tree-based step

        if ((step + 1) % (NUM_STEPS / 10) == 0 || step == NUM_STEPS - 1) {
             double current_omp = omp_get_wtime();
             double elapsed_omp = current_omp - start_omp;
             printf("Step %d / %d completed. Elapsed time: %.2f s\n",
                    step + 1, NUM_STEPS, elapsed_omp);
        }
    }

    double end_omp = omp_get_wtime();
    double total_time = end_omp - start_omp;
    printf("Simulation finished.\n");
    printf("Total simulation time: %.2f seconds\n", total_time);

    printf("Calculating final energy (using O(N^2) potential)...\n");
    double final_energy = calculate_energy(&sys);
    printf("Final System Energy:   %.6e\n", final_energy);

    double energy_diff = final_energy - initial_energy;
    double energy_rel_diff = (initial_energy != 0.0) ? (energy_diff / fabs(initial_energy)) : 0.0;
    printf("Energy Difference (Final - Initial): %.6e\n", energy_diff);
    printf("Relative Energy Difference: %.6e (%.4f%%)\n", energy_rel_diff, energy_rel_diff * 100.0);
    printf("NOTE: Energy difference reflects approximation error from tree method.\n");


    free_system(&sys);
    g_sys_ptr = NULL;

    return 0;
}

// --- Basic System Functions (Allocate, Free, Init) ---
// (allocate_system, free_system, initialize_system remain largely unchanged
// from the OpenMP version, including the rand() warning in initialize_system)
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
        // Proper cleanup of partially allocated memory needed here
        exit(EXIT_FAILURE);
    }
}

void free_system(System *sys) {
    free(sys->posX); free(sys->posY); free(sys->posZ);
    free(sys->velX); free(sys->velY); free(sys->velZ);
    free(sys->mass);
    sys->N = 0;
}

void initialize_system(System *sys, double central_mass, double orbiting_mass,
                       double min_radius, double max_radius, int n_orbiting) {
     if (sys->N != n_orbiting + 1) {
        fprintf(stderr, "Error: System size mismatch during initialization!\n");
        exit(EXIT_FAILURE);
    }
    srand(time(NULL));

    sys->mass[0] = central_mass;
    sys->posX[0] = 0.0; sys->posY[0] = 0.0; sys->posZ[0] = 0.0;
    sys->velX[0] = 0.0; sys->velY[0] = 0.0; sys->velZ[0] = 0.0;

    #pragma omp parallel for
    for (int i = 1; i < sys->N; ++i) {
        // WARNING: rand() called within parallel region - potential issues.
        sys->mass[i] = orbiting_mass;
        double r = min_radius + (max_radius - min_radius) * ((double)rand() / RAND_MAX);
        double angle = 2.0 * M_PI * ((double)rand() / RAND_MAX);
        sys->posX[i] = r * cos(angle);
        sys->posY[i] = r * sin(angle);
        sys->posZ[i] = 0.0;
        if (r > 0) {
            double orbital_speed = sqrt(G * sys->mass[0] / r);
            sys->velX[i] = -orbital_speed * sin(angle);
            sys->velY[i] =  orbital_speed * cos(angle);
            sys->velZ[i] =  0.0;
        } else {
            sys->velX[i] = 0.0; sys->velY[i] = 0.0; sys->velZ[i] = 0.0;
        }
    }
}

// --- O(N^2) Energy Calculation (Unchanged) ---
double calculate_energy(const System *sys) {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;

    #pragma omp parallel for reduction(+:kinetic_energy)
    for (int i = 0; i < sys->N; ++i) {
        double v_sq = sys->velX[i] * sys->velX[i] +
                      sys->velY[i] * sys->velY[i] +
                      sys->velZ[i] * sys->velZ[i];
        kinetic_energy += 0.5 * sys->mass[i] * v_sq;
    }

    #pragma omp parallel for reduction(+:potential_energy)
    for (int i = 0; i < sys->N; ++i) {
        for (int j = i + 1; j < sys->N; ++j) {
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

// --- kD-Tree Functions ---

// Comparison function for qsort based on global dimension g_sort_dim
int compare_sort_particles(const void *a, const void *b) {
    int idxA = ((SortParticle*)a)->index;
    int idxB = ((SortParticle*)b)->index;
    double coordA, coordB;

    switch (g_sort_dim) {
        case 0: // X dimension
            coordA = g_sys_ptr->posX[idxA];
            coordB = g_sys_ptr->posX[idxB];
            break;
        case 1: // Y dimension
            coordA = g_sys_ptr->posY[idxA];
            coordB = g_sys_ptr->posY[idxB];
            break;
        case 2: // Z dimension
        default:
            coordA = g_sys_ptr->posZ[idxA];
            coordB = g_sys_ptr->posZ[idxB];
            break;
    }

    if (coordA < coordB) return -1;
    if (coordA > coordB) return 1;
    return 0;
}


// Compute bounding box for a subset of particles
BoundingBox compute_bounding_box(const System *sys, const int *indices, int count) {
    BoundingBox bbox = { INFINITY, -INFINITY, INFINITY, -INFINITY, INFINITY, -INFINITY };
    if (count == 0) {
         bbox = (BoundingBox){ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // Or handle appropriately
         return bbox;
    }

    for (int i = 0; i < count; ++i) {
        int p_idx = indices[i];
        if (sys->posX[p_idx] < bbox.minX) bbox.minX = sys->posX[p_idx];
        if (sys->posX[p_idx] > bbox.maxX) bbox.maxX = sys->posX[p_idx];
        if (sys->posY[p_idx] < bbox.minY) bbox.minY = sys->posY[p_idx];
        if (sys->posY[p_idx] > bbox.maxY) bbox.maxY = sys->posY[p_idx];
        if (sys->posZ[p_idx] < bbox.minZ) bbox.minZ = sys->posZ[p_idx];
        if (sys->posZ[p_idx] > bbox.maxZ) bbox.maxZ = sys->posZ[p_idx];
    }
    return bbox;
}

// Compute COM, total mass, and size for a node (after children are processed)
void compute_node_properties(KDTreeNode *node, const System *sys) {
    if (node == NULL || node->numParticles == 0) {
        node->totalMass = 0.0;
        node->centerOfMassX = 0.0; node->centerOfMassY = 0.0; node->centerOfMassZ = 0.0;
        node->nodeSize = 0.0;
        return;
    }

    if (node->particleIndex != -1) { // Leaf node with one particle
        int p_idx = node->particleIndex;
        node->totalMass = sys->mass[p_idx];
        node->centerOfMassX = sys->posX[p_idx];
        node->centerOfMassY = sys->posY[p_idx];
        node->centerOfMassZ = sys->posZ[p_idx];
    } else { // Internal node
        double m_left = (node->left) ? node->left->totalMass : 0.0;
        double m_right = (node->right) ? node->right->totalMass : 0.0;
        node->totalMass = m_left + m_right;

        if (node->totalMass > 1e-15) { // Avoid division by zero if mass is tiny
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
             // Handle zero mass case - COM can be midpoint of bbox or (0,0,0)
             node->centerOfMassX = (node->bbox.minX + node->bbox.maxX) / 2.0;
             node->centerOfMassY = (node->bbox.minY + node->bbox.maxY) / 2.0;
             node->centerOfMassZ = (node->bbox.minZ + node->bbox.maxZ) / 2.0;
        }
    }

    // Calculate node size (e.g., max dimension width)
    double sizeX = node->bbox.maxX - node->bbox.minX;
    double sizeY = node->bbox.maxY - node->bbox.minY;
    double sizeZ = node->bbox.maxZ - node->bbox.minZ;
    node->nodeSize = fmax(sizeX, fmax(sizeY, sizeZ));
}


// Recursively build the kD-Tree
KDTreeNode* build_kdtree(const System *sys, int *indices, int count, int depth) {
    KDTreeNode *node = (KDTreeNode*)malloc(sizeof(KDTreeNode));
    if (!node) {
        fprintf(stderr, "Error: Failed to allocate KDTreeNode!\n");
        exit(EXIT_FAILURE);
    }

    node->left = node->right = NULL;
    node->particleIndex = -1; // Assume internal node initially
    node->numParticles = count;
    node->bbox = compute_bounding_box(sys, indices, count);

    // Base cases
    if (count == 0) {
        compute_node_properties(node, sys); // Set mass/COM/size to zero
        return node;
    }
    if (count == 1) {
        node->particleIndex = indices[0];
        compute_node_properties(node, sys); // Set props based on single particle
        return node;
    }

    // Recursive step: Select dimension and partition
    int split_dim = depth % 3;
    g_sort_dim = split_dim; // Set global dimension for qsort comparison

    // Create sortable structs
    SortParticle *sortable_indices = (SortParticle*)malloc(count * sizeof(SortParticle));
     if (!sortable_indices) {
        fprintf(stderr, "Error: Failed to allocate SortParticle array!\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<count; ++i) {
        sortable_indices[i].index = indices[i];
        switch(split_dim) { // Pre-fetch coordinate for sorting
            case 0: sortable_indices[i].coord = sys->posX[indices[i]]; break;
            case 1: sortable_indices[i].coord = sys->posY[indices[i]]; break;
            case 2: sortable_indices[i].coord = sys->posZ[indices[i]]; break;
        }
    }

    // Find median index (integer division intended)
    int median_idx_local = count / 2;
    // Use qsort to partition around the median (modifies sortable_indices)
    qsort(sortable_indices, count, sizeof(SortParticle), compare_sort_particles);

    // Update original index array based on sorted order
    for(int i=0; i<count; ++i) {
        indices[i] = sortable_indices[i].index;
    }
    free(sortable_indices); // Free temporary sortable array

    // Recursively build children
    // Left child: indices[0] to indices[median_idx_local - 1]
    // Right child: indices[median_idx_local] to indices[count - 1]
    node->left = build_kdtree(sys, indices, median_idx_local, depth + 1);
    node->right = build_kdtree(sys, indices + median_idx_local, count - median_idx_local, depth + 1);

    // Compute properties for this internal node based on children
    compute_node_properties(node, sys);

    return node;
}

// Recursively free the tree
void free_kdtree(KDTreeNode *node) {
    if (node == NULL) {
        return;
    }
    free_kdtree(node->left);
    free_kdtree(node->right);
    free(node);
}

// Recursively calculate force on a target particle from tree nodes
void calculate_force_from_node(int targetParticleIdx, KDTreeNode *node, const System *sys, double theta, double *fx, double *fy, double *fz) {
    if (node == NULL || node->numParticles == 0) {
        return; // Empty node, no force
    }

    // Check for self-interaction or interaction with empty leaf
    if (node->particleIndex == targetParticleIdx) {
        return; // Don't interact with self
    }

     // Vector from node's COM to target particle
    double dx = node->centerOfMassX - sys->posX[targetParticleIdx];
    double dy = node->centerOfMassY - sys->posY[targetParticleIdx];
    double dz = node->centerOfMassZ - sys->posZ[targetParticleIdx];

    double distSq = dx * dx + dy * dy + dz * dz;
    double dist = sqrt(distSq + EPSILON_SQ); // Softened distance

    // Barnes-Hut criterion: s/d < theta
    // Or if it's a leaf node containing a single particle (which isn't the target)
    if (node->particleIndex != -1 || (dist > 1e-15 && node->nodeSize / dist < theta)) {
         if (dist > 1e-10) { // Avoid division by zero and ensure mass exists
             double force_mag_over_dist = G * node->totalMass * sys->mass[targetParticleIdx] / (distSq + EPSILON_SQ) / dist;
             // Force acts from target towards node COM
             *fx += force_mag_over_dist * dx;
             *fy += force_mag_over_dist * dy;
             *fz += force_mag_over_dist * dz;
         }
    } else {
        // Node is too close/large, or internal: traverse children
        calculate_force_from_node(targetParticleIdx, node->left, sys, theta, fx, fy, fz);
        calculate_force_from_node(targetParticleIdx, node->right, sys, theta, fx, fy, fz);
    }
}

// Calculate forces for all particles using the kD-Tree (Parallelized)
void calculate_forces_tree(const System *sys, KDTreeNode *root, double theta, double *forceX, double *forceY, double *forceZ) {

    #pragma omp parallel for schedule(dynamic) // Dynamic schedule might be good
    for (int i = 0; i < sys->N; ++i) {
        double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;
        calculate_force_from_node(i, root, sys, theta, &fx_i, &fy_i, &fz_i);
        forceX[i] = fx_i;
        forceY[i] = fy_i;
        forceZ[i] = fz_i;
    }
}


// --- Simulation Step using kD-Tree ---
void simulate_step_tree(System *sys) {
    // --- 1. Build kD-Tree (Sequential) ---
    // Create an array of indices [0, 1, ..., N-1]
    int *particle_indices = (int*)malloc(sys->N * sizeof(int));
    if (!particle_indices) { fprintf(stderr, "Failed to allocate indices array!\n"); exit(EXIT_FAILURE); }
    for (int i = 0; i < sys->N; ++i) { particle_indices[i] = i; }

    KDTreeNode *tree_root = build_kdtree(sys, particle_indices, sys->N, 0);
    free(particle_indices); // Indices array no longer needed after build

    // --- 2. Calculate Forces using Tree (Parallel) ---
    double *forceX = (double*)malloc(sys->N * sizeof(double));
    double *forceY = (double*)malloc(sys->N * sizeof(double));
    double *forceZ = (double*)malloc(sys->N * sizeof(double));
     if (!forceX || !forceY || !forceZ) {
        fprintf(stderr, "Error: Force array allocation failed in simulate_step_tree!\n");
        exit(EXIT_FAILURE);
    }

    calculate_forces_tree(sys, tree_root, THETA, forceX, forceY, forceZ);

    // --- 3. Free kD-Tree ---
    free_kdtree(tree_root);

    // --- 4. Update Velocities (Parallel) ---
    #pragma omp parallel for
    for (int i = 0; i < sys->N; ++i) {
        if (sys->mass[i] > 0) {
            double inv_mass = 1.0 / sys->mass[i];
            sys->velX[i] += forceX[i] * inv_mass * DT;
            sys->velY[i] += forceY[i] * inv_mass * DT;
            sys->velZ[i] += forceZ[i] * inv_mass * DT;
        }
    }

    // --- 5. Update Positions (Parallel) ---
    #pragma omp parallel for
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