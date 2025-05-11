//This response paused because Claude reached its max length for a message. Hit continue to nudge Claude along.  (I did and nothing happened).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <stdbool.h>
#include <immintrin.h>  // For SIMD operations

#define G 6.67430e-11   // Gravitational constant
#define SOLAR_MASS 1.989e30  // Mass of the Sun in kg
#define AU 1.496e11     // Astronomical Unit in meters
#define YEAR 3.154e7    // Year in seconds
#define THETA 0.3       // Barnes-Hut opening angle threshold
#define SOFTENING 1e-9  // Softening parameter to prevent singularities
#define MAX_DEPTH 20    // Maximum octree depth
#define MIN_NODE_SIZE 10 // Minimum number of bodies per node for parallel processing
#define BUCKET_SIZE 16  // Maximum number of bodies in a leaf node

// Structure for bodies with memory alignment for SIMD operations
typedef struct {
    double x, y, z;       // Position
    double vx, vy, vz;    // Velocity
    double ax, ay, az;    // Acceleration
    double mass;          // Mass
    int pad;              // Padding for 64-byte alignment
} __attribute__((aligned(64))) Body;

// Node structure for octree
typedef struct OctNode {
    double x_min, x_max;  // Boundaries
    double y_min, y_max;
    double z_min, z_max;
    
    double center_x, center_y, center_z;  // Center of mass
    double total_mass;                    // Total mass
    double size;                          // Node size (cube side length)
    
    bool is_leaf;                         // Flag for leaf nodes
    int num_bodies;                       // Number of bodies in this node
    int depth;                            // Depth in the tree
    
    int* body_indices;                    // Indices of bodies in this node if it's a leaf
    int capacity;                         // Capacity of body_indices array
    
    struct OctNode* children[8];          // Octants
} OctNode;

// Structure for acceleration cache to avoid recomputation
typedef struct {
    double* ax;
    double* ay;
    double* az;
    bool* calculated;
} AccelerationCache;

// Vector operations with SIMD optimization for groups of 4 doubles
void vector_distance_batch4(const double* x1, const double* y1, const double* z1,
                           const double* x2, const double* y2, const double* z2,
                           double* result) {
    #ifdef __AVX__
    // Use AVX SIMD operations if available
    __m256d vx1 = _mm256_loadu_pd(x1);
    __m256d vy1 = _mm256_loadu_pd(y1);
    __m256d vz1 = _mm256_loadu_pd(z1);
    
    __m256d vx2 = _mm256_loadu_pd(x2);
    __m256d vy2 = _mm256_loadu_pd(y2);
    __m256d vz2 = _mm256_loadu_pd(z2);
    
    __m256d dx = _mm256_sub_pd(vx2, vx1);
    __m256d dy = _mm256_sub_pd(vy2, vy1);
    __m256d dz = _mm256_sub_pd(vz2, vz1);
    
    __m256d dx2 = _mm256_mul_pd(dx, dx);
    __m256d dy2 = _mm256_mul_pd(dy, dy);
    __m256d dz2 = _mm256_mul_pd(dz, dz);
    
    __m256d sum = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
    __m256d dist = _mm256_sqrt_pd(sum);
    
    _mm256_storeu_pd(result, dist);
    #else
    // Fallback for non-AVX systems
    for (int i = 0; i < 4; i++) {
        double dx = x2[i] - x1[i];
        double dy = y2[i] - y1[i];
        double dz = z2[i] - z1[i];
        result[i] = sqrt(dx*dx + dy*dy + dz*dz);
    }
    #endif
}

// Regular vector distance function
double vector_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Initialize the acceleration cache
AccelerationCache* init_acceleration_cache(int n) {
    AccelerationCache* cache = (AccelerationCache*)malloc(sizeof(AccelerationCache));
    if (!cache) {
        fprintf(stderr, "Failed to allocate acceleration cache!\n");
        exit(1);
    }
    
    cache->ax = (double*)calloc(n, sizeof(double));
    cache->ay = (double*)calloc(n, sizeof(double));
    cache->az = (double*)calloc(n, sizeof(double));
    cache->calculated = (bool*)calloc(n, sizeof(bool));
    
    if (!cache->ax || !cache->ay || !cache->az || !cache->calculated) {
        fprintf(stderr, "Failed to allocate acceleration cache arrays!\n");
        exit(1);
    }
    
    return cache;
}

// Free the acceleration cache
void free_acceleration_cache(AccelerationCache* cache) {
    if (cache) {
        free(cache->ax);
        free(cache->ay);
        free(cache->az);
        free(cache->calculated);
        free(cache);
    }
}

// Reset the acceleration cache for the next step
void reset_acceleration_cache(AccelerationCache* cache, int n) {
    memset(cache->ax, 0, n * sizeof(double));
    memset(cache->ay, 0, n * sizeof(double));
    memset(cache->az, 0, n * sizeof(double));
    memset(cache->calculated, 0, n * sizeof(bool));
}

// Create a new node
OctNode* create_node(double x_min, double x_max, double y_min, double y_max, 
                   double z_min, double z_max, int depth) {
    OctNode* node = (OctNode*)malloc(sizeof(OctNode));
    if (!node) {
        fprintf(stderr, "Memory allocation failed for octree node!\n");
        exit(1);
    }
    
    node->x_min = x_min;
    node->x_max = x_max;
    node->y_min = y_min;
    node->y_max = y_max;
    node->z_min = z_min;
    node->z_max = z_max;
    
    node->center_x = 0.0;
    node->center_y = 0.0;
    node->center_z = 0.0;
    node->total_mass = 0.0;
    node->size = x_max - x_min;  // Assuming cubic nodes
    
    node->is_leaf = true;
    node->num_bodies = 0;
    node->depth = depth;
    
    // Allocate space for body indices (with initial capacity)
    node->capacity = BUCKET_SIZE;
    node->body_indices = (int*)malloc(node->capacity * sizeof(int));
    if (!node->body_indices) {
        fprintf(stderr, "Memory allocation failed for body indices!\n");
        free(node);
        exit(1);
    }
    
    for (int i = 0; i < 8; i++) {
        node->children[i] = NULL;
    }
    
    return node;
}

// Determine which octant a body belongs to
int get_octant(OctNode* node, Body* body) {
    int octant = 0;
    double mid_x = (node->x_min + node->x_max) / 2.0;
    double mid_y = (node->y_min + node->y_max) / 2.0;
    double mid_z = (node->z_min + node->z_max) / 2.0;
    
    if (body->x >= mid_x) octant |= 1;
    if (body->y >= mid_y) octant |= 2;
    if (body->z >= mid_z) octant |= 4;
    
    return octant;
}

// Create child nodes for each octant
void subdivide_node(OctNode* node) {
    double mid_x = (node->x_min + node->x_max) / 2.0;
    double mid_y = (node->y_min + node->y_max) / 2.0;
    double mid_z = (node->z_min + node->z_max) / 2.0;
    int next_depth = node->depth + 1;
    
    node->children[0] = create_node(node->x_min, mid_x, node->y_min, mid_y, node->z_min, mid_z, next_depth);
    node->children[1] = create_node(mid_x, node->x_max, node->y_min, mid_y, node->z_min, mid_z, next_depth);
    node->children[2] = create_node(node->x_min, mid_x, mid_y, node->y_max, node->z_min, mid_z, next_depth);
    node->children[3] = create_node(mid_x, node->x_max, mid_y, node->y_max, node->z_min, mid_z, next_depth);
    node->children[4] = create_node(node->x_min, mid_x, node->y_min, mid_y, mid_z, node->z_max, next_depth);
    node->children[5] = create_node(mid_x, node->x_max, node->y_min, mid_y, mid_z, node->z_max, next_depth);
    node->children[6] = create_node(node->x_min, mid_x, mid_y, node->y_max, mid_z, node->z_max, next_depth);
    node->children[7] = create_node(mid_x, node->x_max, mid_y, node->y_max, mid_z, node->z_max, next_depth);
    
    node->is_leaf = false;
}

// Insert a body into the octree
void insert_body(OctNode* node, Body* bodies, int body_index) {
    // Update center of mass and total mass
    double new_mass = node->total_mass + bodies[body_index].mass;
    
    if (new_mass > 0) {
        node->center_x = (node->center_x * node->total_mass + bodies[body_index].x * bodies[body_index].mass) / new_mass;
        node->center_y = (node->center_y * node->total_mass + bodies[body_index].y * bodies[body_index].mass) / new_mass;
        node->center_z = (node->center_z * node->total_mass + bodies[body_index].z * bodies[body_index].mass) / new_mass;
        node->total_mass = new_mass;
    }
    
    // If node is a leaf and has space, add the body to this node
    if (node->is_leaf) {
        if (node->num_bodies < node->capacity) {
            node->body_indices[node->num_bodies++] = body_index;
            return;
        }
        
        // If we've reached max depth, expand the array instead of subdividing
        if (node->depth >= MAX_DEPTH) {
            node->capacity *= 2;
            node->body_indices = (int*)realloc(node->body_indices, node->capacity * sizeof(int));
            if (!node->body_indices) {
                fprintf(stderr, "Memory reallocation failed for body indices!\n");
                exit(1);
            }
            
            node->body_indices[node->num_bodies++] = body_index;
            return;
        }
        
        // Node is full and not at max depth, subdivide it
        subdivide_node(node);
        
        // Redistribute existing bodies to children
        for (int i = 0; i < node->num_bodies; i++) {
            int idx = node->body_indices[i];
            int octant = get_octant(node, &bodies[idx]);
            insert_body(node->children[octant], bodies, idx);
        }
        
        // Free the body indices array as it's no longer needed
        free(node->body_indices);
        node->body_indices = NULL;
        node->num_bodies = 0;  // Reset for internal nodes
    }
    
    // Insert the new body into the appropriate child
    if (!node->is_leaf) {
        int octant = get_octant(node, &bodies[body_index]);
        insert_body(node->children[octant], bodies, body_index);
    }
}

// Free the octree memory
void free_octree(OctNode* node) {
    if (!node) return;
    
    for (int i = 0; i < 8; i++) {
        if (node->children[i]) {
            free_octree(node->children[i]);
        }
    }
    
    if (node->body_indices) {
        free(node->body_indices);
    }
    
    free(node);
}

// Calculate boundaries for the octree with additional margin
void calculate_boundaries(Body* bodies, int n, double* x_min, double* x_max, 
                         double* y_min, double* y_max, double* z_min, double* z_max) {
    *x_min = *y_min = *z_min = INFINITY;
    *x_max = *y_max = *z_max = -INFINITY;
    
    #pragma omp parallel reduction(min:*x_min,*y_min,*z_min) reduction(max:*x_max,*y_max,*z_max)
    {
        #pragma omp for
        for (int i = 0; i < n; i++) {
            if (bodies[i].x < *x_min) *x_min = bodies[i].x;
            if (bodies[i].x > *x_max) *x_max = bodies[i].x;
            if (bodies[i].y < *y_min) *y_min = bodies[i].y;
            if (bodies[i].y > *y_max) *y_max = bodies[i].y;
            if (bodies[i].z < *z_min) *z_min = bodies[i].z;
            if (bodies[i].z > *z_max) *z_max = bodies[i].z;
        }
    }
    
    // Add some padding (10%)
    double padding_x = 0.1 * (*x_max - *x_min);
    double padding_y = 0.1 * (*y_max - *y_min);
    double padding_z = 0.1 * (*z_max - *z_min);
    
    // Ensure minimum padding
    padding_x = fmax(padding_x, 1.0e9);  // 1 billion meters minimum
    padding_y = fmax(padding_y, 1.0e9);
    padding_z = fmax(padding_z, 1.0e9);
    
    *x_min -= padding_x;
    *x_max += padding_x;
    *y_min -= padding_y;
    *y_max += padding_y;
    *z_min -= padding_z;
    *z_max += padding_z;
    
    // Make sure boundaries are cubic for octree
    double max_range = fmax(*x_max - *x_min, fmax(*y_max - *y_min, *z_max - *z_min));
    double x_center = (*x_min + *x_max) / 2.0;
    double y_center = (*y_min + *y_max) / 2.0;
    double z_center = (*z_min + *z_max) / 2.0;
    
    *x_min = x_center - max_range / 2.0;
    *x_max = x_center + max_range / 2.0;
    *y_min = y_center - max_range / 2.0;
    *y_max = y_center + max_range / 2.0;
    *z_min = z_center - max_range / 2.0;
    *z_max = z_center + max_range / 2.0;
}

// Build the octree for a set of bodies using parallel construction when appropriate
OctNode* build_octree(Body* bodies, int n) {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    calculate_boundaries(bodies, n, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);
    
    OctNode* root = create_node(x_min, x_max, y_min, y_max, z_min, z_max, 0);
    
    // Insert bodies into octree
    #pragma omp parallel for schedule(dynamic, 1000) if(n > 100000)
    for (int i = 0; i < n; i++) {
        #pragma omp critical
        {
            insert_body(root, bodies, i);
        }
    }
    
    return root;
}

// Compute acceleration for a single body with fast approximations for distant interactions
void compute_acceleration(OctNode* node, Body* bodies, int body_index, double* ax, double* ay, double* az) {
    // Skip empty nodes
    if (node->num_bodies == 0) return;
    
    // If node is a leaf, calculate direct accelerations for all bodies in it
    if (node->is_leaf) {
        for (int i = 0; i < node->num_bodies; i++) {
            int j = node->body_indices[i];
            if (j == body_index) continue;  // Skip self-interaction
            
            double dx = bodies[j].x - bodies[body_index].x;
            double dy = bodies[j].y - bodies[body_index].y;
            double dz = bodies[j].z - bodies[body_index].z;
            
            double distance_squared = dx*dx + dy*dy + dz*dz + SOFTENING*SOFTENING;
            double distance = sqrt(distance_squared);
            double distance_cubed = distance * distance_squared;
            
            double force = G * bodies[j].mass / distance_cubed;
            
            *ax += force * dx;
            *ay += force * dy;
            *az += force * dz;
        }
        return;
    }
    
    // For non-leaf nodes, check if we can use approximation
    double dx = node->center_x - bodies[body_index].x;
    double dy = node->center_y - bodies[body_index].y;
    double dz = node->center_z - bodies[body_index].z;
    
    double distance_squared = dx*dx + dy*dy + dz*dz;
    
    // Use size/distance < theta as criterion for approximation
    if (node->size * node->size < THETA * THETA * distance_squared) {
        // Use approximation - treat node as a single body at center of mass
        double distance = sqrt(distance_squared + SOFTENING*SOFTENING);
        double distance_cubed = distance * distance_squared;
        
        double force = G * node->total_mass / distance_cubed;
        
        *ax += force * dx;
        *ay += force * dy;
        *az += force * dz;
    } else {
        // Node is too close, recurse into children
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                compute_acceleration(node->children[i], bodies, body_index, ax, ay, az);
            }
        }
    }
}

// Update accelerations of all bodies using the Barnes-Hut approximation
void compute_all_accelerations(Body* bodies, int n, OctNode* root, AccelerationCache* cache) {
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < n; i++) {
        if (!cache->calculated[i]) {
            double ax = 0.0, ay = 0.0, az = 0.0;
            
            // Compute acceleration using Barnes-Hut approximation
            compute_acceleration(root, bodies, i, &ax, &ay, &az);
            
            // Cache the results
            cache->ax[i] = ax;
            cache->ay[i] = ay;
            cache->az[i] = az;
            cache->calculated[i] = true;
        }
        
        // Set accelerations from cache
        bodies[i].ax = cache->ax[i];
        bodies[i].ay = cache->ay[i];
        bodies[i].az = cache->az[i];
    }
}

// Leapfrog integrator first half-step (kick)
void leapfrog_kick(Body* bodies, int n, double dt_half) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        bodies[i].vx += bodies[i].ax * dt_half;
        bodies[i].vy += bodies[i].ay * dt_half;
        bodies[i].vz += bodies[i].az * dt_half;
    }
}

// Update positions (drift)
void leapfrog_drift(Body* bodies, int n, double dt) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

// Leapfrog integrator second half-step (kick)
void leapfrog_kick_final(Body* bodies, int n, double dt_half) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        bodies[i].vx += bodies[i].ax * dt_half;
        bodies[i].vy += bodies[i].ay * dt_half;
        bodies[i].vz += bodies[i].az * dt_half;
    }
}

// Calculate total energy of the system (kinetic + potential) using Barnes-Hut for potential energy
double calculate_energy(Body* bodies, int n, OctNode* root) {
    double total_energy = 0.0;
    
    #pragma omp parallel
    {
        double energy = 0.0;
        
        // Calculate kinetic energy
        #pragma omp for schedule(dynamic, 1000) nowait
        for (int i = 0; i < n; i++) {
            double v_squared = bodies[i].vx * bodies[i].vx + 
                              bodies[i].vy * bodies[i].vy + 
                              bodies[i].vz * bodies[i].vz;
            energy += 0.5 * bodies[i].mass * v_squared;
        }
        
        // Calculate potential energy using direct calculation or Barnes-Hut approximation
        // Direct calculation is more accurate for energy conservation checks
        #pragma omp for schedule(dynamic, 100) reduction(+:energy)
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double distance = vector_distance(
                    bodies[i].x, bodies[i].y, bodies[i].z,
                    bodies[j].x, bodies[j].y, bodies[j].z
                );
                distance = sqrt(distance * distance + SOFTENING * SOFTENING);
                energy -= G * bodies[i].mass * bodies[j].mass / distance;
            }
        }
        
        // Sum up the energy from all threads
        #pragma omp atomic
        total_energy += energy;
    }
    
    return total_energy;
}

// Leapfrog integration with Barnes-Hut approximation
void simulate_step_barnes_hut(Body* bodies, int n, double dt, AccelerationCache* cache) {
    double dt_half = dt * 0.5;
    
    // Reset acceleration cache
    reset_acceleration_cache(cache, n);
    
    // Build octree
    OctNode* root = build_octree(bodies, n);
    
    // Compute initial accelerations
    compute_all_accelerations(bodies, n, root, cache);
    
    // First half of velocity update (kick)
    leapfrog_kick(bodies, n, dt_half);
    
    // Position update (drift)
    leapfrog_drift(bodies, n, dt);
    
    // Free the old octree
    free_octree(root);
    
    // Reset cache for new accelerations
    reset_acceleration_cache(cache, n);
    
    // Rebuild octree with new positions
    root = build_octree(bodies, n);
    
    // Compute accelerations with new positions
    compute_all_accelerations(bodies, n, root, cache);
    
    // Second half of velocity update (kick)
    leapfrog_kick_final(bodies, n, dt_half);
    
    // Free the octree
    free_octree(root);
}

// Initialize system with central body and orbiting bodies
void initialize_system(Body* bodies, int n_bodies, double central_mass, double orbit_radius_min, double orbit_radius_max) {
    if (n_bodies < 1) return;
    
    // Set central body at origin
    bodies[0].x = 0.0;
    bodies[0].y = 0.0;
    bodies[0].z = 0.0;
    bodies[0].vx = 0.0;
    bodies[0].vy = 0.0;
    bodies[0].vz = 0.0;
    bodies[0].ax = 0.0;
    bodies[0].ay = 0.0;
    bodies[0].az = 0.0;
    bodies[0].mass = central_mass;
    
    // Use thread-safe random number generation
    unsigned int seed = time(NULL);
    
    // Pre-compute some values to avoid redundant calculations
    const double orbit_range = orbit_radius_max - orbit_radius_min;
    const double two_pi = 2.0 * M_PI;
    const double half_pi = M_PI / 2.0;
    
    // Create orbiting bodies in parallel
    #pragma omp parallel
    {
        // Each thread needs its own random seed
        unsigned int local_seed = seed + omp_get_thread_num();
        
        #pragma omp for
        for (int i = 1; i < n_bodies; i++) {
            // Thread-safe random number generation with better efficiency
            double rand1 = (double)rand_r(&local_seed) / RAND_MAX;
            double rand2 = (double)rand_r(&local_seed) / RAND_MAX;
            double rand3 = (double)rand_r(&local_seed) / RAND_MAX;
            
            double phi = rand1 * two_pi;
            double theta = rand2 * M_PI - half_pi;
            
            // Distribute radii more uniformly
            double radius = orbit_radius_min + pow(rand3, 0.5) * orbit_range;
            
            // Position on sphere using trig
            double cos_theta = cos(theta);
            bodies[i].x = radius * cos_theta * cos(phi);
            bodies[i].y = radius * cos_theta * sin(phi);
            bodies[i].z = radius * sin(theta);
            
            // Velocity for circular orbit
            double orbit_speed = sqrt(G * central_mass / radius);
            
            // Cross product with up vector to get perpendicular velocity
            // Pre-compute sin/cos for efficiency
            double sin_phi = sin(phi);
            double cos_phi = cos(phi);
            bodies[i].vx = -orbit_speed * sin_phi;
            bodies[i].vy = orbit_speed * cos_phi;
            bodies[i].vz = 0.0;
            
            // Initialize acceleration to zero
            bodies[i].ax = 0.0;
            bodies[i].ay = 0.0;
            bodies[i].az = 0.0;
            
            // Use more realistic mass distribution
            // Small bodies have mass that varies by up to an order of magnitude
            double mass_factor = 0.1 + 0.9 * (double)rand_r(&local_seed) / RAND_MAX;
            bodies[i].mass = central_mass * 1e-6 * mass_factor;
        }
    }
}

// void print_system_stats(Body* bodies, int n) {
//     printf("System statistics:\n");
//     printf("Central body: position (%.2e, %.2e, %.2e), mass: %.2e\n", 
//            bodies[0].x, bodies[0].y, bodies[0].z, bodies[0].mass);
    
//     // Calculate system boundaries and statistics
//     double min_x = bodies[0].x, max_x = bodies[0].x;
//     double min_y = bodies[0].y, max_y = bodies[0].y;
//     double min_z = bodies[0].z, max_z = bodies[0].z;
    
//     double avg_distance = 0.0;
//     double max_velocity = 0.0;
//     double total_mass = bodies[0].mass;
    
//     #pragma omp parallel
//     {
//         double local_min_x = min_x, local_max_x = max_x;
//         double local_min_y = min_y, local_max_y = max_y;
//         double local_min_z = min_z, local_max_z = max_z;
//         double local_avg_distance = 0.0;
//         double local_max_velocity = 0.0;
//         double local_total_mass = 0.0;
        
//         #pragma omp for
//         for (int i = 1; i < n; i++) {
//             if (bodies[i].x < local_min_x) local_min_x = bodies[i].x;
//             if (bodies[i].x > local_max_x) local_max_x = bodies[i].x;
//             if (bodies[i].y < local_min_y) local_min_y = bodies[i].y;
//             if (bodies[i].y > local_max_y) local_max_y = bodies[i].y;
//             if (bodies[i].z < local_min_z) local_min_z = bodies[i].z;
//             if (bodies[i].z > local_max_z) local_max_z = bodies[i].z;
            
//             double distance = sqrt(bodies[i].x * bodies[i].x + 
//                                   bodies[i].y * bodies[i].y + 
//                                   bodies[i].z * bodies[i].z);
//             local_avg_distance += distance;
            
//             double velocity = sqrt(bodies[i].vx * bodies[i].vx + 
//                                   bodies[i].vy * bodies[i].vy + 
//                                   bodies[i].v

void print_system_stats(Body* bodies, int n) {
    printf("System statistics:\n");
    printf("Central body: position (%.2e, %.2e, %.2e), mass: %.2e\n", 
           bodies[0].x, bodies[0].y, bodies[0].z, bodies[0].mass);
    
    // Calculate system boundaries
    double min_x = bodies[0].x, max_x = bodies[0].x;
    double min_y = bodies[0].y, max_y = bodies[0].y;
    double min_z = bodies[0].z, max_z = bodies[0].z;
    
    for (int i = 1; i < n; i++) {
        if (bodies[i].x < min_x) min_x = bodies[i].x;
        if (bodies[i].x > max_x) max_x = bodies[i].x;
        if (bodies[i].y < min_y) min_y = bodies[i].y;
        if (bodies[i].y > max_y) max_y = bodies[i].y;
        if (bodies[i].z < min_z) min_z = bodies[i].z;
        if (bodies[i].z > max_z) max_z = bodies[i].z;
    }
    
    printf("System boundaries:\n");
    printf("  X: [%.2e, %.2e]\n", min_x, max_x);
    printf("  Y: [%.2e, %.2e]\n", min_y, max_y);
    printf("  Z: [%.2e, %.2e]\n", min_z, max_z);
}

int main() {
    // Problem parameters
    const int N = 100000 + 1;  // 1 central + 1 million orbiting bodies
    const int STEPS = 10;
    const double dt = 0.01 * YEAR;  // Timestep of 0.01 years
    
    // OpenMP setup
    int num_threads = omp_get_max_threads();
    printf("Starting N-body simulation with %d bodies for %d steps using %d threads\n", 
           N, STEPS, num_threads);
    printf("Using Barnes-Hut approximation with theta = %.2f\n", THETA);
    printf("Allocating memory...\n");
    
    // Allocate memory for bodies
    Body* bodies = (Body*)malloc(N * sizeof(Body));
    if (!bodies) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    printf("Initializing system...\n");
    
    // Initialize system with Sun-like central body and Earth-like orbits
    initialize_system(bodies, N, SOLAR_MASS, 0.8 * AU, 5.0 * AU);
    
    print_system_stats(bodies, N);
    
    // Calculate initial energy
    printf("Calculating initial energy...\n");
    double initial_energy = calculate_energy(bodies, N);
    printf("Initial energy: %.10e\n", initial_energy);
    
    printf("Starting simulation...\n");
    double start_time = omp_get_wtime();
    
    // Run simulation
    for (int step = 0; step < STEPS; step++) {
        if (step % 100 == 0) {
            printf("Step %d of %d (%.1f%%)\n", step, STEPS, 100.0 * step / STEPS);
        }
        simulate_step_barnes_hut(bodies, N, dt);
    }
    
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("Simulation completed in %.2f seconds\n", elapsed_time);
    
    // Calculate final energy
    printf("Calculating final energy...\n");
    double final_energy = calculate_energy(bodies, N);
    printf("Final energy: %.10e\n", final_energy);
    
    // Calculate energy error
    double energy_error = (final_energy - initial_energy) / initial_energy;
    printf("Relative energy error: %.10e (%.10f%%)\n", energy_error, energy_error * 100.0);
    
    print_system_stats(bodies, N);
    
    // Free memory
    free(bodies);
    
    return 0;
}