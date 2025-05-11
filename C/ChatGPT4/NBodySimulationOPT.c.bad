#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define N_BODIES 100000
#define STEPS 10
#define G 6.67430e-11
#define DT 1.0
#define THETA 0.3

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double mass;
} Body;

typedef struct KDNode {
    double min[3], max[3];
    double cm_x, cm_y, cm_z;
    double mass;
    int start, end;
    struct KDNode *left, *right;
} KDNode;

Body* bodies;

int compare_axis;
int cmp_body(const void* a, const void* b) {
    Body* b1 = (Body*)a;
    Body* b2 = (Body*)b;
    if (compare_axis == 0) return (b1->x > b2->x) - (b1->x < b2->x);
    if (compare_axis == 1) return (b1->y > b2->y) - (b1->y < b2->y);
    return (b1->z > b2->z) - (b1->z < b2->z);
}

void compute_bounds(int start, int end, double min[3], double max[3]) {
    for (int i = 0; i < 3; ++i) min[i] = 1e30, max[i] = -1e30;
    for (int i = start; i < end; ++i) {
        if (bodies[i].x < min[0]) min[0] = bodies[i].x;
        if (bodies[i].y < min[1]) min[1] = bodies[i].y;
        if (bodies[i].z < min[2]) min[2] = bodies[i].z;
        if (bodies[i].x > max[0]) max[0] = bodies[i].x;
        if (bodies[i].y > max[1]) max[1] = bodies[i].y;
        if (bodies[i].z > max[2]) max[2] = bodies[i].z;
    }
}

KDNode* build_kdtree(int start, int end) {
    if (end <= start) return NULL;
    KDNode* node = (KDNode*)malloc(sizeof(KDNode));
    node->start = start;
    node->end = end;
    compute_bounds(start, end, node->min, node->max);

    node->mass = 0.0;
    node->cm_x = node->cm_y = node->cm_z = 0.0;
    for (int i = start; i < end; i++) {
        node->mass += bodies[i].mass;
        node->cm_x += bodies[i].x * bodies[i].mass;
        node->cm_y += bodies[i].y * bodies[i].mass;
        node->cm_z += bodies[i].z * bodies[i].mass;
    }
    node->cm_x /= node->mass;
    node->cm_y /= node->mass;
    node->cm_z /= node->mass;

    int axis = 0;
    double len[3] = {node->max[0] - node->min[0], node->max[1] - node->min[1], node->max[2] - node->min[2]};
    if (len[1] > len[0]) axis = 1;
    if (len[2] > len[axis]) axis = 2;
    compare_axis = axis;

    int mid = (start + end) / 2;
    qsort(&bodies[start], end - start, sizeof(Body), cmp_body);

    #pragma omp task shared(node)
    node->left = build_kdtree(start, mid);
    #pragma omp task shared(node)
    node->right = build_kdtree(mid, end);
    #pragma omp taskwait

    return node;
}

void free_kdtree(KDNode* node) {
    if (!node) return;
    free_kdtree(node->left);
    free_kdtree(node->right);
    free(node);
}

void compute_force_from_node(Body* b, KDNode* node) {
    if (!node || node->start == node->end) return;

    double dx = node->cm_x - b->x;
    double dy = node->cm_y - b->y;
    double dz = node->cm_z - b->z;
    double dist2 = dx*dx + dy*dy + dz*dz + 1e-10;
    double size = fmax(fmax(node->max[0] - node->min[0], node->max[1] - node->min[1]), node->max[2] - node->min[2]);

    if ((node->end - node->start <= 1) || (size / sqrt(dist2) < THETA)) {
        double dist = sqrt(dist2);
        double force = G * node->mass / (dist2 * dist);
        b->ax += dx * force;
        b->ay += dy * force;
        b->az += dz * force;
    } else {
        compute_force_from_node(b, node->left);
        compute_force_from_node(b, node->right);
    }
}

void compute_forces_kdtree(int n) {
    KDNode* root;
    #pragma omp parallel
    {
        #pragma omp single nowait
        root = build_kdtree(0, n+1);
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0; i <= n; i++) {
        bodies[i].ax = bodies[i].ay = bodies[i].az = 0.0;
        compute_force_from_node(&bodies[i], root);
    }

    free_kdtree(root);
}

void initialize_system(int n) {
    bodies = (Body*)malloc((n + 1) * sizeof(Body));

    bodies[0].x = bodies[0].y = bodies[0].z = 0.0;
    bodies[0].vx = bodies[0].vy = bodies[0].vz = 0.0;
    bodies[0].mass = 1e20;

    double radius = 1e7;
    double central_mass = bodies[0].mass;

    double* radii = (double*)malloc((n + 1) * sizeof(double));
    for (int i = 1; i <= n; i++) {
        radii[i] = radius * (1.0 + 0.1 * ((double)rand() / RAND_MAX));
    }

    #pragma omp parallel for schedule(static)
    for (int i = 1; i <= n; i++) {
        double angle = 2.0 * M_PI * i / n;
        double r = radii[i];

        bodies[i].x = r * cos(angle);
        bodies[i].y = r * sin(angle);
        bodies[i].z = 0.0;

        double v = sqrt(G * central_mass / r);
        bodies[i].vx = -v * sin(angle);
        bodies[i].vy = v * cos(angle);
        bodies[i].vz = 0.0;

        bodies[i].mass = 1.0;
    }

    free(radii);
}

double calculate_energy(int n) {
    double kinetic = 0.0, potential = 0.0;

    #pragma omp parallel for reduction(+:kinetic) schedule(static)
    for (int i = 0; i <= n; i++) {
        double v2 = bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy + bodies[i].vz * bodies[i].vz;
        kinetic += 0.5 * bodies[i].mass * v2;
    }

    #pragma omp parallel for reduction(+:potential) schedule(dynamic)
    for (int i = 0; i <= n; i++) {
        for (int j = i + 1; j <= n; j++) {
            double dx = bodies[i].x - bodies[j].x;
            double dy = bodies[i].y - bodies[j].y;
            double dz = bodies[i].z - bodies[j].z;
            double dist = sqrt(dx * dx + dy * dy + dz * dz + 1e-10);
            potential -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }

    return kinetic + potential;
}

void kick_step(int n) {
    compute_forces_kdtree(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i <= n; i++) {
        bodies[i].vx += bodies[i].ax * DT;
        bodies[i].vy += bodies[i].ay * DT;
        bodies[i].vz += bodies[i].az * DT;

        bodies[i].x += bodies[i].vx * DT;
        bodies[i].y += bodies[i].vy * DT;
        bodies[i].z += bodies[i].vz * DT;
    }
}

int main() {
    clock_t start = clock();
    printf("Initializing system...\n");
    initialize_system(N_BODIES);

    printf("Calculating initial energy...\n");
    double initial_energy = calculate_energy(N_BODIES);
    printf("Initial total energy: %.6e\n", initial_energy);

    for (int step = 0; step < STEPS; step++) {
        kick_step(N_BODIES);
        if (step % 100 == 0) printf("Step %d complete\n", step);
    }

    printf("Calculating final energy...\n");
    double final_energy = calculate_energy(N_BODIES);
    printf("Final total energy: %.6e\n", final_energy);
    printf("Energy difference: %.6e\n", fabs(final_energy - initial_energy));

    free(bodies);
    clock_t end = clock();
    printf("Simulation completed in %.2f seconds.\n", (double)(end - start) / CLOCKS_PER_SEC);
    return 0;
}