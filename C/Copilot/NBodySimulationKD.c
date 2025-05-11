#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define G 6.67430e-11 // Gravitational constant
#define DT 3600*24*365*1e-3 // Time step
#define THETA 0.3 // Theta value for approximation

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
} Body;

typedef struct KDNode {
    Body *body;
    struct KDNode *left, *right;
    double min[3], max[3];
} KDNode;

void initialize_bodies(Body *bodies, int n) {
    bodies[0].x = bodies[0].y = bodies[0].z = 0.0;
    bodies[0].vx = bodies[0].vy = bodies[0].vz = 0.0;
    bodies[0].mass = 1e30; // Central body mass

    for (int i = 1; i < n; i++) {
        double angle = 2 * M_PI * i / (n - 1);
        bodies[i].x = cos(angle) * 1e11;
        bodies[i].y = sin(angle) * 1e11;
        bodies[i].z = 0.0;
        bodies[i].vx = -sin(angle) * sqrt(G * bodies[0].mass / 1e11);
        bodies[i].vy = cos(angle) * sqrt(G * bodies[0].mass / 1e11);
        bodies[i].vz = 0.0;
        bodies[i].mass = 1e24 / n; // Small body mass
    }
}

double calculate_energy(Body *bodies, int n) {
    double energy = 0.0;
    #pragma omp parallel for reduction(+:energy)
    for (int i = 0; i < n; i++) {
        double kinetic = 0.5 * bodies[i].mass * (bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy + bodies[i].vz * bodies[i].vz);
        double potential = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dx = bodies[i].x - bodies[j].x;
                double dy = bodies[i].y - bodies[j].y;
                double dz = bodies[i].z - bodies[j].z;
                double distance = sqrt(dx * dx + dy * dy + dz * dz);
                potential -= G * bodies[i].mass * bodies[j].mass / distance;
            }
        }
        energy += kinetic + 0.5 * potential;
    }
    return energy;
}

// double calculate_energy(Body *bodies, int n) {
//     double energy = 0.0;
//     #pragma omp parallel for reduction(+:energy)
//     for (int i = 0; i < n; i++) {
//         double kinetic = 0.5 * bodies[i].mass * (bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy + bodies[i].vz * bodies[i].vz);
//         double potential = 0.0;
//         for (int j = 0; j < n; j++) {
//             if (i != j) {
//                 return energy;
// }

int compare_x(const void *a, const void *b) {
    Body *bodyA = (Body *)a;
    Body *bodyB = (Body *)b;
    return (bodyA->x < bodyB->x) ? -1 : (bodyA->x > bodyB->x);
}

int compare_y(const void *a, const void *b) {
    Body *bodyA = (Body *)a;
    Body *bodyB = (Body *)b;
    return (bodyA->y < bodyB->y) ? -1 : (bodyA->y > bodyB->y);
}

int compare_z(const void *a, const void *b) {
    Body *bodyA = (Body *)a;
    Body *bodyB = (Body *)b;
    return (bodyA->z < bodyB->z) ? -1 : (bodyA->z > bodyB->z);
}

KDNode* build_kdtree(Body *bodies, int n, int depth) {
    if (n <= 0) return NULL;

    int axis = depth % 3;
    int median = n / 2;

    // Sort bodies by the current axis
    qsort(bodies, n, sizeof(Body), (int (*)(const void *, const void *)) (axis == 0 ? compare_x : (axis == 1 ? compare_y : compare_z)));

    KDNode *node = (KDNode *)malloc(sizeof(KDNode));
    node->left =  build_kdtree(bodies, median, depth + 1);
    node->right = build_kdtree(bodies + median + 1, n - median - 1, depth + 1);

    for (int i = 0; i < 3; i++) {
        node->min[i] = node->max[i] = bodies[median].x;
        if (node->left) {
            node->min[i] = fmin(node->min[i], node->left->min[i]);
            node->max[i] = fmax(node->max[i], node->left->max[i]);
        }
        if (node->right) {
            node->min[i] = fmin(node->min[i], node->right->min[i]);
            node->max[i] = fmax(node->max[i], node->right->max[i]);
        }
    }

    return node;
}

void calculate_force(KDNode *node, Body *body, double *ax, double *ay, double *az) {
    if (!node || !node->body) return;

    double dx = node->body->x - body->x;
    double dy = node->body->y - body->y;
    double dz = node->body->z - body->z;
    double distance = sqrt(dx * dx + dy * dy + dz * dz);

    double size = fmax(node->max[0] - node->min[0], fmax(node->max[1] - node->min[1], node->max[2] - node->min[2]));

    if (size / distance < THETA || !node->left && !node->right) {
        double force = G * node->body->mass / (distance * distance * distance);
        *ax += force * dx;
        *ay += force * dy;
        *az += force * dz;
    } else {
        calculate_force(node->left, body, ax, ay, az);
        calculate_force(node->right, body, ax, ay, az);
    }
}

void kick_step(Body *bodies, int n, KDNode *root) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        double ax = 0.0, ay = 0.0, az = 0.0;
        calculate_force(root, &bodies[i], &ax, &ay, &az);
        bodies[i].vx += ax * DT;
        bodies[i].vy += ay * DT;
        bodies[i].vz += az * DT;
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        bodies[i].x += bodies[i].vx * DT;
        bodies[i].y += bodies[i].vy * DT;
        bodies[i].z += bodies[i].vz * DT;
    }
}

int main() {
    int n = 100000; // Number of bodies
    Body *bodies = (Body *)malloc(n * sizeof(Body));
    initialize_bodies(bodies, n);

    double initial_energy = calculate_energy(bodies, n);
    printf("Initial energy: %e\n", initial_energy);

    for (int step = 0; step < 10; step++) {
        KDNode *root = build_kdtree(bodies, n, 0);
        kick_step(bodies, n, root);
        free(root);
    }

    double final_energy = calculate_energy(bodies, n);
    printf("Final energy: %e\n", final_energy);

    free(bodies);
    return 0;
}
