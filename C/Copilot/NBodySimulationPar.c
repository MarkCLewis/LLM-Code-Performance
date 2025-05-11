#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define G 6.67430e-11 // Gravitational constant
#define DT 3600*24*365*1e-3 // Time step

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
} Body;

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

// double calculate_energy(Body *bodies, int n) {
//     double energy = 0.0;
//     #pragma omp parallel for reduction(+:energy)
//     for (int i = 0; i < n; i++) {
//         double kinetic = 0.5 * bodies[i].mass * (bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy + bodies[i].vz * bodies[i].vz);
//         double potential = 0.0;
//         for (int j = 0; j < n; j++) {
//             if (i != j) {
//                 double dx = bodies[i].x - bodies[j].x;
//                 double dy = bodies[i].y - bodies[j].y;
//                 double dz = bodies[i].z - bodies[j].z;
//                 double distance = sqrt(dx * dx + dy * dy + dz * dz);
//                 potential -= G * bodies[i].mass * bodies[j].mass / distance;
//             }
//         }
//         energy += kinetic + 0.5 * potential;
//     }
//     return energy;
// }

void kick_step(Body *bodies, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        double ax = 0.0, ay = 0.0, az = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                double distance = sqrt(dx * dx + dy * dy + dz * dz);
                double force = G * bodies[j].mass / (distance * distance * distance);
                ax += force * dx;
                ay += force * dy;
                az += force * dz;
            }
        }
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
    int n = 10000; // Number of bodies
    Body *bodies = (Body *)malloc(n * sizeof(Body));
    initialize_bodies(bodies, n);

    // double initial_energy = calculate_energy(bodies, n);
    // printf("Initial energy: %e\n", initial_energy);

    for (int step = 0; step < 100; step++) {
        kick_step(bodies, n);
    }

    // double final_energy = calculate_energy(bodies, n);
    // printf("Final energy: %e\n", final_energy);

    free(bodies);
    return 0;
}
