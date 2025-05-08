#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_BODIES 10000
#define STEPS 100
#define G 6.67430e-11
#define DT 1.0

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double mass;
} Body;

Body* bodies;

void initialize_system(int n) {
    bodies = (Body*)malloc((n + 1) * sizeof(Body));

    // Central body at origin
    bodies[0].x = bodies[0].y = bodies[0].z = 0.0;
    bodies[0].vx = bodies[0].vy = bodies[0].vz = 0.0;
    bodies[0].mass = 1e20;

    double radius = 1e7;
    double central_mass = bodies[0].mass;

    for (int i = 1; i <= n; i++) {
        double angle = 2.0 * M_PI * i / n;
        double r = radius * (1.0 + 0.1 * ((double)rand() / RAND_MAX));

        bodies[i].x = r * cos(angle);
        bodies[i].y = r * sin(angle);
        bodies[i].z = 0.0;

        double v = sqrt(G * central_mass / r);
        bodies[i].vx = -v * sin(angle);
        bodies[i].vy = v * cos(angle);
        bodies[i].vz = 0.0;

        bodies[i].mass = 1.0;
    }
}

double calculate_energy(int n) {
    double kinetic = 0.0, potential = 0.0;
    for (int i = 0; i <= n; i++) {
        double v2 = bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy + bodies[i].vz * bodies[i].vz;
        kinetic += 0.5 * bodies[i].mass * v2;

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

void compute_forces(int n) {
    for (int i = 0; i <= n; i++) {
        bodies[i].ax = bodies[i].ay = bodies[i].az = 0.0;
        for (int j = 0; j <= n; j++) {
            if (i == j) continue;
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;
            double dist_sqr = dx * dx + dy * dy + dz * dz + 1e-10;
            double dist = sqrt(dist_sqr);
            double force = G * bodies[j].mass / (dist_sqr * dist);
            bodies[i].ax += dx * force;
            bodies[i].ay += dy * force;
            bodies[i].az += dz * force;
        }
    }
}

void kick_step(int n) {
    compute_forces(n);
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