#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <omp.h>

const double G = 6.67430e-11; // Gravitational constant
const double dt = 1e-3;       // Time step
const int numBodies = 1000000; // Number of small bodies

struct Body {
    double mass;
    double pos[3];
    double vel[3];
    double force[3];
};

void initializeBodies(std::vector<Body>& bodies, int n);
void simulateStep(std::vector<Body>& bodies);
double calculateEnergy(const std::vector<Body>& bodies);

int main() {
    srand(time(0));

    // Initialize bodies
    std::vector<Body> bodies(numBodies + 1);
    initializeBodies(bodies, numBodies);

    // Calculate initial energy
    double initialEnergy = calculateEnergy(bodies);
    std::cout << "Initial Energy: " << initialEnergy << std::endl;

    // Run simulation for 1000 steps
    for (int step = 0; step < 1000; ++step) {
        simulateStep(bodies);
    }

    // Calculate final energy
    double finalEnergy = calculateEnergy(bodies);
    std::cout << "Final Energy: " << finalEnergy << std::endl;

    return 0;
}

void initializeBodies(std::vector<Body>& bodies, int n) {
    // Central body
    bodies[0] = {1e30, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    // Small bodies in circular orbits
    for (int i = 1; i <= n; ++i) {
        double radius = 1e11 * (rand() / (double)RAND_MAX);
        double speed = std::sqrt(G * bodies[0].mass / radius);
        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);

        bodies[i] = {1e24, 
                     {radius * std::cos(angle), radius * std::sin(angle), 0},
                     {-speed * std::sin(angle), speed * std::cos(angle), 0},
                     {0, 0, 0}};
    }
}

void simulateStep(std::vector<Body>& bodies) {
    // Reset forces
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].force[0] = bodies[i].force[1] = bodies[i].force[2] = 0;
    }

    // Calculate forces
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            double dx = bodies[j].pos[0] - bodies[i].pos[0];
            double dy = bodies[j].pos[1] - bodies[i].pos[1];
            double dz = bodies[j].pos[2] - bodies[i].pos[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            double force = G * bodies[i].mass * bodies[j].mass / (dist * dist * dist);

            #pragma omp atomic
            bodies[i].force[0] += force * dx;
            #pragma omp atomic
            bodies[i].force[1] += force * dy;
            #pragma omp atomic
            bodies[i].force[2] += force * dz;

            #pragma omp atomic
            bodies[j].force[0] -= force * dx;
            #pragma omp atomic
            bodies[j].force[1] -= force * dy;
            #pragma omp atomic
            bodies[j].force[2] -= force * dz;
        }
    }

    // Update velocities and positions
    #pragma omp parallel for
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].vel[0] += bodies[i].force[0] / bodies[i].mass * dt;
        bodies[i].vel[1] += bodies[i].force[1] / bodies[i].mass * dt;
        bodies[i].vel[2] += bodies[i].force[2] / bodies[i].mass * dt;

        bodies[i].pos[0] += bodies[i].vel[0] * dt;
        bodies[i].pos[1] += bodies[i].vel[1] * dt;
        bodies[i].pos[2] += bodies[i].vel[2] * dt;
    }
}

double calculateEnergy(const std::vector<Body>& bodies) {
    double kineticEnergy = 0;
    double potentialEnergy = 0;

    #pragma omp parallel for reduction(+:kineticEnergy, potentialEnergy)
    for (size_t i = 0; i < bodies.size(); ++i) {
        kineticEnergy += 0.5 * bodies[i].mass * (bodies[i].vel[0]*bodies[i].vel[0] + bodies[i].vel[1]*bodies[i].vel[1] + bodies[i].vel[2]*bodies[i].vel[2]);

        for (size_t j = i + 1; j < bodies.size(); ++j) {
            double dx = bodies[j].pos[0] - bodies[i].pos[0];
            double dy = bodies[j].pos[1] - bodies[i].pos[1];
            double dz = bodies[j].pos[2] - bodies[i].pos[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            potentialEnergy -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }

    return kineticEnergy + potentialEnergy;
}