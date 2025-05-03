#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

const double G = 6.67430e-11;
const double DT = 1.0;
const double PI = 3.141592653589793;

struct Body {
    double position[3];
    double velocity[3];
    double mass;
};

void compute_energy(const std::vector<Body>& bodies, double& kinetic, double& potential) {
    kinetic = 0.0;
    potential = 0.0;
    
    for (size_t i = 0; i < bodies.size(); ++i) {
        double speed_sq = bodies[i].velocity[0] * bodies[i].velocity[0] +
                          bodies[i].velocity[1] * bodies[i].velocity[1] +
                          bodies[i].velocity[2] * bodies[i].velocity[2];
        kinetic += 0.5 * bodies[i].mass * speed_sq;

        for (size_t j = i + 1; j < bodies.size(); ++j) {
            double dx = bodies[i].position[0] - bodies[j].position[0];
            double dy = bodies[i].position[1] - bodies[j].position[1];
            double dz = bodies[i].position[2] - bodies[j].position[2];
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            if (r > 0.0)
                potential -= G * bodies[i].mass * bodies[j].mass / r;
        }
    }
}

void update_positions(std::vector<Body>& bodies) {
    for (auto& body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.position[i] += body.velocity[i] * DT;
        }
    }
}

void update_velocities(std::vector<Body>& bodies) {
    size_t n = bodies.size();
    std::vector<std::array<double, 3>> accelerations(n, {0.0, 0.0, 0.0});
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                double dx = bodies[j].position[0] - bodies[i].position[0];
                double dy = bodies[j].position[1] - bodies[i].position[1];
                double dz = bodies[j].position[2] - bodies[i].position[2];
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                if (r > 0.0) {
                    double factor = G * bodies[j].mass / (r * r * r);
                    accelerations[i][0] += factor * dx;
                    accelerations[i][1] += factor * dy;
                    accelerations[i][2] += factor * dz;
                }
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        for (int k = 0; k < 3; ++k) {
            bodies[i].velocity[k] += accelerations[i][k] * DT;
        }
    }
}

std::vector<Body> initialize_orbiting_bodies(size_t num_bodies, double central_mass) {
    std::vector<Body> bodies;
    bodies.reserve(num_bodies + 1);
    
    bodies.push_back({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, central_mass});
    double radius = 1.0e9;

    for (size_t i = 0; i < num_bodies; ++i) {
        double angle = 2.0 * PI * (static_cast<double>(i) / num_bodies);
        double pos_x = radius * cos(angle);
        double pos_y = radius * sin(angle);
        double speed = sqrt(G * central_mass / radius);
        double vel_x = -speed * sin(angle);
        double vel_y = speed * cos(angle);
        bodies.push_back({{pos_x, pos_y, 0.0}, {vel_x, vel_y, 0.0}, 1.0});
    }
    return bodies;
}

int main() {
    size_t n_bodies = 1'000'000;
    double central_mass = 1.989e30;
    
    auto bodies = initialize_orbiting_bodies(n_bodies, central_mass);
    
    double initial_kinetic, initial_potential;
    compute_energy(bodies, initial_kinetic, initial_potential);
    std::cout << "Initial Energy: " << initial_kinetic + initial_potential << "\n";
    
    for (int step = 0; step < 1000; ++step) {
        update_velocities(bodies);
        update_positions(bodies);
    }
    
    double final_kinetic, final_potential;
    compute_energy(bodies, final_kinetic, final_potential);
    std::cout << "Final Energy: " << final_kinetic + final_potential << "\n";
    
    return 0;
}
