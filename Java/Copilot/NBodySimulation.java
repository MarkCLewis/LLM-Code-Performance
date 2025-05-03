import java.util.Random;

class Body {
    double mass;
    double[] pos = new double[3];
    double[] vel = new double[3];
    double[] force = new double[3];

    public Body(double mass, double[] pos, double[] vel) {
        this.mass = mass;
        this.pos = pos;
        this.vel = vel;
    }
}

public class NBodySimulation {
    private static final double G = 6.67430e-11; // Gravitational constant
    private static final double dt = 1e-3;       // Time step
    private static final int numBodies = 1000000; // Number of small bodies

    public static void main(String[] args) {
        Random rand = new Random();

        // Initialize bodies
        Body[] bodies = initializeBodies(numBodies, rand);

        // Calculate initial energy
        double initialEnergy = calculateEnergy(bodies);
        System.out.printf("Initial Energy: %e%n", initialEnergy);

        // Run simulation for 1000 steps
        for (int step = 0; step < 1000; step++) {
            simulateStep(bodies);
        }

        // Calculate final energy
        double finalEnergy = calculateEnergy(bodies);
        System.out.printf("Final Energy: %e%n", finalEnergy);
    }

    private static Body[] initializeBodies(int n, Random rand) {
        Body[] bodies = new Body[n + 1];

        // Central body
        bodies[0] = new Body(1e30, new double[]{0, 0, 0}, new double[]{0, 0, 0});

        // Small bodies in circular orbits
        for (int i = 1; i <= n; i++) {
            double radius = 1e11 * rand.nextDouble();
            double speed = Math.sqrt(G * bodies[0].mass / radius);
            double angle = 2 * Math.PI * rand.nextDouble();

            bodies[i] = new Body(1e24,
                    new double[]{radius * Math.cos(angle), radius * Math.sin(angle), 0},
                    new double[]{-speed * Math.sin(angle), speed * Math.cos(angle), 0});
        }

        return bodies;
    }

    private static void simulateStep(Body[] bodies) {
        // Reset forces
        for (Body body : bodies) {
            body.force[0] = body.force[1] = body.force[2] = 0;
        }

        // Calculate forces
        for (int i = 0; i < bodies.length; i++) {
            for (int j = i + 1; j < bodies.length; j++) {
                double dx = bodies[j].pos[0] - bodies[i].pos[0];
                double dy = bodies[j].pos[1] - bodies[i].pos[1];
                double dz = bodies[j].pos[2] - bodies[i].pos[2];
                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                double force = G * bodies[i].mass * bodies[j].mass / (dist * dist * dist);

                bodies[i].force[0] += force * dx;
                bodies[i].force[1] += force * dy;
                bodies[i].force[2] += force * dz;

                bodies[j].force[0] -= force * dx;
                bodies[j].force[1] -= force * dy;
                bodies[j].force[2] -= force * dz;
            }
        }

        // Update velocities and positions
        for (Body body : bodies) {
            body.vel[0] += body.force[0] / body.mass * dt;
            body.vel[1] += body.force[1] / body.mass * dt;
            body.vel[2] += body.force[2] / body.mass * dt;

            body.pos[0] += body.vel[0] * dt;
            body.pos[1] += body.vel[1] * dt;
            body.pos[2] += body.vel[2] * dt;
        }
    }

    private static double calculateEnergy(Body[] bodies) {
        double kineticEnergy = 0;
        double potentialEnergy = 0;

        for (int i = 0; i < bodies.length; i++) {
            kineticEnergy += 0.5 * bodies[i].mass * (bodies[i].vel[0] * bodies[i].vel[0] + bodies[i].vel[1] * bodies[i].vel[1] + bodies[i].vel[2] * bodies[i].vel[2]);

            for (int j = i + 1; j < bodies.length; j++) {
                double dx = bodies[j].pos[0] - bodies[i].pos[0];
                double dy = bodies[j].pos[1] - bodies[i].pos[1];
                double dz = bodies[j].pos[2] - bodies[i].pos[2];
                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                potentialEnergy -= G * bodies[i].mass * bodies[j].mass / dist;
            }
        }

        return kineticEnergy + potentialEnergy;
    }
}