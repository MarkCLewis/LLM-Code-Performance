import java.util.Random;

class Body {
    double[] position;
    double[] velocity;
    double mass;

    public Body(double[] position, double[] velocity, double mass) {
        this.position = position;
        this.velocity = velocity;
        this.mass = mass;
    }
}

public class NBodySimulation {
    private static final double G = 6.67430e-11;
    private static final double DT = 1.0;
    private static final int NUM_BODIES = 1_000_000;
    private static final int STEPS = 1000;

    private static double computeEnergy(Body[] bodies) {
        double energy = 0.0;
        for (int i = 0; i < bodies.length; i++) {
            Body b1 = bodies[i];
            energy += 0.5 * b1.mass * (b1.velocity[0] * b1.velocity[0] +
                                       b1.velocity[1] * b1.velocity[1] +
                                       b1.velocity[2] * b1.velocity[2]);
            for (int j = i + 1; j < bodies.length; j++) {
                Body b2 = bodies[j];
                double dx = b1.position[0] - b2.position[0];
                double dy = b1.position[1] - b2.position[1];
                double dz = b1.position[2] - b2.position[2];
                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                energy -= G * b1.mass * b2.mass / dist;
            }
        }
        return energy;
    }

    private static void updateVelocities(Body[] bodies) {
        for (int i = 0; i < bodies.length; i++) {
            Body b1 = bodies[i];
            double[] force = {0.0, 0.0, 0.0};
            for (int j = 0; j < bodies.length; j++) {
                if (i != j) {
                    Body b2 = bodies[j];
                    double dx = b2.position[0] - b1.position[0];
                    double dy = b2.position[1] - b1.position[1];
                    double dz = b2.position[2] - b1.position[2];
                    double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    double F = G * b1.mass * b2.mass / (dist * dist * dist);
                    force[0] += F * dx;
                    force[1] += F * dy;
                    force[2] += F * dz;
                }
            }
            b1.velocity[0] += force[0] / b1.mass * DT;
            b1.velocity[1] += force[1] / b1.mass * DT;
            b1.velocity[2] += force[2] / b1.mass * DT;
        }
    }

    private static void updatePositions(Body[] bodies) {
        for (Body body : bodies) {
            body.position[0] += body.velocity[0] * DT;
            body.position[1] += body.velocity[1] * DT;
            body.position[2] += body.velocity[2] * DT;
        }
    }

    private static Body[] initializeBodies(int numBodies, double centralMass) {
        Body[] bodies = new Body[numBodies];
        Random random = new Random();
        bodies[0] = new Body(new double[]{0.0, 0.0, 0.0}, new double[]{0.0, 0.0, 0.0}, centralMass);
        double radius = 1.0e9;
        for (int i = 1; i < numBodies; i++) {
            double angle = 2 * Math.PI * random.nextDouble();
            double[] pos = {radius * Math.cos(angle), radius * Math.sin(angle), 0.0};
            double speed = Math.sqrt(G * centralMass / radius);
            double[] vel = {-speed * Math.sin(angle), speed * Math.cos(angle), 0.0};
            bodies[i] = new Body(pos, vel, 1.0);
        }
        return bodies;
    }

    public static void main(String[] args) {
        Body[] bodies = initializeBodies(NUM_BODIES, 1.989e30);

        double initialEnergy = computeEnergy(bodies);
        System.out.printf("Initial Energy: %.6e\n", initialEnergy);

        for (int step = 0; step < STEPS; step++) {
            updateVelocities(bodies);
            updatePositions(bodies);
        }

        double finalEnergy = computeEnergy(bodies);
        System.out.printf("Final Energy: %.6e\n", finalEnergy);
    }
}
