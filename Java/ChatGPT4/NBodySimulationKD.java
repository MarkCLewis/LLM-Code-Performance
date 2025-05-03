import java.util.*;
import java.util.concurrent.*;

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

class KDTree {
    double[] centerOfMass;
    double totalMass;
    double[] boundaryMin;
    double[] boundaryMax;
    List<KDTree> children;
    Body body;

    public KDTree(List<Body> bodies) {
        if (bodies.isEmpty()) return;
        if (bodies.size() == 1) {
            this.body = bodies.get(0);
            this.centerOfMass = Arrays.copyOf(body.position, 3);
            this.totalMass = body.mass;
            this.boundaryMin = Arrays.copyOf(body.position, 3);
            this.boundaryMax = Arrays.copyOf(body.position, 3);
        } else {
            this.totalMass = bodies.stream().mapToDouble(b -> b.mass).sum();
            this.centerOfMass = new double[]{0, 0, 0};
            for (Body b : bodies) {
                for (int i = 0; i < 3; i++) {
                    this.centerOfMass[i] += b.position[i] * b.mass;
                }
            }
            for (int i = 0; i < 3; i++) {
                this.centerOfMass[i] /= this.totalMass;
            }
            this.boundaryMin = new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
            this.boundaryMax = new double[]{Double.MIN_VALUE, Double.MIN_VALUE, Double.MIN_VALUE};
            for (Body b : bodies) {
                for (int i = 0; i < 3; i++) {
                    this.boundaryMin[i] = Math.min(this.boundaryMin[i], b.position[i]);
                    this.boundaryMax[i] = Math.max(this.boundaryMax[i], b.position[i]);
                }
            }
            this.children = new ArrayList<>();
            for (Body b : bodies) {
                this.children.add(new KDTree(Collections.singletonList(b)));
            }
        }
    }
}

public class NBodySimulation {
    private static final double G = 6.67430e-11;
    private static final double DT = 1.0;
    private static final int NUM_BODIES = 1_000_000;
    private static final int STEPS = 1000;
    private static final double THETA = 0.3;

    private static double computeEnergy(List<Body> bodies) {
        return bodies.parallelStream().mapToDouble(b1 -> {
            double energy = 0.5 * b1.mass * (b1.velocity[0] * b1.velocity[0] +
                                             b1.velocity[1] * b1.velocity[1] +
                                             b1.velocity[2] * b1.velocity[2]);
            return energy;
        }).sum();
    }

    private static double[] computeForce(Body body, KDTree tree) {
        if (tree == null || tree.body == body) return new double[]{0.0, 0.0, 0.0};
        double dx = tree.centerOfMass[0] - body.position[0];
        double dy = tree.centerOfMass[1] - body.position[1];
        double dz = tree.centerOfMass[2] - body.position[2];
        double r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        double s = Math.max(tree.boundaryMax[0] - tree.boundaryMin[0],
                            Math.max(tree.boundaryMax[1] - tree.boundaryMin[1],
                                     tree.boundaryMax[2] - tree.boundaryMin[2]));
        if (s / r < THETA || tree.body != null) {
            double F = G * tree.totalMass / (r * r * r);
            return new double[]{F * dx, F * dy, F * dz};
        } else {
            double[] force = {0.0, 0.0, 0.0};
            for (KDTree child : tree.children) {
                double[] childForce = computeForce(body, child);
                force[0] += childForce[0];
                force[1] += childForce[1];
                force[2] += childForce[2];
            }
            return force;
        }
    }

    private static void updateVelocities(List<Body> bodies, KDTree tree) {
        bodies.parallelStream().forEach(body -> {
            double[] force = computeForce(body, tree);
            body.velocity[0] += force[0] / body.mass * DT;
            body.velocity[1] += force[1] / body.mass * DT;
            body.velocity[2] += force[2] / body.mass * DT;
        });
    }

    private static void updatePositions(List<Body> bodies) {
        bodies.parallelStream().forEach(body -> {
            body.position[0] += body.velocity[0] * DT;
            body.position[1] += body.velocity[1] * DT;
            body.position[2] += body.velocity[2] * DT;
        });
    }

    private static List<Body> initializeBodies(int numBodies, double centralMass) {
        List<Body> bodies = new ArrayList<>();
        Random random = new Random();
        bodies.add(new Body(new double[]{0.0, 0.0, 0.0}, new double[]{0.0, 0.0, 0.0}, centralMass));
        double radius = 1.0e9;
        for (int i = 1; i < numBodies; i++) {
            double angle = 2 * Math.PI * random.nextDouble();
            double[] pos = {radius * Math.cos(angle), radius * Math.sin(angle), 0.0};
            double speed = Math.sqrt(G * centralMass / radius);
            double[] vel = {-speed * Math.sin(angle), speed * Math.cos(angle), 0.0};
            bodies.add(new Body(pos, vel, 1.0));
        }
        return bodies;
    }

    public static void main(String[] args) {
        List<Body> bodies = initializeBodies(NUM_BODIES, 1.989e30);
        double initialEnergy = computeEnergy(bodies);
        System.out.printf("Initial Energy: %.6e\n", initialEnergy);

        for (int step = 0; step < STEPS; step++) {
            KDTree tree = new KDTree(bodies);
            updateVelocities(bodies, tree);
            updatePositions(bodies);
        }

        double finalEnergy = computeEnergy(bodies);
        System.out.printf("Final Energy: %.6e\n", finalEnergy);
    }
}