import java.util.Arrays;
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

class Node {
    Body body;
    double[] center = new double[3];
    double size;
    double mass;
    double[] com = new double[3];
    Node[] children = new Node[8];

    public Node() {
        this.mass = 0;
    }
}

public class NBodySimulation {
    private static final double G = 6.67430e-11; // Gravitational constant
    private static final double dt = 1e-3;       // Time step
    private static final int numBodies = 1000000; // Number of small bodies
    private static final double theta = 0.3;     // Barnes-Hut approximation parameter

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
        Arrays.stream(bodies).parallel().forEach(body -> {
            body.force[0] = body.force[1] = body.force[2] = 0;
        });

        // Build kD-tree
        Node root = buildTree(bodies);

        // Calculate forces using kD-tree
        Arrays.stream(bodies).parallel().forEach(body -> {
            calculateForce(body, root);
        });

        // Update velocities and positions
        Arrays.stream(bodies).parallel().forEach(body -> {
            body.vel[0] += body.force[0] / body.mass * dt;
            body.vel[1] += body.force[1] / body.mass * dt;
            body.vel[2] += body.force[2] / body.mass * dt;

            body.pos[0] += body.vel[0] * dt;
            body.pos[1] += body.vel[1] * dt;
            body.pos[2] += body.vel[2] * dt;
        });
    }

    private static Node buildTree(Body[] bodies) {
        Node root = new Node();
        root.center[0] = root.center[1] = root.center[2] = 0;
        root.size = 2 * 1e11;

        for (Body body : bodies) {
            insertBody(root, body);
        }

        return root;
    }

    private static void insertBody(Node node, Body body) {
        if (node.body == null && node.mass == 0) {
            node.body = body;
            node.mass = body.mass;
            System.arraycopy(body.pos, 0, node.com, 0, 3);
            return;
        }

        if (node.body != null) {
            Body oldBody = node.body;
            node.body = null;
            insertBody(node, oldBody);
        }

        node.mass += body.mass;
        for (int i = 0; i < 3; i++) {
            node.com[i] = (node.com[i] * (node.mass - body.mass) + body.pos[i] * body.mass) / node.mass;
        }

        int index = 0;
        for (int i = 0; i < 3; i++) {
            if (body.pos[i] > node.center[i]) {
                index |= 1 << i;
            }
        }

        if (node.children[index] == null) {
            node.children[index] = new Node();
            for (int i = 0; i < 3; i++) {
                node.children[index].center[i] = node.center[i] + (index & (1 << i) != 0 ? 1 : -1) * node.size / 4;
            }
            node.children[index].size = node.size / 2;
        }

        insertBody(node.children[index], body);
    }

    private static void calculateForce(Body body, Node node) {
        if (node.body != null && node.body != body) {
            double dx = node.body.pos[0] - body.pos[0];
            double dy = node.body.pos[1] - body.pos[1];
            double dz = node.body.pos[2] - body.pos[2];
            double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            double force = G * body.mass * node.body.mass / (dist * dist * dist);

            body.force[0] += force * dx;
            body.force[1] += force * dy;
            body.force[2] += force * dz;
            return;
        }

        double dx = node.com[0] - body.pos[0];
        double dy = node.com[1] - body.pos[1];
        double dz = node.com[2] - body.pos[2];
        double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (node.size / dist < theta) {
            double force = G * body.mass * node.mass / (dist * dist * dist);

            body.force[0] += force * dx;
            body.force[1] += force * dy;
            body.force[2] += force * dz;
            return;
        }

        for (Node child : node.children) {
            if (child != null) {
                calculateForce(body, child);
            }
        }
    }

    private static double calculateEnergy(Body[] bodies) {
        return Arrays.stream(bodies).parallel().mapToDouble(body1 -> {
            double kineticEnergy = 0.5 * body1.mass * (body1.vel[0] * body1.vel[0] + body1.vel[1] * body1.vel[1] + body1.vel[2] * body1.vel[2]);
            double potentialEnergy = Arrays.stream(bodies).parallel().filter(body2 -> body1 != body2).mapToDouble(body2 -> {
                double dx = body2.pos[0] - body1.pos[0];
                double dy = body2.pos[1] - body1.pos[1];
                double dz = body2.pos[2] - body1.pos[2];
                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                return -G * body1.mass * body2.mass / dist;
            }).sum();
            return kineticEnergy + potentialEnergy;
        }).sum();
    }
}