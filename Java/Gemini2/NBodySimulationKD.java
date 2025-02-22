import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public class NBodySimulationKD {

    private static final double GRAVITATIONAL_CONSTANT = 6.674e-11;
    private static final double TIME_STEP = 0.001; // Adjust as needed
    private static final int NUM_STEPS = 1000;

    private Body[] bodies;
    private KdTree kdTree;
    private double theta = 0.3; // Barnes-Hut theta value

    public NBodySimulationKD(int numBodies) {
        bodies = new Body[numBodies + 1];
        initializeSystem();
        kdTree = new KdTree(bodies); // Build the kD-tree
    }

    private void initializeSystem() {
        Random random = new Random(0); // Seed for reproducibility

        // Initialize central body (Sun-like)
        bodies[0] = new Body(1.989e30, 0, 0, 0, 0, 0, 0); // Mass, x, y, z, vx, vy, vz

        // Initialize smaller bodies (planet-like) in circular orbits
        for (int i = 1; i < bodies.length; i++) {
            double mass = 1e25 + random.nextDouble() * 1e26; // Varying mass
            double radius = 1e11 + random.nextDouble() * 1e12; // Varying orbital radius
            double speed = Math.sqrt(GRAVITATIONAL_CONSTANT * bodies[0].mass / radius);

            double angle = random.nextDouble() * 2 * Math.PI;
            double x = radius * Math.cos(angle);
            double y = radius * Math.sin(angle);
            double z = 0; // Orbit in x-y plane for simplicity

            double vx = -speed * Math.sin(angle);
            double vy = speed * Math.cos(angle);
            double vz = 0;

            bodies[i] = new Body(mass, x, y, z, vx, vy, vz);
        }
    }

    public void simulate() {
        double initialEnergy = calculateTotalEnergy();
        System.out.println("Initial Energy: " + initialEnergy);

        for (int step = 0; step < NUM_STEPS; step++) {
            // Rebuild the kD-tree each step (important!)
            kdTree = new KdTree(bodies);

            // Kick step (parallelized with kD-tree)
            IntStream.range(0, bodies.length).parallel().forEach(i -> {
                double fx = 0, fy = 0, fz = 0;
                Force force = kdTree.calculateForce(bodies[i], theta); // Use kD-tree
                fx = force.fx;
                fy = force.fy;
                fz = force.fz;


                bodies[i].vx += fx / bodies[i].mass * TIME_STEP;
                bodies[i].vy += fy / bodies[i].mass * TIME_STEP;
                bodies[i].vz += fz / bodies[i].mass * TIME_STEP;
            });

            // Drift step (parallelized)
            IntStream.range(0, bodies.length).parallel().forEach(i -> {
                bodies[i].x += bodies[i].vx * TIME_STEP;
                bodies[i].y += bodies[i].vy * TIME_STEP;
                bodies[i].z += bodies[i].vz * TIME_STEP;
            });
        }

        double finalEnergy = calculateTotalEnergy();
        System.out.println("Final Energy: " + finalEnergy);
        System.out.println("Energy Change: " + (finalEnergy - initialEnergy));
    }

    private double calculateTotalEnergy() {
        double kineticEnergy = Arrays.stream(bodies).parallel()
                .mapToDouble(body -> 0.5 * body.mass * (body.vx * body.vx + body.vy * body.vy + body.vz * body.vz))
                .sum();


        double potentialEnergy = IntStream.range(0, bodies.length).parallel()
                .mapToDouble(i -> {
                    double pe = 0;
                    for (int j = i + 1; j < bodies.length; j++) {
                        double dx = bodies[j].x - bodies[i].x;
                        double dy = bodies[j].y - bodies[i].y;
                        double dz = bodies[j].z - bodies[i].z;
                        double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                        pe -= GRAVITATIONAL_CONSTANT * bodies[i].mass * bodies[j].mass / dist;
                    }
                    return pe;
                }).sum();

        return kineticEnergy + potentialEnergy;
    }

    // KdTree implementation (simplified for 3D)
    private static class KdTree {
        private Node root;

        public KdTree(Body[] bodies) {
            root = buildTree(bodies, 0);
        }

        private Node buildTree(Body[] bodies, int depth) {
            if (bodies.length == 0) return null;

            int axis = depth % 3; // Cycle through x, y, z
            Arrays.sort(bodies, (a, b) -> Double.compare(a.getCoord(axis), b.getCoord(axis)));
            int medianIndex = bodies.length / 2;
            Body medianBody = bodies[medianIndex];

            Node node = new Node(medianBody);
            node.left = buildTree(Arrays.copyOfRange(bodies, 0, medianIndex), depth + 1);
            node.right = buildTree(Arrays.copyOfRange(bodies, medianIndex + 1, bodies.length), depth + 1);
            return node;
        }

        public Force calculateForce(Body body, double theta) {
            Force force = new Force();
            calculateForceRecursive(root, body, force, theta);
            return force;
        }

        private void calculateForceRecursive(Node node, Body body, Force force, double theta) {
            if (node == null) return;

            double dx = node.body.x - body.x;
            double dy = node.body.y - body.y;
            double dz = node.body.z - body.z;
            double distSq = dx * dx + dy * dy + dz * dz;

            if (distSq == 0) return; // Avoid division by zero


            double dist = Math.sqrt(distSq);
            if (isFarEnough(node, body, dist, theta)) {
                double forceMagnitude = GRAVITATIONAL_CONSTANT * body.mass * node.body.mass / distSq;
                force.fx += forceMagnitude * dx / dist;
                force.fy += forceMagnitude * dy / dist;
                force.fz += forceMagnitude * dz / dist;
            } else if (node.left == null && node.right == null) { // Leaf node
                double forceMagnitude = GRAVITATIONAL_CONSTANT * body.mass * node.body.mass / distSq;
                force.fx += forceMagnitude * dx / dist;
                force.fy += forceMagnitude * dy / dist;
                force.fz += forceMagnitude * dz / dist;
            } else {
                calculateForceRecursive(node.left, body, force, theta);
                calculateForceRecursive(node.right, body, force, theta);
            }
        }



        private boolean isFarEnough(Node node, Body body, double dist, double theta) {
            double s = node.size(); // Size of the region represented by the node
            double d = dist;       // Distance from the body to the center of mass
            return (s / d) < theta;
        }

        private static class Node {
            Body body;
            Node left, right;

            public Node(Body body) {
                this.body = body;
            }

            public double size() {
                // Simplified size calculation (bounding box would be more accurate)
                return 1e12; // Example size - replace with actual calculation
            }

        }
    }

    private static class Force {
        double fx, fy, fz;
    }

    private static class Body {
        double mass, x, y, z, vx, vy, vz;

        public Body(double mass, double x, double y, double z, double vx, double vy, double vz) {
            this.mass = mass;
            this.x = x;
            this.y = y;
            this.z = z;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
        }

        // Gemini-2 did not provide this, so there was a compiler error. I had to add it.
        public double getCoord(int axis) {
            switch(axis) {
                case 0: return x;
                case 1: return y;
                case 2: return z;
            }
            throw new IllegalArgumentException("Coordinate "+axis+" is out of range.");
        }
    }

    public static void main(String[] args) {
        int numSmallBodies = 100000;  // One million small bodies
        NBodySimulation simulation = new NBodySimulation(numSmallBodies);
        simulation.simulate();
    }
}
