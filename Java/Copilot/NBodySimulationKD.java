import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

public class NBodySimulation {
    private static final double G = 6.67430e-11; // Gravitational constant
    private static final double DT = 1e-3; // Time step
    private static final int NUM_BODIES = 1000000; // Number of bodies
    private static final double THETA = 0.3; // Theta value for approximation

    static class Body {
        double x, y, z;
        double vx, vy, vz;
        double mass;

        Body(double x, double y, double z, double vx, double vy, double vz, double mass) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
            this.mass = mass;
           }

    static class KDNode {
        Body body;
        KDNode left, right;
        double[] min = new double[3];
        double[] max = new double[3];

        KDNode(Body body) {
            this.body = body;
        }
    }

    private static void initializeBodies(List<Body> bodies) {
        bodies.add(new Body(0, 0, 0, 0, 0, 0, 1e30)); // Central body mass

        IntStream.range(1, NUM_BODIES).parallel().forEach(i -> {
            double angle = 2 * Math.PI * i / (NUM_BODIES - 1);
            bodies.add(new Body(
                Math.cos(angle) * 1e11,
                Math.sin(angle) * 1e11,
                0,
                -Math.sin(angle) * Math.sqrt(G * bodies.get(0).mass / 1e11),
                Math.cos(angle) * Math.sqrt(G * bodies.get(0).mass / 1e11),
                0,
                1e24 // Small body mass
            ));
        });
    }

    private static double calculateEnergy(List<Body> bodies) {
        return bodies.parallelStream().mapToDouble(body -> {
            double kinetic = 0.5 * body.mass * (body.vx * body.vx + body.vy * body.vy + body.vz * body.vz);
            double potential = bodies.stream().filter(other -> body != other).mapToDouble(other -> {
                double dx = body.x - other.x;
                double dy = body.y - other.y;
                double dz = body.z - other.z;
                double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                return -G * body.mass * other.mass / distance;
            }).sum();
            return kinetic + 0.5 * potential;
        }).sum();
    }

    private static KDNode buildKDTree(List<Body> bodies, int depth) {
        if (bodies.isEmpty()) return null;

        int axis = depth % 3;
        bodies.sort((a, b) -> {
            switch (axis) {
                case 0: return Double.compare(a.x, b.x);
                case 1: return Double.compare(a.y, b.y);
                default: return Double.compare(a.z, b.z);
            }
        });

        int median = bodies.size() / 2;
        KDNode node = new KDNode(bodies.get(median));
        node.left = buildKDTree(bodies.subList(0, median), depth + 1);
        node.right = buildKDTree(bodies.subList(median + 1, bodies.size()), depth + 1);

        for (int i = 0; i < 3; i++) {
            node.min[i] = node.max[i] = node.body.x;
            if (node.left != null) {
                node.min[i] = Math.min(node.min[i], node.left.min[i]);
                node.max[i] = Math.max(node.max[i], node.left.max[i]);
            }
            if (node.right != null) {
                node.min[i] = Math.min(node.min[i], node.right.min[i]);
                node.max[i] = Math.max(node.max[i], node.right.max[i]);
            }
        }

        return node;
    }

    private static void calculateForce(KDNode node, Body body, double[] acceleration) {
        if (node == null || node.body == null) return;

        double dx = node.body.x - body.x;
        double dy = node.body.y - body.y;
        double dz = node.body.z - body.z;
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        double size = Math.max(node.max[0] - node.min[0], Math.max(node.max[1] - node.min[1], node.max[2] - node.min[2]));

        if (size / distance < THETA || (node.left == null && node.right == null)) {
            double force = G * node.body.mass / (distance * distance * distance);
            acceleration[0] += force * dx;
            acceleration[1] += force * dy;
            acceleration[2] += force * dz;
        } else {
            calculateForce(node.left, body, acceleration);
            calculateForce(node.right, body, acceleration);
        }
    }

    private static void kickStep(List<Body> bodies, KDNode root) {
        ForkJoinPool.commonPool().invoke(new RecursiveTask<Void>() {
            @Override
            protected Void compute() {
                bodies.parallelStream().forEach(body -> {
                    double[] acceleration = new double[3];
                    calculateForce(root, body, acceleration);
                    body.vx += acceleration[0] * DT;
                    body.vy += acceleration[1] * DT;
                    body.vz += acceleration[2] * DT;
                });

                bodies.parallelStream().forEach(body -> {
                    body.x += body.vx * DT;
                    body.y += body.vy * DT;
                    body.z += body.vz * DT;
                });

                return null;
            }
        });
    }

    public static void main(String[] args) {
        List<Body> bodies = new ArrayList<>(NUM_BODIES);
        initializeBodies(bodies);

        double initialEnergy = calculateEnergy(bodies);
        System.out.printf("Initial energy: %e%n", initialEnergy);

        for (int step = 0; step < 1000; step++) {
            KDNode root = buildKDTree(bodies, 0);
            kickStep(bodies, root);
        }

        double finalEnergy = calculateEnergy(bodies);
        System.out.printf("Final energy: %e%n", finalEnergy);
    }
}
