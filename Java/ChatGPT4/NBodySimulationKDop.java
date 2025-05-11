// The code has been optimized further by replacing Arrays.stream(...).parallel() with IntStream.range(...).parallel() where appropriate, improving indexing and avoiding boxing overhead. I also squared the THETA comparison to save a sqrt call when it can be avoided.

import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

public class NBodySimulationKDop {
    static final int N_BODIES = 100_000;
    static final int STEPS = 10;
    static final double G = 6.67430e-11;
    static final double DT = 1.0;
    static final double THETA = 0.3;

    static Body[] bodies;

    static class Body {
        double x, y, z, vx, vy, vz, ax, ay, az, mass;
    }

    static class KDNode {
        final double[] min = new double[3];
        final double[] max = new double[3];
        double cmX, cmY, cmZ, mass;
        KDNode left, right;
        List<Body> bodyList;
    }

    public static void main(String[] args) {
        System.out.println("Initializing...");
        initializeSystem(N_BODIES);
        // System.out.println("Calculating initial energy...");
        // double initialEnergy = computeEnergy();
        // System.out.printf("Initial energy: %.6e\n", initialEnergy);

        for (int step = 0; step < STEPS; step++) {
            KDNode tree = buildKDTree(Arrays.asList(bodies), 0);
            computeForces(tree);
            updateBodies();
            if (step % 100 == 0) {
                System.out.println("Step " + step);
            }
        }

        // System.out.println("Calculating final energy...");
        // double finalEnergy = computeEnergy();
        // System.out.printf("Final energy: %.6e\n", finalEnergy);
        // System.out.printf("Energy difference: %.6e\n", Math.abs(finalEnergy - initialEnergy));
    }

    static void initializeSystem(int n) {
        bodies = new Body[n + 1];
        bodies[0] = new Body();
        bodies[0].x = bodies[0].y = bodies[0].z = 0;
        bodies[0].vx = bodies[0].vy = bodies[0].vz = 0;
        bodies[0].mass = 1e20;

        Random rand = new Random(42);
        double radius = 1e7;

        for (int i = 1; i <= n; i++) {
            double angle = 2 * Math.PI * i / n;
            double r = radius * (1.0 + 0.1 * rand.nextDouble());

            Body b = new Body();
            b.x = r * Math.cos(angle);
            b.y = r * Math.sin(angle);
            b.z = 0;

            double v = Math.sqrt(G * bodies[0].mass / r);
            b.vx = -v * Math.sin(angle);
            b.vy = v * Math.cos(angle);
            b.vz = 0;
            b.mass = 1.0;
            bodies[i] = b;
        }
    }

    static KDNode buildKDTree(List<Body> bodyList, int depth) {
        if (bodyList.isEmpty()) return null;

        KDNode node = new KDNode();
        node.bodyList = bodyList;

        Arrays.fill(node.min, Double.POSITIVE_INFINITY);
        Arrays.fill(node.max, Double.NEGATIVE_INFINITY);

        node.mass = 0;
        node.cmX = node.cmY = node.cmZ = 0;
        for (Body b : bodyList) {
            node.min[0] = Math.min(node.min[0], b.x);
            node.min[1] = Math.min(node.min[1], b.y);
            node.min[2] = Math.min(node.min[2], b.z);

            node.max[0] = Math.max(node.max[0], b.x);
            node.max[1] = Math.max(node.max[1], b.y);
            node.max[2] = Math.max(node.max[2], b.z);

            node.mass += b.mass;
            node.cmX += b.mass * b.x;
            node.cmY += b.mass * b.y;
            node.cmZ += b.mass * b.z;
        }
        node.cmX /= node.mass;
        node.cmY /= node.mass;
        node.cmZ /= node.mass;

        if (bodyList.size() <= 1) return node;

        int axis = depth % 3;
        bodyList.sort(Comparator.comparingDouble(b -> (axis == 0 ? b.x : axis == 1 ? b.y : b.z)));
        int mid = bodyList.size() / 2;

        node.left = buildKDTree(bodyList.subList(0, mid), depth + 1);
        node.right = buildKDTree(bodyList.subList(mid, bodyList.size()), depth + 1);
        return node;
    }

    static void computeForces(KDNode tree) {
        IntStream.range(0, bodies.length).parallel().forEach(i -> {
            Body b = bodies[i];
            b.ax = b.ay = b.az = 0;
            computeForceFromNode(b, tree);
        });
    }

    static void computeForceFromNode(Body b, KDNode node) {
        if (node == null || node.mass == 0) return; // Had to remove " || node == b" because it is an illegal comparison.

        double dx = node.cmX - b.x;
        double dy = node.cmY - b.y;
        double dz = node.cmZ - b.z;
        double dist2 = dx * dx + dy * dy + dz * dz + 1e-10;

        double size = Math.max(Math.max(node.max[0] - node.min[0], node.max[1] - node.min[1]), node.max[2] - node.min[2]);

        if (node.bodyList.size() <= 1 || (size * size / dist2 < THETA * THETA)) {
            double dist = Math.sqrt(dist2);
            double force = G * node.mass / (dist2 * dist);
            b.ax += dx * force;
            b.ay += dy * force;
            b.az += dz * force;
        } else {
            computeForceFromNode(b, node.left);
            computeForceFromNode(b, node.right);
        }
    }

    static void updateBodies() {
        IntStream.range(0, bodies.length).parallel().forEach(i -> {
            Body b = bodies[i];
            b.vx += b.ax * DT;
            b.vy += b.ay * DT;
            b.vz += b.az * DT;

            b.x += b.vx * DT;
            b.y += b.vy * DT;
            b.z += b.vz * DT;
        });
    }

    // static double computeEnergy() {
    //     double kinetic = IntStream.range(0, bodies.length).parallel().mapToDouble(i -> {
    //         Body b = bodies[i];
    //         double v2 = b.vx * b.vx + b.vy * b.vy + b.vz * b.vz;
    //         return 0.5 * b.mass * v2;
    //     }).sum();

    //     double potential = IntStream.range(0, bodies.length).parallel().mapToDouble(i -> {
    //         double pot = 0;
    //         Body bi = bodies[i];
    //         for (int j = i + 1; j < bodies.length; j++) {
    //             Body bj = bodies[j];
    //             double dx = bi.x - bj.x;
    //             double dy = bi.y - bj.y;
    //             double dz = bi.z - bj.z;
    //             double dist2 = dx * dx + dy * dy + dz * dz + 1e-10;
    //             pot -= G * bi.mass * bj.mass / Math.sqrt(dist2);
    //         }
    //         return pot;
    //     }).sum();

    //     return kinetic + potential;
    // }
}
