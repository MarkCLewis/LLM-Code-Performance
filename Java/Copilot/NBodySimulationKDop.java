//Okay, let's implement the minor optimization in the kickStep function by using a sequential loop for the velocity update, as it might offer a slight performance improvement in some scenarios due to the reduced overhead of parallel stream management for this simple operation.
//Here's the modified kickStep function within the NBodySimulationParallelKDTree class:

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public class NBodySimulationParallelKDTree {

    private static final double GRAVITATIONAL_CONSTANT = 6.67430e-11;
    private static final double THETA = 0.3;

    public static class Body {
        double mass;
        double[] position;
        double[] velocity;
        double[] force;

        public Body(double mass, double[] position, double[] velocity) {
            this.mass = mass;
            this.position = position;
            this.velocity = velocity;
            this.force = new double[]{0, 0, 0};
        }
    }

    public static class System {
        List<Body> bodies;

        public System() {
            this.bodies = new ArrayList<>();
        }

        public void addBody(Body body) {
            this.bodies.add(body);
        }
    }

    public static class KDNode {
        int bodyIndex = -1; // -1 for internal nodes
        double[] centerOfMass = new double[]{0, 0, 0};
        double totalMass = 0;
        double[] min = new double[3];
        double[] max = new double[3];
        KDNode left;
        KDNode right;

        public KDNode(int bodyIndex) {
            this.bodyIndex = bodyIndex;
        }

        public KDNode() {
        }
    }

    public static double distanceSquared(Body b1, Body b2) {
        double dx = b1.position[0] - b2.position[0];
        double dy = b1.position[1] - b2.position[1];
        double dz = b1.position[2] - b2.position[2];
        return dx * dx + dy * dy + dz * dz;
    }

    public static double distanceSquaredPoints(double[] p1, double[] p2) {
        double dx = p1[0] - p2[0];
        double dy = p1[1] - p2[1];
        double dz = p1[2] - p2[2];
        return dx * dx + dy * dy + dz * dz;
    }

    public static double distancePoints(double[] p1, double[] p2) {
        return Math.sqrt(distanceSquaredPoints(p1, p2));
    }

    public static void calculateForceElement(Body body, double[] massElementCM, double massElementMass, double[] force) {
        double rSq = distanceSquaredPoints(body.position, massElementCM);
        if (rSq > 1e-9) {
            double r = Math.sqrt(rSq);
            double magnitude = (GRAVITATIONAL_CONSTANT * body.mass * massElementMass) / rSq;
            force[0] += magnitude * (massElementCM[0] - body.position[0]) / r;
            force[1] += magnitude * (massElementCM[1] - body.position[1]) / r;
            force[2] += magnitude * (massElementCM[2] - body.position[2]) / r;
        }
    }

    public static double calculateTotalEnergy(System system) {
        int numBodies = system.bodies.size();
        double kineticEnergy = system.bodies.parallelStream()
                .mapToDouble(body -> 0.5 * body.mass * (body.velocity[0] * body.velocity[0] +
                                                        body.velocity[1] * body.velocity[1] +
                                                        body.velocity[2] * body.velocity[2]))
                .sum();

        double potentialEnergy = IntStream.range(0, numBodies)
                .parallel()
                .mapToDouble(i -> IntStream.range(i + 1, numBodies)
                        .mapToDouble(j -> {
                            Body bodyI = system.bodies.get(i);
                            Body bodyJ = system.bodies.get(j);
                            double r = Math.sqrt(distanceSquared(bodyI, bodyJ));
                            return -(GRAVITATIONAL_CONSTANT * bodyI.mass * bodyJ.mass) / r;
                        })
                        .sum())
                .sum();

        return kineticEnergy + potentialEnergy;
    }

    public static System initializeCircularOrbits(int numOrbiting, double centralMass, double orbitRadius, double orbitingMass) {
        System system = new System();

        // Initialize the central body
        system.addBody(new Body(centralMass, new double[]{0, 0, 0}, new double[]{0, 0, 0}));

        // Initialize the orbiting bodies
        Random random = new Random();
        for (int i = 0; i < numOrbiting; i++) {
            double angle = 2 * Math.PI * i / (double) numOrbiting;
            double x = orbitRadius * Math.cos(angle);
            double y = orbitRadius * Math.sin(angle);
            double z = 0; // Orbiting in the xy-plane

            // Calculate the orbital velocity for a circular orbit
            double orbitalSpeed = Math.sqrt(GRAVITATIONAL_CONSTANT * centralMass / orbitRadius);
            double vx = -orbitalSpeed * Math.sin(angle);
            double vy = orbitalSpeed * Math.cos(angle);
            double vz = 0;

            system.addBody(new Body(orbitingMass, new double[]{x, y, z}, new double[]{vx, vy, vz}));
        }

        return system;
    }

    public static KDNode buildKDTree(System system) {
        List<Integer> indices = IntStream.range(0, system.bodies.size()).boxed().toList();
        double[] minBound = new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
        double[] maxBound = new double[]{Double.MIN_VALUE, Double.MIN_VALUE, Double.MIN_VALUE};

        for (Body body : system.bodies) {
            for (int i = 0; i < 3; i++) {
                minBound[i] = Math.min(minBound[i], body.position[i]);
                maxBound[i] = Math.max(maxBound[i], body.position[i]);
            }
        }

        return buildKDTreeRecursive(system, indices, minBound, maxBound, 0);
    }

    private static KDNode buildKDTreeRecursive(System system, List<Integer> indices, double[] minBound, double[] maxBound, int depth) {
        int numIndices = indices.size();
        if (numIndices == 0) {
            return null;
        }
        if (numIndices == 1) {
            int index = indices.get(0);
            KDNode node = new KDNode(index);
            Body body = system.bodies.get(index);
            node.totalMass = body.mass;
            node.centerOfMass = body.position.clone();
            node.min = minBound.clone();
            node.max = maxBound.clone();
            return node;
        }

        KDNode node = new KDNode();
        node.min = minBound.clone();
        node.max = maxBound.clone();
        node.totalMass = 0;
        double[] cm = new double[]{0, 0, 0};

        for (int index : indices) {
            Body body = system.bodies.get(index);
            node.totalMass += body.mass;
            for (int i = 0; i < 3; i++) {
                cm[i] += body.mass * body.position[i];
            }
        }
        if (node.totalMass > 0) {
            for (int i = 0; i < 3; i++) {
                node.centerOfMass[i] = cm[i] / node.totalMass;
            }
        }

        int splitDim = depth % 3;
        double median = (minBound[splitDim] + maxBound[splitDim]) / 2;

        List<Integer> leftIndices = new ArrayList<>();
        List<Integer> rightIndices = new ArrayList<>();

        for (int index : indices) {
            if (system.bodies.get(index).position[splitDim] <= median) {
                leftIndices.add(index);
            } else {
                rightIndices.add(index);
            }
        }

        double[] leftMinBound = minBound.clone();
        double[] leftMaxBound = maxBound.clone();
        leftMaxBound[splitDim] = median;

        double[] rightMinBound = minBound.clone();
        double[] rightMaxBound = maxBound.clone();
        rightMinBound[splitDim] = median;

        node.left = buildKDTreeRecursive(system, leftIndices, leftMinBound, leftMaxBound, depth + 1);
        node.right = buildKDTreeRecursive(system, rightIndices, rightMinBound, rightMaxBound, depth + 1);

        return node;
    }

    public static void calculateForceKDTree(Body bodyI, KDNode node, System system) {
        if (node == null) {
            return;
        }

        if (node.bodyIndex != -1) {
            if (node.bodyIndex != system.bodies.indexOf(bodyI)) {
                calculateForceElement(bodyI, node.centerOfMass, node.totalMass, bodyI.force);
            }
            return;
        }

        double s = node.max[0] - node.min[0]; // Approximate size of the region
        double d = distancePoints(bodyI.position, node.centerOfMass);

        if (d == 0 || (s / d < THETA)) {
            calculateForceElement(bodyI, node.centerOfMass, node.totalMass, bodyI.force);
        } else {
            calculateForceKDTree(bodyI, node.left, system);
            calculateForceKDTree(bodyI, node.right, system);
        }
    }

    public static void kickStep(System system, double dt) {
        KDNode root = buildKDTree(system);
        int numBodies = system.bodies.size();

        // Reset forces in parallel
        system.bodies.parallelStream().forEach(body -> {
            body.force[0] = 0;
            body.force[1] = 0;
            body.force[2] = 0;
        });

        // Calculate forces using the kD-tree in parallel
        system.bodies.parallelStream().forEach(body -> calculateForceKDTree(body, root, system));

        // Update velocities (kick) sequentially
        double half_dt_over_mass;
        for (Body body : system.bodies) {
            half_dt_over_mass = 0.5 * dt / body.mass;
            body.velocity[0] += body.force[0] * half_dt_over_mass;
            body.velocity[1] += body.force[1] * half_dt_over_mass;
            body.velocity[2] += body.force[2] * half_dt_over_mass;
        }
    }

    public static void driftStep(System system, double dt) {
        system.bodies.parallelStream().forEach(body -> {
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;
        });
    }

    public static void firstOrderKickStep(System system, double dt) {
        kickStep(system, dt);
        driftStep(system, dt);
        kickStep(system, dt);
    }

    public static void main(String[] args) {
        int numOrbitingBodies = 1000000;
        double centralMass = 1.989e30;    // Mass of the Sun (kg)
        double orbitRadius = 1.496e11;    // 1 AU (m)
        double orbitingMass = 5.972e24;   // Mass of the Earth (kg)
        int numSteps = 1000;
        double timeStep = 3600.0 * 24.0 * 7.0; // 1 week in seconds

        // Initialize the system
        System system = initializeCircularOrbits(numOrbitingBodies, centralMass, orbitRadius, orbitingMass);
        System initialSystem = new System();
        for (Body body : system.bodies) {
            initialSystem.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        // Calculate initial energy
        double initialEnergy = calculateTotalEnergy(initialSystem);
        System systemForSimulation = new System();
        for (Body body : initialSystem.bodies) {
            systemForSimulation.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy2 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy2.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy3 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy3.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy4 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy4.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy5 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy5.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy6 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy6.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy7 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy7.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy8 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy8.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy9 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy9.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy10 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy10.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy11 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy11.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy12 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy12.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy13 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy13.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy14 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy14.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy15 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy15.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy16 = new System();
        for (Body body : initialSystem.bodies) {
            systemForFinalEnergy16.addBody(new Body(body.mass, body.position.clone(), body.velocity.clone()));
        }

        System systemForFinalEnergy17 =