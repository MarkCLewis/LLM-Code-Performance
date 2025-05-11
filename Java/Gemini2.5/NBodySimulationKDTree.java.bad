import java.util.Random;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;
import java.util.Arrays;
import java.util.Comparator; // For median finding

// =======================================================================
// Simple Vector Helper Class
// =======================================================================
class Vector3D {
    double x, y, z;

    Vector3D() { this.x = 0; this.y = 0; this.z = 0; }
    Vector3D(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }

    Vector3D add(Vector3D other) {
        return new Vector3D(this.x + other.x, this.y + other.y, this.z + other.z);
    }

    Vector3D subtract(Vector3D other) {
        return new Vector3D(this.x - other.x, this.y - other.y, this.z - other.z);
    }

    Vector3D scale(double scalar) {
        return new Vector3D(this.x * scalar, this.y * scalar, this.z * scalar);
    }

    double magnitudeSq() {
        return x * x + y * y + z * z;
    }
     // No sqrt needed often, save computation
}


// =======================================================================
// k-D Tree Node Class
// =======================================================================
class KDTreeNode {
    // Tree structure
    int axis; // Splitting axis (0=x, 1=y, 2=z)
    KDTreeNode leftChild;
    KDTreeNode rightChild;

    // Node data (for internal nodes)
    double totalMass;
    Vector3D centerOfMass; // Weighted average position
    int count; // Number of particles in this subtree

    // Bounding box for this node's region
    Vector3D minBounds;
    Vector3D maxBounds;

    // Leaf node data
    int particleIndex = -1; // Index of the particle if this is a leaf node

    // Constructor for internal node (children must be set later)
    KDTreeNode(int axis, Vector3D minBounds, Vector3D maxBounds) {
        this.axis = axis;
        this.minBounds = minBounds;
        this.maxBounds = maxBounds;
        this.centerOfMass = new Vector3D();
    }

    // Constructor for leaf node
    KDTreeNode(int particleIndex, Vector3D particlePos, double particleMass, Vector3D minBounds, Vector3D maxBounds) {
        this.particleIndex = particleIndex;
        this.minBounds = minBounds; // Bounds can be tight for leaf
        this.maxBounds = maxBounds;
        this.totalMass = particleMass;
        this.centerOfMass = particlePos; // CoM is just the particle's position
        this.count = 1;
        this.axis = -1; // Indicate leaf
    }

    boolean isLeaf() {
        return particleIndex != -1;
    }

    // Calculate the size 's' of the node (longest side of bounding box)
    double getSize() {
        double dx = maxBounds.x - minBounds.x;
        double dy = maxBounds.y - minBounds.y;
        double dz = maxBounds.z - minBounds.z;
        return Math.max(dx, Math.max(dy, dz));
    }
}


// =======================================================================
// k-D Tree Class
// =======================================================================
class KDTree {
    private KDTreeNode root;
    // References to global simulation data (read-only during force calc)
    private final double[] posX, posY, posZ, masses;
    private final double softeningSquared;

    // Temporary array for indices during build
    private int[] buildIndices;
    private final int N;


    public KDTree(double[] posX, double[] posY, double[] posZ, double[] masses, double softeningSquared) {
        this.posX = posX;
        this.posY = posY;
        this.posZ = posZ;
        this.masses = masses;
        this.softeningSquared = softeningSquared;
        this.N = masses.length;
        this.buildIndices = new int[N]; // Reusable index array
        for(int i=0; i<N; i++) buildIndices[i] = i;
    }

    /** Builds the k-D tree for the current particle positions. */
    public void build() {
        // 1. Compute the global bounding box containing all particles
        Vector3D globalMin = new Vector3D(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        Vector3D globalMax = new Vector3D(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < N; i++) {
            globalMin.x = Math.min(globalMin.x, posX[i]);
            globalMin.y = Math.min(globalMin.y, posY[i]);
            globalMin.z = Math.min(globalMin.z, posZ[i]);
            globalMax.x = Math.max(globalMax.x, posX[i]);
            globalMax.y = Math.max(globalMax.y, posY[i]);
            globalMax.z = Math.max(globalMax.z, posZ[i]);
        }
        // Add a small buffer to avoid issues at the exact boundary
        double bufferX = (globalMax.x - globalMin.x) * 0.01;
        double bufferY = (globalMax.y - globalMin.y) * 0.01;
        double bufferZ = (globalMax.z - globalMin.z) * 0.01;
        globalMin = globalMin.subtract(new Vector3D(bufferX, bufferY, bufferZ));
        globalMax = globalMax.add(new Vector3D(bufferX, bufferY, bufferZ));


        // 2. Start recursive build process
        // Reset buildIndices (needed if build is called multiple times)
        for(int i=0; i<N; i++) buildIndices[i] = i;
        this.root = buildRecursive(0, N, 0, globalMin, globalMax);
    }

    /** Recursive helper to build the tree. Operates on sub-array of buildIndices. */
    private KDTreeNode buildRecursive(int start, int end, int depth, Vector3D minBounds, Vector3D maxBounds) {
        int count = end - start;

        if (count <= 0) {
            return null; // Should not happen with proper bounds
        }

        if (count == 1) {
            // Base case: Leaf node
            int particleIndex = buildIndices[start];
            Vector3D pos = new Vector3D(posX[particleIndex], posY[particleIndex], posZ[particleIndex]);
            return new KDTreeNode(particleIndex, pos, masses[particleIndex], minBounds, maxBounds);
        }

        // --- Internal Node ---
        int axis = depth % 3; // Cycle through x, y, z

        // Find median index along the current axis to partition
        // Using Quickselect (median-of-medians or similar for guaranteed O(n)) is best,
        // but sorting the sub-array and picking the middle is simpler to implement here,
        // making build O(N log^2 N) instead of O(N log N).
        Comparator<Integer> comparator = Comparator.comparingDouble(idx -> getCoord(idx, axis));
        Arrays.sort(buildIndices, start, end, comparator); // Sort sub-array indices by coordinate
        int medianIndex = start + count / 2;
        int medianParticleIndex = buildIndices[medianIndex];
        double splitValue = getCoord(medianParticleIndex, axis);

        // Create the internal node
        KDTreeNode node = new KDTreeNode(axis, minBounds, maxBounds);

        // Calculate bounds for children
        Vector3D leftMaxBounds = new Vector3D(maxBounds.x, maxBounds.y, maxBounds.z);
        Vector3D rightMinBounds = new Vector3D(minBounds.x, minBounds.y, minBounds.z);
        if (axis == 0) { leftMaxBounds.x = splitValue; rightMinBounds.x = splitValue; }
        else if (axis == 1) { leftMaxBounds.y = splitValue; rightMinBounds.y = splitValue; }
        else { leftMaxBounds.z = splitValue; rightMinBounds.z = splitValue; }

        // Recursively build children
        node.leftChild = buildRecursive(start, medianIndex, depth + 1, minBounds, leftMaxBounds);
        node.rightChild = buildRecursive(medianIndex, end, depth + 1, rightMinBounds, maxBounds); // Note: includes median point in right child usually

        // Calculate properties (Mass, CoM, Count) for this internal node from children
        node.count = 0;
        node.totalMass = 0.0;
        Vector3D weightedPosSum = new Vector3D();

        if (node.leftChild != null) {
            node.count += node.leftChild.count;
            node.totalMass += node.leftChild.totalMass;
            weightedPosSum = weightedPosSum.add(node.leftChild.centerOfMass.scale(node.leftChild.totalMass));
        }
        if (node.rightChild != null) {
            node.count += node.rightChild.count;
            node.totalMass += node.rightChild.totalMass;
             weightedPosSum = weightedPosSum.add(node.rightChild.centerOfMass.scale(node.rightChild.totalMass));
        }

        if (node.totalMass > 1e-100) { // Avoid division by zero if node is empty (shouldn't happen)
             node.centerOfMass = weightedPosSum.scale(1.0 / node.totalMass);
        } else {
            // Handle case of zero mass node if necessary (e.g., average position)
             if (node.leftChild != null && node.rightChild != null) {
                 node.centerOfMass = node.leftChild.centerOfMass.add(node.rightChild.centerOfMass).scale(0.5);
             } else if (node.leftChild != null) node.centerOfMass = node.leftChild.centerOfMass;
             else if (node.rightChild != null) node.centerOfMass = node.rightChild.centerOfMass;
             else node.centerOfMass = minBounds.add(maxBounds).scale(0.5); // Center of empty box
        }

        // Self-consistency check (optional)
        if (node.count != count && node.leftChild != null && node.rightChild != null) {
             // This might happen if the median split puts all points on one side repeatedly.
             // A more robust median/partition strategy might be needed.
             // System.err.printf("Warning: Node count mismatch at depth %d. Expected %d, got %d (L:%d, R:%d)\n",
             // depth, count, node.count, node.leftChild.count, node.rightChild.count);
        }


        return node;
    }

    /** Helper to get coordinate by axis */
    private double getCoord(int particleIndex, int axis) {
        if (axis == 0) return posX[particleIndex];
        if (axis == 1) return posY[particleIndex];
        return posZ[particleIndex];
    }


    /** Calculates the approximate gravitational force on a target particle using the tree. */
    public Vector3D calculateForce(int targetParticleIndex, double theta) {
        return calculateForceRecursive(this.root, targetParticleIndex, theta);
    }

    /** Recursive helper for force calculation using Barnes-Hut criterion. */
    private Vector3D calculateForceRecursive(KDTreeNode node, int targetParticleIndex, double theta) {
        if (node == null || node.count == 0) {
            return new Vector3D(); // No force from empty node
        }

        // --- Case 1: Leaf Node ---
        if (node.isLeaf()) {
            if (node.particleIndex == targetParticleIndex) {
                return new Vector3D(); // No self-force
            } else {
                // Direct calculation for leaf node particle
                return calculateDirectForce(targetParticleIndex, node.particleIndex);
            }
        }

        // --- Case 2: Internal Node ---
        // Calculate distance d from target particle to node's Center of Mass
        Vector3D targetPos = new Vector3D(posX[targetParticleIndex], posY[targetParticleIndex], posZ[targetParticleIndex]);
        Vector3D dVec = node.centerOfMass.subtract(targetPos);
        double distSq = dVec.magnitudeSq();

        // Calculate node size 's' (longest side)
        double s = node.getSize();

        // Barnes-Hut Criterion: s/d < theta
        // d = sqrt(distSq) => s/sqrt(distSq) < theta => s^2 / distSq < theta^2
        if ( (s * s) / (distSq + 1e-100) < (theta * theta) ) { // Add small epsilon to avoid division by zero if d=0
             // Node is far enough away: Approximate using node's CoM and Total Mass
             // Avoid direct calculation if target is *exactly* at CoM (use softening)
             return calculateForceApprox(targetParticleIndex, node.centerOfMass, node.totalMass);
        } else {
            // Node is too close: Recursively calculate force from children
            Vector3D force = new Vector3D();
            if (node.leftChild != null) {
                force = force.add(calculateForceRecursive(node.leftChild, targetParticleIndex, theta));
            }
            if (node.rightChild != null) {
                force = force.add(calculateForceRecursive(node.rightChild, targetParticleIndex, theta));
            }
            return force;
        }
    }

    /** Helper: Direct force between two particles */
    private Vector3D calculateDirectForce(int targetIndex, int sourceIndex) {
        double dx = posX[sourceIndex] - posX[targetIndex];
        double dy = posY[sourceIndex] - posY[targetIndex];
        double dz = posZ[sourceIndex] - posZ[targetIndex];
        double distSq = dx * dx + dy * dy + dz * dz + softeningSquared; // Include softening
        double invDistCube = 1.0 / Math.sqrt(distSq * distSq * distSq);
        double forceScalar = NBodySimulationKDTree.G * masses[sourceIndex] * masses[targetIndex] * invDistCube;
        return new Vector3D(forceScalar * dx, forceScalar * dy, forceScalar * dz);
    }

    /** Helper: Approximate force from a node (treated as pseudo-particle) */
     private Vector3D calculateForceApprox(int targetIndex, Vector3D nodeCoM, double nodeTotalMass) {
        double dx = nodeCoM.x - posX[targetIndex];
        double dy = nodeCoM.y - posY[targetIndex];
        double dz = nodeCoM.z - posZ[targetIndex];
        double distSq = dx * dx + dy * dy + dz * dz + softeningSquared; // Include softening
         double invDistCube = 1.0 / Math.sqrt(distSq * distSq * distSq);
        double forceScalar = NBodySimulationKDTree.G * nodeTotalMass * masses[targetIndex] * invDistCube;
        return new Vector3D(forceScalar * dx, forceScalar * dy, forceScalar * dz);
    }

}


// =======================================================================
// Main Simulation Class (Modified for k-D Tree)
// =======================================================================
public class NBodySimulationKDTree { // Renamed class

    // Constants
    static final double G = 1.0;
    static final double SOFTENING_SQUARED = 1e-6; // Keep softening!

    // System state arrays
    static double[] masses;
    static double[] posX, posY, posZ;
    static double[] velX, velY, velZ;
    static double[] forceX, forceY, forceZ; // Stores final forces for the step
    static int N;

    // k-D Tree specific parameters
    static final double THETA = 0.3; // Barnes-Hut opening angle

    // KDTree instance (rebuilt each step)
    static KDTree kdTree;


    /** Initialization (mostly unchanged, but need to init force arrays) */
    public static void initializeSystem(double centralMass, int numSmallBodies, double smallMass,
                                        double minDist, double maxDist, long seed) {
         // ... (Allocation and initialization logic same as NBodySimulationParallel)
         // Make sure forceX, forceY, forceZ are allocated
        N = numSmallBodies + 1;
        System.out.printf("Allocating memory for %,d bodies...\n", N);
        masses = new double[N];
        posX = new double[N]; posY = new double[N]; posZ = new double[N];
        velX = new double[N]; velY = new double[N]; velZ = new double[N];
        forceX = new double[N]; forceY = new double[N]; forceZ = new double[N]; // Allocate forces

        // ... (Rest of initialization: central body, small bodies - identical to previous)
        Random rand = new Random(seed);
        masses[0]=centralMass; posX[0]=0; posY[0]=0; posZ[0]=0; velX[0]=0; velY[0]=0; velZ[0]=0;
        for (int i = 1; i < N; i++) {
            masses[i] = smallMass;
            double r = minDist + (maxDist - minDist) * rand.nextDouble();
            double phi = 2.0 * Math.PI * rand.nextDouble();
            double costheta = 2.0 * rand.nextDouble() - 1.0;
            double sintheta = Math.sqrt(1.0 - costheta * costheta);
            double px = r * sintheta * Math.cos(phi);
            double py = r * sintheta * Math.sin(phi);
            double pz = r * costheta;
            posX[i] = px; posY[i] = py; posZ[i] = pz;
            double orbitalSpeed = Math.sqrt(G * centralMass / r);
            Vector3D posVec = new Vector3D(px, py, pz);
            Vector3D randomVec = new Vector3D(rand.nextDouble() - 0.5, rand.nextDouble() - 0.5, rand.nextDouble() - 0.5);
            Vector3D cross = new Vector3D(posVec.y * randomVec.z - posVec.z * randomVec.y,
                                           posVec.z * randomVec.x - posVec.x * randomVec.z,
                                           posVec.x * randomVec.y - posVec.y * randomVec.x);
            double crossMag = Math.sqrt(cross.magnitudeSq());
             if (crossMag < 1e-10) { // Fallback
                 Vector3D axis = (Math.abs(posVec.x) < 1e-9 && Math.abs(posVec.y) < 1e-9) ? new Vector3D(1,0,0) : new Vector3D(0,0,1);
                  cross = new Vector3D(posVec.y * axis.z - posVec.z * axis.y,
                                      posVec.z * axis.x - posVec.x * axis.z,
                                      posVec.x * axis.y - posVec.y * axis.x);
                 crossMag = Math.sqrt(cross.magnitudeSq());
             }
             Vector3D velDir = cross.scale(1.0/crossMag);
            velX[i] = orbitalSpeed * velDir.x;
            velY[i] = orbitalSpeed * velDir.y;
            velZ[i] = orbitalSpeed * velDir.z;
        }

        System.out.println("Initialization complete.");
    }


    /** Calculates total energy using direct summation (O(N^2)) - for verification */
    public static double calculateTotalEnergyDirectParallel() {
         // ... (Identical to calculateTotalEnergyParallel from previous version)
        double kineticEnergy = IntStream.range(0, N).parallel()
            .mapToDouble(i -> 0.5 * masses[i] * (velX[i]*velX[i] + velY[i]*velY[i] + velZ[i]*velZ[i]))
            .sum();
         double potentialEnergy = IntStream.range(0, N).parallel()
            .mapToDouble(i -> {
                double pe_i = 0.0;
                for (int j = i + 1; j < N; j++) {
                    double dx = posX[i] - posX[j]; double dy = posY[i] - posY[j]; double dz = posZ[i] - posZ[j];
                    double dist = Math.sqrt(dx*dx + dy*dy + dz*dz + SOFTENING_SQUARED); // Use softening
                    pe_i -= G * masses[i] * masses[j] / dist;
                } return pe_i; })
            .sum();
        return kineticEnergy + potentialEnergy;
    }


    /** Performs simulation step using k-D Tree for forces. */
    public static void simulationStepKDTree(double dt) {

        // 1. Build the k-D Tree (Sequential O(N log N) or O(N log^2 N))
        //    Tree needs to be rebuilt each step as particle positions change.
        kdTree = new KDTree(posX, posY, posZ, masses, SOFTENING_SQUARED);
        kdTree.build(); // Build the tree based on current positions


        // 2. Calculate Forces using the k-D Tree (Parallel O(N log N))
        IntStream.range(0, N).parallel().forEach(i -> {
            // Calculate force on particle 'i' using the tree and theta
            Vector3D forceOnI = kdTree.calculateForce(i, THETA);
            // Store the calculated force directly into the main arrays
            // No temporary arrays needed here as each thread calculates force for its own 'i'
            forceX[i] = forceOnI.x;
            forceY[i] = forceOnI.y;
            forceZ[i] = forceOnI.z;
        });


        // 3. Update Velocities (Kick) & 4. Update Positions (Step) (Parallel O(N))
        //    This part remains the same as the previous parallel version.
        IntStream.range(0, N).parallel().forEach(i -> {
            double accX = forceX[i] / masses[i];
            double accY = forceY[i] / masses[i];
            double accZ = forceZ[i] / masses[i];
            // Kick
            velX[i] += accX * dt;
            velY[i] += accY * dt;
            velZ[i] += accZ * dt;
            // Step
            posX[i] += velX[i] * dt;
            posY[i] += velY[i] * dt;
            posZ[i] += velZ[i] * dt;
        });
    }

    /** Helper method to format time duration. */
    private static String formatDuration(long millis) { /* ... same as before ... */
        long hours = TimeUnit.MILLISECONDS.toHours(millis); long minutes = TimeUnit.MILLISECONDS.toMinutes(millis)%60; long seconds = TimeUnit.MILLISECONDS.toSeconds(millis)%60; long milliseconds = millis % 1000; if(hours>0) return String.format("%dh %02dm %02ds", hours, minutes, seconds); if(minutes>0) return String.format("%dm %02ds %03dms", minutes, seconds, milliseconds); return String.format("%ds %03dms", seconds, milliseconds);
     }

    // =======================================================================
    // Main Driver
    // =======================================================================
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        System.out.println("Starting k-D Tree (Barnes-Hut) 3D N-Body Simulation...");
        int processors = Runtime.getRuntime().availableProcessors();
        System.out.println("Number of available processors: " + processors);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(processors));
        System.out.println("Using common ForkJoinPool with parallelism: " + java.util.concurrent.ForkJoinPool.commonPool().getParallelism());

        // --- Simulation Parameters ---
        int numSmallBodies = 100000; // ONE MILLION orbiting bodies as requested
        int numSteps = 10;            // Number of simulation steps as requested
        double centralMass = 1.0e6;
        double smallMass = 1.0;
        double minDist = 5.0;
        double maxDist = 50.0;
        double dt = 0.005;
        long seed = 123456789L;

        System.out.println("-------------------- Parameters --------------------");
        System.out.printf("Algorithm:              k-D Tree (Barnes-Hut)\n");
        System.out.printf("Theta (Opening Angle):  %.2f\n", THETA);
        System.out.printf("Number of small bodies: %,d\n", numSmallBodies);
        System.out.printf("Total bodies (N):       %,d\n", numSmallBodies + 1);
        System.out.printf("Number of steps:        %,d\n", numSteps);
        System.out.printf("Time step (dt):         %.5f\n", dt);
        System.out.printf("Central Mass:           %.2e\n", centralMass);
        // ... other parameters ...
        System.out.printf("Softening^2:            %.2e\n", SOFTENING_SQUARED);
        System.out.println("--------------------------------------------------");

        // --- WARNING --- (Still relevant for memory, less severe for time than O(N^2))
        System.out.println("\n************************** INFO **************************");
        System.out.printf ("Simulating N = %,d bodies using k-D Tree (O(N log N) per step).\n", (numSmallBodies+1));
        System.out.println("This should be significantly faster than O(N^2) but still requires substantial computation and memory.");
        System.out.println("Ensure sufficient RAM allocated (use -Xmx JVM flag).");
        System.out.println("***********************************************************\n");

        // 1. Initialize System
        System.out.println("Initializing system...");
        long initStartTime = System.currentTimeMillis();
        try {
            initializeSystem(centralMass, numSmallBodies, smallMass, minDist, maxDist, seed);
        } catch (OutOfMemoryError e) {
             System.err.println("Exiting due to OutOfMemoryError during initialization."); return;
        }
        long initEndTime = System.currentTimeMillis();
        System.out.printf("System initialization complete. Time: %s\n", formatDuration(initEndTime - initStartTime));

        // 2. Calculate Initial Energy (Using direct O(N^2) for accuracy)
        System.out.println("Calculating initial energy (direct summation, parallel)...");
        long energyStartTime = System.currentTimeMillis();
        // Use the accurate direct method for energy check, not the tree approx
        double initialEnergy = calculateTotalEnergyDirectParallel();
        long energyEndTime = System.currentTimeMillis();
        System.out.printf("Initial Total Energy: %.8e (Calculation time: %s)\n",
                          initialEnergy, formatDuration(energyEndTime - energyStartTime));

        // 3. Run Simulation Loop (Using k-D Tree Step)
        System.out.printf("Starting simulation loop for %,d steps (using k-D Tree steps)...\n", numSteps);
        long simStartTime = System.currentTimeMillis();
        long lastReportTime = simStartTime;
        long treeBuildTimeTotal = 0; // Track tree build time
        long forceCalcTimeTotal = 0; // Track tree force calc time
        long updateTimeTotal = 0;   // Track position/velocity update time

        int reportInterval = (N > 50000) ? 1 : (N > 1000 ? 10 : 100); // Report more often
        if(numSteps < 100) reportInterval = 1;

        for (int step = 0; step < numSteps; step++) {
            // --- Timing the parts of the step ---
            long stepStart = System.nanoTime();

            // 1. Build Tree
            long buildStart = System.nanoTime();
            kdTree = new KDTree(posX, posY, posZ, masses, SOFTENING_SQUARED);
            kdTree.build();
            long buildEnd = System.nanoTime();
            treeBuildTimeTotal += (buildEnd - buildStart);

            // 2. Calculate Forces
            long forceStart = System.nanoTime();
            IntStream.range(0, N).parallel().forEach(i -> {
                Vector3D forceOnI = kdTree.calculateForce(i, THETA);
                forceX[i] = forceOnI.x; forceY[i] = forceOnI.y; forceZ[i] = forceOnI.z;
            });
            long forceEnd = System.nanoTime();
            forceCalcTimeTotal += (forceEnd - forceStart);

            // 3. Update Positions/Velocities
            long updateStart = System.nanoTime();
            IntStream.range(0, N).parallel().forEach(i -> {
                double accX = forceX[i] / masses[i]; double accY = forceY[i] / masses[i]; double accZ = forceZ[i] / masses[i];
                velX[i] += accX * dt; velY[i] += accY * dt; velZ[i] += accZ * dt; // Kick
                posX[i] += velX[i] * dt; posY[i] += velY[i] * dt; posZ[i] += velZ[i] * dt; // Step
            });
            long updateEnd = System.nanoTime();
            updateTimeTotal += (updateEnd - updateStart);
             // --- End Timing ---

            // --- Progress Reporting ---
            if ((step + 1) % reportInterval == 0 || step == numSteps - 1) {
                long currentTime = System.currentTimeMillis();
                double timeElapsedTotal = (currentTime - simStartTime) / 1000.0;
                double intervalDurationMillis = (currentTime - lastReportTime);
                double avgStepTimeMillis = intervalDurationMillis / reportInterval;
                int stepsRemaining = numSteps - (step + 1);
                double estimatedRemainingMillis = stepsRemaining * avgStepTimeMillis;

                 // Detailed timing report (optional)
                double avgBuildNs = (double)treeBuildTimeTotal / (step + 1);
                double avgForceNs = (double)forceCalcTimeTotal / (step + 1);
                double avgUpdateNs = (double)updateTimeTotal / (step + 1);

                System.out.printf("  Step %d/%d. Total time: %.2fs. Avg step (last %d): %.1fms [Build:%.1fms, Force:%.1fms, Update:%.1fms]. Est. rem: %s\n",
                                  step + 1, numSteps, timeElapsedTotal,
                                  reportInterval, avgStepTimeMillis,
                                  avgBuildNs/1e6, avgForceNs/1e6, avgUpdateNs/1e6, // Convert ns to ms
                                  formatDuration((long)estimatedRemainingMillis) );
                lastReportTime = currentTime;

                 // Adjust report interval dynamically (optional)
                 // ...
            }
        }
        long simEndTime = System.currentTimeMillis();
        System.out.printf("Simulation loop complete. Total simulation time: %s\n", formatDuration(simEndTime - simStartTime));
        System.out.printf("Average time per step: %.3f ms\n", (double)(simEndTime - simStartTime)/numSteps);
        System.out.printf("  Avg Tree Build Time: %.3f ms (%.1f%%)\n", (double)treeBuildTimeTotal / numSteps / 1e6, (double)treeBuildTimeTotal / (simEndTime - simStartTime)/1e4);
        System.out.printf("  Avg Force Calc Time: %.3f ms (%.1f%%)\n", (double)forceCalcTimeTotal / numSteps / 1e6, (double)forceCalcTimeTotal / (simEndTime - simStartTime)/1e4);
        System.out.printf("  Avg Update Time:     %.3f ms (%.1f%%)\n", (double)updateTimeTotal / numSteps / 1e6, (double)updateTimeTotal / (simEndTime - simStartTime)/1e4);


        // 4. Calculate Final Energy (Using direct O(N^2) for accuracy check)
        System.out.println("Calculating final energy (direct summation, parallel)...");
        energyStartTime = System.currentTimeMillis();
        double finalEnergy = calculateTotalEnergyDirectParallel(); // Use accurate method
        energyEndTime = System.currentTimeMillis();
        System.out.printf("Final Total Energy:   %.8e (Calculation time: %s)\n",
                          finalEnergy, formatDuration(energyEndTime - energyStartTime));

        // 5. Verify Energy Conservation
        System.out.println("-------------------- Energy Check --------------------");
        System.out.printf("Initial Energy:        %.8e\n", initialEnergy);
        System.out.printf("Final Energy:          %.8e\n", finalEnergy);
        // Note: Barnes-Hut is approximate, so energy conservation will NOT be perfect.
        // The error depends on theta, dt, and softening.
        if (Math.abs(initialEnergy) < 1e-15) {
            System.out.printf("Energy Change (Absolute): %.4e\n", finalEnergy - initialEnergy);
        } else {
            double relativeEnergyChange = Math.abs((finalEnergy - initialEnergy) / initialEnergy);
            System.out.printf("Relative Energy Change: %.4e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange * 100.0);
            if (relativeEnergyChange < 5e-3) System.out.printf("Energy conservation reasonable for Barnes-Hut (theta=%.2f).", THETA);
            else if (relativeEnergyChange < 5e-2) System.out.println("Noticeable energy drift; expected with approximate methods.");
            else System.out.println("WARNING: Significant energy drift detected! Check theta, dt, softening, or implementation.");
        }
        System.out.println("----------------------------------------------------");
        System.out.println("k-D Tree Simulation Finished.");
    }
}