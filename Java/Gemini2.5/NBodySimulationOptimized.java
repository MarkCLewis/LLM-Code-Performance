import java.util.Random;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;
// Removed Arrays and Comparator imports as we replace sort with Quickselect

// =======================================================================
// Simple Vector Helper Class (Now MUTABLE for optimization)
// =======================================================================
class Vector3D {
    double x, y, z;

    Vector3D() { this.x = 0; this.y = 0; this.z = 0; }
    Vector3D(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }

    // --- MUTABLE OPERATIONS ---
    // Add 'other' to this vector
    void addInPlace(Vector3D other) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;
    }
     // Add scaled 'other' to this vector
    void addScaledInPlace(Vector3D other, double scalar) {
        this.x += other.x * scalar;
        this.y += other.y * scalar;
        this.z += other.z * scalar;
    }
     // Add components directly
    void addInPlace(double dx, double dy, double dz) {
        this.x += dx;
        this.y += dy;
        this.z += dz;
    }
     // Set this vector to be 'other' scaled
    void setToScale(Vector3D other, double scalar){
        this.x = other.x * scalar;
        this.y = other.y * scalar;
        this.z = other.z * scalar;
    }
     // Reset to zero
    void zero() {
        this.x = 0; this.y = 0; this.z = 0;
    }

    // --- IMMUTABLE OPERATIONS (Create new vector) ---
    Vector3D subtract(Vector3D other) { // Often needed for differences
        return new Vector3D(this.x - other.x, this.y - other.y, this.z - other.z);
    }
    /** Gemini didn't provide this one method. */
    Vector3D add(Vector3D other) { // Often needed for differences
        return new Vector3D(this.x + other.x, this.y + other.y, this.z + other.z);
    }
    Vector3D scale(double scalar) { // Needed for CoM calc sometimes
        return new Vector3D(this.x * scalar, this.y * scalar, this.z * scalar);
    }

    double magnitudeSq() {
        return x * x + y * y + z * z;
    }
}


// =======================================================================
// k-D Tree Node Class (Mostly Unchanged)
// =======================================================================
class KDTreeNode {
    // Tree structure
    int axis;
    KDTreeNode leftChild;
    KDTreeNode rightChild;

    // Node data
    double totalMass;
    Vector3D centerOfMass; // Still immutable Vector3D here is fine for node state
    int count;

    // Bounding box
    Vector3D minBounds; // Immutable is fine
    Vector3D maxBounds; // Immutable is fine

    // Leaf data
    int particleIndex = -1;

    // Constructor for internal node
    KDTreeNode(int axis, Vector3D minBounds, Vector3D maxBounds) {
        this.axis = axis; this.minBounds = minBounds; this.maxBounds = maxBounds;
        this.centerOfMass = new Vector3D(); // Initialize CoM
    }
    // Constructor for leaf node
    KDTreeNode(int particleIndex, Vector3D particlePos, double particleMass, Vector3D minBounds, Vector3D maxBounds) {
        this.particleIndex = particleIndex; this.minBounds = minBounds; this.maxBounds = maxBounds;
        this.totalMass = particleMass; this.centerOfMass = particlePos; // CoM is particle pos
        this.count = 1; this.axis = -1;
    }
    boolean isLeaf() { return particleIndex != -1; }
    double getSize() { /* ... same as before ... */
        double dx = maxBounds.x - minBounds.x; double dy = maxBounds.y - minBounds.y; double dz = maxBounds.z - minBounds.z; return Math.max(dx, Math.max(dy, dz));
    }
}


// =======================================================================
// k-D Tree Class (Optimized Build and Force Calculation)
// =======================================================================
class KDTree {
    private KDTreeNode root;
    private final double[] posX, posY, posZ, masses;
    private final double softeningSquared;
    private final int N;
    private int[] buildIndices; // indices used during build

    // Pre-allocated reusable vectors for force calculation to reduce GC
    private final Vector3D forceCalc_targetPos = new Vector3D();
    private final Vector3D forceCalc_dVec = new Vector3D();


    public KDTree(double[] posX, double[] posY, double[] posZ, double[] masses, double softeningSquared) {
        this.posX = posX; this.posY = posY; this.posZ = posZ; this.masses = masses;
        this.softeningSquared = softeningSquared; this.N = masses.length;
        this.buildIndices = new int[N];
    }

    /** Builds the k-D tree using Quickselect for O(N log N) average time. */
    public void build() {
        // Reset buildIndices array
        for(int i=0; i<N; i++) buildIndices[i] = i;

        // Compute global bounding box
        Vector3D globalMin = new Vector3D(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        Vector3D globalMax = new Vector3D(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < N; i++) {
             globalMin.x = Math.min(globalMin.x, posX[i]); globalMin.y = Math.min(globalMin.y, posY[i]); globalMin.z = Math.min(globalMin.z, posZ[i]);
             globalMax.x = Math.max(globalMax.x, posX[i]); globalMax.y = Math.max(globalMax.y, posY[i]); globalMax.z = Math.max(globalMax.z, posZ[i]);
        }
        double bufferX = (globalMax.x - globalMin.x) * 0.01 + 1e-9; // Add tiny absolute buffer too
        double bufferY = (globalMax.y - globalMin.y) * 0.01 + 1e-9;
        double bufferZ = (globalMax.z - globalMin.z) * 0.01 + 1e-9;
        globalMin = globalMin.subtract(new Vector3D(bufferX, bufferY, bufferZ));
        globalMax = globalMax.add(new Vector3D(bufferX, bufferY, bufferZ));

        this.root = buildRecursive(0, N, 0, globalMin, globalMax);
    }

     /** Recursive build using Quickselect */
    private KDTreeNode buildRecursive(int start, int end, int depth, Vector3D minBounds, Vector3D maxBounds) {
        int count = end - start;
        if (count <= 0) return null;
        if (count == 1) {
            int pIdx = buildIndices[start];
            Vector3D pos = new Vector3D(posX[pIdx], posY[pIdx], posZ[pIdx]);
            return new KDTreeNode(pIdx, pos, masses[pIdx], minBounds, maxBounds);
        }

        int axis = depth % 3;
        int medianIndex = start + count / 2; // Target rank for median

        // --- Quickselect to find the median element ---
        quickSelect(medianIndex, start, end - 1, axis);
        // Now buildIndices[medianIndex] holds the index of the particle at the median coordinate

        int medianParticleIndex = buildIndices[medianIndex];
        double splitValue = getCoord(medianParticleIndex, axis);

        // Create internal node
        KDTreeNode node = new KDTreeNode(axis, minBounds, maxBounds);

        // Calculate bounds for children
        Vector3D leftMaxBounds = new Vector3D(maxBounds.x, maxBounds.y, maxBounds.z);
        Vector3D rightMinBounds = new Vector3D(minBounds.x, minBounds.y, minBounds.z);
        if (axis == 0) { leftMaxBounds.x = splitValue; rightMinBounds.x = splitValue; }
        else if (axis == 1) { leftMaxBounds.y = splitValue; rightMinBounds.y = splitValue; }
        else { leftMaxBounds.z = splitValue; rightMinBounds.z = splitValue; }

        // Recursively build (pass median to the right child)
        node.leftChild = buildRecursive(start, medianIndex, depth + 1, minBounds, leftMaxBounds);
        // NOTE: The element *at* medianIndex is passed to the right subtree here.
        // If we wanted strict < median on left, >= median on right, the partition logic inside
        // quickSelect/partition needs to handle equal elements carefully, and the recursive
        // calls adjusted. This simpler approach usually works well enough.
        node.rightChild = buildRecursive(medianIndex, end, depth + 1, rightMinBounds, maxBounds);

        // Calculate node properties (Mass, CoM, Count)
        node.count = 0; node.totalMass = 0.0;
        Vector3D weightedPosSum = new Vector3D(); // Use temp Vector3D
        if (node.leftChild != null) {
            node.count += node.leftChild.count;
            node.totalMass += node.leftChild.totalMass;
            weightedPosSum.addScaledInPlace(node.leftChild.centerOfMass, node.leftChild.totalMass);
        }
        if (node.rightChild != null) {
            node.count += node.rightChild.count;
            node.totalMass += node.rightChild.totalMass;
            weightedPosSum.addScaledInPlace(node.rightChild.centerOfMass, node.rightChild.totalMass);
        }
        if (node.totalMass > 1e-100) {
             node.centerOfMass = weightedPosSum.scale(1.0 / node.totalMass); // Need immutable scale here for CoM
        } else if (node.count > 0) {
             // Handle zero mass case if needed (average position might be okay)
             node.centerOfMass = minBounds.add(maxBounds).scale(0.5); // Geometric center
        }
        // Ensure count matches if needed for debugging, usually it should
         if(node.count != count && node.leftChild != null && node.rightChild != null && node.leftChild.count + node.rightChild.count != count) {
            // This indicates an issue in partitioning or recursive bounds
            // System.err.printf("Warning: Node count mismatch at depth %d. Start %d, End %d, Count %d, NodeCount %d (L:%d, R:%d)\n",
            //    depth, start, end, count, node.count, node.leftChild.count, node.rightChild.count);
         }


        return node;
    }

    // --- Quickselect Implementation ---
    private void quickSelect(int k, int start, int end, int axis) {
        if (start == end) return;
        int pivotIndex = start + (end - start) / 2; // Simple pivot choice
        pivotIndex = partition(start, end, pivotIndex, axis);
        if (k == pivotIndex) {
            return;
        } else if (k < pivotIndex) {
            quickSelect(k, start, pivotIndex - 1, axis);
        } else {
            quickSelect(k, pivotIndex + 1, end, axis);
        }
    }

    private int partition(int start, int end, int pivotIndex, int axis) {
        double pivotValue = getCoord(buildIndices[pivotIndex], axis);
        swap(pivotIndex, end); // Move pivot to end
        int storeIndex = start;
        for (int i = start; i < end; i++) {
            if (getCoord(buildIndices[i], axis) < pivotValue) {
                swap(storeIndex, i);
                storeIndex++;
            }
        }
        swap(storeIndex, end); // Move pivot to its final place
        return storeIndex;
    }

    private void swap(int i, int j) {
        int temp = buildIndices[i];
        buildIndices[i] = buildIndices[j];
        buildIndices[j] = temp;
    }
    // --- End Quickselect ---


    /** Helper to get coordinate */
    private double getCoord(int particleIndex, int axis) { /* ... same ... */
         if (axis == 0) return posX[particleIndex]; if (axis == 1) return posY[particleIndex]; return posZ[particleIndex];
    }

    /** Calculates force, adding result to the mutable resultForce vector. */
    public void calculateForce(int targetParticleIndex, double theta, Vector3D resultForce) {
        // Initialize / reset the result vector for this particle
        resultForce.zero();
        // Set the reusable target position vector
        forceCalc_targetPos.x = posX[targetParticleIndex];
        forceCalc_targetPos.y = posY[targetParticleIndex];
        forceCalc_targetPos.z = posZ[targetParticleIndex];

        calculateForceRecursive(this.root, targetParticleIndex, theta, resultForce);
    }

    /** Recursive force calculation, adds force to resultForce. */
    private void calculateForceRecursive(KDTreeNode node, int targetParticleIndex, double theta, Vector3D resultForce) {
        if (node == null || node.count == 0) return;

        // --- Leaf Node ---
        if (node.isLeaf()) {
            if (node.particleIndex != targetParticleIndex) {
                // Calculate direct force and add it to resultForce
                addDirectForce(targetParticleIndex, node.particleIndex, resultForce);
            }
            return; // Done with leaf
        }

        // --- Internal Node ---
        // Distance from target particle to node's CoM (use pre-allocated dVec)
        forceCalc_dVec.x = node.centerOfMass.x - forceCalc_targetPos.x;
        forceCalc_dVec.y = node.centerOfMass.y - forceCalc_targetPos.y;
        forceCalc_dVec.z = node.centerOfMass.z - forceCalc_targetPos.z;
        double distSq = forceCalc_dVec.x * forceCalc_dVec.x +
                        forceCalc_dVec.y * forceCalc_dVec.y +
                        forceCalc_dVec.z * forceCalc_dVec.z;

        double s = node.getSize(); // Node size

        // Barnes-Hut Criterion: s^2 / distSq < theta^2
        if ( (s * s) / (distSq + 1e-100) < (theta * theta) ) {
            // Approximate force from node's CoM (add to resultForce)
            addForceApprox(targetParticleIndex, node.centerOfMass, node.totalMass, resultForce);
        } else {
            // Recurse on children
            calculateForceRecursive(node.leftChild, targetParticleIndex, theta, resultForce);
            calculateForceRecursive(node.rightChild, targetParticleIndex, theta, resultForce);
        }
    }

    /** Helper: Adds direct force between two particles to resultForce. */
    private void addDirectForce(int targetIndex, int sourceIndex, Vector3D resultForce) {
        double dx = posX[sourceIndex] - posX[targetIndex]; // Use target pos from member vec? No, need targetIndex here.
        double dy = posY[sourceIndex] - posY[targetIndex];
        double dz = posZ[sourceIndex] - posZ[targetIndex];
        double distSq = dx * dx + dy * dy + dz * dz + softeningSquared;
        // Check for zero distance even with softening (should be rare)
        if (distSq < 1e-100) return; // Avoid NaN/Infinity
        double invDist = 1.0 / Math.sqrt(distSq); // 1/dist
        double invDistCube = invDist * invDist * invDist; // 1/dist^3
        double forceScalar = NBodySimulationOptimized.G * masses[sourceIndex] * masses[targetIndex] * invDistCube;
        resultForce.addInPlace(forceScalar * dx, forceScalar * dy, forceScalar * dz);
    }

    /** Helper: Adds approximate force from a node to resultForce. */
     private void addForceApprox(int targetIndex, Vector3D nodeCoM, double nodeTotalMass, Vector3D resultForce) {
        double dx = nodeCoM.x - posX[targetIndex]; // Use target pos from member vec? No, need targetIndex here.
        double dy = nodeCoM.y - posY[targetIndex];
        double dz = nodeCoM.z - posZ[targetIndex];
        double distSq = dx * dx + dy * dy + dz * dz + softeningSquared;
        if (distSq < 1e-100) return;
        double invDist = 1.0 / Math.sqrt(distSq);
        double invDistCube = invDist * invDist * invDist;
        double forceScalar = NBodySimulationOptimized.G * nodeTotalMass * masses[targetIndex] * invDistCube;
        resultForce.addInPlace(forceScalar * dx, forceScalar * dy, forceScalar * dz);
    }
}


// =======================================================================
// Main Simulation Class (Using Optimized k-D Tree)
// =======================================================================
public class NBodySimulationOptimized { // Renamed class

    static final double G = 1.0;
    static final double SOFTENING_SQUARED = 1e-6;
    static final double THETA = 0.3;

    static double[] masses, invMasses; // Added inverse mass array
    static double[] posX, posY, posZ;
    static double[] velX, velY, velZ;
    // Force arrays are technically not needed if updates happen right after calc,
    // but keep them for clarity and potential future modifications.
    static double[] forceX, forceY, forceZ;
    static int N;
    static KDTree kdTree;

    // Thread-local storage for reusable force vectors to avoid contention in parallel stream
    // Each thread working on the parallel force calculation gets its own Vector3D instance.
    private static final ThreadLocal<Vector3D> threadLocalForceVector =
        ThreadLocal.withInitial(Vector3D::new);


    public static void initializeSystem(double centralMass, int numSmallBodies, double smallMass,
                                        double minDist, double maxDist, long seed) {
        N = numSmallBodies + 1;
        System.out.printf("Allocating memory for %,d bodies...\n", N);
        masses = new double[N]; invMasses = new double[N]; // Allocate invMass
        posX = new double[N]; posY = new double[N]; posZ = new double[N];
        velX = new double[N]; velY = new double[N]; velZ = new double[N];
        forceX = new double[N]; forceY = new double[N]; forceZ = new double[N];

        // ... Initialization logic identical to previous version ...
        Random rand = new Random(seed);
        masses[0]=centralMass; invMasses[0] = 1.0/centralMass; // Store invMass
        posX[0]=0; posY[0]=0; posZ[0]=0; velX[0]=0; velY[0]=0; velZ[0]=0;
        for (int i = 1; i < N; i++) {
            masses[i] = smallMass; invMasses[i] = 1.0/smallMass; // Store invMass
             // ... rest of position/velocity init ...
            double r = minDist + (maxDist - minDist) * rand.nextDouble(); double phi = 2.0 * Math.PI * rand.nextDouble(); double costheta = 2.0 * rand.nextDouble() - 1.0; double sintheta = Math.sqrt(1.0 - costheta*costheta);
            double px = r * sintheta * Math.cos(phi); double py = r * sintheta * Math.sin(phi); double pz = r * costheta; posX[i]=px; posY[i]=py; posZ[i]=pz;
            double orbitalSpeed = Math.sqrt(G * centralMass / r); Vector3D posVec = new Vector3D(px, py, pz); Vector3D randomVec = new Vector3D(rand.nextDouble()-0.5, rand.nextDouble()-0.5, rand.nextDouble()-0.5);
            Vector3D cross = new Vector3D(posVec.y * randomVec.z - posVec.z * randomVec.y, posVec.z * randomVec.x - posVec.x * randomVec.z, posVec.x * randomVec.y - posVec.y * randomVec.x); double crossMag = Math.sqrt(cross.magnitudeSq());
            if(crossMag < 1e-10){ Vector3D axis = (Math.abs(posVec.x)<1e-9 && Math.abs(posVec.y)<1e-9) ? new Vector3D(1,0,0): new Vector3D(0,0,1); cross = new Vector3D(posVec.y*axis.z-posVec.z*axis.y, posVec.z*axis.x-posVec.x*axis.z, posVec.x*axis.y-posVec.y*axis.x); crossMag=Math.sqrt(cross.magnitudeSq()); }
            Vector3D velDir = cross.scale(1.0/crossMag); velX[i]=orbitalSpeed*velDir.x; velY[i]=orbitalSpeed*velDir.y; velZ[i]=orbitalSpeed*velDir.z;
        }
        System.out.println("Initialization complete.");
    }

    /** Accurate Energy Calc (Unchanged) */
    public static double calculateTotalEnergyDirectParallel() { /* ... same as before ... */
         double kineticEnergy = IntStream.range(0, N).parallel().mapToDouble(i -> 0.5 * masses[i] * (velX[i]*velX[i] + velY[i]*velY[i] + velZ[i]*velZ[i])).sum();
         double potentialEnergy = IntStream.range(0, N).parallel().mapToDouble(i -> { double pe_i=0.0; for(int j=i+1; j<N; j++){ double dx=posX[i]-posX[j]; double dy=posY[i]-posY[j]; double dz=posZ[i]-posZ[j]; double dist = Math.sqrt(dx*dx+dy*dy+dz*dz + SOFTENING_SQUARED); pe_i -= G*masses[i]*masses[j]/dist; } return pe_i; }).sum();
         return kineticEnergy+potentialEnergy;
     }


    /** Optimized Simulation Step */
    public static void simulationStepOptimized(double dt) {

        // 1. Build k-D Tree (Now O(N log N) average)
        kdTree = new KDTree(posX, posY, posZ, masses, SOFTENING_SQUARED);
        kdTree.build();


        // 2. Calculate Forces using k-D Tree (Parallel O(N log N), Reduced GC)
        IntStream.range(0, N).parallel().forEach(i -> {
            // Get thread-local vector to store the result, avoids creating new Vector3D per particle
            Vector3D forceOnI = threadLocalForceVector.get();
            // Calculate force, result is placed into forceOnI
            kdTree.calculateForce(i, THETA, forceOnI);
            // Store result in main arrays
            forceX[i] = forceOnI.x;
            forceY[i] = forceOnI.y;
            forceZ[i] = forceOnI.z;
        });


        // 3. Update Velocities & Positions (Parallel O(N), uses invMass)
        IntStream.range(0, N).parallel().forEach(i -> {
            // Use pre-calculated inverse mass
            double accX = forceX[i] * invMasses[i];
            double accY = forceY[i] * invMasses[i];
            double accZ = forceZ[i] * invMasses[i];
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

    /** Format duration helper */
    private static String formatDuration(long millis) { /* ... same ... */
        long h=TimeUnit.MILLISECONDS.toHours(millis), m=TimeUnit.MILLISECONDS.toMinutes(millis)%60, s=TimeUnit.MILLISECONDS.toSeconds(millis)%60, ms=millis%1000; if(h>0) return String.format("%dh %02dm %02ds",h,m,s); if(m>0) return String.format("%dm %02ds %03dms",m,s,ms); return String.format("%ds %03dms",s,ms);
    }


    // =======================================================================
    // Main Driver
    // =======================================================================
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        System.out.println("Starting OPTIMIZED k-D Tree (Barnes-Hut) 3D N-Body Simulation...");
        // ... (Processor info printout same) ...
        int processors = Runtime.getRuntime().availableProcessors(); System.out.println("Available processors: " + processors); System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(processors)); System.out.println("Using common ForkJoinPool parallelism: " + java.util.concurrent.ForkJoinPool.commonPool().getParallelism());

        // --- Parameters ---
        int numSmallBodies = 100000; // ONE MILLION orbiting bodies as requested
        int numSteps = 10;            // Number of simulation steps as requested
        double centralMass = 1.0e6; double smallMass = 1.0;
        double minDist = 5.0; double maxDist = 50.0;
        double dt = 0.005; long seed = 123456789L;

        // ... (Parameter printout same, mention "Optimized") ...
        System.out.println("-------------------- Parameters --------------------");
        System.out.printf("Algorithm:              k-D Tree (Barnes-Hut, Optimized)\n");
        System.out.printf("Theta (Opening Angle):  %.2f\n", THETA);
        System.out.printf("Number of small bodies: %,d\n", numSmallBodies);
        // ... rest of params ...
        System.out.println("--------------------------------------------------");
        System.out.println("\n************************** INFO **************************"); // Same info message
        System.out.printf ("Simulating N = %,d bodies using optimized k-D Tree (O(N log N) per step).\n", (numSmallBodies+1));
        System.out.println("Should be faster than previous k-D tree version due to O(N log N) build and reduced object allocation.");
        // ... rest of info ...
        System.out.println("***********************************************************\n");


        // 1. Initialize System
        System.out.println("Initializing system...");
        long initStartTime = System.currentTimeMillis();
        try { initializeSystem(centralMass, numSmallBodies, smallMass, minDist, maxDist, seed); }
        catch (OutOfMemoryError e) { System.err.println("Exiting: OutOfMemoryError."); return; }
        long initEndTime = System.currentTimeMillis();
        System.out.printf("Initialization complete. Time: %s\n", formatDuration(initEndTime - initStartTime));


        // 2. Calculate Initial Energy (Direct, Parallel)
        System.out.println("Calculating initial energy (direct summation, parallel)...");
        long energyStartTime = System.currentTimeMillis();
        double initialEnergy = calculateTotalEnergyDirectParallel();
        long energyEndTime = System.currentTimeMillis();
        System.out.printf("Initial Total Energy: %.8e (Calculation time: %s)\n", initialEnergy, formatDuration(energyEndTime - energyStartTime));


        // 3. Run Simulation Loop (Optimized Step)
        System.out.printf("Starting simulation loop for %,d steps (using Optimized k-D Tree steps)...\n", numSteps);
        long simStartTime = System.currentTimeMillis(); long lastReportTime = simStartTime;
        long treeBuildTimeTotal = 0; long forceCalcTimeTotal = 0; long updateTimeTotal = 0;
        int reportInterval = (N > 50000) ? 1 : (N > 1000 ? 10 : 100); if(numSteps < 100) reportInterval = 1;

        for (int step = 0; step < numSteps; step++) {
             long stepStart = System.nanoTime();

             // --- Build Tree ---
             long buildStart = System.nanoTime();
             kdTree = new KDTree(posX, posY, posZ, masses, SOFTENING_SQUARED);
             kdTree.build(); // O(N log N) avg
             long buildEnd = System.nanoTime();
             treeBuildTimeTotal += (buildEnd - buildStart);

             // --- Calculate Forces ---
             long forceStart = System.nanoTime();
             IntStream.range(0, N).parallel().forEach(i -> {
                 Vector3D forceOnI = threadLocalForceVector.get(); // Get thread-local vector
                 kdTree.calculateForce(i, THETA, forceOnI); // Calculate into forceOnI
                 forceX[i] = forceOnI.x; forceY[i] = forceOnI.y; forceZ[i] = forceOnI.z;
             });
             long forceEnd = System.nanoTime();
             forceCalcTimeTotal += (forceEnd - forceStart);

             // --- Update Positions/Velocities ---
             long updateStart = System.nanoTime();
             IntStream.range(0, N).parallel().forEach(i -> {
                 double accX = forceX[i] * invMasses[i]; double accY = forceY[i] * invMasses[i]; double accZ = forceZ[i] * invMasses[i];
                 velX[i] += accX * dt; velY[i] += accY * dt; velZ[i] += accZ * dt; // Kick
                 posX[i] += velX[i] * dt; posY[i] += velY[i] * dt; posZ[i] += velZ[i] * dt; // Step
             });
             long updateEnd = System.nanoTime();
             updateTimeTotal += (updateEnd - updateStart);

            // --- Progress Reporting --- (same as previous version)
            if ((step + 1) % reportInterval == 0 || step == numSteps - 1) {
                 long currentTime = System.currentTimeMillis(); double timeElapsedTotal = (currentTime - simStartTime)/1000.0;
                 double intervalDurationMillis = (currentTime - lastReportTime); double avgStepTimeMillis = intervalDurationMillis / reportInterval;
                 int stepsRemaining = numSteps - (step + 1); double estimatedRemainingMillis = stepsRemaining * avgStepTimeMillis;
                 double avgBuildNs = (double)treeBuildTimeTotal / (step + 1); double avgForceNs = (double)forceCalcTimeTotal / (step + 1); double avgUpdateNs = (double)updateTimeTotal / (step + 1);
                 System.out.printf("  Step %d/%d. Total: %.2fs. Avg step(last %d): %.1fms [B:%.1f, F:%.1f, U:%.1f]ms. Est. rem: %s\n",
                    step + 1, numSteps, timeElapsedTotal, reportInterval, avgStepTimeMillis, avgBuildNs/1e6, avgForceNs/1e6, avgUpdateNs/1e6, formatDuration((long)estimatedRemainingMillis));
                 lastReportTime = currentTime;
            }
        }
        // ... (Final timing printout same) ...
        long simEndTime = System.currentTimeMillis(); System.out.printf("Simulation loop complete. Total time: %s\n", formatDuration(simEndTime-simStartTime)); System.out.printf("Avg time/step: %.3f ms\n", (double)(simEndTime-simStartTime)/numSteps); System.out.printf("  Avg Build: %.3f ms (%.1f%%)\n", (double)treeBuildTimeTotal/numSteps/1e6, (double)treeBuildTimeTotal/(simEndTime-simStartTime)/1e4); System.out.printf("  Avg Force: %.3f ms (%.1f%%)\n", (double)forceCalcTimeTotal/numSteps/1e6, (double)forceCalcTimeTotal/(simEndTime-simStartTime)/1e4); System.out.printf("  Avg Update: %.3f ms (%.1f%%)\n", (double)updateTimeTotal/numSteps/1e6, (double)updateTimeTotal/(simEndTime-simStartTime)/1e4);


        // 4. Calculate Final Energy (Direct, Parallel)
        System.out.println("Calculating final energy (direct summation, parallel)...");
        energyStartTime = System.currentTimeMillis();
        double finalEnergy = calculateTotalEnergyDirectParallel();
        energyEndTime = System.currentTimeMillis();
        System.out.printf("Final Total Energy:   %.8e (Calculation time: %s)\n", finalEnergy, formatDuration(energyEndTime - energyStartTime));

        // 5. Verify Energy Conservation (same as before)
        // ... (Energy check printout same) ...
        System.out.println("-------------------- Energy Check --------------------"); System.out.printf("Initial Energy:        %.8e\n", initialEnergy); System.out.printf("Final Energy:          %.8e\n", finalEnergy);
        if(Math.abs(initialEnergy)<1e-15) { System.out.printf("Energy Change (Absolute): %.4e\n", finalEnergy-initialEnergy); }
        else { double relativeEnergyChange = Math.abs((finalEnergy-initialEnergy)/initialEnergy); System.out.printf("Relative Energy Change: %.4e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange*100.0);
            if(relativeEnergyChange < 5e-3) System.out.printf("Energy conservation reasonable for Barnes-Hut (theta=%.2f).\n", THETA);
            else if(relativeEnergyChange < 5e-2) System.out.println("Noticeable energy drift; expected with approximate methods.");
            else System.out.println("WARNING: Significant energy drift detected!");
        }
        System.out.println("----------------------------------------------------");

        System.out.println("Optimized k-D Tree Simulation Finished.");
    }
}