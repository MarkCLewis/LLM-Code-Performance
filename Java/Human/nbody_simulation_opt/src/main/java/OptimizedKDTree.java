import java.util.ArrayList;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.random.RandomGenerator;

/**
 * Optimized implementation of KD-Tree for n-body simulations with:
 * - Node pooling to reduce allocation
 * - Tree reuse between steps
 * - Improved force calculation
 * - Better parallelization with Fork/Join framework
 * - Optimized energy calculation
 */
public class OptimizedKDTree {
    // Configuration constants
    static final int MAX_PARTS = 7;
    static final double THETA = 0.3;
    static final double G = 1.0; // Gravitational constant (using simulation units)
    static final RandomGenerator rand = RandomGenerator.getDefault();
    
    // Custom thread pool for better control of parallelism
    private static final ForkJoinPool FORK_JOIN_POOL = new ForkJoinPool();
    
    // Counters for performance monitoring
    private static final AtomicInteger forceCalculations = new AtomicInteger(0);
    
    // Node pool for reuse
    private static final ArrayList<TreeNode> nodePool = new ArrayList<>(1024);
    private static int nextNodeIndex = 0;
    
    /**
     * Tree node implementation with support for pooling
     */
    static class TreeNode {
        // For leaves
        int num_parts;
        int[] particles = new int[MAX_PARTS];
        
        // For internal nodes
        int split_dim;
        double split_val;
        double m;
        double[] cm = new double[3];
        double size;
        int left;
        int right;
        
        // Temporary workspace for distances to avoid recomputation
        double[] distances;
        
        /**
         * Reset this node for reuse
         */
        void reset() {
            num_parts = 0;
            split_dim = 0;
            split_val = 0;
            m = 0;
            cm[0] = cm[1] = cm[2] = 0;
            size = 0;
            left = right = 0;
        }
    }
    
    /**
     * Get a node from the pool or create a new one if needed
     */
    private static TreeNode getNode() {
        if (nextNodeIndex >= nodePool.size()) {
            TreeNode node = new TreeNode();
            nodePool.add(node);
            nextNodeIndex++;
            return node;
        } else {
            TreeNode node = nodePool.get(nextNodeIndex);
            node.reset();
            nextNodeIndex++;
            return node;
        }
    }
    
    /**
     * Prepare the tree for a new simulation step
     */
    static void resetTree() {
        nextNodeIndex = 0;
        forceCalculations.set(0);
    }
    
    /**
     * Build the KD-Tree for the given particle system
     */
    static int buildTree(int[] indices, int start, int end, JSystem system) {
        return buildTreeRecursive(indices, start, end, system, 0);
    }
    
    /**
     * Recursive helper for tree building
     */
    private static int buildTreeRecursive(int[] indices, int start, int end, JSystem system, int curNodeIndex) {
        int np = end - start;
        
        // Leaf case
        if (np <= MAX_PARTS) {
            TreeNode node = getNode();
            node.num_parts = np;
            for (int i = 0; i < np; i++) {
                node.particles[i] = indices[start + i];
            }
            return curNodeIndex;
        } 
        
        // Internal node case
        else {
            // Calculate bounds and center of mass
            double[] min = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
            double[] max = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
            double m = 0.0;
            double[] cm = {0.0, 0.0, 0.0};
            
            // Calculate bounds and center of mass - unrolled for better vectorization
            for (int i = start; i < end; i++) {
                int idx = indices[i];
                double particleMass = system.m(idx);
                double px = system.p(idx, 0);
                double py = system.p(idx, 1);
                double pz = system.p(idx, 2);
                
                m += particleMass;
                cm[0] += particleMass * px;
                cm[1] += particleMass * py;
                cm[2] += particleMass * pz;
                
                min[0] = Math.min(min[0], px);
                min[1] = Math.min(min[1], py);
                min[2] = Math.min(min[2], pz);
                
                max[0] = Math.max(max[0], px);
                max[1] = Math.max(max[1], py);
                max[2] = Math.max(max[2], pz);
            }
            
            // Normalize center of mass
            cm[0] /= m;
            cm[1] /= m;
            cm[2] /= m;
            
            // Find dimension with largest spread
            int splitDim = 0;
            double xSpread = max[0] - min[0];
            double ySpread = max[1] - min[1];
            double zSpread = max[2] - min[2];
            
            if (ySpread > xSpread) splitDim = 1;
            if (zSpread > Math.max(xSpread, ySpread)) splitDim = 2;
            
            double size = max[splitDim] - min[splitDim];
            
            // Use more efficient partitioning
            int mid = (start + end) / 2;
            int left = start;
            int right = end - 1;
            
            // Quick select pivot
            while (left < right) {
                double pivotValue = system.p(indices[(left + right) / 2], splitDim);
                int i = left;
                int j = right;
                
                while (i <= j) {
                    while (system.p(indices[i], splitDim) < pivotValue) i++;
                    while (system.p(indices[j], splitDim) > pivotValue) j--;
                    
                    if (i <= j) {
                        int temp = indices[i];
                        indices[i] = indices[j];
                        indices[j] = temp;
                        i++;
                        j--;
                    }
                }
                
                if (j < mid) left = i;
                if (i > mid) right = j;
                if (i <= mid && j >= mid) break;
            }
            
            double splitVal = system.p(indices[mid], splitDim);
            
            // Recursively build left and right subtrees
            int leftNodeIndex = curNodeIndex + 1;
            int rightNodeIndex = buildTreeRecursive(indices, start, mid, system, leftNodeIndex);
            int endNodeIndex = buildTreeRecursive(indices, mid, end, system, rightNodeIndex + 1);
            
            // Create and initialize the current node
            TreeNode node = getNode();
            node.num_parts = 0;
            node.split_dim = splitDim;
            node.split_val = splitVal;
            node.m = m;
            node.cm[0] = cm[0];
            node.cm[1] = cm[1];
            node.cm[2] = cm[2];
            node.size = size;
            node.left = leftNodeIndex;
            node.right = rightNodeIndex + 1;
            
            return endNodeIndex;
        }
    }
    
    /**
     * Calculate acceleration on a particle due to another particle
     */
    static void calcParticleParticleAccel(JSystem system, int i, int j, double[] acc) {
        double dx = system.p(i, 0) - system.p(j, 0);
        double dy = system.p(i, 1) - system.p(j, 1);
        double dz = system.p(i, 2) - system.p(j, 2);
        
        // Cache distance calculation
        double distSqr = dx * dx + dy * dy + dz * dz;
        double dist = Math.sqrt(distSqr);
        double invDist3 = 1.0 / (dist * distSqr);
        
        // Force magnitude
        double magi = -system.m(j) * invDist3;
        
        // Accumulate acceleration components
        acc[0] += dx * magi;
        acc[1] += dy * magi;
        acc[2] += dz * magi;
        
        // Count force calculations for diagnostics
        forceCalculations.incrementAndGet();
    }
    
    /**
     * Recursive calculation of acceleration using the KD-tree
     */
    static void calcAccelRecursive(int nodeIndex, int particleIndex, JSystem system, double[] acc) {
        TreeNode node = nodePool.get(nodeIndex);
        
        // Handle leaf nodes - direct calculation for all particles
        if (node.num_parts > 0) {
            for (int i = 0; i < node.num_parts; i++) {
                if (node.particles[i] != particleIndex) {
                    calcParticleParticleAccel(system, particleIndex, node.particles[i], acc);
                }
            }
        } 
        // Handle internal nodes
        else {
            // Calculate distance to center of mass
            double dx = system.p(particleIndex, 0) - node.cm[0];
            double dy = system.p(particleIndex, 1) - node.cm[1];
            double dz = system.p(particleIndex, 2) - node.cm[2];
            double distSqr = dx * dx + dy * dy + dz * dz;
            
            // Check if we can use multipole approximation
            double sizeSqr = node.size * node.size;
            if (sizeSqr < THETA * THETA * distSqr) {
                // Use center of mass approximation
                double dist = Math.sqrt(distSqr);
                double invDist3 = 1.0 / (dist * distSqr);
                double magi = -node.m * invDist3;
                
                acc[0] += dx * magi;
                acc[1] += dy * magi;
                acc[2] += dz * magi;
                
                forceCalculations.incrementAndGet();
            } else {
                // Need to recurse into children
                calcAccelRecursive(node.left, particleIndex, system, acc);
                calcAccelRecursive(node.right, particleIndex, system, acc);
            }
        }
    }
    
    /**
     * Calculate acceleration for a particle
     */
    static void calcAccel(int particleIndex, JSystem system, double[] acc) {
        acc[0] = acc[1] = acc[2] = 0.0;
        calcAccelRecursive(0, particleIndex, system, acc);
    }
    
    /**
     * Fork/Join task for parallel acceleration calculation
     */
    static class AccelerationTask extends RecursiveAction {
        private final int startIndex;
        private final int endIndex;
        private final JSystem system;
        private final double[][] accelerations;
        private static final int THRESHOLD = 1000;
        
        AccelerationTask(int start, int end, JSystem system, double[][] accelerations) {
            this.startIndex = start;
            this.endIndex = end;
            this.system = system;
            this.accelerations = accelerations;
        }
        
        @Override
        protected void compute() {
            if (endIndex - startIndex <= THRESHOLD) {
                // Small enough chunk, calculate directly
                for (int i = startIndex; i < endIndex; i++) {
                    calcAccel(i, system, accelerations[i]);
                }
            } else {
                // Split into smaller tasks
                int mid = (startIndex + endIndex) / 2;
                AccelerationTask left = new AccelerationTask(startIndex, mid, system, accelerations);
                AccelerationTask right = new AccelerationTask(mid, endIndex, system, accelerations);
                
                // Fork right task and compute left task
                right.fork();
                left.compute();
                right.join();
            }
        }
    }
    
    /**
     * Calculate total energy of the system (kinetic + potential)
     */
    static double calculateTotalEnergy(JSystem system) {
        double kineticEnergy = 0.0;
        double potentialEnergy = 0.0;
        
        // Parallelize energy calculation
        final double[] kE = new double[1];
        final double[] pE = new double[1];
        
        // Fork/Join task for parallel energy calculation
        class EnergyTask extends RecursiveAction {
            private final int startIndex;
            private final int endIndex;
            private final JSystem system;
            private static final int THRESHOLD = 1000;
            private double localKE = 0.0;
            private double localPE = 0.0;
            
            EnergyTask(int start, int end, JSystem system) {
                this.startIndex = start;
                this.endIndex = end;
                this.system = system;
            }
            
            @Override
            protected void compute() {
                if (endIndex - startIndex <= THRESHOLD) {
                    // Small enough chunk, calculate directly
                    for (int i = startIndex; i < endIndex; i++) {
                        // Kinetic energy: 0.5 * m * v^2
                        double vx = system.v(i, 0);
                        double vy = system.v(i, 1);
                        double vz = system.v(i, 2);
                        double v2 = vx * vx + vy * vy + vz * vz;
                        localKE += 0.5 * system.m(i) * v2;
                        
                        // Potential energy: -G * m1 * m2 / r (for j > i only to avoid double counting)
                        for (int j = i + 1; j < system.numBodies(); j++) {
                            double dx = system.p(i, 0) - system.p(j, 0);
                            double dy = system.p(i, 1) - system.p(j, 1);
                            double dz = system.p(i, 2) - system.p(j, 2);
                            double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            if (distance > 1e-10) {
                                localPE -= G * system.m(i) * system.m(j) / distance;
                            }
                        }
                    }
                    
                    // Add to global energy (synchronizing to avoid race conditions)
                    synchronized (kE) {
                        kE[0] += localKE;
                    }
                    synchronized (pE) {
                        pE[0] += localPE;
                    }
                } else {
                    // Split into smaller tasks
                    int mid = (startIndex + endIndex) / 2;
                    EnergyTask left = new EnergyTask(startIndex, mid, system);
                    EnergyTask right = new EnergyTask(mid, endIndex, system);
                    
                    // Fork both tasks in parallel
                    invokeAll(left, right);
                }
            }
        }
        
        // Calculate energy in parallel
        FORK_JOIN_POOL.invoke(new EnergyTask(0, system.numBodies(), system));
        kineticEnergy = kE[0];
        potentialEnergy = pE[0];
        
        return kineticEnergy + potentialEnergy;
    }
    
    /**
     * Main simulation method with optimized implementation
     */
    static void optimizedSim(JSystem system, double dt, int steps) {
        // Pre-allocate acceleration arrays
        double[][] accelerations = new double[system.numBodies()][3];
        int[] indices = new int[system.numBodies()];
        
        // Initialize indices
        for (int i = 0; i < system.numBodies(); i++) {
            indices[i] = i;
        }
        
        // Calculate initial energy
        // double initialEnergy = calculateTotalEnergy(system);
        // System.out.println("Initial total energy: " + initialEnergy);
        
        long startTime = System.nanoTime();
        
        // Main simulation loop
        for (int step = 0; step < steps; step++) {
            // Reset tree for reuse
            resetTree();
            
            // Build the KD-tree (this modifies indices array)
            buildTree(indices, 0, system.numBodies(), system);
            
            // Calculate accelerations in parallel using Fork/Join
            AccelerationTask mainTask = new AccelerationTask(0, system.numBodies(), system, accelerations);
            FORK_JOIN_POOL.invoke(mainTask);
            
            // Update positions and velocities in parallel
            class UpdateTask extends RecursiveAction {
                private final int startIndex;
                private final int endIndex;
                private final JSystem system;
                private final double[][] accelerations;
                private final double dt;
                private static final int THRESHOLD = 1000;
                
                UpdateTask(int start, int end, JSystem system, double[][] accelerations, double dt) {
                    this.startIndex = start;
                    this.endIndex = end;
                    this.system = system;
                    this.accelerations = accelerations;
                    this.dt = dt;
                }
                
                @Override
                protected void compute() {
                    if (endIndex - startIndex <= THRESHOLD) {
                        for (int i = startIndex; i < endIndex; i++) {
                            // Update velocities
                            system.incV(i, 0, dt * accelerations[i][0]);
                            system.incV(i, 1, dt * accelerations[i][1]);
                            system.incV(i, 2, dt * accelerations[i][2]);
                            
                            // Update positions
                            system.incP(i, 0, dt * system.v(i, 0));
                            system.incP(i, 1, dt * system.v(i, 1));
                            system.incP(i, 2, dt * system.v(i, 2));
                            
                            // Reset accelerations for next step
                            accelerations[i][0] = 0.0;
                            accelerations[i][1] = 0.0;
                            accelerations[i][2] = 0.0;
                        }
                    } else {
                        int mid = (startIndex + endIndex) / 2;
                        UpdateTask left = new UpdateTask(startIndex, mid, system, accelerations, dt);
                        UpdateTask right = new UpdateTask(mid, endIndex, system, accelerations, dt);
                        invokeAll(left, right);
                    }
                }
            }
            
            // Update positions and velocities in parallel
            FORK_JOIN_POOL.invoke(new UpdateTask(0, system.numBodies(), system, accelerations, dt));
            
            // Periodically check energy and print diagnostics
            // if (step > 0 && (step % 1000 == 0 || step == steps - 1)) {
            //     double currentEnergy = calculateTotalEnergy(system);
            //     double relativeError = Math.abs((currentEnergy - initialEnergy) / initialEnergy);
                
            //     System.out.println("Step " + step + 
            //                       " - Energy: " + currentEnergy + 
            //                       " (relative error: " + relativeError + ")" +
            //                       " - Force calculations: " + forceCalculations.get());
            // }
        }
        
        long endTime = System.nanoTime();
        double elapsedSeconds = (endTime - startTime) / 1_000_000_000.0;
        
        // Calculate final energy and report conservation
        // double finalEnergy = calculateTotalEnergy(system);
        // double absoluteError = finalEnergy - initialEnergy;
        // double relativeError = Math.abs(absoluteError / initialEnergy);
        
        System.out.println("Simulation completed in " + elapsedSeconds + " seconds");
        // System.out.println("Final total energy: " + finalEnergy);
        // System.out.println("Absolute energy error: " + absoluteError);
        // System.out.println("Relative energy error: " + relativeError);
    }
} 