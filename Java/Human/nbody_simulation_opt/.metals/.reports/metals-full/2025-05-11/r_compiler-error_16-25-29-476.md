file://<WORKSPACE>/src/main/java/PerformanceTest.java
### java.util.NoSuchElementException: next on empty iterator

occurred in the presentation compiler.

presentation compiler configuration:


action parameters:
uri: file://<WORKSPACE>/src/main/java/PerformanceTest.java
text:
```scala
/**
 * Performance and accuracy test for the n-body simulation implementations.
 * This class compares the original implementation with the optimized one.
 */
public class PerformanceTest {
    public static void main(String[] args) {
        // Parameters for the test
        int numParticles = 10000;
        int numSteps = 10;
        double dt = 1e-3;
        
        // Run each implementation and measure performance and accuracy
        runTest("original", numParticles, numSteps, dt);
        runTest("optimized", numParticles, numSteps, dt);
    }
    
    private static void runTest(String implementation, int numParticles, int numSteps, double dt) {
        System.out.println("\n=== Testing " + implementation + " implementation with " + 
                          numParticles + " particles for " + numSteps + " steps ===");
        
        JSystem system;
        long startTime = System.nanoTime();
        
        if (implementation.equals("original")) {
            // Create system with original implementation
            system = new JClassSystem(numParticles + 1);
            initializeSystem(system, numParticles);
            
            // Measure initial energy
            double initialEnergy = JKDTree.calculateTotalEnergy(system);
            System.out.println("Initial energy: " + initialEnergy);
            
            // Run simulation
            JKDTree.simple_sim(system, dt, numSteps);
            
            // Measure final energy
            double finalEnergy = JKDTree.calculateTotalEnergy(system);
            double relativeError = Math.abs((finalEnergy - initialEnergy) / initialEnergy);
            System.out.println("Final energy: " + finalEnergy);
            System.out.println("Relative energy error: " + relativeError);
        } else {
            // Create system with optimized implementation
            system = new JOptimizedArraySystem(numParticles + 1);
            initializeSystem(system, numParticles);
            
            // Measure initial energy
            double initialEnergy = OptimizedKDTree.calculateTotalEnergy(system);
            System.out.println("Initial energy: " + initialEnergy);
            
            // Run simulation
            OptimizedKDTree.optimizedSim(system, dt, numSteps);
            
            // Measure final energy is done inside the optimizedSim method
        }
        
        long endTime = System.nanoTime();
        double elapsedSeconds = (endTime - startTime) / 1_000_000_000.0;
        System.out.println("Total execution time: " + elapsedSeconds + " seconds");
    }
    
    /**
     * Initialize the system with the circular orbit test case
     */
    private static void initializeSystem(JSystem system, int n) {
        var rand = java.util.random.RandomGenerator.getDefault();
        
        // Central body
        system.init(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0);
        
        // Orbiting bodies
        for (int i = 0; i < n; ++i) {
            double d = 0.1 + (i * 5.0 / n);
            double v = Math.sqrt(1.0 / d);
            double theta = rand.nextDouble(0.0, 6.28);
            double x = d * Math.cos(theta);
            double y = d * Math.sin(theta);
            double vx = -v * Math.sin(theta);
            double vy = v * Math.cos(theta);
            system.init(i+1, 
                x, y, 0.0,
                vx, vy, 0.0,
                1e-14,
                1e-7
            );
        }
    }
} 
```



#### Error stacktrace:

```
scala.collection.Iterator$$anon$19.next(Iterator.scala:973)
	scala.collection.Iterator$$anon$19.next(Iterator.scala:971)
	scala.collection.mutable.MutationTracker$CheckedIterator.next(MutationTracker.scala:76)
	scala.collection.IterableOps.head(Iterable.scala:222)
	scala.collection.IterableOps.head$(Iterable.scala:222)
	scala.collection.AbstractIterable.head(Iterable.scala:935)
	dotty.tools.dotc.interactive.InteractiveDriver.run(InteractiveDriver.scala:164)
	dotty.tools.pc.CachingDriver.run(CachingDriver.scala:45)
	dotty.tools.pc.WithCompilationUnit.<init>(WithCompilationUnit.scala:31)
	dotty.tools.pc.SimpleCollector.<init>(PcCollector.scala:351)
	dotty.tools.pc.PcSemanticTokensProvider$Collector$.<init>(PcSemanticTokensProvider.scala:63)
	dotty.tools.pc.PcSemanticTokensProvider.Collector$lzyINIT1(PcSemanticTokensProvider.scala:63)
	dotty.tools.pc.PcSemanticTokensProvider.Collector(PcSemanticTokensProvider.scala:63)
	dotty.tools.pc.PcSemanticTokensProvider.provide(PcSemanticTokensProvider.scala:88)
	dotty.tools.pc.ScalaPresentationCompiler.semanticTokens$$anonfun$1(ScalaPresentationCompiler.scala:111)
```
#### Short summary: 

java.util.NoSuchElementException: next on empty iterator