import java.util.random.RandomGenerator;

public class JavaMain {
  public static void main(String[] args) {
    if (args.length < 2) {
      System.out.println("Specify a number of particles and number of steps.");
      return;
    }
    System.out.println("Running optimized n-body simulation.");
  
    int steps = Integer.parseInt(args[0]);
    int n = Integer.parseInt(args[1]);
    
    // Parse optional type argument
    String systemType = args.length > 2 ? args[2] : "optimized";
  
    double dt = 1e-3; // * 2.0 * std::f64::consts::PI;
  
    long startTime = System.nanoTime();
    
    // Create the appropriate system type based on command line argument
    JSystem system;
    if ("classic".equals(systemType)) {
      system = circular_orbits_classic(n);
      JKDTree.simple_sim(system, dt, steps);
    } else if ("array".equals(systemType)) {
      system = circular_orbits_array(n);
      JKDTree.simple_sim(system, dt, steps);
    } else {
      // Default to optimized implementation
      system = circular_orbits_optimized(n);
      OptimizedKDTree.optimizedSim(system, dt, steps);
    }
    
    long endTime = System.nanoTime();
    double elapsedSeconds = (endTime - startTime) / 1_000_000_000.0;
    System.out.println("Total execution time: " + elapsedSeconds + " seconds");
  }

  static JSystem circular_orbits_classic(int n) {
    var rand = RandomGenerator.getDefault();
    var system = new JClassSystem(n+1);
    system.init(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0 );
  
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
    return system;
  }
  
  static JSystem circular_orbits_array(int n) {
    var rand = RandomGenerator.getDefault();
    var system = new JArraySystem(n+1);
    system.init(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0 );
  
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
    return system;
  }
  
  static JSystem circular_orbits_optimized(int n) {
    var rand = RandomGenerator.getDefault();
    var system = new JOptimizedArraySystem(n+1);
    system.init(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0 );
  
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
    return system;
  }
}
