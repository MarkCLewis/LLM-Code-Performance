use nbody_simulation_kd::array_kd_tree;
use nbody_simulation_kd::array_particle;
use std::time::Instant;

fn main() {
    println!("N-Body Simulation Benchmark");
    println!("---------------------------");
    
    // Test parameters
    let particle_counts = [100, 1000, 10000];
    let steps = 10;
    let dt = 1e-3;
    
    println!("Running {} steps with dt={}", steps, dt);
    println!("\nParticle Count | Implementation | Runtime (s) | Avg Step Time (s) | Energy Change (%)");
    println!("-------------|----------------|------------|------------------|---------------");
    
    for &n in &particle_counts {
        // Test original implementation
        let mut bodies = array_particle::circular_orbits(n);
        let start = Instant::now();
        array_kd_tree::simple_sim(&mut bodies, dt, steps);
        let elapsed = start.elapsed();
        let runtime = elapsed.as_nanos() as f64 / 1e9;
        let avg_step_time = runtime / steps as f64;
        
        // Get the energy change from the last line of output
        // (We're relying on the output format from simple_sim)

        println!("{:13} | {:14} | {:10.4} | {:18.6} | See above", 
                 n, "Original AoS", runtime, avg_step_time);
        
        // Test optimized implementation
        let mut system = array_particle::circular_orbits_soa(n);
        let start = Instant::now();
        array_kd_tree::simple_sim_soa(&mut system, dt, steps);
        let elapsed = start.elapsed();
        let runtime = elapsed.as_nanos() as f64 / 1e9;
        let avg_step_time = runtime / steps as f64;
        
        println!("{:13} | {:14} | {:10.4} | {:18.6} | See above", 
                 n, "Optimized SoA", runtime, avg_step_time);
        
        println!("-------------|----------------|------------|------------------|---------------");
    }
} 