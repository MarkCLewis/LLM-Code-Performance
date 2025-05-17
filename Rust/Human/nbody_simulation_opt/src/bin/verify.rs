use nbody_simulation_kd::array_kd_tree;
use nbody_simulation_kd::array_particle::{self, ParticleSystem, Particle};
use std::time::Instant;

// Function to compare positions between original and optimized versions
fn verify_positions(original: &[Particle], optimized: &ParticleSystem, tolerance: f64) -> bool {
    if original.len() != optimized.count {
        println!("ERROR: Particle count mismatch: {} vs {}", original.len(), optimized.count);
        return false;
    }
    
    let mut max_diff = 0.0;
    let mut avg_diff = 0.0;
    let mut diff_count = 0;
    
    for i in 0..original.len() {
        for d in 0..3 {
            let diff = (original[i].p[d] - optimized.positions[i][d]).abs();
            max_diff = f64::max(max_diff, diff);
            avg_diff += diff;
            
            if diff > tolerance {
                diff_count += 1;
                println!("Large position difference at particle {}, dimension {}: Original={}, Optimized={}, Diff={}",
                    i, d, original[i].p[d], optimized.positions[i][d], diff);
                
                if diff_count >= 10 {
                    println!("Too many differences, stopping comparison");
                    return false;
                }
            }
        }
    }
    
    avg_diff /= (original.len() * 3) as f64;
    
    println!("Position verification results:");
    println!("  - Maximum difference: {}", max_diff);
    println!("  - Average difference: {}", avg_diff);
    println!("  - Differences above tolerance: {}", diff_count);
    
    diff_count == 0
}

fn verify_energy(original_energy: f64, optimized_energy: f64, tolerance: f64) -> bool {
    let diff = (original_energy - optimized_energy).abs();
    let rel_diff = diff / f64::max(original_energy.abs(), optimized_energy.abs());
    
    println!("Energy verification results:");
    println!("  - Original energy: {}", original_energy);
    println!("  - Optimized energy: {}", optimized_energy);
    println!("  - Absolute difference: {}", diff);
    println!("  - Relative difference: {:.6}%", rel_diff * 100.0);
    
    rel_diff <= tolerance
}

fn main() {
    println!("N-Body Simulation Verification");
    println!("------------------------------");
    
    // Test parameters
    let n = 1000;
    let steps = 10;
    let dt = 1e-3;
    let position_tolerance = 1e-10;
    let energy_tolerance = 1e-6;
    
    println!("Running verification with {} particles for {} steps", n, steps);
    println!("Position tolerance: {}", position_tolerance);
    println!("Energy tolerance: {}%", energy_tolerance * 100.0);
    
    // Create initial particles with fixed seed for reproducibility
    fastrand::seed(12345);
    let original_bodies = array_particle::circular_orbits(n);
    
    // Reset seed to get identical positions in SoA version
    fastrand::seed(12345);
    let mut system_optimized = array_particle::circular_orbits_soa(n);
    
    // Make a copy of the original bodies for simulation
    let mut bodies_original = Vec::new();
    for p in &original_bodies {
        bodies_original.push(Particle {
            p: p.p,
            v: p.v,
            r: p.r,
            m: p.m
        });
    }
    
    let original_initial_energy = array_particle::calc_total_energy(&bodies_original);
    let optimized_initial_energy = array_particle::calc_total_energy_soa(&system_optimized);
    
    // Verify initial conditions match
    println!("\nVerifying initial conditions...");
    let initial_positions_ok = verify_positions(&bodies_original, &system_optimized, 1e-15);
    let initial_energy_ok = verify_energy(original_initial_energy, optimized_initial_energy, 1e-15);
    
    if !initial_positions_ok || !initial_energy_ok {
        println!("ERROR: Initial conditions don't match! Cannot continue verification.");
        return;
    }
    
    // Run original implementation
    let start = Instant::now();
    array_kd_tree::simple_sim(&mut bodies_original, dt, steps);
    let original_runtime = start.elapsed().as_nanos() as f64 / 1e9;
    let original_final_energy = array_particle::calc_total_energy(&bodies_original);
    
    // Run optimized implementation
    let start = Instant::now();
    array_kd_tree::simple_sim_soa(&mut system_optimized, dt, steps);
    let optimized_runtime = start.elapsed().as_nanos() as f64 / 1e9;
    let optimized_final_energy = array_particle::calc_total_energy_soa(&system_optimized);
    
    // Verify final positions
    println!("\nVerifying final positions...");
    let positions_ok = verify_positions(&bodies_original, &system_optimized, position_tolerance);
    
    // Verify energy conservation
    println!("\nVerifying energy conservation...");
    let original_energy_change = (original_final_energy - original_initial_energy) / original_initial_energy.abs();
    let optimized_energy_change = (optimized_final_energy - optimized_initial_energy) / optimized_initial_energy.abs();
    
    println!("Original energy change: {:.6}%", original_energy_change * 100.0);
    println!("Optimized energy change: {:.6}%", optimized_energy_change * 100.0);
    
    let energy_ok = verify_energy(original_final_energy, optimized_final_energy, energy_tolerance);
    
    // Performance comparison
    println!("\nPerformance comparison:");
    println!("  - Original runtime: {:.6} seconds", original_runtime);
    println!("  - Optimized runtime: {:.6} seconds", optimized_runtime);
    println!("  - Speedup: {:.2}x", original_runtime / optimized_runtime);
    
    // Final verdict
    if positions_ok && energy_ok {
        println!("\nVERIFICATION PASSED: Optimized version produces correct results!");
    } else {
        println!("\nVERIFICATION FAILED: Optimizations may have introduced errors.");
        
        if !positions_ok {
            println!("  - Position differences exceed tolerance");
        }
        
        if !energy_ok {
            println!("  - Energy differences exceed tolerance");
        }
    }
} 