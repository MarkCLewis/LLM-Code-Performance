mod array_particle;
mod array_kd_tree;
mod quickstat;

use clap::Parser;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of particules to generate.
    #[arg(short, long)]
    number: usize,

    /// Number of steps to run the simulation.
    #[arg(short, long, default_value_t = 1)]
    steps: i64,

    /// Use optimized SoA implementation.
    #[arg(short, long, default_value_t = false)]
    optimized: bool,
}

fn main() {
    let args = Args::parse();
    let dt = 1e-3;

    if args.optimized {
        // Run optimized SoA implementation
        println!("Running optimized SoA implementation with {} particles for {} steps", args.number, args.steps);
        let mut system = array_particle::circular_orbits_soa(args.number);
        
        let start = Instant::now();
        array_kd_tree::simple_sim_soa(&mut system, dt, args.steps);
        let elapsed = start.elapsed();
        
        println!("Execution time: {} seconds", elapsed.as_nanos() as f64 / 1e9);
        println!("Average step time: {} seconds", elapsed.as_nanos() as f64 / (1e9 * args.steps as f64));
    } else {
        // Run original implementation
        println!("Running original AoS implementation with {} particles for {} steps", args.number, args.steps);
        let mut bodies = array_particle::circular_orbits(args.number);
        
        let start = Instant::now();
        array_kd_tree::simple_sim(&mut bodies, dt, args.steps);
        let elapsed = start.elapsed();
        
        println!("Execution time: {} seconds", elapsed.as_nanos() as f64 / 1e9);
        println!("Average step time: {} seconds", elapsed.as_nanos() as f64 / (1e9 * args.steps as f64));
    }
}
