use std::f64::consts::PI;
use std::time::Instant;

const GRAVITATIONAL_CONSTANT: f64 = 6.67430e-11;

#[derive(Clone, Debug)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    fn zero() -> Self {
        Vec3::new(0.0, 0.0, 0.0)
    }

    fn distance_squared(&self, other: &Self) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)
    }

    fn distance(&self, other: &Self) -> f64 {
        self.distance_squared(other).sqrt()
    }

    fn sub(&self, other: &Self) -> Self {
        Vec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    fn add(&self, other: &Self) -> Self {
        Vec3::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    fn mul(&self, scalar: f64) -> Self {
        Vec3::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }

    fn div(&self, scalar: f64) -> Self {
        Vec3::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }

    fn norm(&self) -> f64 {
        self.distance(&Vec3::zero())
    }
}

#[derive(Clone, Debug)]
struct Body {
    mass: f64,
    position: Vec3,
    velocity: Vec3,
}

impl Body {
    fn new(mass: f64, position: Vec3, velocity: Vec3) -> Self {
        Body { mass, position, velocity }
    }
}

struct System {
    bodies: Vec<Body>,
}

impl System {
    fn new() -> Self {
        System { bodies: Vec::new() }
    }

    fn add_body(&mut self, body: Body) {
        self.bodies.push(body);
    }

    fn calculate_force(&self, body_i_index: usize, body_j_index: usize) -> Vec3 {
        let body_i = &self.bodies[body_i_index];
        let body_j = &self.bodies[body_j_index];
        let r_vec = body_j.position.sub(&body_i.position);
        let r_sq = r_vec.distance_squared(&Vec3::zero());
        if r_sq > 1e-9 {
            let r = r_sq.sqrt();
            let magnitude = (GRAVITATIONAL_CONSTANT * body_i.mass * body_j.mass) / r_sq;
            r_vec.mul(magnitude / r)
        } else {
            Vec3::zero()
        }
    }

    fn calculate_total_energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;
        let num_bodies = self.bodies.len();

        for i in 0..num_bodies {
            let body_i = &self.bodies[i];
            let v_sq = body_i.velocity.distance_squared(&Vec3::zero());
            kinetic_energy += 0.5 * body_i.mass * v_sq;

            for j in i + 1..num_bodies {
                let body_j = &self.bodies[j];
                let r = body_i.position.distance(&body_j.position);
                potential_energy -= (GRAVITATIONAL_CONSTANT * body_i.mass * body_j.mass) / r;
            }
        }
        kinetic_energy + potential_energy
    }

    fn initialize_circular_orbits(
        num_orbiting: usize,
        central_mass: f64,
        orbit_radius: f64,
        orbiting_mass: f64,
    ) -> Self {
        let mut system = System::new();

        // Initialize the central body
        system.add_body(Body::new(central_mass, Vec3::zero(), Vec3::zero()));

        // Initialize the orbiting bodies
        for i in 0..num_orbiting {
            let angle = 2.0 * PI * (i as f64) / (num_orbiting as f64);
            let x = orbit_radius * angle.cos();
            let y = orbit_radius * angle.sin();
            let z = 0.0;
            let position = Vec3::new(x, y, z);

            // Calculate the orbital velocity for a circular orbit
            let orbital_speed = (GRAVITATIONAL_CONSTANT * central_mass / orbit_radius).sqrt();
            let vx = -orbital_speed * angle.sin();
            let vy = orbital_speed * angle.cos();
            let vz = 0.0;
            let velocity = Vec3::new(vx, vy, vz);

            system.add_body(Body::new(orbiting_mass, position, velocity));
        }
        system
    }

    fn kick_step(&mut self, dt: f64) {
        let num_bodies = self.bodies.len();
        let mut forces = vec![Vec3::zero(); num_bodies];

        // Calculate forces
        for i in 0..num_bodies {
            for j in 0..num_bodies {
                if i != j {
                    let force = self.calculate_force(i, j);
                    forces[i].x += force.x;
                    forces[i].y += force.y;
                    forces[i].z += force.z;
                }
            }
        }

        // Update velocities (kick)
        for i in 0..num_bodies {
            self.bodies[i].velocity.x += (forces[i].x / self.bodies[i].mass) * (dt / 2.0);
            self.bodies[i].velocity.y += (forces[i].y / self.bodies[i].mass) * (dt / 2.0);
            self.bodies[i].velocity.z += (forces[i].z / self.bodies[i].mass) * (dt / 2.0);
        }
    }

    fn drift_step(&mut self, dt: f64) {
        for body in &mut self.bodies {
            body.position.x += body.velocity.x * dt;
            body.position.y += body.velocity.y * dt;
            body.position.z += body.velocity.z * dt;
        }
    }

    fn first_order_kick_step(&mut self, dt: f64) {
        self.kick_step(dt);
        self.drift_step(dt);
        self.kick_step(dt);
    }
}

fn main() {
    let num_orbiting_bodies = 1_000_000;
    let central_mass = 1.989e30; // Mass of the Sun (kg)
    let orbit_radius = 1.496e11; // 1 AU (m)
    let orbiting_mass = 5.972e24; // Mass of the Earth (kg)
    let num_steps = 1000;
    let time_step = 3600.0 * 24.0 * 7.0; // 1 week in seconds

    // Initialize the system
    let mut system =
        System::initialize_circular_orbits(num_orbiting_bodies, central_mass, orbit_radius, orbiting_mass);
    let initial_system = system.clone();

    println!("Initial number of bodies: {}", system.bodies.len());

    // Calculate initial energy
    let initial_energy = initial_system.calculate_total_energy();
    println!("Initial total energy: {:.e} J", initial_energy);

    println!("Running simulation for {} steps...", num_steps);
    let start_time = Instant::now();
    for step in 0..num_steps {
        system.first_order_kick_step(time_step);
        if (step + 1) % 100 == 0 {
            println!("Step {} completed.", step + 1);
        }
    }
    let elapsed_time = start_time.elapsed();
    println!("Simulation finished in {:.2} seconds.", elapsed_time.as_secs_f64());

    // Calculate final energy
    let final_energy = system.calculate_total_energy();
    println!("Final total energy: {:.e} J", final_energy);

    // Calculate the energy difference
    let energy_difference = (final_energy - initial_energy).abs();
    let relative_energy_difference = energy_difference / initial_energy.abs();
    println!("Absolute energy difference: {:.e} J", energy_difference);
    println!("Relative energy difference: {:.e}", relative_energy_difference);
}