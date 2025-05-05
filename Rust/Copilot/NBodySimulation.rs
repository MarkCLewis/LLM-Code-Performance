use std::f64::consts::PI;

const G: f64 = 6.67430e-11; // Gravitational constant
const DT: f64 = 1e-3; // Time step
const NUM_BODIES: usize = 1000000; // Number of bodies

#[derive(Clone, Copy)]
struct Body {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    mass: f64,
}

fn initialize_bodies() -> Vec<Body> {
    let mut bodies = Vec::with_capacity(NUM_BODIES);
    bodies.push(Body {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: 0.0,
        vz: 0.0,
        mass: 1e30, // Central body mass
    });

    for i in 1..NUM_BODIES {
        let angle = 2.0 * PI * i as f64 / (NUM_BODIES - 1) as f64;
        bodies.push(Body {
            x: angle.cos() * 1e11,
            y: angle.sin() * 1e11,
            z: 0.0,
            vx: -angle.sin() * (G * bodies[0].mass / 1e11).sqrt(),
            vy: angle.cos() * (G * bodies[0].mass / 1e11).sqrt(),
            vz: 0.0,
            mass: 1e24, // Small body mass
        });
    }
    bodies
}

fn calculate_energy(bodies: &[Body]) -> f64 {
    let mut energy = 0.0;
    for body in bodies {
        let kinetic = 0.5 * body.mass * (body.vx.powi(2) + body.vy.powi(2) + body.vz.powi(2));
        let mut potential = 0.0;
        for other in bodies {
            if body != other {
                let dx = body.x - other.x;
                let dy = body.y - other.y;
                let dz = body.z - other.z;
                let distance = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
                potential -= G * body.mass * other.mass / distance;
            }
        }
        energy += kinetic + 0.5 * potential;
    }
    energy
}

fn kick_step(bodies: &mut [Body]) {
    for i in 0..bodies.len() {
        let mut ax = 0.0;
        let mut ay = 0.0;
        let mut az = 0.0;
        for j in 0..bodies.len() {
            if i != j {
                let dx = bodies[j].x - bodies[i].x;
                let dy = bodies[j].y - bodies[i].y;
                let dz = bodies[j].z - bodies[i].z;
                let distance = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
                let force = G * bodies[j].mass / distance.powi(3);
                ax += force * dx;
                ay += force * dy;
                az += force * dz;
            }
        }
        bodies[i].vx += ax * DT;
        bodies[i].vy += ay * DT;
        bodies[i].vz += az * DT;
    }

    for body in bodies.iter_mut() {
        body.x += body.vx * DT;
        body.y += body.vy * DT;
        body.z += body.vz * DT;
    }
}

fn main() {
    let mut bodies = initialize_bodies();

    let initial_energy = calculate_energy(&bodies);
    println!("Initial energy: {}", initial_energy);

    for _ in 0..1000 {
        kick_step(&mut bodies);
    }

    let final_energy = calculate_energy(&bodies);
    println!("Final energy: {}", final_energy);
}
