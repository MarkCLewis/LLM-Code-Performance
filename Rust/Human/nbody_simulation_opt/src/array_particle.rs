#[derive(Clone)]
pub struct Particle {
  pub p: [f64; 3],
  pub v: [f64; 3],
  pub r: f64,
  pub m: f64
}

pub struct ParticleSystem {
  pub positions: Vec<[f64; 3]>,
  pub velocities: Vec<[f64; 3]>,
  pub radii: Vec<f64>,
  pub masses: Vec<f64>,
  pub count: usize,
}

impl ParticleSystem {
  pub fn new(count: usize) -> Self {
    ParticleSystem {
      positions: vec![[0.0, 0.0, 0.0]; count],
      velocities: vec![[0.0, 0.0, 0.0]; count],
      radii: vec![0.0; count],
      masses: vec![0.0; count],
      count,
    }
  }

  pub fn from_particles(particles: Vec<Particle>) -> Self {
    let count = particles.len();
    let mut system = ParticleSystem::new(count);
    
    for (i, p) in particles.into_iter().enumerate() {
      system.positions[i] = p.p;
      system.velocities[i] = p.v;
      system.radii[i] = p.r;
      system.masses[i] = p.m;
    }
    
    system
  }

  pub fn to_particles(&self) -> Vec<Particle> {
    let mut particles = Vec::with_capacity(self.count);
    for i in 0..self.count {
      particles.push(Particle {
        p: self.positions[i],
        v: self.velocities[i],
        r: self.radii[i],
        m: self.masses[i],
      });
    }
    particles
  }
}

pub fn two_bodies() -> Vec<Particle> {
  let mut bodies = Vec::new();
  bodies.push(Particle { p: [0.0, 0.0, 0.0], 
                         v: [0.0, 0.0, 0.0], r: 1.0, m: 1.0 });
  bodies.push(Particle { p: [1.0, 0.0, 0.0], 
                         v: [0.0, 1.0, 0.0], r: 1e-4, m: 1e-20 });
  bodies
}

pub fn circular_orbits(n: usize) -> Vec<Particle> {
  let mut particle_buf = vec![];
    particle_buf.push(Particle {
      p: [0.0, 0.0, 0.0],
      v: [0.0, 0.0, 0.0],
      r: 0.00465047,
      m: 1.0,
    });

    for i in 0..n {
        let d = 0.1 + ((i as f64) * 5.0 / (n as f64));
        let v = f64::sqrt(1.0 / d);
        let theta = fastrand::f64() * 6.28;
        let x = d * f64::cos(theta);
        let y = d * f64::sin(theta);
        let vx = -v * f64::sin(theta);
        let vy = v * f64::cos(theta);
        particle_buf.push(Particle {
            p: [x, y, 0.0],
            v: [vx, vy, 0.0],
            m: 1e-14,
            r: 1e-7,
        });
    }
    particle_buf
}

pub fn two_bodies_soa() -> ParticleSystem {
  let mut system = ParticleSystem::new(2);
  system.positions[0] = [0.0, 0.0, 0.0];
  system.velocities[0] = [0.0, 0.0, 0.0];
  system.radii[0] = 1.0;
  system.masses[0] = 1.0;
  
  system.positions[1] = [1.0, 0.0, 0.0];
  system.velocities[1] = [0.0, 1.0, 0.0];
  system.radii[1] = 1e-4;
  system.masses[1] = 1e-20;
  
  system
}

pub fn circular_orbits_soa(n: usize) -> ParticleSystem {
  let mut system = ParticleSystem::new(n + 1);
  
  // Central body
  system.positions[0] = [0.0, 0.0, 0.0];
  system.velocities[0] = [0.0, 0.0, 0.0];
  system.radii[0] = 0.00465047;
  system.masses[0] = 1.0;

  for i in 0..n {
    let d = 0.1 + ((i as f64) * 5.0 / (n as f64));
    let v = f64::sqrt(1.0 / d);
    let theta = fastrand::f64() * 6.28;
    let x = d * f64::cos(theta);
    let y = d * f64::sin(theta);
    let vx = -v * f64::sin(theta);
    let vy = v * f64::cos(theta);
    
    system.positions[i + 1] = [x, y, 0.0];
    system.velocities[i + 1] = [vx, vy, 0.0];
    system.radii[i + 1] = 1e-7;
    system.masses[i + 1] = 1e-14;
  }
  
  system
}

pub fn distance_sqr(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
  let dx = x1[0] - x2[0];
  let dy = x1[1] - x2[1];
  let dz = x1[2] - x2[2];
  dx*dx + dy*dy + dz*dz
}

pub fn distance(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
  f64::sqrt(distance_sqr(x1, x2))
}

pub fn calc_pp_accel(pi: &Particle, pj: &Particle) -> [f64; 3] {
  let dx = pi.p[0] - pj.p[0];
  let dy = pi.p[1] - pj.p[1];
  let dz = pi.p[2] - pj.p[2];
  let dp2 = dx * dx + dy * dy + dz * dz;
  let dist = f64::sqrt(dp2);
  let magi = -pj.m / (dist*dist*dist);
  [magi*dx, magi*dy, magi*dz]
}

pub fn calc_cm_accel(pi: &Particle, m: f64, cm: [f64; 3]) -> [f64; 3] {
  let dx = pi.p[0] - cm[0];
  let dy = pi.p[1] - cm[1];
  let dz = pi.p[2] - cm[2];
  let dp2 = dx * dx + dy * dy + dz * dz;
  let dist = f64::sqrt(dp2);
  let magi = -m / (dist*dp2);
  [dx * magi, dy * magi, dz * magi]
}

pub fn calc_pp_accel_soa(i: usize, j: usize, system: &ParticleSystem) -> [f64; 3] {
  let p_i = &system.positions[i];
  let p_j = &system.positions[j];
  
  let dx = p_i[0] - p_j[0];
  let dy = p_i[1] - p_j[1];
  let dz = p_i[2] - p_j[2];
  let dp2 = dx * dx + dy * dy + dz * dz;
  let dist = f64::sqrt(dp2);
  let magi = -system.masses[j] / (dist*dist*dist);
  
  [magi*dx, magi*dy, magi*dz]
}

pub fn calc_cm_accel_soa(i: usize, m: f64, cm: [f64; 3], system: &ParticleSystem) -> [f64; 3] {
  let p_i = &system.positions[i];
  
  let dx = p_i[0] - cm[0];
  let dy = p_i[1] - cm[1];
  let dz = p_i[2] - cm[2];
  let dp2 = dx * dx + dy * dy + dz * dz;
  let dist = f64::sqrt(dp2);
  let magi = -m / (dist*dp2);
  
  [dx * magi, dy * magi, dz * magi]
}

pub fn calc_kinetic_energy(particles: &[Particle]) -> f64 {
    particles.iter().fold(0.0, |ke, p| {
        let v2 = p.v[0] * p.v[0] + p.v[1] * p.v[1] + p.v[2] * p.v[2];
        ke + 0.5 * p.m * v2
    })
}

pub fn calc_potential_energy(particles: &[Particle]) -> f64 {
    let mut pe = 0.0;
    for i in 0..particles.len() {
        for j in (i + 1)..particles.len() {
            let dist = distance(&particles[i].p, &particles[j].p);
            pe -= particles[i].m * particles[j].m / dist;
        }
    }
    pe
}

pub fn calc_total_energy(particles: &[Particle]) -> f64 {
    calc_kinetic_energy(particles) + calc_potential_energy(particles)
}

pub fn calc_kinetic_energy_soa(system: &ParticleSystem) -> f64 {
    let mut ke = 0.0;
    for i in 0..system.count {
        let v = &system.velocities[i];
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        ke += 0.5 * system.masses[i] * v2;
    }
    ke
}

pub fn calc_potential_energy_soa(system: &ParticleSystem) -> f64 {
    let mut pe = 0.0;
    for i in 0..system.count {
        for j in (i + 1)..system.count {
            let dist = distance(&system.positions[i], &system.positions[j]);
            pe -= system.masses[i] * system.masses[j] / dist;
        }
    }
    pe
}

pub fn calc_total_energy_soa(system: &ParticleSystem) -> f64 {
    calc_kinetic_energy_soa(system) + calc_potential_energy_soa(system)
}

// Add SIMD implementations when the feature is available
#[cfg(feature = "portable_simd")]
use std::simd::*;

#[cfg(feature = "portable_simd")]
pub fn distance_simd(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
    let p1 = f64x4::from_array([x1[0], x1[1], x1[2], 0.0]);
    let p2 = f64x4::from_array([x2[0], x2[1], x2[2], 0.0]);
    let diff = p1 - p2;
    let squared = diff * diff;
    f64::sqrt(squared[0] + squared[1] + squared[2])
}

#[cfg(feature = "portable_simd")]
pub fn calc_pp_accel_simd(p_idx: usize, j_idx: usize, system: &ParticleSystem) -> [f64; 3] {
    let p_i = &system.positions[p_idx];
    let p_j = &system.positions[j_idx];
    
    // Load positions into SIMD registers
    let pos_i = f64x4::from_array([p_i[0], p_i[1], p_i[2], 0.0]);
    let pos_j = f64x4::from_array([p_j[0], p_j[1], p_j[2], 0.0]);
    
    // Calculate displacement vector
    let disp = pos_i - pos_j;
    
    // Calculate distance squared
    let disp_squared = disp * disp;
    let dist_squared = disp_squared[0] + disp_squared[1] + disp_squared[2];
    let dist = f64::sqrt(dist_squared);
    
    // Calculate magnitude
    let dist_cubed = dist * dist_squared;
    let magi = -system.masses[j_idx] / dist_cubed;
    
    // Apply magnitude to displacement vector
    let acc = disp * f64x4::splat(magi);
    
    [acc[0], acc[1], acc[2]]
}

// Memory-aligned version for better SIMD performance
#[repr(align(32))]
pub struct AlignedParticleSystem {
    pub positions: Vec<[f64; 4]>,  // Use 4-element arrays for better SIMD alignment
    pub velocities: Vec<[f64; 4]>, // Use 4-element arrays for better SIMD alignment
    pub radii: Vec<f64>,
    pub masses: Vec<f64>,
    pub count: usize,
}

impl AlignedParticleSystem {
    pub fn new(count: usize) -> Self {
        AlignedParticleSystem {
            positions: vec![[0.0, 0.0, 0.0, 0.0]; count],
            velocities: vec![[0.0, 0.0, 0.0, 0.0]; count],
            radii: vec![0.0; count],
            masses: vec![0.0; count],
            count,
        }
    }
    
    pub fn from_system(system: &ParticleSystem) -> Self {
        let mut aligned = AlignedParticleSystem::new(system.count);
        
        for i in 0..system.count {
            aligned.positions[i] = [
                system.positions[i][0],
                system.positions[i][1],
                system.positions[i][2],
                0.0
            ];
            
            aligned.velocities[i] = [
                system.velocities[i][0],
                system.velocities[i][1],
                system.velocities[i][2],
                0.0
            ];
            
            aligned.radii[i] = system.radii[i];
            aligned.masses[i] = system.masses[i];
        }
        
        aligned
    }
}