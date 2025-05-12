package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

// Constants
const (
	G              = 6.67430e-11 // Gravitational constant (m^3 kg^-1 s^-2)
	dt             = 86400.0     // Time step in seconds (1 day)
	softening      = 1.0e7       // Softening parameter to avoid singularities
	centralMass    = 1.989e30    // Mass of central body (kg) - roughly solar mass
	planetMass     = 1.0e24      // Mass of smaller bodies (kg)
	numSteps       = 100         // Number of simulation steps
	numSmallBodies = 10000       // Number of small bodies
)

// Vector3 represents a 3D vector
type Vector3 struct {
	X, Y, Z float64
}

// Body represents a celestial body
type Body struct {
	Position Vector3
	Velocity Vector3
	Mass     float64
}

// Add returns the sum of two vectors
func (v Vector3) Add(other Vector3) Vector3 {
	return Vector3{
		X: v.X + other.X,
		Y: v.Y + other.Y,
		Z: v.Z + other.Z,
	}
}

// Sub returns the difference of two vectors
func (v Vector3) Sub(other Vector3) Vector3 {
	return Vector3{
		X: v.X - other.X,
		Y: v.Y - other.Y,
		Z: v.Z - other.Z,
	}
}

// Mul returns the product of a vector and a scalar
func (v Vector3) Mul(scalar float64) Vector3 {
	return Vector3{
		X: v.X * scalar,
		Y: v.Y * scalar,
		Z: v.Z * scalar,
	}
}

// SqrMagnitude returns the squared magnitude of a vector
func (v Vector3) SqrMagnitude() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// Magnitude returns the magnitude of a vector
func (v Vector3) Magnitude() float64 {
	return math.Sqrt(v.SqrMagnitude())
}

// Normalized returns a normalized version of the vector
func (v Vector3) Normalized() Vector3 {
	mag := v.Magnitude()
	if mag == 0 {
		return Vector3{}
	}
	return Vector3{
		X: v.X / mag,
		Y: v.Y / mag,
		Z: v.Z / mag,
	}
}

// CalculateGravitationalForce calculates the gravitational force between two bodies
func CalculateGravitationalForce(a, b Body) Vector3 {
	// Vector from a to b
	r := b.Position.Sub(a.Position)

	// Squared distance with softening parameter
	distSqr := r.SqrMagnitude() + softening*softening

	// Force magnitude: F = G*m1*m2/r^2
	forceMag := G * a.Mass * b.Mass / distSqr

	// Force direction (normalized)
	direction := r.Normalized()

	// Return the force vector
	return direction.Mul(forceMag)
}

// CalculateAcceleration calculates the acceleration of a body due to gravitational forces
func CalculateAcceleration(bodies []Body, index int) Vector3 {
	body := bodies[index]
	totalForce := Vector3{}

	// Sum forces from all other bodies
	for i, otherBody := range bodies {
		if i != index {
			force := CalculateGravitationalForce(body, otherBody)
			totalForce = totalForce.Add(force)
		}
	}

	// a = F/m
	return Vector3{
		X: totalForce.X / body.Mass,
		Y: totalForce.Y / body.Mass,
		Z: totalForce.Z / body.Mass,
	}
}

// KickStep performs a single kick-step (leapfrog) integration
func KickStep(bodies []Body) {
	// Temporary array to store accelerations
	accelerations := make([]Vector3, len(bodies))

	// Calculate accelerations for all bodies
	for i := range bodies {
		accelerations[i] = CalculateAcceleration(bodies, i)
	}

	// Update positions and velocities
	for i := range bodies {
		// Update velocity (kick)
		bodies[i].Velocity = bodies[i].Velocity.Add(accelerations[i].Mul(dt))

		// Update position (drift)
		bodies[i].Position = bodies[i].Position.Add(bodies[i].Velocity.Mul(dt))
	}
}

// CalculateTotalEnergy calculates the total energy of the system (kinetic + potential)
// func CalculateTotalEnergy(bodies []Body) float64 {
// 	kineticEnergy := 0.0
// 	potentialEnergy := 0.0

// 	// Calculate kinetic energy: KE = 0.5 * m * v^2
// 	for _, body := range bodies {
// 		kineticEnergy += 0.5 * body.Mass * body.Velocity.SqrMagnitude()
// 	}

// 	// Calculate potential energy: PE = -G * m1 * m2 / r
// 	for i := 0; i < len(bodies); i++ {
// 		for j := i + 1; j < len(bodies); j++ {
// 			r := bodies[i].Position.Sub(bodies[j].Position)
// 			distance := r.Magnitude()

// 			// Avoid division by zero
// 			if distance > 0 {
// 				potentialEnergy -= G * bodies[i].Mass * bodies[j].Mass / distance
// 			}
// 		}
// 	}

// 	return kineticEnergy + potentialEnergy
// }

// InitializeSystem creates a system with a central body and smaller bodies on circular orbits
func InitializeSystem(numSmallBodies int) []Body {
	bodies := make([]Body, numSmallBodies+1)

	// Set up central body at the origin
	bodies[0] = Body{
		Position: Vector3{},
		Velocity: Vector3{},
		Mass:     centralMass,
	}

	rand.Seed(time.Now().UnixNano())

	// Set up smaller bodies in circular orbits around the central body
	for i := 1; i <= numSmallBodies; i++ {
		// Random orbital radius between 1 AU and 5 AU (in meters)
		// 1 AU ≈ 1.496e11 meters
		radius := (1.0 + 4.0*rand.Float64()) * 1.496e11

		// Random angle for position in the orbit
		theta := rand.Float64() * 2 * math.Pi
		phi := (rand.Float64() - 0.5) * 0.1 * math.Pi // Near-planar with some variation

		// Position in spherical coordinates converted to Cartesian
		x := radius * math.Cos(theta) * math.Cos(phi)
		y := radius * math.Sin(theta) * math.Cos(phi)
		z := radius * math.Sin(phi)

		position := Vector3{X: x, Y: y, Z: z}

		// Calculate orbital velocity for a circular orbit
		// v = sqrt(G*M/r)
		orbitalSpeed := math.Sqrt(G * centralMass / radius)

		// Velocity vector is perpendicular to the position vector
		// For a circular orbit in the xy-plane, it's (-y, x, 0) normalized * speed
		velocity := Vector3{
			X: -position.Y,
			Y: position.X,
			Z: 0,
		}.Normalized().Mul(orbitalSpeed)

		bodies[i] = Body{
			Position: position,
			Velocity: velocity,
			Mass:     planetMass,
		}
	}

	return bodies
}

func main() {
	fmt.Println("Initializing N-body simulation...")
	startTime := time.Now()

	// Initialize system with one central body and numSmallBodies smaller bodies
	bodies := InitializeSystem(numSmallBodies)

	fmt.Printf("System initialized with %d bodies\n", len(bodies))
	fmt.Printf("Central body mass: %.3e kg\n", bodies[0].Mass)
	fmt.Printf("Small body mass: %.3e kg\n", bodies[1].Mass)

	// Calculate initial energy
	// initialEnergy := CalculateTotalEnergy(bodies)
	// fmt.Printf("Initial total energy: %.6e J\n", initialEnergy)

	// Run simulation for numSteps steps
	fmt.Printf("Running simulation for %d steps...\n", numSteps)
	for step := 1; step <= numSteps; step++ {
		if step%(numSteps/10) == 0 {
			fmt.Printf("Progress: %d%%\n", step*100/numSteps)
		}
		KickStep(bodies)
	}

	// Calculate final energy
	// finalEnergy := CalculateTotalEnergy(bodies)
	// fmt.Printf("Final total energy: %.6e J\n", finalEnergy)

	// // Calculate energy conservation
	// energyDiff := (finalEnergy - initialEnergy) / initialEnergy
	// fmt.Printf("Energy difference: %.6e (%.6f%%)\n", finalEnergy-initialEnergy, energyDiff*100)
	fmt.Printf("Body[0] %e %e %e", bodies[0].Position.X, bodies[0].Position.Y, bodies[0].Position.Z)

	// Print execution time
	elapsedTime := time.Since(startTime)
	fmt.Printf("Simulation completed in %v\n", elapsedTime)
}
