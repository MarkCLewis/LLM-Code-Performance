package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

// --- Constants ---
const (
	G                 = 1.0   // Gravitational constant
	dt                = 0.001 // Time step
	epsilon           = 0.01  // Softening factor to avoid division by zero/large forces
	numSteps          = 100   // Number of simulation steps
	numOrbitingBodies = 10000 // Number of orbiting bodies (WARNING: Very large for sequential!)
	// numOrbitingBodies = 1000 // More reasonable number for testing
)

// --- Data Structures ---

// Vector3D represents a point or vector in 3D space
type Vector3D struct {
	X, Y, Z float64
}

// Body represents a celestial body with mass, position, and velocity
type Body struct {
	Mass     float64
	Position Vector3D
	Velocity Vector3D
}

// --- Vector Operations ---

// Add returns the sum of two vectors
func (v Vector3D) Add(other Vector3D) Vector3D {
	return Vector3D{v.X + other.X, v.Y + other.Y, v.Z + other.Z}
}

// Subtract returns the difference between two vectors (v - other)
func (v Vector3D) Subtract(other Vector3D) Vector3D {
	return Vector3D{v.X - other.X, v.Y - other.Y, v.Z - other.Z}
}

// Scale returns the vector scaled by a scalar factor
func (v Vector3D) Scale(scalar float64) Vector3D {
	return Vector3D{v.X * scalar, v.Y * scalar, v.Z * scalar}
}

// MagnitudeSquared returns the squared magnitude (length) of the vector
func (v Vector3D) MagnitudeSquared() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// Magnitude returns the magnitude (length) of the vector
func (v Vector3D) Magnitude() float64 {
	return math.Sqrt(v.MagnitudeSquared())
}

// Normalize returns a unit vector in the same direction
func (v Vector3D) Normalize() Vector3D {
	mag := v.Magnitude()
	if mag == 0 {
		return Vector3D{0, 0, 0} // Avoid division by zero
	}
	return v.Scale(1.0 / mag)
}

// --- Simulation Logic ---

// calculateTotalEnergy calculates the total kinetic and potential energy of the system
// func calculateTotalEnergy(bodies []Body) float64 {
// 	kineticEnergy := 0.0
// 	potentialEnergy := 0.0

// 	// Kinetic Energy: Sum(0.5 * m * v^2)
// 	for i := range bodies {
// 		kineticEnergy += 0.5 * bodies[i].Mass * bodies[i].Velocity.MagnitudeSquared()
// 	}

// 	// Potential Energy: Sum(-G * m_i * m_j / |r_i - r_j|) for i < j
// 	for i := 0; i < len(bodies); i++ {
// 		for j := i + 1; j < len(bodies); j++ {
// 			deltaPos := bodies[i].Position.Subtract(bodies[j].Position)
// 			distance := deltaPos.Magnitude()
// 			if distance > 0 { // Avoid self-interaction potential if distance is zero
// 				potentialEnergy -= G * bodies[i].Mass * bodies[j].Mass / distance
// 			}
// 		}
// 	}

// 	return kineticEnergy + potentialEnergy
// }

// simulateStep performs one step of the N-body simulation using the kick-step method
func simulateStep(bodies []Body, timeStep, gravitationalConstant, softeningFactor float64) {
	n := len(bodies)
	accelerations := make([]Vector3D, n) // Store accelerations calculated at the start of the step

	// 1. Calculate forces and accelerations for all bodies based on current positions
	// Iterate over unique pairs (i, j) where i < j
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			deltaPos := bodies[i].Position.Subtract(bodies[j].Position)
			distSq := deltaPos.MagnitudeSquared()
			invDistCubed := 1.0 / math.Pow(distSq+softeningFactor*softeningFactor, 1.5) // Include softening

			forceDirection := deltaPos // Direction from j to i
			forceMagnitude := gravitationalConstant * bodies[i].Mass * bodies[j].Mass * invDistCubed
			forceVector := forceDirection.Scale(forceMagnitude)

			// Newton's Third Law: F_ij = -F_ji
			// Accumulate acceleration (a = F/m)
			accelerations[i] = accelerations[i].Subtract(forceVector.Scale(1.0 / bodies[i].Mass)) // Force ON i BY j is in -deltaPos direction
			accelerations[j] = accelerations[j].Add(forceVector.Scale(1.0 / bodies[j].Mass))      // Force ON j BY i is in +deltaPos direction
		}
	}

	// 2. Kick: Update velocities using calculated accelerations
	for i := range bodies {
		bodies[i].Velocity = bodies[i].Velocity.Add(accelerations[i].Scale(timeStep))
	}

	// 3. Step: Update positions using the *new* velocities
	for i := range bodies {
		bodies[i].Position = bodies[i].Position.Add(bodies[i].Velocity.Scale(timeStep))
	}
}

// --- Initialization ---

// initializeCircularOrbitSystem creates a central body and N orbiting bodies
// Orbiting bodies are placed randomly between minRadius and maxRadius in 3D space
// with initial velocities for circular orbits around the central body.
func initializeCircularOrbitSystem(numOrbiting int, centralMass, orbitingMass, minRadius, maxRadius float64) []Body {
	bodies := make([]Body, numOrbiting+1)
	rand.Seed(time.Now().UnixNano()) // Seed random number generator

	// Central Body (at origin, initially stationary)
	bodies[0] = Body{
		Mass:     centralMass,
		Position: Vector3D{0, 0, 0},
		Velocity: Vector3D{0, 0, 0},
	}

	// Orbiting Bodies
	for i := 1; i <= numOrbiting; i++ {
		// Random radius
		radius := minRadius + rand.Float64()*(maxRadius-minRadius)

		// Random direction for position (uniform on sphere surface)
		theta := 2 * math.Pi * rand.Float64()  // Azimuthal angle
		phi := math.Acos(1 - 2*rand.Float64()) // Polar angle (using correct sampling)
		posX := radius * math.Sin(phi) * math.Cos(theta)
		posY := radius * math.Sin(phi) * math.Sin(theta)
		posZ := radius * math.Cos(phi)
		positionVec := Vector3D{posX, posY, posZ}

		// Calculate circular orbit speed: v = sqrt(G * M_central / r)
		orbitSpeed := math.Sqrt(G * centralMass / radius)

		// Find a velocity vector perpendicular to the position vector.
		// Create a random vector not parallel to position, take cross product.
		// Exception: If position is along Z-axis, use (0,1,0) as the helper.
		helperVec := Vector3D{rand.Float64()*2 - 1, rand.Float64()*2 - 1, rand.Float64()*2 - 1}
		// Ensure helper is not parallel to position (unlikely but possible)
		if math.Abs(positionVec.Normalize().X-helperVec.Normalize().X) < 1e-6 &&
			math.Abs(positionVec.Normalize().Y-helperVec.Normalize().Y) < 1e-6 &&
			math.Abs(positionVec.Normalize().Z-helperVec.Normalize().Z) < 1e-6 {
			helperVec = Vector3D{0, 1, 0} // Use a default if parallel
		}

		// Calculate cross product: pos x helper gives a vector perpendicular to both
		velDirection := Vector3D{
			positionVec.Y*helperVec.Z - positionVec.Z*helperVec.Y,
			positionVec.Z*helperVec.X - positionVec.X*helperVec.Z,
			positionVec.X*helperVec.Y - positionVec.Y*helperVec.X,
		}.Normalize() // Normalize to get unit direction vector

		velocityVec := velDirection.Scale(orbitSpeed)

		bodies[i] = Body{
			Mass:     orbitingMass,
			Position: positionVec,
			Velocity: velocityVec,
		}
	}

	return bodies
}

// --- Main Function ---

func main() {
	fmt.Printf("Starting N-body simulation...\n")
	fmt.Printf("Number of bodies: %d (+1 central body)\n", numOrbitingBodies)
	fmt.Printf("Number of steps: %d\n", numSteps)
	fmt.Printf("Time step (dt): %.4f\n", dt)
	fmt.Printf("WARNING: N=%d is very large for a sequential O(N^2) simulation. This may take a very long time.\n", numOrbitingBodies+1)

	// --- Initialization ---
	startTime := time.Now()
	centralMass := 10000.0 // Mass of the central body
	orbitingMass := 0.01   // Mass of the smaller orbiting bodies
	minRadius := 5.0
	maxRadius := 50.0

	bodies := initializeCircularOrbitSystem(numOrbitingBodies, centralMass, orbitingMass, minRadius, maxRadius)
	initializationTime := time.Since(startTime)
	fmt.Printf("System initialized in %v\n", initializationTime)

	// --- Initial Energy Calculation ---
	// initialEnergy := calculateTotalEnergy(bodies)
	// fmt.Printf("Initial Total Energy: %.5e\n", initialEnergy)

	// --- Simulation Loop ---
	fmt.Println("Running simulation...")
	simulationStartTime := time.Now()

	for step := 0; step < numSteps; step++ {
		simulateStep(bodies, dt, G, epsilon)
		if (step+1)%10 == 0 { // Print progress occasionally
			fmt.Printf("Step %d/%d completed\n", step+1, numSteps)
		}
	}

	simulationTime := time.Since(simulationStartTime)
	fmt.Printf("Simulation completed in %v\n", simulationTime)

	// --- Final Energy Calculation ---
	// finalEnergy := calculateTotalEnergy(bodies)
	// fmt.Printf("Final Total Energy:   %.5e\n", finalEnergy)

	// // --- Energy Verification ---
	// energyChange := finalEnergy - initialEnergy
	// relativeEnergyChange := math.Abs(energyChange / initialEnergy)
	// fmt.Printf("Absolute Energy Change: %.5e\n", energyChange)
	// fmt.Printf("Relative Energy Change: %.5e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange*100.0)
	fmt.Printf("Body[0] %e %e %e", bodies[0].Position.X, bodies[0].Position.Y, bodies[0].Position.Z)

	totalTime := time.Since(startTime)
	fmt.Printf("Total execution time: %v\n", totalTime)

	// Optional: Print final state of a few bodies
	// fmt.Println("Final state of first few bodies:")
	// for i := 0; i < min(5, len(bodies)); i++ {
	// 	fmt.Printf(" Body %d: Pos=(%.2f, %.2f, %.2f), Vel=(%.2f, %.2f, %.2f)\n",
	// 		i, bodies[i].Position.X, bodies[i].Position.Y, bodies[i].Position.Z,
	// 		bodies[i].Velocity.X, bodies[i].Velocity.Y, bodies[i].Velocity.Z)
	// }
}

// Helper function min
// func min(a, b int) int {
// 	if a < b {
// 		return a
// 	}
// 	return b
// }
