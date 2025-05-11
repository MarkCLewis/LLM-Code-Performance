package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime" // Added for NumCPU
	"sync"    // Added for WaitGroup
	"time"
)

// --- Constants ---
const (
	G                 = 1.0   // Gravitational constant
	dt                = 0.001 // Time step
	epsilon           = 0.01  // Softening factor to avoid division by zero/large forces
	numSteps          = 100   // Number of simulation steps
	numOrbitingBodies = 10000 // Number of orbiting bodies (WARNING: Very large even for parallel!)
	// numOrbitingBodies = 1000 // More reasonable number for testing
)

// --- Data Structures ---

// Vector3D represents a point or vector in 3D space
type Vector3D struct {
	X, Y, Z float64
}

// Body represents a celestial body with mass, position, and velocity
// Using pointers for Position and Velocity might be slightly more idiomatic
// when modifying them in place frequently, but structs work fine here too.
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
// (Remains sequential as it's mainly for verification before/after)
func calculateTotalEnergy(bodies []Body) float64 {
	kineticEnergy := 0.0
	potentialEnergy := 0.0

	// Kinetic Energy: Sum(0.5 * m * v^2)
	for i := range bodies {
		kineticEnergy += 0.5 * bodies[i].Mass * bodies[i].Velocity.MagnitudeSquared()
	}

	// Potential Energy: Sum(-G * m_i * m_j / |r_i - r_j|) for i < j
	// This O(N^2) part is only run twice, so parallelizing it isn't critical path
	for i := 0; i < len(bodies); i++ {
		for j := i + 1; j < len(bodies); j++ {
			deltaPos := bodies[i].Position.Subtract(bodies[j].Position)
			distance := deltaPos.Magnitude()
			if distance > 0 { // Avoid self-interaction potential if distance is zero
				potentialEnergy -= G * bodies[i].Mass * bodies[j].Mass / distance
			}
		}
	}

	return kineticEnergy + potentialEnergy
}

// --- Parallel Simulation Step ---

// calculateAccelerationChunk calculates acceleration for a subset of bodies
func calculateAccelerationChunk(bodies []Body, accelerations []Vector3D, start, end int, wg *sync.WaitGroup, gravitationalConstant, softeningFactor float64) {
	defer wg.Done() // Signal completion when function exits
	n := len(bodies)
	softeningSq := softeningFactor * softeningFactor

	for i := start; i < end; i++ {
		acc := Vector3D{0, 0, 0} // Initialize acceleration for body i
		for j := 0; j < n; j++ {
			if i == j {
				continue // Don't interact with self
			}
			deltaPos := bodies[j].Position.Subtract(bodies[i].Position) // Vector from i to j
			distSq := deltaPos.MagnitudeSquared()
			invDistCubed := 1.0 / math.Pow(distSq+softeningSq, 1.5) // Include softening

			// Force magnitude G * m_i * m_j / (r^2 + eps^2)^(3/2) * deltaPos
			// Acc = F/m_i = G * m_j / (r^2 + eps^2)^(3/2) * deltaPos
			accScale := gravitationalConstant * bodies[j].Mass * invDistCubed
			acc = acc.Add(deltaPos.Scale(accScale))
		}
		accelerations[i] = acc // Store calculated acceleration
	}
}

// updateVelocityChunk updates velocity for a subset of bodies
func updateVelocityChunk(bodies []Body, accelerations []Vector3D, start, end int, wg *sync.WaitGroup, timeStep float64) {
	defer wg.Done()
	for i := start; i < end; i++ {
		// Use pointer to modify the slice element directly if using []*Body
		// bodies[i].Velocity = bodies[i].Velocity.Add(accelerations[i].Scale(timeStep))
		// If using []Body, we need to assign back to the slice index
		bodies[i].Velocity = bodies[i].Velocity.Add(accelerations[i].Scale(timeStep))
	}
}

// updatePositionChunk updates position for a subset of bodies
func updatePositionChunk(bodies []Body, start, end int, wg *sync.WaitGroup, timeStep float64) {
	defer wg.Done()
	for i := start; i < end; i++ {
		// Use pointer to modify the slice element directly if using []*Body
		// bodies[i].Position = bodies[i].Position.Add(bodies[i].Velocity.Scale(timeStep))
		// If using []Body, we need to assign back to the slice index
		bodies[i].Position = bodies[i].Position.Add(bodies[i].Velocity.Scale(timeStep))
	}
}

// simulateStepParallel performs one step using multiple goroutines
func simulateStepParallel(bodies []Body, accelerations []Vector3D, numWorkers int, timeStep, gravitationalConstant, softeningFactor float64) {
	n := len(bodies)
	var wg sync.WaitGroup

	// --- 1. Parallel Acceleration Calculation ---
	wg.Add(numWorkers)
	chunkSize := (n + numWorkers - 1) / numWorkers // Ceiling division
	for i := 0; i < numWorkers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if end > n {
			end = n
		}
		if start >= end { // Handle case where numWorkers > n
			wg.Done() // Decrement counter if worker has no work
			continue
		}
		// Pass necessary variables to the goroutine
		go calculateAccelerationChunk(bodies, accelerations, start, end, &wg, gravitationalConstant, softeningFactor)
	}
	wg.Wait() // Wait for all acceleration calculations to finish

	// --- 2. Parallel Velocity Update (Kick) ---
	wg.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if end > n {
			end = n
		}
		if start >= end {
			wg.Done()
			continue
		}
		go updateVelocityChunk(bodies, accelerations, start, end, &wg, timeStep)
	}
	wg.Wait() // Wait for all velocity updates to finish

	// --- 3. Parallel Position Update (Step) ---
	wg.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if end > n {
			end = n
		}
		if start >= end {
			wg.Done()
			continue
		}
		go updatePositionChunk(bodies, start, end, &wg, timeStep)
	}
	wg.Wait() // Wait for all position updates to finish
}

// --- Initialization --- (Same as before)

// initializeCircularOrbitSystem creates a central body and N orbiting bodies
func initializeCircularOrbitSystem(numOrbiting int, centralMass, orbitingMass, minRadius, maxRadius float64) []Body {
	bodies := make([]Body, numOrbiting+1)
	// Use a fixed seed for reproducibility during testing if needed
	// rand.Seed(12345)
	rand.Seed(time.Now().UnixNano()) // Seed random number generator

	// Central Body (at origin, initially stationary)
	bodies[0] = Body{
		Mass:     centralMass,
		Position: Vector3D{0, 0, 0},
		Velocity: Vector3D{0, 0, 0},
	}

	// Orbiting Bodies
	for i := 1; i <= numOrbiting; i++ {
		radius := minRadius + rand.Float64()*(maxRadius-minRadius)
		theta := 2 * math.Pi * rand.Float64()
		phi := math.Acos(1 - 2*rand.Float64())
		posX := radius * math.Sin(phi) * math.Cos(theta)
		posY := radius * math.Sin(phi) * math.Sin(theta)
		posZ := radius * math.Cos(phi)
		positionVec := Vector3D{posX, posY, posZ}

		orbitSpeed := math.Sqrt(G * centralMass / radius)

		// Simplified perpendicular velocity calculation (assuming Z is up)
		// For a position vector (x, y, z), a perpendicular vector in the XY plane projection's
		// perpendicular direction is (-y, x, 0). Normalize and scale.
		// This prefers orbits somewhat aligned with the XY plane initially.
		// A more robust method (like the cross-product one) ensures better 3D distribution.
		// Let's stick to the cross-product method for better generality:
		helperVec := Vector3D{rand.Float64()*2 - 1, rand.Float64()*2 - 1, rand.Float64()*2 - 1}.Normalize()
		// Avoid parallel vectors (cross product would be zero)
		dotProd := positionVec.X*helperVec.X + positionVec.Y*helperVec.Y + positionVec.Z*helperVec.Z
		if math.Abs(dotProd) > 0.99 { // If vectors are nearly parallel/anti-parallel
			helperVec = Vector3D{0, 1, 0}       // Use a default perpendicular seed vector
			if math.Abs(positionVec.Y) > 0.99 { // If position is along Y axis, use X axis
				helperVec = Vector3D{1, 0, 0}
			}
		}

		velDirection := Vector3D{
			positionVec.Y*helperVec.Z - positionVec.Z*helperVec.Y,
			positionVec.Z*helperVec.X - positionVec.X*helperVec.Z,
			positionVec.X*helperVec.Y - positionVec.Y*helperVec.X,
		}.Normalize()

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
	numWorkers := runtime.NumCPU()
	runtime.GOMAXPROCS(numWorkers) // Ensure Go uses all available cores

	fmt.Printf("Starting N-body simulation (Parallel)...\n")
	fmt.Printf("Using %d CPU cores/workers.\n", numWorkers)
	fmt.Printf("Number of bodies: %d (+1 central body)\n", numOrbitingBodies)
	fmt.Printf("Number of steps: %d\n", numSteps)
	fmt.Printf("Time step (dt): %.4f\n", dt)
	fmt.Printf("WARNING: N=%d is very large. Even parallelized, this O(N^2) simulation may take an extremely long time.\n", numOrbitingBodies+1)

	// --- Initialization ---
	startTime := time.Now()
	centralMass := 10000.0 // Mass of the central body
	orbitingMass := 0.01   // Mass of the smaller orbiting bodies
	minRadius := 5.0
	maxRadius := 50.0

	bodies := initializeCircularOrbitSystem(numOrbitingBodies, centralMass, orbitingMass, minRadius, maxRadius)
	initializationTime := time.Since(startTime)
	fmt.Printf("System initialized in %v\n", initializationTime)

	// Pre-allocate acceleration slice (avoids allocation in loop)
	accelerations := make([]Vector3D, len(bodies))

	// --- Initial Energy Calculation ---
	initialEnergy := calculateTotalEnergy(bodies)
	fmt.Printf("Initial Total Energy: %.5e\n", initialEnergy)

	// --- Simulation Loop ---
	fmt.Println("Running simulation...")
	simulationStartTime := time.Now()

	for step := 0; step < numSteps; step++ {
		// Pass the pre-allocated accelerations slice
		simulateStepParallel(bodies, accelerations, numWorkers, dt, G, epsilon)
		if (step+1)%10 == 0 || step == 0 { // Print progress more often initially/periodically
			fmt.Printf("Step %d/%d completed (Elapsed: %v)\n", step+1, numSteps, time.Since(simulationStartTime))
		}
	}

	simulationTime := time.Since(simulationStartTime)
	fmt.Printf("Simulation completed in %v\n", simulationTime)

	// --- Final Energy Calculation ---
	finalEnergy := calculateTotalEnergy(bodies)
	fmt.Printf("Final Total Energy:   %.5e\n", finalEnergy)

	// --- Energy Verification ---
	energyChange := finalEnergy - initialEnergy
	relativeEnergyChange := math.Abs(energyChange / initialEnergy)
	fmt.Printf("Absolute Energy Change: %.5e\n", energyChange)
	fmt.Printf("Relative Energy Change: %.5e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange*100.0)

	totalTime := time.Since(startTime)
	fmt.Printf("Total execution time: %v\n", totalTime)

}
