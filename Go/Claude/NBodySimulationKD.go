package main

// Import sort package
import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"
	"time"
)

// Constants
const (
	G              = 6.67430e-11 // Gravitational constant (m^3 kg^-1 s^-2)
	dt             = 86400.0     // Time step in seconds (1 day)
	softening      = 1.0e7       // Softening parameter to avoid singularities
	centralMass    = 1.989e30    // Mass of central body (kg) - roughly solar mass
	planetMass     = 1.0e24      // Mass of smaller bodies (kg)
	numSteps       = 10          // Number of simulation steps
	numSmallBodies = 100000      // Number of small bodies
	theta          = 0.3         // Barnes-Hut opening angle parameter
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

// KDNode represents a node in the k-d tree
type KDNode struct {
	// Spatial information
	Min, Max  Vector3 // Bounding box
	Center    Vector3 // Center of mass
	Size      float64 // Size of the node (maximum dimension of bounding box)
	TotalMass float64 // Total mass of bodies in this node

	// Tree structure
	Left, Right *KDNode // Child nodes
	Bodies      []int   // Indices of bodies (only for leaf nodes)
	IsLeaf      bool    // Whether this is a leaf node
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

// BuildKDTree constructs a k-d tree from a set of bodies
func BuildKDTree(bodies []Body, indices []int, depth int) *KDNode {
	if len(indices) == 0 {
		return nil
	}

	// Create a new node
	node := &KDNode{
		Min:       Vector3{X: math.MaxFloat64, Y: math.MaxFloat64, Z: math.MaxFloat64},
		Max:       Vector3{X: -math.MaxFloat64, Y: -math.MaxFloat64, Z: -math.MaxFloat64},
		TotalMass: 0,
		Center:    Vector3{},
	}

	// Find the bounding box and calculate center of mass
	weightedCenter := Vector3{}
	for _, idx := range indices {
		body := bodies[idx]

		// Update bounding box
		node.Min.X = math.Min(node.Min.X, body.Position.X)
		node.Min.Y = math.Min(node.Min.Y, body.Position.Y)
		node.Min.Z = math.Min(node.Min.Z, body.Position.Z)

		node.Max.X = math.Max(node.Max.X, body.Position.X)
		node.Max.Y = math.Max(node.Max.Y, body.Position.Y)
		node.Max.Z = math.Max(node.Max.Z, body.Position.Z)

		// Update center of mass
		weightedCenter.X += body.Position.X * body.Mass
		weightedCenter.Y += body.Position.Y * body.Mass
		weightedCenter.Z += body.Position.Z * body.Mass

		node.TotalMass += body.Mass
	}

	// Calculate final center of mass
	if node.TotalMass > 0 {
		node.Center = Vector3{
			X: weightedCenter.X / node.TotalMass,
			Y: weightedCenter.Y / node.TotalMass,
			Z: weightedCenter.Z / node.TotalMass,
		}
	}

	// Calculate node size (largest dimension of bounding box)
	sizeX := math.Abs(node.Max.X - node.Min.X)
	sizeY := math.Abs(node.Max.Y - node.Min.Y)
	sizeZ := math.Abs(node.Max.Z - node.Min.Z)
	node.Size = math.Max(sizeX, math.Max(sizeY, sizeZ))

	// Base case: if there are few enough bodies, make this a leaf node
	if len(indices) <= 16 { // threshold for leaf nodes
		node.Bodies = indices
		node.IsLeaf = true
		return node
	}

	// Determine split axis (cycle through X, Y, Z based on depth)
	axis := depth % 3

	// Sort indices based on position along the split axis
	sortIndices := make([]int, len(indices))
	copy(sortIndices, indices)

	switch axis {
	case 0: // X-axis
		sort.Slice(sortIndices, func(i, j int) bool {
			return bodies[sortIndices[i]].Position.X < bodies[sortIndices[j]].Position.X
		})
	case 1: // Y-axis
		sort.Slice(sortIndices, func(i, j int) bool {
			return bodies[sortIndices[i]].Position.Y < bodies[sortIndices[j]].Position.Y
		})
	case 2: // Z-axis
		sort.Slice(sortIndices, func(i, j int) bool {
			return bodies[sortIndices[i]].Position.Z < bodies[sortIndices[j]].Position.Z
		})
	}

	// Find median index for splitting
	median := len(sortIndices) / 2

	// Recursively build left and right subtrees
	leftIndices := sortIndices[:median]
	rightIndices := sortIndices[median:]

	// Create child nodes
	node.Left = BuildKDTree(bodies, leftIndices, depth+1)
	node.Right = BuildKDTree(bodies, rightIndices, depth+1)
	node.IsLeaf = false

	return node
}

// CalculateForceFromNode computes the gravitational force on a body from a k-d tree node
func CalculateForceFromNode(body Body, node *KDNode, bodies []Body) Vector3 {
	if node == nil {
		return Vector3{}
	}

	// Calculate distance from body to node's center of mass
	r := node.Center.Sub(body.Position)
	distSqr := r.SqrMagnitude()

	// If this is a leaf node or the node is sufficiently far away (Barnes-Hut approximation)
	// s/d < θ, where s is the node size, d is the distance, and θ is the opening angle parameter
	if node.IsLeaf || (node.Size*node.Size < theta*theta*distSqr) {
		if node.IsLeaf {
			// For leaf nodes, calculate force directly from each body
			force := Vector3{}
			for _, idx := range node.Bodies {
				if &bodies[idx] != &body { // Skip self-interaction
					// Vector from body to other body
					r := bodies[idx].Position.Sub(body.Position)
					distSqr := r.SqrMagnitude() + softening*softening

					// Force magnitude: F = G*m1*m2/r^2
					forceMag := G * body.Mass * bodies[idx].Mass / distSqr

					// Force direction (normalized)
					direction := r.Normalized()

					// Add to total force
					force = force.Add(direction.Mul(forceMag))
				}
			}
			return force
		} else {
			// For distant nodes, use center of mass approximation
			// Avoid division by zero with softening
			distSqr += softening * softening

			// Force magnitude: F = G*m1*m2/r^2
			forceMag := G * body.Mass * node.TotalMass / distSqr

			// Force direction (normalized)
			direction := r.Normalized()

			return direction.Mul(forceMag)
		}
	}

	// Node is too close for approximation, traverse children
	leftForce := CalculateForceFromNode(body, node.Left, bodies)
	rightForce := CalculateForceFromNode(body, node.Right, bodies)

	// Combine forces from both children
	return leftForce.Add(rightForce)
}

// CalculateAcceleration calculates the acceleration of a body due to gravitational forces using k-d tree
func CalculateAcceleration(body Body, root *KDNode, bodies []Body) Vector3 {
	// Calculate the total gravitational force on the body from the k-d tree
	totalForce := CalculateForceFromNode(body, root, bodies)

	// a = F/m
	return Vector3{
		X: totalForce.X / body.Mass,
		Y: totalForce.Y / body.Mass,
		Z: totalForce.Z / body.Mass,
	}
}

// KickStep performs a single kick-step (leapfrog) integration using goroutines for parallelization
// and a k-d tree for force calculations
func KickStep(bodies []Body, root *KDNode) {
	// Temporary array to store accelerations
	accelerations := make([]Vector3, len(bodies))

	// Determine the number of goroutines based on CPU cores
	numCPU := runtime.NumCPU()
	numGoroutines := numCPU * 2 // Use 2x the number of CPU cores for better utilization

	// Ensure we don't create more goroutines than bodies
	if numGoroutines > len(bodies) {
		numGoroutines = len(bodies)
	}

	// Calculate the chunk size for each goroutine
	chunkSize := len(bodies) / numGoroutines

	// Wait group to synchronize goroutines
	var wg sync.WaitGroup
	wg.Add(numGoroutines)

	// Launch goroutines to calculate accelerations
	for g := 0; g < numGoroutines; g++ {
		// Calculate start and end indices for this goroutine
		start := g * chunkSize
		end := start + chunkSize
		if g == numGoroutines-1 {
			end = len(bodies) // Make sure the last goroutine handles any remainder
		}

		go func(start, end int) {
			defer wg.Done()

			// Calculate accelerations for assigned bodies using k-d tree
			for i := start; i < end; i++ {
				accelerations[i] = CalculateAcceleration(bodies[i], root, bodies)
			}
		}(start, end)
	}

	// Wait for all acceleration calculations to complete
	wg.Wait()

	// Reset wait group for position/velocity updates
	wg.Add(numGoroutines)

	// Launch goroutines to update positions and velocities
	for g := 0; g < numGoroutines; g++ {
		// Calculate start and end indices for this goroutine
		start := g * chunkSize
		end := start + chunkSize
		if g == numGoroutines-1 {
			end = len(bodies) // Make sure the last goroutine handles any remainder
		}

		go func(start, end int) {
			defer wg.Done()

			// Update positions and velocities for assigned bodies
			for i := start; i < end; i++ {
				// Update velocity (kick)
				bodies[i].Velocity = bodies[i].Velocity.Add(accelerations[i].Mul(dt))

				// Update position (drift)
				bodies[i].Position = bodies[i].Position.Add(bodies[i].Velocity.Mul(dt))
			}
		}(start, end)
	}

	// Wait for all position/velocity updates to complete
	wg.Wait()
}

// CalculateTotalEnergy calculates the total energy of the system (kinetic + potential) using parallel computation
// func CalculateTotalEnergy(bodies []Body) float64 {
// 	// Determine the number of goroutines based on CPU cores
// 	numCPU := runtime.NumCPU()
// 	numGoroutines := numCPU * 2 // Use 2x the number of CPU cores for better utilization

// 	// Ensure we don't create more goroutines than bodies
// 	if numGoroutines > len(bodies) {
// 		numGoroutines = len(bodies)
// 	}

// 	// Calculate kinetic energy in parallel
// 	var keMutex sync.Mutex
// 	kineticEnergy := 0.0

// 	// Wait group for synchronization
// 	var wg sync.WaitGroup
// 	wg.Add(numGoroutines)

// 	// Calculate the chunk size for each goroutine
// 	chunkSize := len(bodies) / numGoroutines

// 	// Launch goroutines to calculate kinetic energy
// 	for g := 0; g < numGoroutines; g++ {
// 		// Calculate start and end indices for this goroutine
// 		start := g * chunkSize
// 		end := start + chunkSize
// 		if g == numGoroutines-1 {
// 			end = len(bodies) // Make sure the last goroutine handles any remainder
// 		}

// 		go func(start, end int) {
// 			defer wg.Done()

// 			localKE := 0.0
// 			// Calculate kinetic energy: KE = 0.5 * m * v^2
// 			for i := start; i < end; i++ {
// 				localKE += 0.5 * bodies[i].Mass * bodies[i].Velocity.SqrMagnitude()
// 			}

// 			// Update total kinetic energy with mutex protection
// 			keMutex.Lock()
// 			kineticEnergy += localKE
// 			keMutex.Unlock()
// 		}(start, end)
// 	}

// 	// Wait for kinetic energy calculations to complete
// 	wg.Wait()

// 	// Calculate potential energy in parallel
// 	// For potential energy, we need to handle the double-loop differently
// 	potentialEnergy := 0.0
// 	var peMutex sync.Mutex

// 	// Determine how to divide the work for potential energy calculation
// 	// We'll assign each goroutine a chunk of the first loop
// 	wg.Add(numGoroutines)

// 	for g := 0; g < numGoroutines; g++ {
// 		// Calculate start and end indices for this goroutine
// 		start := g * chunkSize
// 		end := start + chunkSize
// 		if g == numGoroutines-1 {
// 			end = len(bodies) // Make sure the last goroutine handles any remainder
// 		}

// 		go func(start, end int) {
// 			defer wg.Done()

// 			localPE := 0.0
// 			// Calculate potential energy: PE = -G * m1 * m2 / r
// 			for i := start; i < end; i++ {
// 				for j := i + 1; j < len(bodies); j++ {
// 					r := bodies[i].Position.Sub(bodies[j].Position)
// 					distance := r.Magnitude()

// 					// Avoid division by zero
// 					if distance > 0 {
// 						localPE -= G * bodies[i].Mass * bodies[j].Mass / distance
// 					}
// 				}
// 			}

// 			// Update total potential energy with mutex protection
// 			peMutex.Lock()
// 			potentialEnergy += localPE
// 			peMutex.Unlock()
// 		}(start, end)
// 	}

// 	// Wait for potential energy calculations to complete
// 	wg.Wait()

// 	return kineticEnergy + potentialEnergy
// }

// InitializeSystem creates a system with a central body and smaller bodies on circular orbits
// Uses parallel processing for generating the bodies
func InitializeSystem(numSmallBodies int) []Body {
	bodies := make([]Body, numSmallBodies+1)

	// Set up central body at the origin
	bodies[0] = Body{
		Position: Vector3{},
		Velocity: Vector3{},
		Mass:     centralMass,
	}

	// Initialize random seed
	rand.Seed(time.Now().UnixNano())

	// Determine the number of goroutines based on CPU cores
	numCPU := runtime.NumCPU()
	numGoroutines := numCPU * 2 // Use 2x the number of CPU cores for better utilization

	// Ensure we don't create more goroutines than bodies
	if numGoroutines > numSmallBodies {
		numGoroutines = numSmallBodies
	}

	// Calculate chunk size for each goroutine
	chunkSize := numSmallBodies / numGoroutines

	// Use a wait group to synchronize goroutines
	var wg sync.WaitGroup
	wg.Add(numGoroutines)

	// Launch goroutines to initialize bodies in parallel
	for g := 0; g < numGoroutines; g++ {
		// Calculate start and end indices for this goroutine
		start := g*chunkSize + 1 // +1 to skip the central body
		end := start + chunkSize
		if g == numGoroutines-1 {
			end = numSmallBodies + 1 // Make sure the last goroutine handles any remainder
		}

		go func(start, end int) {
			defer wg.Done()

			// Create a local RNG for this goroutine to avoid contention
			localRand := rand.New(rand.NewSource(time.Now().UnixNano() + int64(start)))

			// Set up smaller bodies in circular orbits around the central body
			for i := start; i < end; i++ {
				// Random orbital radius between 1 AU and 5 AU (in meters)
				// 1 AU ≈ 1.496e11 meters
				radius := (1.0 + 4.0*localRand.Float64()) * 1.496e11

				// Random angle for position in the orbit
				theta := localRand.Float64() * 2 * math.Pi
				phi := (localRand.Float64() - 0.5) * 0.1 * math.Pi // Near-planar with some variation

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
		}(start, end)
	}

	// Wait for all goroutines to complete
	wg.Wait()

	return bodies
}

func main() {
	fmt.Println("Initializing N-body simulation with Barnes-Hut algorithm...")
	startTime := time.Now()

	// Set maximum number of CPUs to use
	numCPU := runtime.NumCPU()
	runtime.GOMAXPROCS(numCPU)
	fmt.Printf("Using %d CPU cores\n", numCPU)

	// Initialize system with one central body and numSmallBodies smaller bodies
	fmt.Printf("Creating system with %d bodies...\n", numSmallBodies+1)
	bodies := InitializeSystem(numSmallBodies)

	fmt.Printf("System initialized with %d bodies\n", len(bodies))
	fmt.Printf("Central body mass: %.3e kg\n", bodies[0].Mass)
	fmt.Printf("Small body mass: %.3e kg\n", bodies[1].Mass)

	// Calculate initial energy
	// fmt.Println("Calculating initial energy...")
	// initialEnergy := CalculateTotalEnergy(bodies)
	// fmt.Printf("Initial total energy: %.6e J\n", initialEnergy)

	// Run simulation for numSteps steps
	fmt.Printf("Running parallel simulation with Barnes-Hut algorithm (theta=%.2f) for %d steps...\n", theta, numSteps)
	simulationStartTime := time.Now()

	// Create a channel for progress updates
	progressChan := make(chan int, numSteps/10)

	// Start a goroutine to display progress
	go func() {
		for progress := range progressChan {
			fmt.Printf("Progress: %d%%\n", progress)
		}
	}()

	// Create indices for the initial k-d tree
	indices := make([]int, len(bodies))
	for i := range indices {
		indices[i] = i
	}

	// Run the simulation
	for step := 1; step <= numSteps; step++ {
		// Rebuild the k-d tree every step
		root := BuildKDTree(bodies, indices, 0)

		// Perform integration step using the k-d tree
		KickStep(bodies, root)

		// Send progress updates at 10% intervals
		if step%(numSteps/10) == 0 {
			progressChan <- step * 100 / numSteps
		}
	}

	// Close the progress channel
	close(progressChan)

	simulationTime := time.Since(simulationStartTime)
	fmt.Printf("Simulation steps completed in %v\n", simulationTime)

	// Calculate final energy
	// fmt.Println("Calculating final energy...")
	// finalEnergy := CalculateTotalEnergy(bodies)
	// fmt.Printf("Final total energy: %.6e J\n", finalEnergy)

	// // Calculate energy conservation
	// energyDiff := (finalEnergy - initialEnergy) / initialEnergy
	// fmt.Printf("Energy difference: %.6e (%.6f%%)\n", finalEnergy-initialEnergy, energyDiff*100)
	fmt.Printf("Body[0] %e %e %e", bodies[0].Position.X, bodies[0].Position.Y, bodies[0].Position.Z)

	// Performance metrics
	totalElapsedTime := time.Since(startTime)
	bodiesToMillion := float64(numSmallBodies+1) / 1000000.0
	stepsPerSecond := float64(numSteps) / simulationTime.Seconds()
	bodyStepsPerSecond := stepsPerSecond * float64(numSmallBodies+1)

	fmt.Printf("Total simulation time: %v\n", totalElapsedTime)
	fmt.Printf("Performance: %.2f steps/sec\n", stepsPerSecond)
	fmt.Printf("Body-steps per second: %.2e (%.2f million body-steps/sec)\n",
		bodyStepsPerSecond, bodyStepsPerSecond/1000000.0)
	fmt.Printf("Simulation of %.2f million bodies for %d steps\n",
		bodiesToMillion, numSteps)
}
