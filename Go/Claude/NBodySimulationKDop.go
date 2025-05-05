package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"
	"time"
	"sync/atomic"
)

// Constants
const (
	G              = 6.67430e-11 // Gravitational constant (m^3 kg^-1 s^-2)
	dt             = 86400.0     // Time step in seconds (1 day)
	softening      = 1.0e7       // Softening parameter to avoid singularities
	centralMass    = 1.989e30    // Mass of central body (kg) - roughly solar mass
	planetMass     = 1.0e24      // Mass of smaller bodies (kg)
	numSteps       = 1000        // Number of simulation steps
	numSmallBodies = 1000000     // Number of small bodies
	theta          = 0.3         // Barnes-Hut opening angle parameter
	leafMaxBodies  = 16          // Maximum bodies in a leaf node
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
	Min, Max   Vector3    // Bounding box
	Center     Vector3    // Center of mass
	Size       float64    // Size of the node (maximum dimension of bounding box)
	TotalMass  float64    // Total mass of bodies in this node
	
	// Tree structure
	Left, Right *KDNode   // Child nodes
	Bodies      []int     // Indices of bodies (only for leaf nodes)
	IsLeaf      bool      // Whether this is a leaf node
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

// Pre-allocate these vectors to reduce allocations in hot paths
var (
	nodePool      sync.Pool
	indicesPool   = sync.Pool{New: func() interface{} { return make([]int, 0, 1024) }}
	bodyIndices   []int // Reuse this slice across simulation steps
	accelerations []Vector3
	thetaSquared  = theta * theta
)

// InitNodePool initializes the node pool
func InitNodePool() {
	nodePool = sync.Pool{
		New: func() interface{} {
			return &KDNode{
				Min:       Vector3{X: math.MaxFloat64, Y: math.MaxFloat64, Z: math.MaxFloat64},
				Max:       Vector3{X: -math.MaxFloat64, Y: -math.MaxFloat64, Z: -math.MaxFloat64},
				Bodies:    make([]int, 0, leafMaxBodies),
			}
		},
	}
}

// GetNode retrieves a node from the pool
func GetNode() *KDNode {
	node := nodePool.Get().(*KDNode)
	// Reset node properties
	node.Min = Vector3{X: math.MaxFloat64, Y: math.MaxFloat64, Z: math.MaxFloat64}
	node.Max = Vector3{X: -math.MaxFloat64, Y: -math.MaxFloat64, Z: -math.MaxFloat64}
	node.Center = Vector3{}
	node.Size = 0
	node.TotalMass = 0
	node.Left = nil
	node.Right = nil
	node.Bodies = node.Bodies[:0]
	node.IsLeaf = false
	return node
}

// PutNode returns a node to the pool
func PutNode(node *KDNode) {
	if node == nil {
		return
	}
	PutNode(node.Left)
	PutNode(node.Right)
	nodePool.Put(node)
}

// BuildKDTree constructs a k-d tree from a set of bodies
func BuildKDTree(bodies []Body, indices []int, depth int) *KDNode {
	if len(indices) == 0 {
		return nil
	}
	
	// Create a new node from the pool
	node := GetNode()
	
	// Find the bounding box and calculate center of mass
	weightedCenter := Vector3{}
	for _, idx := range indices {
		body := &bodies[idx]
		
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
	sizeX := node.Max.X - node.Min.X
	sizeY := node.Max.Y - node.Min.Y
	sizeZ := node.Max.Z - node.Min.Z
	node.Size = math.Max(sizeX, math.Max(sizeY, sizeZ))
	
	// Base case: if there are few enough bodies, make this a leaf node
	if len(indices) <= leafMaxBodies {
		node.Bodies = append(node.Bodies, indices...)
		node.IsLeaf = true
		return node
	}
	
	// Determine split axis (cycle through X, Y, Z based on depth)
	axis := depth % 3
	
	// Sort indices based on position along the split axis
	if axis == 0 { // X-axis
		sort.Slice(indices, func(i, j int) bool {
			return bodies[indices[i]].Position.X < bodies[indices[j]].Position.X
		})
	} else if axis == 1 { // Y-axis
		sort.Slice(indices, func(i, j int) bool {
			return bodies[indices[i]].Position.Y < bodies[indices[j]].Position.Y
		})
	} else { // Z-axis
		sort.Slice(indices, func(i, j int) bool {
			return bodies[indices[i]].Position.Z < bodies[indices[j]].Position.Z
		})
	}
	
	// Find median index for splitting
	median := len(indices) / 2
	
	// Recursively build left and right subtrees
	node.Left = BuildKDTree(bodies, indices[:median], depth+1)
	node.Right = BuildKDTree(bodies, indices[median:], depth+1)
	node.IsLeaf = false
	
	return node
}

// CalculateForceFromNode computes the gravitational force on a body from a k-d tree node
func CalculateForceFromNode(body *Body, node *KDNode, bodies []Body) Vector3 {
	if node == nil {
		return Vector3{}
	}
	
	// Calculate distance from body to node's center of mass
	dx := node.Center.X - body.Position.X
	dy := node.Center.Y - body.Position.Y
	dz := node.Center.Z - body.Position.Z
	distSqr := dx*dx + dy*dy + dz*dz
	
	// If this is a leaf node or the node is sufficiently far away (Barnes-Hut approximation)
	// s/d < θ, where s is the node size, d is the distance, and θ is the opening angle parameter
	if node.IsLeaf || (node.Size*node.Size < thetaSquared*distSqr) {
		if node.IsLeaf {
			// For leaf nodes, calculate force directly from each body
			var forceX, forceY, forceZ float64
			for _, idx := range node.Bodies {
				// Skip self-interaction
				if &bodies[idx] == body {
					continue
				}
				
				// Vector from body to other body
				otherBody := &bodies[idx]
				rx := otherBody.Position.X - body.Position.X
				ry := otherBody.Position.Y - body.Position.Y
				rz := otherBody.Position.Z - body.Position.Z
				distSqr := rx*rx + ry*ry + rz*rz + softening*softening
				
				// Precalculate terms for efficiency
				invDist := 1.0 / math.Sqrt(distSqr)
				invDistCube := invDist * invDist * invDist
				
				// Force magnitude: F = G*m1*m2/r^2
				// Force direction: F * (r/|r|)
				factor := G * body.Mass * otherBody.Mass * invDistCube
				
				// Add to total force components
				forceX += rx * factor
				forceY += ry * factor
				forceZ += rz * factor
			}
			return Vector3{X: forceX, Y: forceY, Z: forceZ}
		} else {
			// For distant nodes, use center of mass approximation
			// Avoid division by zero with softening
			distSqr += softening*softening
			
			// Precalculate terms for efficiency
			invDist := 1.0 / math.Sqrt(distSqr)
			invDistCube := invDist * invDist * invDist
			
			// Force magnitude: F = G*m1*m2/r^2
			// Force direction: F * (r/|r|)
			factor := G * body.Mass * node.TotalMass * invDistCube
			
			return Vector3{
				X: dx * factor,
				Y: dy * factor,
				Z: dz * factor,
			}
		}
	}
	
	// Node is too close for approximation, traverse children
	leftForce := CalculateForceFromNode(body, node.Left, bodies)
	rightForce := CalculateForceFromNode(body, node.Right, bodies)
	
	// Combine forces from both children
	return Vector3{
		X: leftForce.X + rightForce.X,
		Y: leftForce.Y + rightForce.Y,
		Z: leftForce.Z + rightForce.Z,
	}
}

// CalculateAcceleration calculates the acceleration of a body due to gravitational forces using k-d tree
func CalculateAcceleration(body *Body, root *KDNode, bodies []Body) Vector3 {
	// Calculate the total gravitational force on the body from the k-d tree
	totalForce := CalculateForceFromNode(body, root, bodies)
	
	// a = F/m
	invMass := 1.0 / body.Mass
	return Vector3{
		X: totalForce.X * invMass,
		Y: totalForce.Y * invMass,
		Z: totalForce.Z * invMass,
	}
}

// WorkerPool represents a pool of workers for parallel computation
type WorkerPool struct {
	numWorkers   int
	tasks        chan func()
	wg           sync.WaitGroup
	initialized  bool
}

// NewWorkerPool creates and initializes a new worker pool
func NewWorkerPool(numWorkers int) *WorkerPool {
	pool := &WorkerPool{
		numWorkers:  numWorkers,
		tasks:       make(chan func(), numWorkers*2), // Buffer channel for better performance
		initialized: false,
	}
	
	return pool
}

// Initialize starts the workers
func (p *WorkerPool) Initialize() {
	if p.initialized {
		return
	}
	
	// Start workers
	for i := 0; i < p.numWorkers; i++ {
		go func() {
			for task := range p.tasks {
				task()
				p.wg.Done()
			}
		}()
	}
	
	p.initialized = true
}

// Submit adds a task to the worker pool
func (p *WorkerPool) Submit(task func()) {
	p.wg.Add(1)
	p.tasks <- task
}

// Wait waits for all tasks to complete
func (p *WorkerPool) Wait() {
	p.wg.Wait()
}

// Close stops all workers
func (p *WorkerPool) Close() {
	close(p.tasks)
	p.initialized = false
}

// Global worker pool
var workerPool *WorkerPool

// KickStep performs a single kick-step (leapfrog) integration using a worker pool
// and a k-d tree for force calculations
func KickStep(bodies []Body, root *KDNode) {
	// Calculate accelerations for all bodies
	for i := 0; i < len(bodies); i++ {
		bodyIndex := i // Capture loop variable
		workerPool.Submit(func() {
			accelerations[bodyIndex] = CalculateAcceleration(&bodies[bodyIndex], root, bodies)
		})
	}
	
	// Wait for all acceleration calculations to complete
	workerPool.Wait()
	
	// Update positions and velocities for all bodies
	for i := 0; i < len(bodies); i++ {
		bodyIndex := i // Capture loop variable
		workerPool.Submit(func() {
			// Update velocity (kick)
			bodies[bodyIndex].Velocity.X += accelerations[bodyIndex].X * dt
			bodies[bodyIndex].Velocity.Y += accelerations[bodyIndex].Y * dt
			bodies[bodyIndex].Velocity.Z += accelerations[bodyIndex].Z * dt
			
			// Update position (drift)
			bodies[bodyIndex].Position.X += bodies[bodyIndex].Velocity.X * dt
			bodies[bodyIndex].Position.Y += bodies[bodyIndex].Velocity.Y * dt
			bodies[bodyIndex].Position.Z += bodies[bodyIndex].Velocity.Z * dt
		})
	}
	
	// Wait for all position/velocity updates to complete
	workerPool.Wait()
}

// CalculateTotalEnergy calculates the total energy of the system (kinetic + potential) using worker pool
func CalculateTotalEnergy(bodies []Body) float64 {
	// Use atomic operations for thread-safe increments
	var kineticEnergyBits uint64
	var potentialEnergyBits uint64
	
	// Calculate kinetic energy in parallel
	for i := 0; i < len(bodies); i++ {
		bodyIndex := i // Capture loop variable
		workerPool.Submit(func() {
			// Calculate kinetic energy: KE = 0.5 * m * v^2
			body := &bodies[bodyIndex]
			speed2 := body.Velocity.X*body.Velocity.X + 
			          body.Velocity.Y*body.Velocity.Y + 
			          body.Velocity.Z*body.Velocity.Z
			ke := 0.5 * body.Mass * speed2
			
			// Atomically add to total kinetic energy
			atomic.AddUint64(&kineticEnergyBits, math.Float64bits(ke))
		})
	}
	
	// Calculate potential energy in parallel
	// Using a more efficient distribution of work
	totalBodies := len(bodies)
	pairsPerWorker := totalBodies / workerPool.numWorkers
	
	for w := 0; w < workerPool.numWorkers; w++ {
		startIdx := w * pairsPerWorker
		endIdx := (w + 1) * pairsPerWorker
		if w == workerPool.numWorkers-1 {
			endIdx = totalBodies // Make sure the last worker handles any remainder
		}
		
		workerPool.Submit(func() {
			localPE := 0.0
			
			// Calculate potential energy: PE = -G * m1 * m2 / r
			for i := startIdx; i < endIdx; i++ {
				for j := i + 1; j < totalBodies; j++ {
					body1 := &bodies[i]
					body2 := &bodies[j]
					
					// Calculate distance
					dx := body1.Position.X - body2.Position.X
					dy := body1.Position.Y - body2.Position.Y
					dz := body1.Position.Z - body2.Position.Z
					distSquared := dx*dx + dy*dy + dz*dz
					
					// Avoid division by zero
					if distSquared > 0 {
						dist := math.Sqrt(distSquared)
						localPE -= G * body1.Mass * body2.Mass / dist
					}
				}
			}
			
			// Atomically add to total potential energy
			atomic.AddUint64(&potentialEnergyBits, math.Float64bits(localPE))
		})
	}
	
	// Wait for all energy calculations to complete
	workerPool.Wait()
	
	// Convert bits back to float64
	kineticEnergy := math.Float64frombits(kineticEnergyBits)
	potentialEnergy := math.Float64frombits(potentialEnergyBits)
	
	return kineticEnergy + potentialEnergy
}

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
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	
	// Create a slice of local RNGs for worker threads
	localRNGs := make([]*rand.Rand, workerPool.numWorkers)
	for i := range localRNGs {
		localRNGs[i] = rand.New(rand.NewSource(rng.Int63()))
	}
	
	// Calculate chunk size for each worker
	chunkSize := numSmallBodies / workerPool.numWorkers
	
	// Initialize bodies in parallel using our worker pool
	for w := 0; w < workerPool.numWorkers; w++ {
		workerID := w
		start := w*chunkSize + 1 // +1 to skip the central body
		end := start + chunkSize
		if w == workerPool.numWorkers-1 {
			end = numSmallBodies + 1 // Make sure the last worker handles any remainder
		}
		
		workerPool.Submit(func() {
			// Use this worker's local RNG
			localRand := localRNGs[workerID]
			
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
				
				// Calculate orbital velocity for a circular orbit
				// v = sqrt(G*M/r)
				orbitalSpeed := math.Sqrt(G * centralMass / radius)
				
				// Velocity vector is perpendicular to the position vector
				// For a circular orbit in the xy-plane, it's (-y, x, 0) normalized * speed
				vx := -y
				vy := x
				vz := 0.0
				vmag := math.Sqrt(vx*vx + vy*vy + vz*vz)
				
				// Normalized velocity vector * orbital speed
				vx = (vx / vmag) * orbitalSpeed
				vy = (vy / vmag) * orbitalSpeed
				vz = (vz / vmag) * orbitalSpeed
				
				bodies[i] = Body{
					Position: Vector3{X: x, Y: y, Z: z},
					Velocity: Vector3{X: vx, Y: vy, Z: vz},
					Mass:     planetMass,
				}
			}
		})
	}
	
	// Wait for all initialization to complete
	workerPool.Wait()
	
	return bodies
}

func main() {
	fmt.Println("Initializing optimized N-body simulation with Barnes-Hut algorithm...")
	startTime := time.Now()
	
	// Set maximum number of CPUs to use
	numCPU := runtime.NumCPU()
	runtime.GOMAXPROCS(numCPU)
	fmt.Printf("Using %d CPU cores\n", numCPU)
	
	// Initialize node pool
	InitNodePool()
	
	// Initialize worker pool with twice as many workers as CPU cores
	workerPool = NewWorkerPool(numCPU * 2)
	workerPool.Initialize()
	defer workerPool.Close()
	
	// Pre-allocate global accelerations array
	accelerations = make([]Vector3, numSmallBodies+1)
	
	// Initialize bodyIndices for k-d tree
	bodyIndices = make([]int, numSmallBodies+1)
	for i := range bodyIndices {
		bodyIndices[i] = i
	}
	
	// Initialize system with one central body and numSmallBodies smaller bodies
	fmt.Printf("Creating system with %d bodies...\n", numSmallBodies+1)
	bodies := InitializeSystem(numSmallBodies)
	
	fmt.Printf("System initialized with %d bodies\n", len(bodies))
	fmt.Printf("Central body mass: %.3e kg\n", bodies[0].Mass)
	fmt.Printf("Small body mass: %.3e kg\n", bodies[1].Mass)
	
	// Calculate initial energy
	fmt.Println("Calculating initial energy...")
	initialEnergy := CalculateTotalEnergy(bodies)
	fmt.Printf("Initial total energy: %.6e J\n", initialEnergy)
	
	// Run simulation for numSteps steps
	fmt.Printf("Running optimized parallel simulation with Barnes-Hut algorithm (theta=%.2f) for %d steps...\n", theta, numSteps)
	simulationStartTime := time.Now()
	
	// Create atomic counter for progress tracking
	var completedSteps int32
	
	// Start a goroutine to display progress
	go func() {
		lastProgress := -1
		for {
			current := atomic.LoadInt32(&completedSteps)
			progress := int(current * 100 / numSteps)
			
			if progress != lastProgress && progress%10 == 0 {
				fmt.Printf("Progress: %d%%\n", progress)
				lastProgress = progress
			}
			
			if current >= numSteps {
				break
			}
			
			time.Sleep(100 * time.Millisecond)
		}
	}()
	
	// Run the simulation
	for step := 1; step <= numSteps; step++ {
		// Rebuild the k-d tree every step
		root := BuildKDTree(bodies, bodyIndices, 0)
		
		// Perform integration step using the k-d tree
		KickStep(bodies, root)
		
		// Return the tree nodes to the pool
		PutNode(root)
		
		// Update progress counter
		atomic.AddInt32(&completedSteps, 1)
	}
	
	simulationTime := time.Since(simulationStartTime)
	fmt.Printf("Simulation steps completed in %v\n", simulationTime)
	
	// Calculate final energy
	fmt.Println("Calculating final energy...")
	finalEnergy := CalculateTotalEnergy(bodies)
	fmt.Printf("Final total energy: %.6e J\n", finalEnergy)
	
	// Calculate energy conservation
	energyDiff := (finalEnergy - initialEnergy) / initialEnergy
	fmt.Printf("Energy difference: %.6e (%.6f%%)\n", finalEnergy-initialEnergy, energyDiff*100)
	
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