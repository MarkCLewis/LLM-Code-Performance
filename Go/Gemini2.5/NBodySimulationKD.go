package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort" // Added for sorting during k-D tree build
	"sync"
	"time"
)

// --- Constants ---
const (
	G                 = 1.0    // Gravitational constant
	dt                = 0.001  // Time step
	epsilon           = 0.01   // Softening factor
	numSteps          = 10     // Number of simulation steps
	numOrbitingBodies = 100000 // Number of orbiting bodies (Feasible with N log N)
	// numOrbitingBodies = 5000 // Smaller number for quicker testing
	theta = 0.3 // Barnes-Hut opening angle parameter
)

// --- Data Structures ---

// Vector3D represents a point or vector in 3D space
type Vector3D struct {
	X, Y, Z float64
}

// Body represents a celestial body
type Body struct {
	Mass     float64
	Position Vector3D
	Velocity Vector3D
	// Acceleration stored temporarily during step calculation
	Acceleration Vector3D // Store acceleration here instead of separate slice
}

// BoundingBox defines an axis-aligned bounding box
type BoundingBox struct {
	Min, Max Vector3D
}

// KDNode represents a node in the k-D tree
type KDNode struct {
	Axis          int         // Splitting axis (0=X, 1=Y, 2=Z)
	SplitValue    float64     // Value at which space is split
	Bounds        BoundingBox // Spatial bounds of this node
	CenterOfMass  Vector3D    // Center of mass of particles in this node
	TotalMass     float64     // Total mass of particles in this node
	ParticleIndex int         // Index of the particle if this is a leaf node (-1 otherwise)
	LeftChild     *KDNode     // Pointer to the left/below child
	RightChild    *KDNode     // Pointer to the right/above child
	NodeSize      float64     // Characteristic size of the node (e.g., width along split axis)
}

// --- Vector Operations --- (Same as before)
// Add, Subtract, Scale, MagnitudeSquared, Magnitude, Normalize
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

// Component returns the vector component along a given axis (0=X, 1=Y, 2=Z)
func (v Vector3D) Component(axis int) float64 {
	switch axis {
	case 0:
		return v.X
	case 1:
		return v.Y
	case 2:
		return v.Z
	default:
		panic(fmt.Sprintf("Invalid axis: %d", axis))
	}
}

// --- Bounding Box Operations ---

// Contains checks if a point is within the bounding box
func (bb BoundingBox) Contains(p Vector3D) bool {
	return p.X >= bb.Min.X && p.X <= bb.Max.X &&
		p.Y >= bb.Min.Y && p.Y <= bb.Max.Y &&
		p.Z >= bb.Min.Z && p.Z <= bb.Max.Z
}

// Expand includes another point or bounding box
func (bb *BoundingBox) Expand(other BoundingBox) {
	bb.Min.X = math.Min(bb.Min.X, other.Min.X)
	bb.Min.Y = math.Min(bb.Min.Y, other.Min.Y)
	bb.Min.Z = math.Min(bb.Min.Z, other.Min.Z)
	bb.Max.X = math.Max(bb.Max.X, other.Max.X)
	bb.Max.Y = math.Max(bb.Max.Y, other.Max.Y)
	bb.Max.Z = math.Max(bb.Max.Z, other.Max.Z)
}

// ExpandWithPoint includes a single point
func (bb *BoundingBox) ExpandWithPoint(p Vector3D) {
	bb.Min.X = math.Min(bb.Min.X, p.X)
	bb.Min.Y = math.Min(bb.Min.Y, p.Y)
	bb.Min.Z = math.Min(bb.Min.Z, p.Z)
	bb.Max.X = math.Max(bb.Max.X, p.X)
	bb.Max.Y = math.Max(bb.Max.Y, p.Y)
	bb.Max.Z = math.Max(bb.Max.Z, p.Z)
}

// Calculate initial bounding box for all bodies
func calculateOverallBounds(bodies []Body) BoundingBox {
	if len(bodies) == 0 {
		return BoundingBox{} // Or handle error
	}
	bounds := BoundingBox{Min: bodies[0].Position, Max: bodies[0].Position}
	for i := 1; i < len(bodies); i++ {
		bounds.ExpandWithPoint(bodies[i].Position)
	}
	// Add a small margin to avoid issues with particles exactly on the boundary
	margin := 1e-5
	bounds.Min = bounds.Min.Subtract(Vector3D{margin, margin, margin})
	bounds.Max = bounds.Max.Add(Vector3D{margin, margin, margin})
	return bounds
}

// --- k-D Tree Construction ---

// Interface for sorting particle indices based on position along an axis
type bodySorter struct {
	indices []int
	bodies  []Body
	axis    int
}

func (s *bodySorter) Len() int      { return len(s.indices) }
func (s *bodySorter) Swap(i, j int) { s.indices[i], s.indices[j] = s.indices[j], s.indices[i] }
func (s *bodySorter) Less(i, j int) bool {
	return s.bodies[s.indices[i]].Position.Component(s.axis) < s.bodies[s.indices[j]].Position.Component(s.axis)
}

// buildKDTree recursively builds the tree
func buildKDTree(particleIndices []int, bodies []Body, depth int, currentBounds BoundingBox) *KDNode {
	n := len(particleIndices)
	if n == 0 {
		return nil
	}

	node := &KDNode{
		Bounds:        currentBounds, // Start with the parent's bounds
		ParticleIndex: -1,            // Assume internal node initially
	}

	// Base Case: Leaf node
	if n == 1 {
		idx := particleIndices[0]
		node.ParticleIndex = idx
		node.TotalMass = bodies[idx].Mass
		node.CenterOfMass = bodies[idx].Position
		// Leaf bounds are tight around the particle (using parent bounds is also ok)
		node.Bounds = BoundingBox{Min: bodies[idx].Position, Max: bodies[idx].Position}
		node.NodeSize = 0 // Size is irrelevant for leaf force calc
		return node
	}

	// Recursive Step: Internal node
	node.Axis = depth % 3 // Cycle through X, Y, Z

	// Sort indices based on position along the current axis to find the median
	sorter := &bodySorter{indices: particleIndices, bodies: bodies, axis: node.Axis}
	sort.Sort(sorter) // O(N log N) sort

	medianIndex := n / 2
	medianParticleIdx := particleIndices[medianIndex]
	node.SplitValue = bodies[medianParticleIdx].Position.Component(node.Axis)

	// --- Calculate Bounds for Children ---
	leftBounds := currentBounds
	leftBounds.Max = leftBounds.Max.withComponent(node.Axis, node.SplitValue)

	rightBounds := currentBounds
	rightBounds.Min = rightBounds.Min.withComponent(node.Axis, node.SplitValue)

	// Handle potential empty partitions if multiple particles have the exact same coordinate
	leftIndices := particleIndices[:medianIndex]
	rightIndices := particleIndices[medianIndex:] // Include median in right side typically

	// Build children recursively
	// Use separate goroutines? Maybe later, adds complexity with COM calculation sync
	node.LeftChild = buildKDTree(leftIndices, bodies, depth+1, leftBounds)
	node.RightChild = buildKDTree(rightIndices, bodies, depth+1, rightBounds)

	// --- Calculate COM and Total Mass for this internal node ---
	node.TotalMass = 0
	comNumerator := Vector3D{0, 0, 0}
	if node.LeftChild != nil {
		node.TotalMass += node.LeftChild.TotalMass
		comNumerator = comNumerator.Add(node.LeftChild.CenterOfMass.Scale(node.LeftChild.TotalMass))
		node.Bounds.Expand(node.LeftChild.Bounds) // Expand bounds based on children
	}
	if node.RightChild != nil {
		node.TotalMass += node.RightChild.TotalMass
		comNumerator = comNumerator.Add(node.RightChild.CenterOfMass.Scale(node.RightChild.TotalMass))
		node.Bounds.Expand(node.RightChild.Bounds) // Expand bounds based on children
	}

	if node.TotalMass > 1e-12 { // Avoid division by zero if node is effectively empty
		node.CenterOfMass = comNumerator.Scale(1.0 / node.TotalMass)
	} else {
		// Handle empty node case - COM can be arbitrary, maybe center of bounds
		node.CenterOfMass = node.Bounds.Min.Add(node.Bounds.Max).Scale(0.5)
	}

	// Calculate node size (width along the splitting axis)
	node.NodeSize = node.Bounds.Max.Component(node.Axis) - node.Bounds.Min.Component(node.Axis)

	return node
}

// Helper to set a component of a vector (returns a new vector)
func (v Vector3D) withComponent(axis int, value float64) Vector3D {
	switch axis {
	case 0:
		return Vector3D{value, v.Y, v.Z}
	case 1:
		return Vector3D{v.X, value, v.Z}
	case 2:
		return Vector3D{v.X, v.Y, value}
	default:
		panic(fmt.Sprintf("Invalid axis: %d", axis))
	}
}

// --- Force Calculation using k-D Tree ---

// calculateForceRecursive traverses the tree to calculate force on targetParticleIdx
func calculateForceRecursive(node *KDNode, targetParticleIdx int, bodies []Body, thetaValue, gravitationalConstant, softeningFactor float64) Vector3D {
	if node == nil {
		return Vector3D{0, 0, 0} // Empty node contributes no force
	}

	// If the node is a leaf
	if node.ParticleIndex != -1 {
		// Don't calculate force of particle on itself
		if node.ParticleIndex == targetParticleIdx {
			return Vector3D{0, 0, 0}
		}
		// Calculate direct force between targetParticle and the leaf particle
		bodyI := bodies[targetParticleIdx]
		bodyJ := bodies[node.ParticleIndex]
		deltaPos := bodyJ.Position.Subtract(bodyI.Position) // Vector from i to j
		distSq := deltaPos.MagnitudeSquared()
		distSq += softeningFactor * softeningFactor // Add softening squared
		invDistCubed := 1.0 / math.Pow(distSq, 1.5)

		forceScale := gravitationalConstant * bodyJ.Mass * invDistCubed // Force = G * m_j / dist^3 * deltaPos (scaled by m_i later)
		return deltaPos.Scale(forceScale)
	}

	// If the node is internal, apply Barnes-Hut criterion
	deltaPosCOM := node.CenterOfMass.Subtract(bodies[targetParticleIdx].Position) // Vector from particle to COM
	distSqCOM := deltaPosCOM.MagnitudeSquared()
	distCOM := math.Sqrt(distSqCOM)

	// Criterion: s / d < theta. If true, use approximation.
	// Avoid division by zero if particle is exactly at COM (unlikely with softening)
	if distCOM > 1e-9 && (node.NodeSize/distCOM) < thetaValue {
		// Node is far enough away, use its Center of Mass approximation
		distSqCOM += softeningFactor * softeningFactor // Add softening squared
		invDistCubed := 1.0 / math.Pow(distSqCOM, 1.5)
		forceScale := gravitationalConstant * node.TotalMass * invDistCubed // Force = G * M_node / dist^3 * deltaPos (scaled by m_i later)
		return deltaPosCOM.Scale(forceScale)
	} else {
		// Node is too close, recurse into children
		forceLeft := calculateForceRecursive(node.LeftChild, targetParticleIdx, bodies, thetaValue, gravitationalConstant, softeningFactor)
		forceRight := calculateForceRecursive(node.RightChild, targetParticleIdx, bodies, thetaValue, gravitationalConstant, softeningFactor)
		return forceLeft.Add(forceRight)
	}
}

// --- Simulation Step ---

// calculateForceKDTreeChunk calculates acceleration for a subset of bodies using the k-D tree
func calculateForceKDTreeChunk(treeRoot *KDNode, bodies []Body, start, end int, wg *sync.WaitGroup, thetaValue, gravitationalConstant, softeningFactor float64) {
	defer wg.Done()
	for i := start; i < end; i++ {
		// Calculate the force sum using the tree traversal
		// The result is Force / m_i
		forceSumFactor := calculateForceRecursive(treeRoot, i, bodies, thetaValue, gravitationalConstant, softeningFactor)
		// Actual acceleration a = F / m_i. The recursive func returns F/m_i * m_node or F/m_i * m_leaf
		// We need to multiply by target particle mass m_i here.
		// NO - check the force calculation: acc = G * m_j / dist^3 * deltaPos. Mass m_i is not needed here.
		bodies[i].Acceleration = forceSumFactor // Store acceleration directly in body struct
	}
}

// updateVelocityChunk updates velocity for a subset of bodies
func updateVelocityChunk(bodies []Body, start, end int, wg *sync.WaitGroup, timeStep float64) {
	defer wg.Done()
	for i := start; i < end; i++ {
		bodies[i].Velocity = bodies[i].Velocity.Add(bodies[i].Acceleration.Scale(timeStep))
	}
}

// updatePositionChunk updates position for a subset of bodies
func updatePositionChunk(bodies []Body, start, end int, wg *sync.WaitGroup, timeStep float64) {
	defer wg.Done()
	for i := start; i < end; i++ {
		bodies[i].Position = bodies[i].Position.Add(bodies[i].Velocity.Scale(timeStep))
	}
}

// simulateStepParallel performs one step using k-D tree and multiple goroutines
func simulateStepParallelKDTree(bodies []Body, numWorkers int, timeStep, gravitationalConstant, softeningFactor, thetaValue float64) *KDNode { // Return tree root for potential reuse/analysis
	n := len(bodies)
	var wg sync.WaitGroup

	// --- 0. Build k-D Tree (Sequential for simplicity, could be parallelized) ---
	// buildStart := time.Now()
	initialBounds := calculateOverallBounds(bodies)
	particleIndices := make([]int, n)
	for i := range particleIndices {
		particleIndices[i] = i
	}
	treeRoot := buildKDTree(particleIndices, bodies, 0, initialBounds)
	// buildDuration := time.Since(buildStart)

	// --- 1. Parallel Force/Acceleration Calculation using k-D Tree ---
	// forceStart := time.Now()
	wg.Add(numWorkers)
	chunkSize := (n + numWorkers - 1) / numWorkers
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
		go calculateForceKDTreeChunk(treeRoot, bodies, start, end, &wg, thetaValue, gravitationalConstant, softeningFactor)
	}
	wg.Wait()
	// forceDuration := time.Since(forceStart)

	// --- 2. Parallel Velocity Update (Kick) ---
	// velStart := time.Now()
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
		go updateVelocityChunk(bodies, start, end, &wg, timeStep)
	}
	wg.Wait()
	// velDuration := time.Since(velStart)

	// --- 3. Parallel Position Update (Step) ---
	// posStart := time.Now()
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
	wg.Wait()
	// posDuration := time.Since(posStart)

	// Optionally return durations for profiling
	// fmt.Printf(" Step Timings: Build=%.3fs Force=%.3fs Vel=%.3fs Pos=%.3fs\n", buildDuration.Seconds(), forceDuration.Seconds(), velDuration.Seconds(), posDuration.Seconds())

	return treeRoot // Might be useful, but typically discarded next step
}

// --- Initialization --- (Same as before)
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
		radius := minRadius + rand.Float64()*(maxRadius-minRadius)
		theta := 2 * math.Pi * rand.Float64()
		phi := math.Acos(1 - 2*rand.Float64())
		posX := radius * math.Sin(phi) * math.Cos(theta)
		posY := radius * math.Sin(phi) * math.Sin(theta)
		posZ := radius * math.Cos(phi)
		positionVec := Vector3D{posX, posY, posZ}

		orbitSpeed := math.Sqrt(G * centralMass / radius)

		helperVec := Vector3D{rand.Float64()*2 - 1, rand.Float64()*2 - 1, rand.Float64()*2 - 1}.Normalize()
		dotProd := positionVec.Normalize().Dot(helperVec) // Use Dot product helper if defined
		if math.Abs(dotProd) > 0.99 {
			helperVec = Vector3D{0, 1, 0}
			if math.Abs(positionVec.Normalize().Y) > 0.99 {
				helperVec = Vector3D{1, 0, 0}
			}
		}

		velDirection := positionVec.Cross(helperVec).Normalize() // Use Cross product helper if defined
		if velDirection.MagnitudeSquared() < 1e-12 {             // Handle cases where cross product is zero
			// If position is along Z axis, velocity can be in XY plane
			if math.Abs(positionVec.X) < 1e-6 && math.Abs(positionVec.Y) < 1e-6 {
				velDirection = Vector3D{1, 0, 0}
			} else { // Otherwise, try another perpendicular direction, e.g., cross with Z axis
				zAxis := Vector3D{0, 0, 1}
				velDirection = positionVec.Cross(zAxis).Normalize()
				if velDirection.MagnitudeSquared() < 1e-12 { // If position is along Z axis
					velDirection = Vector3D{1, 0, 0} // Default to X direction
				}
			}
		}

		velocityVec := velDirection.Scale(orbitSpeed)

		bodies[i] = Body{
			Mass:     orbitingMass,
			Position: positionVec,
			Velocity: velocityVec,
		}
	}
	return bodies
}

// Add Dot and Cross product helpers to Vector3D
func (v Vector3D) Dot(other Vector3D) float64 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

func (v Vector3D) Cross(other Vector3D) Vector3D {
	return Vector3D{
		v.Y*other.Z - v.Z*other.Y,
		v.Z*other.X - v.X*other.Z,
		v.X*other.Y - v.Y*other.X,
	}
}

// --- Energy Calculation --- (Same as before, but less accurate due to approximation)
func calculateTotalEnergy(bodies []Body) float64 {
	// ... (implementation remains the same) ...
	kineticEnergy := 0.0
	potentialEnergy := 0.0

	// Kinetic Energy: Sum(0.5 * m * v^2)
	for i := range bodies {
		kineticEnergy += 0.5 * bodies[i].Mass * bodies[i].Velocity.MagnitudeSquared()
	}

	// Potential Energy: Sum(-G * m_i * m_j / |r_i - r_j|) for i < j
	// NOTE: This calculates the *exact* potential energy. The simulation uses an approximation.
	// So, the total energy calculated here might not appear well-conserved by the simulation step.
	for i := 0; i < len(bodies); i++ {
		for j := i + 1; j < len(bodies); j++ {
			deltaPos := bodies[i].Position.Subtract(bodies[j].Position)
			distance := deltaPos.Magnitude()
			if distance > 0 {
				potentialEnergy -= G * bodies[i].Mass * bodies[j].Mass / distance
			}
		}
	}
	return kineticEnergy + potentialEnergy

}

// --- Main Function ---
func main() {
	numWorkers := runtime.NumCPU()
	runtime.GOMAXPROCS(numWorkers)

	fmt.Printf("Starting N-body simulation (Parallel k-D Tree, theta=%.2f)...\n", theta)
	fmt.Printf("Using %d CPU cores/workers.\n", numWorkers)
	fmt.Printf("Number of bodies: %d (+1 central body)\n", numOrbitingBodies)
	fmt.Printf("Number of steps: %d\n", numSteps)
	fmt.Printf("Time step (dt): %.4f\n", dt)
	if numOrbitingBodies >= 100000 {
		fmt.Printf("N=%d is large, tree building and simulation will take time, but faster than O(N^2).\n", numOrbitingBodies+1)
	}

	// --- Initialization ---
	startTime := time.Now()
	centralMass := 10000.0
	orbitingMass := 0.01
	minRadius := 5.0
	maxRadius := 50.0

	bodies := initializeCircularOrbitSystem(numOrbitingBodies, centralMass, orbitingMass, minRadius, maxRadius)
	initializationTime := time.Since(startTime)
	fmt.Printf("System initialized in %v\n", initializationTime)

	// --- Initial Energy Calculation ---
	initialEnergy := calculateTotalEnergy(bodies)
	fmt.Printf("Initial Total Energy (Exact): %.5e\n", initialEnergy)
	fmt.Println("Note: Simulation uses approximated forces (theta > 0), so exact energy conservation is not expected.")

	// --- Simulation Loop ---
	fmt.Println("Running simulation...")
	simulationStartTime := time.Now()

	for step := 0; step < numSteps; step++ {
		stepStart := time.Now()
		_ = simulateStepParallelKDTree(bodies, numWorkers, dt, G, epsilon, theta) // Ignore returned tree root for now
		stepDuration := time.Since(stepStart)

		if (step+1)%10 == 0 || step == 0 || step == numSteps-1 { // Print progress
			fmt.Printf("Step %d/%d completed (Step duration: %v, Elapsed: %v)\n",
				step+1, numSteps, stepDuration, time.Since(simulationStartTime))
		}
	}

	simulationTime := time.Since(simulationStartTime)
	fmt.Printf("Simulation completed in %v\n", simulationTime)

	// --- Final Energy Calculation ---
	finalEnergy := calculateTotalEnergy(bodies)
	fmt.Printf("Final Total Energy (Exact):   %.5e\n", finalEnergy)

	// --- Energy Verification ---
	energyChange := finalEnergy - initialEnergy
	// Avoid division by zero if initial energy is near zero
	relativeEnergyChange := 0.0
	if math.Abs(initialEnergy) > 1e-12 {
		relativeEnergyChange = math.Abs(energyChange / initialEnergy)
	}

	fmt.Printf("Absolute Energy Change (Exact): %.5e\n", energyChange)
	fmt.Printf("Relative Energy Change (Exact): %.5e (%.4f%%)\n", relativeEnergyChange, relativeEnergyChange*100.0)
	fmt.Println("Energy drift is expected due to force approximation (theta > 0).")

	totalTime := time.Since(startTime)
	fmt.Printf("Total execution time: %v\n", totalTime)
}
