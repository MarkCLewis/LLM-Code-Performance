package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort" // Still needed for sort.Interface, now using sort.Select
	"sync"
	"time"
)

// --- Constants --- (Same as before)
const (
	G                 = 1.0
	dt                = 0.001
	epsilon           = 0.01
	numSteps          = 10
	numOrbitingBodies = 100000
	// numOrbitingBodies = 5000 // Smaller number for quicker testing
	theta = 0.3
)

var thetaSquared = theta * theta // Pre-calculate theta^2

// --- Data Structures ---

// Vector3D and Body structs remain the same
type Vector3D struct{ X, Y, Z float64 }
type Body struct {
	Mass         float64
	Position     Vector3D
	Velocity     Vector3D
	Acceleration Vector3D // Keep storing acceleration here
}
type BoundingBox struct{ Min, Max Vector3D }

// KDNode struct remains the same
type KDNode struct {
	Axis          int
	SplitValue    float64
	Bounds        BoundingBox
	CenterOfMass  Vector3D
	TotalMass     float64
	ParticleIndex int
	LeftChild     *KDNode // Pointer into nodePool
	RightChild    *KDNode // Pointer into nodePool
	NodeSize      float64
}

// --- Node Pool ---
var (
	nodePool      []KDNode   // Pre-allocated pool of nodes
	nextNodeIndex int        // Index for the next available node in the pool
	poolMutex     sync.Mutex // Mutex to protect nextNodeIndex if build were parallel
)

// initNodePool initializes or resizes the node pool
func initNodePool(numBodies int) {
	// Estimate required nodes: Max depth logN, roughly 2N nodes total is safe upper bound
	requiredNodes := numBodies * 2
	if cap(nodePool) < requiredNodes {
		nodePool = make([]KDNode, requiredNodes)
		fmt.Printf("Resized node pool to %d nodes\n", requiredNodes)
	}
	nextNodeIndex = 0 // Reset index for the new simulation step
}

// getNode retrieves a node from the pool
func getNode() *KDNode {
	// In a parallel build, this index increment would need protection (mutex or atomic)
	// poolMutex.Lock() // Not needed for sequential build
	if nextNodeIndex >= len(nodePool) {
		// This should not happen if initNodePool is sized correctly
		panic(fmt.Sprintf("Node pool exhausted! Index %d, Size %d. Increase pool size.", nextNodeIndex, len(nodePool)))
	}
	node := &nodePool[nextNodeIndex]
	nextNodeIndex++
	// poolMutex.Unlock() // Not needed for sequential build

	// Reset fields (important!)
	node.LeftChild = nil
	node.RightChild = nil
	node.ParticleIndex = -1
	// Other fields (Axis, SplitValue, Bounds, COM, TotalMass, NodeSize) will be overwritten by buildKDTree

	return node
}

// --- Vector Operations --- (Same as before)
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

// Implementations omitted for brevity - they are the same as before

// --- Bounding Box Operations --- (Same as before)
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

// Implementations omitted for brevity - they are the same as before

// --- k-D Tree Construction (Optimized) ---

// bodySorter remains the same (needed for sort.Interface for sort.Select)
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

// Implementations omitted for brevity - they are the same as before

// buildKDTree recursively builds the tree using node pool and sort.Select
func buildKDTree(particleIndices []int, bodies []Body, depth int, currentBounds BoundingBox) *KDNode {
	n := len(particleIndices)
	if n == 0 {
		return nil
	}

	// Get node from pool instead of allocating
	node := getNode()
	node.Bounds = currentBounds
	node.ParticleIndex = -1

	// Base Case: Leaf node
	if n == 1 {
		idx := particleIndices[0]
		node.ParticleIndex = idx
		node.TotalMass = bodies[idx].Mass
		node.CenterOfMass = bodies[idx].Position
		node.Bounds = BoundingBox{Min: bodies[idx].Position, Max: bodies[idx].Position}
		node.NodeSize = 0
		return node
	}

	// Recursive Step: Internal node
	node.Axis = depth % 3

	// --- Use sort.Select to find median and partition (O(N) average) ---
	medianIndex := n / 2
	sorter := &bodySorter{indices: particleIndices, bodies: bodies, axis: node.Axis}
	// sort.Select partitions such that element at medianIndex is the true median,
	// and elements <= median are before it, elements >= median are after it.
	// Note: This is partial sorting, not full sorting.
	quickSelect(sorter, 0, n-1, medianIndex) // Use helper for clarity

	medianParticleIdx := particleIndices[medianIndex]
	node.SplitValue = bodies[medianParticleIdx].Position.Component(node.Axis)
	// --- End sort.Select ---

	// --- Calculate Bounds for Children (Same as before) ---
	leftBounds := currentBounds
	leftBounds.Max = leftBounds.Max.withComponent(node.Axis, node.SplitValue)
	rightBounds := currentBounds
	rightBounds.Min = rightBounds.Min.withComponent(node.Axis, node.SplitValue)

	// Build children recursively (still sequential build)
	// Partition indices based on the median element found by Select
	node.LeftChild = buildKDTree(particleIndices[:medianIndex], bodies, depth+1, leftBounds)
	node.RightChild = buildKDTree(particleIndices[medianIndex:], bodies, depth+1, rightBounds) // NOTE: Includes median in right child

	// --- Calculate COM and Total Mass (Same as before) ---
	node.TotalMass = 0
	comNumerator := Vector3D{0, 0, 0}
	nodeBounds := BoundingBox{} // Initialize bounds for expansion
	hasLeft, hasRight := false, false

	if node.LeftChild != nil {
		node.TotalMass += node.LeftChild.TotalMass
		comNumerator = comNumerator.Add(node.LeftChild.CenterOfMass.Scale(node.LeftChild.TotalMass))
		nodeBounds = node.LeftChild.Bounds // Initialize with left bounds
		hasLeft = true
	}
	if node.RightChild != nil {
		node.TotalMass += node.RightChild.TotalMass
		comNumerator = comNumerator.Add(node.RightChild.CenterOfMass.Scale(node.RightChild.TotalMass))
		if hasLeft {
			nodeBounds.Expand(node.RightChild.Bounds) // Expand with right bounds
		} else {
			nodeBounds = node.RightChild.Bounds // Initialize with right bounds
		}
		hasRight = true
	}

	// If node only had one child type or was created from empty partitions (shouldn't happen with len > 1 check)
	if !hasLeft && !hasRight {
		// Fallback if children were nil, use original bounds? Or from single particle?
		// This path ideally shouldn't be hit if n > 1 initially.
		// If hit, COM/Mass calculation below handles TotalMass=0
		nodeBounds = currentBounds // Or tighter bounds?
	}
	node.Bounds = nodeBounds // Assign the calculated combined bounds

	if node.TotalMass > 1e-12 {
		node.CenterOfMass = comNumerator.Scale(1.0 / node.TotalMass)
	} else {
		node.CenterOfMass = node.Bounds.Min.Add(node.Bounds.Max).Scale(0.5)
	}

	// Calculate node size (width along the splitting axis) using the potentially expanded bounds
	node.NodeSize = node.Bounds.Max.Component(node.Axis) - node.Bounds.Min.Component(node.Axis)

	return node
}

// --- QuickSelect Helper (using sort.Interface) ---
// This implements the core logic of finding the k-th element
// It modifies the underlying slice (particleIndices) in place.

func quickSelect(sorter sort.Interface, left, right, k int) {
	// Standard recursive QuickSelect implementation
	if left == right {
		return
	}

	pivotIndex := partition(sorter, left, right)

	if k == pivotIndex {
		return
	} else if k < pivotIndex {
		quickSelect(sorter, left, pivotIndex-1, k)
	} else {
		quickSelect(sorter, pivotIndex+1, right, k)
	}
}

// partition helper for quickSelect (Lomuto partition scheme example)
func partition(sorter sort.Interface, left, right int) int {
	// Choose a pivot (e.g., the rightmost element)
	pivotIndex := right
	storeIndex := left
	// Move pivot to end (optional, but common in some Lomuto variants)
	// sorter.Swap(pivotIndex, right)

	for i := left; i < right; i++ {
		if sorter.Less(i, pivotIndex) { // Compare elements to pivot value
			sorter.Swap(storeIndex, i)
			storeIndex++
		}
	}
	// Move pivot to its final place
	sorter.Swap(storeIndex, right)
	return storeIndex
}

// --- Force Calculation using k-D Tree (Optimized) ---

// calculateForceRecursive traverses the tree (optimized MAC)
func calculateForceRecursive(node *KDNode, targetParticleIdx int, bodies []Body, thetaSq, gravitationalConstant, softeningFactor float64) Vector3D {
	if node == nil {
		return Vector3D{0, 0, 0}
	}

	// Leaf Node Case
	if node.ParticleIndex != -1 {
		if node.ParticleIndex == targetParticleIdx {
			return Vector3D{0, 0, 0}
		} // Skip self

		bodyI := bodies[targetParticleIdx]
		bodyJ := bodies[node.ParticleIndex]
		deltaPos := bodyJ.Position.Subtract(bodyI.Position)
		distSq := deltaPos.MagnitudeSquared() + softeningFactor*softeningFactor // Add softening squared
		if distSq < 1e-18 {
			return Vector3D{0, 0, 0}
		} // Avoid division by zero if softening is tiny/zero

		// Use optimized Pow: math.Pow(distSq, -1.5) == 1 / (distSq * sqrt(distSq))
		invDistCubed := math.Pow(distSq, -1.5) // Optimized inverse cube
		forceScale := gravitationalConstant * bodyJ.Mass * invDistCubed
		return deltaPos.Scale(forceScale)
	}

	// Internal Node Case - Apply Optimized Barnes-Hut Criterion
	deltaPosCOM := node.CenterOfMass.Subtract(bodies[targetParticleIdx].Position)
	distSqCOM := deltaPosCOM.MagnitudeSquared()

	// Optimized Criterion: s*s < theta*theta * d*d  (avoids sqrt)
	// Check distSqCOM > ~0 to avoid issues if particle is exactly at COM
	if distSqCOM > 1e-18 && (node.NodeSize*node.NodeSize) < (thetaSq*distSqCOM) {
		// Node is far enough away, use its Center of Mass approximation
		distSqCOM += softeningFactor * softeningFactor // Add softening squared
		if distSqCOM < 1e-18 {
			return Vector3D{0, 0, 0}
		} // Check again after adding softening

		invDistCubed := math.Pow(distSqCOM, -1.5) // Optimized inverse cube
		forceScale := gravitationalConstant * node.TotalMass * invDistCubed
		return deltaPosCOM.Scale(forceScale)
	} else {
		// Node is too close or large, recurse into children
		forceLeft := calculateForceRecursive(node.LeftChild, targetParticleIdx, bodies, thetaSq, gravitationalConstant, softeningFactor)
		forceRight := calculateForceRecursive(node.RightChild, targetParticleIdx, bodies, thetaSq, gravitationalConstant, softeningFactor)
		return forceLeft.Add(forceRight)
	}
}

// --- Simulation Step (Using Optimized Functions) ---

// calculateForceKDTreeChunk remains structurally similar, calls optimized recursive func
func calculateForceKDTreeChunk(treeRoot *KDNode, bodies []Body, start, end int, wg *sync.WaitGroup, thetaSq, gravitationalConstant, softeningFactor float64) {
	defer wg.Done()
	for i := start; i < end; i++ {
		// Pass thetaSquared instead of theta
		forceSumFactor := calculateForceRecursive(treeRoot, i, bodies, thetaSq, gravitationalConstant, softeningFactor)
		bodies[i].Acceleration = forceSumFactor
	}
}

// updateVelocityChunk and updatePositionChunk remain the same
func updateVelocityChunk(bodies []Body, start, end int, wg *sync.WaitGroup, timeStep float64) { /* ... */
}
func updatePositionChunk(bodies []Body, start, end int, wg *sync.WaitGroup, timeStep float64) { /* ... */
}

// Implementations omitted for brevity - they are the same as before

// simulateStepParallelKDTree orchestrates the optimized step
func simulateStepParallelKDTree(bodies []Body, numWorkers int, timeStep, gravitationalConstant, softeningFactor, thetaSq float64) /* *KDNode - removed return */ {
	n := len(bodies)
	var wg sync.WaitGroup

	// --- 0. Reset node pool index and Build k-D Tree ---
	// In a real application, you might init the pool once outside the loop
	// and just reset the index here. initNodePool includes the reset.
	initNodePool(n) // Ensure pool is ready and reset index
	// buildStart := time.Now()
	initialBounds := calculateOverallBounds(bodies)
	particleIndices := make([]int, n) // Still need indices slice
	for i := range particleIndices {
		particleIndices[i] = i
	}

	treeRoot := buildKDTree(particleIndices, bodies, 0, initialBounds)
	// buildDuration := time.Since(buildStart) // Optional profiling

	// --- 1. Parallel Force/Acceleration Calculation using k-D Tree ---
	// forceStart := time.Now() // Optional profiling
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
		// Pass thetaSquared to the chunk worker
		go calculateForceKDTreeChunk(treeRoot, bodies, start, end, &wg, thetaSq, gravitationalConstant, softeningFactor)
	}
	wg.Wait()
	// forceDuration := time.Since(forceStart) // Optional profiling

	// --- 2. Parallel Velocity Update (Kick) ---
	// velStart := time.Now() // Optional profiling
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
	// velDuration := time.Since(velStart) // Optional profiling

	// --- 3. Parallel Position Update (Step) ---
	// posStart := time.Now() // Optional profiling
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
	// posDuration := time.Since(posStart) // Optional profiling

	// Optional Step Timing Print
	// if n < 10000 { // Avoid spamming console for large N
	//     fmt.Printf("    Step Timings: Build=%.4fs Force=%.4fs Vel=%.4fs Pos=%.4fs\n", buildDuration.Seconds(), forceDuration.Seconds(), velDuration.Seconds(), posDuration.Seconds())
	// }
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

// Implementation omitted for brevity - it is the same as before

// --- Energy Calculation --- (Same as before)
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

// Implementation omitted for brevity - it is the same as before

// --- Main Function --- (Updated to pass thetaSquared)
func main() {
	numWorkers := runtime.NumCPU()
	runtime.GOMAXPROCS(numWorkers) // Good practice

	// Initialize node pool once before starting
	initNodePool(numOrbitingBodies + 1)

	fmt.Printf("Starting N-body simulation (Optimized Parallel k-D Tree, theta=%.2f)...\n", theta)
	fmt.Printf("Using %d CPU cores/workers.\n", numWorkers)
	fmt.Printf("Number of bodies: %d (+1 central body)\n", numOrbitingBodies)
	// ... other print statements ...

	// --- Initialization ---
	startTime := time.Now()
	// ... body initialization ...
	bodies := initializeCircularOrbitSystem(numOrbitingBodies, 10000.0, 0.01, 5.0, 50.0)
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
		// Pass thetaSquared to the simulation step function
		simulateStepParallelKDTree(bodies, numWorkers, dt, G, epsilon, thetaSquared)
		stepDuration := time.Since(stepStart)

		// Progress reporting (same as before)
		if (step+1)%10 == 0 || step == 0 || step == numSteps-1 {
			fmt.Printf("Step %d/%d completed (Step duration: %v, Elapsed: %v)\n",
				step+1, numSteps, stepDuration, time.Since(simulationStartTime))
		}
	}

	simulationTime := time.Since(simulationStartTime)
	fmt.Printf("Simulation completed in %v\n", simulationTime)

	// --- Final Energy Calculation & Reporting --- (Same as before)
	finalEnergy := calculateTotalEnergy(bodies)
	// ... energy reporting ...
	fmt.Printf("Final Total Energy (Exact):   %.5e\n", finalEnergy)
	// ... rest of energy reporting ...

	totalTime := time.Since(startTime)
	fmt.Printf("Total execution time: %v\n", totalTime)
}

// Ensure all helper functions (vector ops, bounding box ops) are included correctly
// ... (Implementations for Add, Subtract, Scale, etc., here if not shown above) ...
// Example:
// func (v Vector3D) Add(other Vector3D) Vector3D {
// 	return Vector3D{v.X + other.X, v.Y + other.Y, v.Z + other.Z}
// }

// ... Add other vector/bbox methods back in ...
