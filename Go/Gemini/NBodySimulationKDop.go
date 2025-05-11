// Okay, let's implement the kickStepKDTreeOptimized function with chunking for the parallel force calculation. This is a relatively straightforward optimization to try and can potentially improve performance by reducing the overhead of launching a goroutine for each individual body.
/* Changes Made:

kickStepKDTree Modified:
The function now takes an additional chunkSize parameter.
The loop that launches goroutines now iterates with a step of chunkSize.
Each goroutine processes a slice of the bodies and forces arrays.
firstOrderKickStepKDTree Modified:
It now calls kickStepKDTree with the chunkSize.
main Function Modified:
A chunkSize variable is introduced (set to 1000 as an example). You should experiment with different values for this based on your system's core count and performance.
firstOrderKickStepKDTree is called with the chunkSize.
How Chunking Can Help:

Reduced Goroutine Overhead: Launching a single goroutine to process a batch of bodies is generally more efficient than launching a goroutine for each individual body, especially when the number of bodies is very large. The overhead of goroutine creation and scheduling can become significant.
Improved Cache Locality: Processing a contiguous chunk of bodies might lead to better cache utilization as the data for these bodies is likely to be closer in memory.
**Next Steps for Further Optimization (Beyond this immediate change):
*/

package main

import (
	"fmt"
	"math"
	"math/rand"
	"sync"
	"time"

	"github.com/rudty/go-parallel"
)

// Body represents a celestial body with mass, position, and velocity.
type Body struct {
	Mass     float64
	Position [3]float64
	Velocity [3]float64
}

// System represents a collection of celestial bodies.
type System struct {
	Bodies []Body
}

// KDNode represents a node in the kD-tree.
type KDNode struct {
	bodyIndex    int // Index of the body in this leaf node (-1 if internal)
	centerOfMass [3]float64
	totalMass    float64
	min          [3]float64 // Minimum bounds of the node
	max          [3]float64 // Maximum bounds of the node
	left         *KDNode
	right        *KDNode
}

// GravitationalConstant is the universal gravitational constant.
const GravitationalConstant = 6.67430e-11
const theta = 0.3 // Barnes-Hut opening angle

// distanceSquared calculates the squared Euclidean distance between two points.
func distanceSquaredPoints(p1, p2 [3]float64) float64 {
	dx := p1[0] - p2[0]
	dy := p1[1] - p2[1]
	dz := p1[2] - p2[2]
	return dx*dx + dy*dy + dz*dz
}

// distancePoints calculates the Euclidean distance between two points.
func distancePoints(p1, p2 [3]float64) float64 {
	return math.Sqrt(distanceSquaredPoints(p1, p2))
}

// calculateForceElement calculates the force on a body due to a mass element (node or body).
func calculateForceElement(body *Body, massElementCM [3]float64, massElementMass float64, force *[3]float64) {
	rSq := distanceSquaredPoints(body.Position, massElementCM)
	if rSq > 1e-9 {
		r := math.Sqrt(rSq)
		magnitude := (GravitationalConstant * body.Mass * massElementMass) / rSq
		force[0] += magnitude * (massElementCM[0] - body.Position[0]) / r
		force[1] += magnitude * (massElementCM[1] - body.Position[1]) / r
		force[2] += magnitude * (massElementCM[2] - body.Position[2]) / r
	}
}

// calculateTotalEnergy calculates the total energy of the system.
func calculateTotalEnergy(s *System) float64 {
	kineticEnergy := 0.0
	potentialEnergy := 0.0
	numBodies := len(s.Bodies)

	var wg sync.WaitGroup
	var kineticMutex sync.Mutex
	var potentialMutex sync.Mutex

	wg.Add(numBodies)
	for i := 0; i < numBodies; i++ {
		go func(i int) {
			defer wg.Done()
			vSq := s.Bodies[i].Velocity[0]*s.Bodies[i].Velocity[0] +
				s.Bodies[i].Velocity[1]*s.Bodies[i].Velocity[1] +
				s.Bodies[i].Velocity[2]*s.Bodies[i].Velocity[2]
			kineticMutex.Lock()
			kineticEnergy += 0.5 * s.Bodies[i].Mass * vSq
			kineticMutex.Unlock()

			for j := i + 1; j < numBodies; j++ {
				r := math.Sqrt(distanceSquared(&s.Bodies[i], &s.Bodies[j]))
				potentialMutex.Lock()
				potentialEnergy -= (GravitationalConstant * s.Bodies[i].Mass * s.Bodies[j].Mass) / r
				potentialMutex.Unlock()
			}
		}(i)
	}
	wg.Wait()

	return kineticEnergy + potentialEnergy
}

// initializeCircularOrbits initializes a system with a central body and orbiting smaller bodies.
func initializeCircularOrbits(numOrbiting int, centralMass, orbitRadius, orbitingMass float64) *System {
	system := &System{
		Bodies: make([]Body, numOrbiting+1),
	}

	// Initialize the central body
	system.Bodies[0] = Body{
		Mass:     centralMass,
		Position: [3]float64{0, 0, 0},
		Velocity: [3]float64{0, 0, 0},
	}

	// Initialize the orbiting bodies
	for i := 0; i < numOrbiting; i++ {
		angle := 2 * math.Pi * float64(i) / float64(numOrbiting)
		system.Bodies[i+1] = Body{
			Mass: orbitingMass,
			Position: [3]float64{
				orbitRadius * math.Cos(angle),
				orbitRadius * math.Sin(angle),
				0,
			},
		}

		// Calculate the orbital velocity for a circular orbit
		orbitalSpeed := math.Sqrt(GravitationalConstant * centralMass / orbitRadius)
		system.Bodies[i+1].Velocity = [3]float64{
			-orbitalSpeed * math.Sin(angle),
			orbitalSpeed * math.Cos(angle),
			0,
		}
	}

	return system
}

// calculateBounds calculates the bounding box of a set of bodies.
func calculateBounds(bodies []Body) ([3]float64, [3]float64) {
	if len(bodies) == 0 {
		return [3]float64{}, [3]float64{}
	}
	min := bodies[0].Position
	max := bodies[0].Position
	for _, b := range bodies {
		for i := 0; i < 3; i++ {
			if b.Position[i] < min[i] {
				min[i] = b.Position[i]
			}
			if b.Position[i] > max[i] {
				max[i] = b.Position[i]
			}
		}
	}
	return min, max
}

// newKDNode creates a new kD-tree node.
func newKDNode(bodyIndex int) *KDNode {
	return &KDNode{
		bodyIndex: bodyIndex,
		totalMass: 0,
		left:      nil,
		right:     nil,
	}
}

// buildKDTreeRecursive recursively builds the kD-tree.
func buildKDTreeRecursive(bodies []Body, indices []int, minBound, maxBound [3]float64, depth int) *KDNode {
	numIndices := len(indices)
	if numIndices == 0 {
		return nil
	}

	if numIndices == 1 {
		node := newKDNode(indices[0])
		node.totalMass = bodies[indices[0]].Mass
		node.centerOfMass = bodies[indices[0]].Position
		node.min = minBound
		node.max = maxBound
		return node
	}

	node := newKDNode(-1) // Internal node
	node.min = minBound
	node.max = maxBound
	node.totalMass = 0
	centerOfMass := [3]float64{0, 0, 0}
	for _, index := range indices {
		node.totalMass += bodies[index].Mass
		for i := 0; i < 3; i++ {
			centerOfMass[i] += bodies[index].Mass * bodies[index].Position[i]
		}
	}
	if node.totalMass > 0 {
		for i := 0; i < 3; i++ {
			node.centerOfMass[i] = centerOfMass[i] / node.totalMass
		}
	}

	splitDim := depth % 3
	median := (minBound[splitDim] + maxBound[splitDim]) / 2

	leftIndices := make([]int, 0)
	rightIndices := make([]int, 0)

	for _, index := range indices {
		if bodies[index].Position[splitDim] <= median {
			leftIndices = append(leftIndices, index)
		} else {
			rightIndices = append(rightIndices, index)
		}
	}

	leftMinBound := minBound
	leftMaxBound := maxBound
	leftMaxBound[splitDim] = median

	rightMinBound := minBound
	rightMaxBound := maxBound
	rightMinBound[splitDim] = median

	node.left = buildKDTreeRecursive(bodies, leftIndices, leftMinBound, leftMaxBound, depth+1)
	node.right = buildKDTreeRecursive(bodies, rightIndices, rightMinBound, rightMaxBound, depth+1)

	return node
}

// buildKDTree builds the kD-tree for the given system.
func buildKDTree(s *System) *KDNode {
	indices := make([]int, len(s.Bodies))
	for i := range s.Bodies {
		indices[i] = i
	}
	minBound, maxBound := calculateBounds(s.Bodies)
	return buildKDTreeRecursive(s.Bodies, indices, minBound, maxBound, 0)
}

// calculateForceKDTreeRecursive recursively calculates the force on a body using the kD-tree.
func calculateForceKDTreeRecursive(body *Body, node *KDNode, s *System, force *[3]float64) {
	if node == nil {
		return
	}

	if node.bodyIndex != -1 { // Leaf node
		if node.bodyIndex != (body - &s.Bodies[0]) {
			calculateForceElement(body, node.centerOfMass, node.totalMass, force)
		}
		return
	}

	// Internal node
	sNode := [3]float64{}
	for i := 0; i < 3; i++ {
		sNode[i] = node.max[i] - node.min[i]
	}
	d := distancePoints(body.Position, node.centerOfMass)

	if d == 0 || (sNode[0]/d < theta && sNode[1]/d < theta && sNode[2]/d < theta) {
		calculateForceElement(body, node.centerOfMass, node.totalMass, force)
	} else {
		calculateForceKDTreeRecursive(body, node.left, s, force)
		calculateForceKDTreeRecursive(body, node.right, s, force)
	}
}

// calculateForceKDTree calculates the force on a body using the kD-tree.
func calculateForceKDTree(body *Body, root *KDNode, s *System, force *[3]float64) {
	force[0] = 0
	force[1] = 0
	force[2] = 0
	calculateForceKDTreeRecursive(body, root, s, force)
}

// kickStepKDTree performs the "kick" step using the kD-tree for force calculation with chunking.
func kickStepKDTree(s *System, dt float64, chunkSize int) {
	numBodies := len(s.Bodies)
	forces := make([][3]float64, numBodies)
	root := buildKDTree(s)

	// pool := parallel.NewGroup()
	// pool.SetMaxWorkers(runtime.NumCPU())

	for i := 0; i < numBodies; i += chunkSize {
		start := i
		end := i + chunkSize
		if end > numBodies {
			end = numBodies
		}

		// pool.Go(func() error {
		parallel.For(0, numBodies, func(i int) {
			for j := start; j < end; j++ {
				calculateForceKDTree(&s.Bodies[j], root, s, &forces[j])
			}
			// return nil
		})
	}
	// pool.Wait()

	// Update velocities
	for i := 0; i < numBodies; i++ {
		s.Bodies[i].Velocity[0] += (forces[i][0] / s.Bodies[i].Mass) * dt * 0.5
		s.Bodies[i].Velocity[1] += (forces[i][1] / s.Bodies[i].Mass) * dt * 0.5
		s.Bodies[i].Velocity[2] += (forces[i][2] / s.Bodies[i].Mass) * dt * 0.5
	}
}

// driftStep performs the "drift" part of the first-order kick-step method (parallelized).
func driftStep(s *System, dt float64) {
	numBodies := len(s.Bodies)
	// pool := parallel.NewGroup()
	// pool.SetMaxWorkers(runtime.NumCPU())

	for i := 0; i < numBodies; i++ {
		// i := i
		// pool.Go(func() error {
		parallel.For(0, numBodies, func(i int) {
			s.Bodies[i].Position[0] += s.Bodies[i].Velocity[0] * dt
			s.Bodies[i].Position[1] += s.Bodies[i].Velocity[1] * dt
			s.Bodies[i].Position[2] += s.Bodies[i].Velocity[2] * dt
			// return nil
		})
	}
	// pool.Wait()
}

// firstOrderKickStepKDTree performs one full first-order kick-step using the kD-tree with chunking.
func firstOrderKickStepKDTree(s *System, dt float64, chunkSize int) {
	kickStepKDTree(s, dt, chunkSize)
	driftStep(s, dt)
	kickStepKDTree(s, dt, chunkSize)
}

func main() {
	numOrbitingBodies := 100000
	centralMass := 1.989e30  // Mass of the Sun (kg)
	orbitRadius := 1.496e11  // 1 AU (m)
	orbitingMass := 5.972e24 // Mass of the Earth (kg)
	numSteps := 10
	timeStep := 3600.0 * 24.0 * 7.0 // 1 week in seconds
	chunkSize := 1000               // Experiment with different chunk sizes

	rand.Seed(time.Now().UnixNano())

	// Initialize the system
	system := initializeCircularOrbits(numOrbitingBodies, centralMass, orbitRadius, orbitingMass)
	fmt.Printf("Initial number of bodies: %d\n", len(system.Bodies))

	// Calculate initial energy
	initialEnergy := calculateTotalEnergy(system)
	fmt.Printf("Initial total energy: %e J\n", initialEnergy)

	// Run the simulation
	fmt.Printf("Running simulation for %d steps...\n", numSteps)
	startTime := time.Now()
	for step := 0; step < numSteps; step++ {
		firstOrderKickStepKDTree(system, timeStep, chunkSize)
		if (step+1)%100 == 0 {
			fmt.Printf("Step %d completed.\n", step+1)
		}
	}
	elapsedTime := time.Since(startTime)
	fmt.Printf("Simulation finished in %s.\n", elapsedTime)

	// Calculate final energy
	finalEnergy := calculateTotalEnergy(system)
	fmt.Printf("Final total energy: %e J\n", finalEnergy)

	// Calculate the energy difference
	energyDifference := math.Abs(finalEnergy - initialEnergy)
	relativeEnergyDifference := energyDifference / math.Abs(initialEnergy)
	fmt.Printf("Absolute energy difference: %e J\n", energyDifference)
	fmt.Printf("Relative energy difference: %e\n", relativeEnergyDifference)
}
