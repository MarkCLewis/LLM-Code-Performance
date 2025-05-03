package main

import (
	"fmt"
	"math"
	"math/rand"
	"sync"
	"time"
)

const (
	G          = 6.67430e-11 // Gravitational constant
	DT         = 1.0         // Time step
	NUM_BODIES = 1000000     // Number of bodies
	STEPS      = 1000        // Simulation steps
	THETA      = 0.3         // Barnes-Hut approximation threshold
)

type Body struct {
	Position [3]float64
	Velocity [3]float64
	Mass     float64
}

type KDTree struct {
	centerOfMass [3]float64
	totalMass    float64
	boundMin     [3]float64
	boundMax     [3]float64
	children     []*KDTree
	body         *Body
}

func buildKDTree(bodies []Body, boundMin, boundMax [3]float64) *KDTree {
	if len(bodies) == 0 {
		return nil
	}

	node := &KDTree{}
	if len(bodies) == 1 {
		node.body = &bodies[0]
		node.totalMass = bodies[0].Mass
		node.centerOfMass = bodies[0].Position
		return node
	}

	node.children = make([]*KDTree, 8)
	for _, body := range bodies {
		node.totalMass += body.Mass
		for i := 0; i < 3; i++ {
			node.centerOfMass[i] += body.Position[i] * body.Mass
		}
	}
	for i := 0; i < 3; i++ {
		node.centerOfMass[i] /= node.totalMass
	}

	// Calculate the bounding box for the KD-tree node
	node.boundMin = boundMin
	node.boundMax = boundMax

	// Distribute bodies into 8 children nodes based on their position
	for _, body := range bodies {
		idx := 0
		for i := 0; i < 3; i++ {
			if body.Position[i] > (boundMin[i]+boundMax[i])/2 {
				idx |= 1 << i
			}
		}
		node.children[idx] = buildKDTree([]Body{body}, boundMin, boundMax)
	}
	return node
}

func computeForce(body *Body, tree *KDTree, theta float64) (ax, ay, az float64) {
	if tree == nil || tree.body == body {
		return 0, 0, 0
	}

	// Calculate distance to the center of mass of the node
	dx := tree.centerOfMass[0] - body.Position[0]
	dy := tree.centerOfMass[1] - body.Position[1]
	dz := tree.centerOfMass[2] - body.Position[2]
	dist := math.Sqrt(dx*dx + dy*dy + dz*dz)
	if dist == 0 {
		return 0, 0, 0
	}

	// Calculate the size of the node in comparison to the distance
	s := tree.boundMax[0] - tree.boundMin[0]
	if s/dist < theta || tree.body != nil {
		factor := G * tree.totalMass / (dist * dist * dist)
		ax += factor * dx
		ay += factor * dy
		az += factor * dz
	} else {
		for _, child := range tree.children {
			cx, cy, cz := computeForce(body, child, theta)
			ax += cx
			ay += cy
			az += cz
		}
	}
	return ax, ay, az
}

func updateVelocities(bodies []Body, tree *KDTree) {
	var wg sync.WaitGroup
	n := len(bodies)
	wg.Add(n)
	for i := 0; i < n; i++ {
		go func(i int) {
			defer wg.Done()
			ax, ay, az := computeForce(&bodies[i], tree, THETA)
			bodies[i].Velocity[0] += ax * DT
			bodies[i].Velocity[1] += ay * DT
			bodies[i].Velocity[2] += az * DT
		}(i)
	}
	wg.Wait()
}

func updatePositions(bodies []Body) {
	for i := 0; i < len(bodies); i++ {
		bodies[i].Position[0] += bodies[i].Velocity[0] * DT
		bodies[i].Position[1] += bodies[i].Velocity[1] * DT
		bodies[i].Position[2] += bodies[i].Velocity[2] * DT
	}
}

func initializeOrbitingBodies(numBodies int, centralMass float64) []Body {
	bodies := make([]Body, numBodies)
	bodies[0] = Body{Mass: centralMass}
	radius := 1.0e9
	for i := 1; i < numBodies; i++ {
		angle := 2.0 * math.Pi * rand.Float64()
		posX := radius * math.Cos(angle)
		posY := radius * math.Sin(angle)
		speed := math.Sqrt(G * centralMass / radius)
		velX := -speed * math.Sin(angle)
		velY := speed * math.Cos(angle)
		bodies[i] = Body{
			Position: [3]float64{posX, posY, 0},
			Velocity: [3]float64{velX, velY, 0},
			Mass:     1.0,
		}
	}
	return bodies
}

func main() {
	rand.Seed(time.Now().UnixNano())
	bodies := initializeOrbitingBodies(NUM_BODIES, 1.989e30)

	// Set initial bounding box for KD-tree (assuming a large enough region)
	boundMin := [3]float64{-1.0e10, -1.0e10, -1.0e10}
	boundMax := [3]float64{1.0e10, 1.0e10, 1.0e10}

	// Build the KD-tree once
	initialTree := buildKDTree(bodies, boundMin, boundMax)

	// Main simulation loop
	for step := 0; step < STEPS; step++ {
		// Update the velocities based on the forces calculated from the tree
		updateVelocities(bodies, initialTree)

		// Update the positions of the bodies based on the updated velocities
		updatePositions(bodies)

		// Optionally, you could rebuild the tree at intervals, or you can update it incrementally
	}

	fmt.Println("Simulation complete.")
}