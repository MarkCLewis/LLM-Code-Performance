package main

import (
	"fmt"
	"math"
	"sort"
	"sync"
)

const (
	G         = 6.67430e-11            // Gravitational constant
	DT        = 1e-3 * 3600 * 24 * 365 // Time step
	NumBodies = 100000                 // Number of bodies
	THETA     = 0.3                    // Theta value for approximation
)

type Body struct {
	x, y, z    float64
	vx, vy, vz float64
	mass       float64
}

type KDNode struct {
	body        *Body
	left, right *KDNode
	min, max    [3]float64
}

func initializeBodies(bodies []Body) {
	bodies[0] = Body{0, 0, 0, 0, 0, 0, 1e30} // Central body mass

	for i := 1; i < len(bodies); i++ {
		angle := 2 * math.Pi * float64(i) / float64(len(bodies)-1)
		bodies[i] = Body{
			x:    math.Cos(angle) * 1e11,
			y:    math.Sin(angle) * 1e11,
			z:    0,
			vx:   -math.Sin(angle) * math.Sqrt(G*bodies[0].mass/1e11),
			vy:   math.Cos(angle) * math.Sqrt(G*bodies[0].mass/1e11),
			vz:   0,
			mass: 1e24 / float64(len(bodies)), // Small body mass
		}
	}
}

func calculateEnergy(bodies []Body) float64 {
	energy := 0.0
	var wg sync.WaitGroup
	var mu sync.Mutex

	for i := range bodies {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			kinetic := 0.5 * bodies[i].mass * (bodies[i].vx*bodies[i].vx + bodies[i].vy*bodies[i].vy + bodies[i].vz*bodies[i].vz)
			potential := 0.0
			for j := range bodies {
				if i != j {
					dx := bodies[i].x - bodies[j].x
					dy := bodies[i].y - bodies[j].y
					dz := bodies[i].z - bodies[j].z
					distance := math.Sqrt(dx*dx + dy*dy + dz*dz)
					potential -= G * bodies[i].mass * bodies[j].mass / distance
				}
			}
			mu.Lock()
			energy += kinetic + 0.5*potential
			mu.Unlock()
		}(i)
	}
	wg.Wait()
	return energy
}

func buildKDTree(bodies []Body, depth int) *KDNode {
	if len(bodies) == 0 {
		return nil
	}

	axis := depth % 3
	sort.Slice(bodies, func(i, j int) bool {
		switch axis {
		case 0:
			return bodies[i].x < bodies[j].x
		case 1:
			return bodies[i].y < bodies[j].y
		default:
			return bodies[i].z < bodies[j].z
		}
	})

	median := len(bodies) / 2
	node := &KDNode{
		body:  &bodies[median],
		left:  buildKDTree(bodies[:median], depth+1),
		right: buildKDTree(bodies[median+1:], depth+1),
	}

	for i := 0; i < 3; i++ {
		node.min[i] = node.body.x
		node.max[i] = node.body.x
		if node.left != nil {
			node.min[i] = math.Min(node.min[i], node.left.min[i])
			node.max[i] = math.Max(node.max[i], node.left.max[i])
		}
		if node.right != nil {
			node.min[i] = math.Min(node.min[i], node.right.min[i])
			node.max[i] = math.Max(node.max[i], node.right.max[i])
		}
	}

	return node
}

func calculateForce(node *KDNode, body *Body, ax, ay, az *float64) {
	if node == nil || node.body == nil {
		return
	}

	dx := node.body.x - body.x
	dy := node.body.y - body.y
	dz := node.body.z - body.z
	distance := math.Sqrt(dx*dx + dy*dy + dz*dz)

	size := math.Max(node.max[0]-node.min[0], math.Max(node.max[1]-node.min[1], node.max[2]-node.min[2]))

	if size/distance < THETA || (node.left == nil && node.right == nil) {
		force := G * node.body.mass / (distance * distance * distance)
		*ax += force * dx
		*ay += force * dy
		*az += force * dz
	} else {
		calculateForce(node.left, body, ax, ay, az)
		calculateForce(node.right, body, ax, ay, az)
	}
}

func kickStep(bodies []Body, root *KDNode) {
	var wg sync.WaitGroup

	for i := range bodies {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			ax, ay, az := 0.0, 0.0, 0.0
			calculateForce(root, &bodies[i], &ax, &ay, &az)
			bodies[i].vx += ax * DT
			bodies[i].vy += ay * DT
			bodies[i].vz += az * DT
		}(i)
	}
	wg.Wait()

	for i := range bodies {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			bodies[i].x += bodies[i].vx * DT
			bodies[i].y += bodies[i].vy * DT
			bodies[i].z += bodies[i].vz * DT
		}(i)
	}
	wg.Wait()
}

func main() {
	bodies := make([]Body, NumBodies)
	initializeBodies(bodies)

	initialEnergy := calculateEnergy(bodies)
	fmt.Printf("Initial energy: %e\n", initialEnergy)

	for step := 0; step < 10; step++ {
		root := buildKDTree(bodies, 0)
		kickStep(bodies, root)
	}

	finalEnergy := calculateEnergy(bodies)
	fmt.Printf("Final energy: %e\n", finalEnergy)
}
