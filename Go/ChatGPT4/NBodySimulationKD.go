package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"
	"time"
)

const (
	G       = 6.67430e-11
	DT      = 3600.0 * 24.0
	N       = 100000
	STEPS   = 10
	THETA   = 0.3
	CENTERM = 1e20
)

type Body struct {
	X, Y, Z    float64
	VX, VY, VZ float64
	AX, AY, AZ float64
	Mass       float64
}

type KDNode struct {
	Min, Max     [3]float64
	CenterOfMass [3]float64
	TotalMass    float64
	Body         *Body
	Left, Right  *KDNode
}

var bodies []Body

func initializeSystem(n int) {
	bodies = make([]Body, n+1)

	// Central massive body
	bodies[0] = Body{0, 0, 0, 0, 0, 0, 0, 0, 0, CENTERM}

	radius := 1e7
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i := 1; i <= n; i++ {
		angle := 2 * math.Pi * float64(i) / float64(n)
		r := radius * (1.0 + 0.1*rng.Float64())
		x := r * math.Cos(angle)
		y := r * math.Sin(angle)
		v := math.Sqrt(G * CENTERM / r)
		bodies[i] = Body{
			X:    x,
			Y:    y,
			Z:    0,
			VX:   -v * math.Sin(angle),
			VY:   v * math.Cos(angle),
			VZ:   0,
			Mass: 1.0,
		}
	}
}

func buildKDTree(indices []int, depth int) *KDNode {
	if len(indices) == 0 {
		return nil
	}
	axis := depth % 3
	sortIndices := func(i, j int) bool {
		switch axis {
		case 0:
			return bodies[indices[i]].X < bodies[indices[j]].X
		case 1:
			return bodies[indices[i]].Y < bodies[indices[j]].Y
		default:
			return bodies[indices[i]].Z < bodies[indices[j]].Z
		}
	}
	sort.Slice(indices, sortIndices)
	mid := len(indices) / 2
	node := &KDNode{}
	node.Body = &bodies[indices[mid]]

	if len(indices) == 1 {
		b := node.Body
		node.CenterOfMass = [3]float64{b.X, b.Y, b.Z}
		node.TotalMass = b.Mass
		node.Min = [3]float64{b.X, b.Y, b.Z}
		node.Max = [3]float64{b.X, b.Y, b.Z}
		return node
	}

	node.Left = buildKDTree(indices[:mid], depth+1)
	node.Right = buildKDTree(indices[mid+1:], depth+1)

	var m, cx, cy, cz float64
	min := [3]float64{math.Inf(1), math.Inf(1), math.Inf(1)}
	max := [3]float64{math.Inf(-1), math.Inf(-1), math.Inf(-1)}

	nodes := []*KDNode{node.Left, node.Right}
	for _, n := range nodes {
		if n == nil {
			continue
		}
		m += n.TotalMass
		cx += n.CenterOfMass[0] * n.TotalMass
		cy += n.CenterOfMass[1] * n.TotalMass
		cz += n.CenterOfMass[2] * n.TotalMass
		for i := 0; i < 3; i++ {
			if n.Min[i] < min[i] {
				min[i] = n.Min[i]
			}
			if n.Max[i] > max[i] {
				max[i] = n.Max[i]
			}
		}
	}
	node.TotalMass = m
	node.CenterOfMass = [3]float64{cx / m, cy / m, cz / m}
	node.Min = min
	node.Max = max
	return node
}

func computeForceKD(b *Body, node *KDNode, theta float64) (fx, fy, fz float64) {
	if node == nil || (node.Body == b && node.Left == nil && node.Right == nil) {
		return
	}

	dx := node.CenterOfMass[0] - b.X
	dy := node.CenterOfMass[1] - b.Y
	dz := node.CenterOfMass[2] - b.Z
	distSq := dx*dx + dy*dy + dz*dz + 1e-10
	dist := math.Sqrt(distSq)

	s := math.Max(math.Max(node.Max[0]-node.Min[0], node.Max[1]-node.Min[1]), node.Max[2]-node.Min[2])
	if s/dist < theta || (node.Left == nil && node.Right == nil) {
		force := G * node.TotalMass / (distSq * dist)
		return dx * force, dy * force, dz * force
	}

	lx, ly, lz := computeForceKD(b, node.Left, theta)
	rx, ry, rz := computeForceKD(b, node.Right, theta)
	return lx + rx, ly + ry, lz + rz
}

func computeForcesKD(tree *KDNode) {
	numWorkers := runtime.NumCPU()
	wg := sync.WaitGroup{}
	chunk := (len(bodies) + numWorkers - 1) / numWorkers
	wg.Add(numWorkers)

	for w := 0; w < numWorkers; w++ {
		start := w * chunk
		end := (w + 1) * chunk
		if end > len(bodies) {
			end = len(bodies)
		}
		go func(start, end int) {
			defer wg.Done()
			for i := start; i < end; i++ {
				fx, fy, fz := computeForceKD(&bodies[i], tree, THETA)
				bodies[i].AX = fx
				bodies[i].AY = fy
				bodies[i].AZ = fz
			}
		}(start, end)
	}
	wg.Wait()
}

func kickStepKD(tree *KDNode) {
	computeForcesKD(tree)
	numWorkers := runtime.NumCPU()
	wg := sync.WaitGroup{}
	chunk := (len(bodies) + numWorkers - 1) / numWorkers
	wg.Add(numWorkers)

	for w := 0; w < numWorkers; w++ {
		start := w * chunk
		end := (w + 1) * chunk
		if end > len(bodies) {
			end = len(bodies)
		}
		go func(start, end int) {
			defer wg.Done()
			for i := start; i < end; i++ {
				b := &bodies[i]
				b.VX += b.AX * DT
				b.VY += b.AY * DT
				b.VZ += b.AZ * DT
				b.X += b.VX * DT
				b.Y += b.VY * DT
				b.Z += b.VZ * DT
			}
		}(start, end)
	}
	wg.Wait()
}

func calculateEnergy() float64 {
	kinetic := 0.0
	for i := range bodies {
		v2 := bodies[i].VX*bodies[i].VX + bodies[i].VY*bodies[i].VY + bodies[i].VZ*bodies[i].VZ
		kinetic += 0.5 * bodies[i].Mass * v2
	}
	return kinetic // Approximate, omits potential energy for brevity
}

func main() {
	start := time.Now()
	fmt.Println("Initializing system...")
	initializeSystem(N)

	fmt.Println("Building KD-tree...")
	indices := make([]int, len(bodies))
	for i := range bodies {
		indices[i] = i
	}
	tree := buildKDTree(indices, 0)

	fmt.Println("Calculating initial energy...")
	initialEnergy := calculateEnergy()
	fmt.Printf("Initial kinetic energy: %.6e\n", initialEnergy)

	for step := 0; step < STEPS; step++ {
		tree = buildKDTree(indices, 0)
		kickStepKD(tree)
		if step%10 == 0 {
			fmt.Printf("Step %d complete\n", step)
		}
	}

	fmt.Println("Calculating final energy...")
	finalEnergy := calculateEnergy()
	fmt.Printf("Final kinetic energy: %.6e\n", finalEnergy)
	fmt.Printf("Energy difference:    %.6e\n", math.Abs(finalEnergy-initialEnergy))
	fmt.Printf("Simulation time:      %.2f seconds\n", time.Since(start).Seconds())
}
