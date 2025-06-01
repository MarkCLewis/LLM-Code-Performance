package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"runtime"
	"sync"

	"github.com/dgravesa/go-parallel/parallel"
)

const (
	maxParts int     = 7
	THETA    float64 = 0.3
)

// KDTree structure with optimized memory layout
type KDTree struct {
	// For leaves
	num_parts int
	particles [maxParts]int

	// For internal nodes
	split_dim int
	split_val float64
	m         float64
	cm        [3]float64
	size      float64
	left      int
	right     int
}

// ObjectPool implements a simple object pool for KDTree nodes
type TreeNodePool struct {
	pool [][]KDTree
	mu   sync.Mutex
}

// Global node pool
var nodePool = &TreeNodePool{
	pool: make([][]KDTree, 0),
}

// Get a node slice from the pool or create a new one
func (p *TreeNodePool) Get(size int) []KDTree {
	p.mu.Lock()
	defer p.mu.Unlock()

	for i, nodes := range p.pool {
		if len(nodes) >= size {
			// Found a suitable slice in the pool
			p.pool = append(p.pool[:i], p.pool[i+1:]...)
			return nodes
		}
	}

	// Create a new slice if none available in pool
	return make([]KDTree, size)
}

// Return a node slice to the pool
func (p *TreeNodePool) Put(nodes []KDTree) {
	p.mu.Lock()
	defer p.mu.Unlock()

	// Clear node data to avoid memory leaks
	for i := range nodes {
		nodes[i] = empty_leaf()
	}

	p.pool = append(p.pool, nodes)
}

func empty_leaf() KDTree {
	return KDTree{
		0,
		[maxParts]int{},
		-1,
		0.0,
		0.0,
		[3]float64{0.0, 0.0, 0.0},
		0.0,
		-1,
		-1,
	}
}

func allocate_node_vec(num_parts int) []KDTree {
	num_nodes := 20 * (num_parts/(maxParts-1) + 1)
	// Get nodes from the pool instead of allocating new ones
	return nodePool.Get(num_nodes)
}

// Returns the index of the last Node used in the construction.
func build_tree(
	indices []int,
	start int,
	end int,
	particles []Particle,
	cur_node int,
	nodes *[]KDTree,
) int {
	// println!("start = {} end = {} cur_node = {}", start, end, cur_node);
	np := end - start
	// println!("s = {}, e = {}, cn = {}", start, end, cur_node);
	if np <= maxParts {
		for cur_node >= len(*nodes) {
			*nodes = append(*nodes, empty_leaf())
		}
		(*nodes)[cur_node].num_parts = np
		for i := 0; i < np; i++ {
			(*nodes)[cur_node].particles[i] = indices[start+i]
		}
		return cur_node
	} else {
		// Pick split dim and value
		min := [3]float64{1e100, 1e100, 1e100}
		max := [3]float64{-1e100, -1e100, -1e100}
		m := 0.0
		cm := [3]float64{0.0, 0.0, 0.0}

		// Batch process particles for better cache locality
		for i := start; i < end; i++ {
			p := &particles[indices[i]]
			m += p.m
			cm[0] += p.m * p.p[0]
			cm[1] += p.m * p.p[1]
			cm[2] += p.m * p.p[2]

			// Calculate min/max in the same loop
			if p.p[0] < min[0] {
				min[0] = p.p[0]
			}
			if p.p[1] < min[1] {
				min[1] = p.p[1]
			}
			if p.p[2] < min[2] {
				min[2] = p.p[2]
			}
			if p.p[0] > max[0] {
				max[0] = p.p[0]
			}
			if p.p[1] > max[1] {
				max[1] = p.p[1]
			}
			if p.p[2] > max[2] {
				max[2] = p.p[2]
			}
		}

		// Calculate center of mass
		if m > 0 {
			invM := 1.0 / m
			cm[0] *= invM
			cm[1] *= invM
			cm[2] *= invM
		}

		// Find dimension with greatest spread
		split_dim := 0
		max_spread := max[0] - min[0]
		for dim := 1; dim < 3; dim++ {
			spread := max[dim] - min[dim]
			if spread > max_spread {
				max_spread = spread
				split_dim = dim
			}
		}
		size := max_spread

		// Optimized partitioning for better cache behavior
		mid := (start + end) / 2

		// Use quickselect algorithm for better performance
		quickSelect(indices, start, end-1, mid, split_dim, particles)
		split_val := particles[indices[mid]].p[split_dim]

		// Recurse on children and build this node.
		left := build_tree(indices, start, mid, particles, cur_node+1, nodes)
		right := build_tree(indices, mid, end, particles, left+1, nodes)

		for cur_node >= len(*nodes) {
			*nodes = append(*nodes, empty_leaf())
		}
		(*nodes)[cur_node].num_parts = 0
		(*nodes)[cur_node].split_dim = split_dim
		(*nodes)[cur_node].split_val = split_val
		(*nodes)[cur_node].m = m
		(*nodes)[cur_node].cm = cm
		(*nodes)[cur_node].size = size
		(*nodes)[cur_node].left = cur_node + 1
		(*nodes)[cur_node].right = left + 1

		return right
	}
}

// QuickSelect algorithm for faster median finding - O(n) average case
func quickSelect(indices []int, left, right, k, split_dim int, particles []Particle) {
	if left == right {
		return
	}

	pivotIndex := left + rand.Intn(right-left+1)
	pivotIndex = partition(indices, left, right, pivotIndex, split_dim, particles)

	if k == pivotIndex {
		return
	} else if k < pivotIndex {
		quickSelect(indices, left, pivotIndex-1, k, split_dim, particles)
	} else {
		quickSelect(indices, pivotIndex+1, right, k, split_dim, particles)
	}
}

func partition(indices []int, left, right, pivotIndex, split_dim int, particles []Particle) int {
	pivotValue := particles[indices[pivotIndex]].p[split_dim]

	// Move pivot to end
	indices[pivotIndex], indices[right] = indices[right], indices[pivotIndex]

	// Move all elements smaller than pivot to the left
	storeIndex := left
	for i := left; i < right; i++ {
		if particles[indices[i]].p[split_dim] < pivotValue {
			indices[storeIndex], indices[i] = indices[i], indices[storeIndex]
			storeIndex++
		}
	}

	// Move pivot to its final place
	indices[right], indices[storeIndex] = indices[storeIndex], indices[right]

	return storeIndex
}

// Optimized force calculation with iterative approach to reduce recursion overhead
func accel_recur(cur_node int, p int, particles []Particle, nodes []KDTree) [3]float64 {
	acc := [3]float64{0.0, 0.0, 0.0}

	// Use a stack to simulate recursion
	type StackItem struct {
		nodeIndex int
	}

	stack := make([]StackItem, 0, 64) // Preallocate stack to avoid resizing
	stack = append(stack, StackItem{cur_node})

	for len(stack) > 0 {
		// Pop from stack
		n := len(stack) - 1
		item := stack[n]
		stack = stack[:n]

		node := nodes[item.nodeIndex]

		if node.num_parts > 0 {
			// Leaf node - direct calculation with all particles
			for i := 0; i < node.num_parts; i++ {
				if node.particles[i] != p {
					pp_acc := Calc_pp_accel(&particles[p], &particles[node.particles[i]])
					acc[0] += pp_acc[0]
					acc[1] += pp_acc[1]
					acc[2] += pp_acc[2]
				}
			}
		} else {
			// Calculate distance to center of mass
			dx := particles[p].p[0] - node.cm[0]
			dy := particles[p].p[1] - node.cm[1]
			dz := particles[p].p[2] - node.cm[2]
			dist_sqr := dx*dx + dy*dy + dz*dz

			// Use Barnes-Hut approximation criterion
			if node.size*node.size < THETA*THETA*dist_sqr {
				// Far enough to use approximation
				dist := math.Sqrt(dist_sqr)
				magi := -node.m / (dist_sqr * dist)
				acc[0] += dx * magi
				acc[1] += dy * magi
				acc[2] += dz * magi
			} else {
				// Too close, need to check children
				// Push right child first so left is processed first (depth-first traversal)
				stack = append(stack, StackItem{node.right})
				stack = append(stack, StackItem{node.left})
			}
		}
	}

	return acc
}

// Original recursive implementation kept for reference
func accel_recur_original(cur_node int, p int, particles []Particle, nodes []KDTree) [3]float64 {
	if nodes[cur_node].num_parts > 0 {
		acc := [3]float64{0.0, 0.0, 0.0}
		for i := 0; i < nodes[cur_node].num_parts; i++ {
			if nodes[cur_node].particles[i] != p {
				pp_acc := Calc_pp_accel(&particles[p], &particles[nodes[cur_node].particles[i]])
				acc[0] += pp_acc[0]
				acc[1] += pp_acc[1]
				acc[2] += pp_acc[2]
			}
		}
		return acc
	} else {
		dx := particles[p].p[0] - nodes[cur_node].cm[0]
		dy := particles[p].p[1] - nodes[cur_node].cm[1]
		dz := particles[p].p[2] - nodes[cur_node].cm[2]
		dist_sqr := dx*dx + dy*dy + dz*dz
		// println!("dist = {}, size = {}", dist, nodes[cur_node].size);
		if nodes[cur_node].size*nodes[cur_node].size < THETA*THETA*dist_sqr {
			dist := math.Sqrt(dist_sqr)
			magi := -nodes[cur_node].m / (dist_sqr * dist)
			return [3]float64{dx * magi, dy * magi, dz * magi}
		} else {
			left_acc := accel_recur_original(nodes[cur_node].left, p, particles, nodes)
			right_acc := accel_recur_original(nodes[cur_node].right, p, particles, nodes)
			return [3]float64{left_acc[0] + right_acc[0], left_acc[1] + right_acc[1], left_acc[2] + right_acc[2]}
		}
	}
}

func calc_accel(p int, particles []Particle, nodes []KDTree) [3]float64 {
	return accel_recur(0, p, particles, nodes)
}

// Optimized simulation with better memory management and parallelism
func Simple_sim(bodies []Particle, dt float64, steps int, p int) {
	// Determine optimal parallelism level if not specified
	if p <= 0 {
		p = runtime.NumCPU()
	}

	// Pre-allocate all memory needed for the simulation
	n := len(bodies)
	acc := make([][3]float64, n)
	tree := allocate_node_vec(n)
	indices := make([]int, n)

	// Pre-initialize indices to avoid repeated initialization
	for i := 0; i < n; i++ {
		indices[i] = i
	}

	// Calculate initial energy for verification
	// initialEnergy := CalculateSystemEnergy(bodies)
	// fmt.Printf("Initial system energy: %e\n", initialEnergy)

	for step := 0; step < steps; step++ {
		// Build tree - no need to reinitialize indices every step
		build_tree(indices, 0, n, bodies, 0, &tree)

		// Calculate accelerations in parallel
		parallel.WithNumGoroutines(p).For(n, func(i, _ int) {
			acc[i] = calc_accel(i, bodies, tree)
		})

		// Update positions in parallel
		parallel.WithNumGoroutines(p).For(n, func(i, _ int) {
			// Update velocity
			bodies[i].v[0] += dt * acc[i][0]
			bodies[i].v[1] += dt * acc[i][1]
			bodies[i].v[2] += dt * acc[i][2]

			// Update position
			bodies[i].p[0] += dt * bodies[i].v[0]
			bodies[i].p[1] += dt * bodies[i].v[1]
			bodies[i].p[2] += dt * bodies[i].v[2]

			// Clear acceleration for next step
			acc[i][0] = 0.0
			acc[i][1] = 0.0
			acc[i][2] = 0.0
		})
	}

	// Return tree to the pool
	nodePool.Put(tree)

	// Calculate final energy and compare with initial energy
	// finalEnergy := CalculateSystemEnergy(bodies)
	// fmt.Printf("Final system energy: %e\n", finalEnergy)
	// energyDiff := math.Abs(finalEnergy - initialEnergy)
	// relativeError := energyDiff / math.Abs(initialEnergy)
	// fmt.Printf("Energy difference: %e (relative error: %e)\n", energyDiff, relativeError)
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func print_tree(step int, tree []KDTree, particles []Particle) {
	fname := fmt.Sprintf("tree%d.txt", step)
	file, err := os.Create(fname)
	if err != nil {
		panic(err)
	}

	len_line := fmt.Sprintf("%d\n", len(tree))
	file.WriteString(len_line)
	for _, n := range tree {
		if n.num_parts > 0 {
			line := fmt.Sprintf("L %d\n", n.num_parts)
			file.WriteString(line)
			for i := 0; i < n.num_parts; i++ {
				p := n.particles[i]
				node_line := fmt.Sprintf(
					"%e %e %e\n",
					particles[p].p[0], particles[p].p[1], particles[p].p[2],
				)
				file.WriteString(node_line)
			}
		} else {
			leaf_line := fmt.Sprintf("I %d %e %d %d\n", n.split_dim, n.split_val, n.left, n.right)
			file.WriteString(leaf_line)
		}
	}
}
