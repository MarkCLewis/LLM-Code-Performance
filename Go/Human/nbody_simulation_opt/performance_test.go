package main

import (
	"fmt"
	"math/rand"
	"sync"
	"testing"
)

// Benchmark the KD-tree construction
func BenchmarkBuildTree(b *testing.B) {
	particleCounts := []int{1000, 10000, 100000}
	
	for _, count := range particleCounts {
		b.Run(fmt.Sprintf("Particles-%d", count), func(b *testing.B) {
			particles := circular_orbits(count)
			indices := make([]int, count)
			for i := 0; i < count; i++ {
				indices[i] = i
			}
			
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				tree := allocate_node_vec(count)
				build_tree(indices, 0, count, particles, 0, &tree)
				// Return tree to pool when done
				nodePool.Put(tree)
			}
		})
	}
}

// Benchmark the force calculation
func BenchmarkForceCalculation(b *testing.B) {
	particleCounts := []int{1000, 10000}
	
	for _, count := range particleCounts {
		b.Run(fmt.Sprintf("Particles-%d", count), func(b *testing.B) {
			particles := circular_orbits(count)
			indices := make([]int, count)
			for i := 0; i < count; i++ {
				indices[i] = i
			}
			
			// Build tree once
			tree := allocate_node_vec(count)
			build_tree(indices, 0, count, particles, 0, &tree)
			
			// Select a random particle to calculate forces on
			p := rand.Intn(count)
			
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				calc_accel(p, particles, tree)
			}
			
			// Return tree to pool when done
			nodePool.Put(tree)
		})
	}
}

// Benchmark complete simulation step
func BenchmarkSimulationStep(b *testing.B) {
	particleCounts := []int{1000, 10000}
	threadCounts := []int{1, 4, 8}
	
	for _, count := range particleCounts {
		for _, threads := range threadCounts {
			b.Run(fmt.Sprintf("Particles-%d-Threads-%d", count, threads), func(b *testing.B) {
				particles := circular_orbits(count)
				acc := make([][3]float64, count)
				indices := make([]int, count)
				for i := 0; i < count; i++ {
					indices[i] = i
				}
				
				dt := 1e-3
				tree := allocate_node_vec(count)
				
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					// Simulate one step
					build_tree(indices, 0, count, particles, 0, &tree)
					
					// Calculate accelerations in parallel
					workloadSize := (count + threads - 1) / threads
					var wg sync.WaitGroup
					for t := 0; t < threads; t++ {
						wg.Add(1)
						go func(threadID int) {
							defer wg.Done()
							start := threadID * workloadSize
							end := min((threadID+1)*workloadSize, count)
							
							for i := start; i < end; i++ {
								acc[i] = calc_accel(i, particles, tree)
							}
						}(t)
					}
					wg.Wait()
					
					// Update positions
					for t := 0; t < threads; t++ {
						wg.Add(1)
						go func(threadID int) {
							defer wg.Done()
							start := threadID * workloadSize
							end := min((threadID+1)*workloadSize, count)
							
							for i := start; i < end; i++ {
								particles[i].v[0] += dt * acc[i][0]
								particles[i].v[1] += dt * acc[i][1]
								particles[i].v[2] += dt * acc[i][2]
								
								particles[i].p[0] += dt * particles[i].v[0]
								particles[i].p[1] += dt * particles[i].v[1]
								particles[i].p[2] += dt * particles[i].v[2]
								
								acc[i][0] = 0.0
								acc[i][1] = 0.0
								acc[i][2] = 0.0
							}
						}(t)
					}
					wg.Wait()
				}
				
				// Return tree to pool when done
				nodePool.Put(tree)
			})
		}
	}
} 