//To make the Go N-body simulation multithreaded, we can use goroutines with worker pools and channels to parallelize:
//Force computation
//Energy calculation
//(Optionally) integration steps (kick-step)
//Go doesn’t have a built-in "go-parallel" library like OpenMP, but goroutines and channels can be combined to achieve effective parallelism.
//Updated: Multithreaded N-body Simulation in Go
//This version:
//Uses a worker pool to parallelize computeForces() and calculateEnergy().
//Keeps the rest of the logic similar for clarity.

package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sync"
	"time"
)

const (
	G       = 6.67430e-11
	DT      = 3600.0 * 24.0
	N       = 10000
	STEPS   = 100
	CENTERM = 1e20
)

type Body struct {
	X, Y, Z    float64
	VX, VY, VZ float64
	AX, AY, AZ float64
	Mass       float64
}

var bodies []Body

func initializeSystem(n int) {
	bodies = make([]Body, n+1)

	// Central body
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

func computeForces(n int) {
	numWorkers := runtime.NumCPU()
	wg := sync.WaitGroup{}
	wg.Add(numWorkers)
	chunk := (n + 1 + numWorkers - 1) / numWorkers

	for w := 0; w < numWorkers; w++ {
		start := w * chunk
		end := (w + 1) * chunk
		if end > n+1 {
			end = n + 1
		}
		go func(start, end int) {
			defer wg.Done()
			for i := start; i < end; i++ {
				ax, ay, az := 0.0, 0.0, 0.0
				for j := 0; j <= n; j++ {
					if i == j {
						continue
					}
					dx := bodies[j].X - bodies[i].X
					dy := bodies[j].Y - bodies[i].Y
					dz := bodies[j].Z - bodies[i].Z
					distSq := dx*dx + dy*dy + dz*dz + 1e-10
					dist := math.Sqrt(distSq)
					force := G * bodies[j].Mass / (distSq * dist)
					ax += dx * force
					ay += dy * force
					az += dz * force
				}
				bodies[i].AX = ax
				bodies[i].AY = ay
				bodies[i].AZ = az
			}
		}(start, end)
	}
	wg.Wait()
}

func kickStep(n int) {
	computeForces(n)
	numWorkers := runtime.NumCPU()
	wg := sync.WaitGroup{}
	wg.Add(numWorkers)
	chunk := (n + 1 + numWorkers - 1) / numWorkers

	for w := 0; w < numWorkers; w++ {
		start := w * chunk
		end := (w + 1) * chunk
		if end > n+1 {
			end = n + 1
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

// func calculateEnergy(n int) float64 {
// 	kinetic := 0.0
// 	numWorkers := runtime.NumCPU()
// 	wg := sync.WaitGroup{}
// 	kinChan := make(chan float64, numWorkers)

// 	// Parallel kinetic energy
// 	wg.Add(numWorkers)
// 	chunk := (n + 1 + numWorkers - 1) / numWorkers
// 	for w := 0; w < numWorkers; w++ {
// 		start := w * chunk
// 		end := (w + 1) * chunk
// 		if end > n+1 {
// 			end = n + 1
// 		}
// 		go func(start, end int) {
// 			defer wg.Done()
// 			sum := 0.0
// 			for i := start; i < end; i++ {
// 				v2 := bodies[i].VX*bodies[i].VX + bodies[i].VY*bodies[i].VY + bodies[i].VZ*bodies[i].VZ
// 				sum += 0.5 * bodies[i].Mass * v2
// 			}
// 			kinChan <- sum
// 		}(start, end)
// 	}
// 	wg.Wait()
// 	close(kinChan)
// 	for k := range kinChan {
// 		kinetic += k
// 	}

// 	// Potential energy (still O(N^2), parallelize outer loop)
// 	potential := 0.0
// 	wg = sync.WaitGroup{}
// 	potChan := make(chan float64, numWorkers)
// 	wg.Add(numWorkers)
// 	for w := 0; w < numWorkers; w++ {
// 		start := w * chunk
// 		end := (w + 1) * chunk
// 		if end > n+1 {
// 			end = n + 1
// 		}
// 		go func(start, end int) {
// 			defer wg.Done()
// 			sum := 0.0
// 			for i := start; i < end; i++ {
// 				for j := i + 1; j <= n; j++ {
// 					dx := bodies[i].X - bodies[j].X
// 					dy := bodies[i].Y - bodies[j].Y
// 					dz := bodies[i].Z - bodies[j].Z
// 					dist := math.Sqrt(dx*dx + dy*dy + dz*dz + 1e-10)
// 					sum -= G * bodies[i].Mass * bodies[j].Mass / dist
// 				}
// 			}
// 			potChan <- sum
// 		}(start, end)
// 	}
// 	wg.Wait()
// 	close(potChan)
// 	for p := range potChan {
// 		potential += p
// 	}

// 	return kinetic + potential
// }

func main() {
	start := time.Now()
	fmt.Println("Initializing system...")
	initializeSystem(N)

	// fmt.Println("Calculating initial energy...")
	// initialEnergy := calculateEnergy(N)
	// fmt.Printf("Initial total energy: %.6e\n", initialEnergy)

	for step := 0; step < STEPS; step++ {
		kickStep(N)
		if step%100 == 0 {
			fmt.Printf("Step %d complete\n", step)
		}
	}

	// fmt.Println("Calculating final energy...")
	// finalEnergy := calculateEnergy(N)
	// fmt.Printf("Final total energy: %.6e\n", finalEnergy)
	// fmt.Printf("Energy difference:  %.6e\n", math.Abs(finalEnergy-initialEnergy))
	fmt.Printf("Body[0] %e %e %e", bodies[0].X, bodies[0].Y, bodies[0].Z)

	fmt.Printf("Simulation completed in %.2f seconds.\n", time.Since(start).Seconds())
}
