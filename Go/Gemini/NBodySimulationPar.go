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
)

type Body struct {
	Position [3]float64
	Velocity [3]float64
	Mass     float64
}

func computeEnergy(bodies []Body) (kinetic, potential float64) {
	var wg sync.WaitGroup
	var mu sync.Mutex
	n := len(bodies)
	wg.Add(n)

	for i := 0; i < n; i++ {
		go func(i int) {
			defer wg.Done()
			speedSq := bodies[i].Velocity[0]*bodies[i].Velocity[0] +
				bodies[i].Velocity[1]*bodies[i].Velocity[1] +
				bodies[i].Velocity[2]*bodies[i].Velocity[2]
			kineticPart := 0.5 * bodies[i].Mass * speedSq
			potentialPart := 0.0
			for j := i + 1; j < n; j++ {
				dx := bodies[i].Position[0] - bodies[j].Position[0]
				dy := bodies[i].Position[1] - bodies[j].Position[1]
				dz := bodies[i].Position[2] - bodies[j].Position[2]
				r := math.Sqrt(dx*dx + dy*dy + dz*dz)
				if r > 0 {
					potentialPart -= G * bodies[i].Mass * bodies[j].Mass / r
				}
			}
			mu.Lock()
			kinetic += kineticPart
			potential += potentialPart
			mu.Unlock()
		}(i)
	}
	wg.Wait()
	return
}

func updatePositions(bodies []Body) {
	var wg sync.WaitGroup
	n := len(bodies)
	wg.Add(n)
	for i := 0; i < n; i++ {
		go func(i int) {
			defer wg.Done()
			for j := 0; j < 3; j++ {
				bodies[i].Position[j] += bodies[i].Velocity[j] * DT
			}
		}(i)
	}
	wg.Wait()
}

func updateVelocities(bodies []Body) {
	var wg sync.WaitGroup
	n := len(bodies)
	wg.Add(n)
	for i := 0; i < n; i++ {
		go func(i int) {
			defer wg.Done()
			ax, ay, az := 0.0, 0.0, 0.0
			for j := 0; j < n; j++ {
				if i != j {
					dx := bodies[j].Position[0] - bodies[i].Position[0]
					dy := bodies[j].Position[1] - bodies[i].Position[1]
					dz := bodies[j].Position[2] - bodies[i].Position[2]
					r := math.Sqrt(dx*dx + dy*dy + dz*dz)
					if r > 0 {
						factor := G * bodies[j].Mass / (r * r * r)
						ax += factor * dx
						ay += factor * dy
						az += factor * dz
					}
				}
			}
			bodies[i].Velocity[0] += ax * DT
			bodies[i].Velocity[1] += ay * DT
			bodies[i].Velocity[2] += az * DT
		}(i)
	}
	wg.Wait()
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

	initialKinetic, initialPotential := computeEnergy(bodies)
	fmt.Printf("Initial Energy: %e\n", initialKinetic+initialPotential)

	for step := 0; step < STEPS; step++ {
		updateVelocities(bodies)
		updatePositions(bodies)
	}

	finalKinetic, finalPotential := computeEnergy(bodies)
	fmt.Printf("Final Energy: %e\n", finalKinetic+finalPotential)
}
