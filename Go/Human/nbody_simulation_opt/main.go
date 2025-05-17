package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"time"
)

func main() {
	steps, err := strconv.Atoi(os.Args[1])
	if err != nil {
		fmt.Println("You need to provide a number of steps.")
		return
	}
	n, err := strconv.Atoi(os.Args[2])
	if err != nil {
		fmt.Println("You need to provide a number of particles and number of steps.")
		return
	}
	p, err := strconv.Atoi(os.Args[3])
	if err != nil {
		fmt.Println("You need to provide a number of threads.")
		return
	}

	fmt.Printf("Running %d steps with %d particles using %d threads\n", steps, n, p)
	
	dt := 1e-3 // * 2.0 * std::f64::consts::PI;

	// Generate particles
	particles := circular_orbits(n)
	
	// Measure execution time
	start := time.Now()
	
	// Run simulation
	Simple_sim(particles, dt, steps, p)
	
	elapsed := time.Since(start)
	fmt.Printf("Simulation completed in %v\n", elapsed)
	fmt.Printf("Performance: %.2f particle-steps/second\n", float64(n*steps)/elapsed.Seconds())
}

func circular_orbits(n int) []Particle {
	particle_buf := make([]Particle, n)
	particle_buf[0] = Particle{
		[3]float64{0.0, 0.0, 0.0},
		[3]float64{0.0, 0.0, 0.0},
		0.00465047,
		1.0,
	}

	for i := 1; i < n; i++ {
		d := float64(0.1) + (float64(i) * 5.0 / float64(n))
		v := math.Sqrt(1.0 / d)
		theta := rand.Float64() * 6.28
		x := d * math.Cos(theta)
		y := d * math.Sin(theta)
		vx := -v * math.Sin(theta)
		vy := v * math.Cos(theta)
		particle_buf[i] = Particle{
			[3]float64{x, y, 0.0},
			[3]float64{vx, vy, 0.0},
			1e-14,
			1e-7,
		}
	}
	return particle_buf
}
