package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
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

// GravitationalConstant is the universal gravitational constant.
const GravitationalConstant = 6.67430e-11

// distanceSquared calculates the squared Euclidean distance between two bodies.
func distanceSquared(b1, b2 *Body) float64 {
	dx := b1.Position[0] - b2.Position[0]
	dy := b1.Position[1] - b2.Position[1]
	dz := b1.Position[2] - b2.Position[2]
	return dx*dx + dy*dy + dz*dz
}

// calculateForce calculates the gravitational force on body i due to body j.
func calculateForce(bodyI, bodyJ *Body, force *[3]float64) {
	rSq := distanceSquared(bodyI, bodyJ)
	if rSq > 1e-9 { // Avoid division by zero for very close bodies
		r := math.Sqrt(rSq)
		magnitude := (GravitationalConstant * bodyI.Mass * bodyJ.Mass) / rSq
		force[0] += magnitude * (bodyJ.Position[0] - bodyI.Position[0]) / r
		force[1] += magnitude * (bodyJ.Position[1] - bodyI.Position[1]) / r
		force[2] += magnitude * (bodyJ.Position[2] - bodyI.Position[2]) / r
	}
}

// calculateTotalEnergy calculates the total energy of the system.
// func calculateTotalEnergy(s *System) float64 {
// 	kineticEnergy := 0.0
// 	potentialEnergy := 0.0

// 	for i := range s.Bodies {
// 		// Kinetic energy
// 		vSq := s.Bodies[i].Velocity[0]*s.Bodies[i].Velocity[0] +
// 			s.Bodies[i].Velocity[1]*s.Bodies[i].Velocity[1] +
// 			s.Bodies[i].Velocity[2]*s.Bodies[i].Velocity[2]
// 		kineticEnergy += 0.5 * s.Bodies[i].Mass * vSq

// 		// Potential energy
// 		for j := i + 1; j < len(s.Bodies); j++ {
// 			r := math.Sqrt(distanceSquared(&s.Bodies[i], &s.Bodies[j]))
// 			potentialEnergy -= (GravitationalConstant * s.Bodies[i].Mass * s.Bodies[j].Mass) / r
// 		}
// 	}

// 	return kineticEnergy + potentialEnergy
// }

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

// kickStep performs the "kick" part of the first-order kick-step method.
func kickStep(s *System, dt float64) {
	numBodies := len(s.Bodies)
	forces := make([][3]float64, numBodies)

	// Calculate forces
	for i := 0; i < numBodies; i++ {
		for j := 0; j < numBodies; j++ {
			if i != j {
				calculateForce(&s.Bodies[i], &s.Bodies[j], &forces[i])
			}
		}
	}

	// Update velocities
	for i := 0; i < numBodies; i++ {
		s.Bodies[i].Velocity[0] += (forces[i][0] / s.Bodies[i].Mass) * dt * 0.5
		s.Bodies[i].Velocity[1] += (forces[i][1] / s.Bodies[i].Mass) * dt * 0.5
		s.Bodies[i].Velocity[2] += (forces[i][2] / s.Bodies[i].Mass) * dt * 0.5
	}
}

// driftStep performs the "drift" part of the first-order kick-step method.
func driftStep(s *System, dt float64) {
	for i := range s.Bodies {
		s.Bodies[i].Position[0] += s.Bodies[i].Velocity[0] * dt
		s.Bodies[i].Position[1] += s.Bodies[i].Velocity[1] * dt
		s.Bodies[i].Position[2] += s.Bodies[i].Velocity[2] * dt
	}
}

// firstOrderKickStep performs one full first-order kick-step.
func firstOrderKickStep(s *System, dt float64) {
	kickStep(s, dt)
	driftStep(s, dt)
	kickStep(s, dt)
}

func main() {
	numOrbitingBodies := 10000
	centralMass := 1.989e30  // Mass of the Sun (kg)
	orbitRadius := 1.496e11  // 1 AU (m)
	orbitingMass := 5.972e24 // Mass of the Earth (kg)
	numSteps := 100
	timeStep := 3600.0 * 24.0 * 7.0 // 1 week in seconds

	rand.Seed(time.Now().UnixNano())

	// Initialize the system
	system := initializeCircularOrbits(numOrbitingBodies, centralMass, orbitRadius, orbitingMass)
	fmt.Printf("Initial number of bodies: %d\n", len(system.Bodies))

	// Calculate initial energy
	// initialEnergy := calculateTotalEnergy(system)
	// fmt.Printf("Initial total energy: %e J\n", initialEnergy)

	// Run the simulation
	fmt.Printf("Running simulation for %d steps...\n", numSteps)
	startTime := time.Now()
	for step := 0; step < numSteps; step++ {
		firstOrderKickStep(system, timeStep)
		if (step+1)%100 == 0 {
			fmt.Printf("Step %d completed.\n", step+1)
		}
	}
	elapsedTime := time.Since(startTime)
	fmt.Printf("Simulation finished in %s.\n", elapsedTime)

	// Calculate final energy
	// finalEnergy := calculateTotalEnergy(system)
	// fmt.Printf("Final total energy: %e J\n", finalEnergy)

	// // Calculate the energy difference
	// energyDifference := math.Abs(finalEnergy - initialEnergy)
	// relativeEnergyDifference := energyDifference / math.Abs(initialEnergy)
	// fmt.Printf("Absolute energy difference: %e J\n", energyDifference)
	// fmt.Printf("Relative energy difference: %e\n", relativeEnergyDifference)
	fmt.Printf("Body[0] %e %e %e", system.Bodies[0].Position[0], system.Bodies[0].Position[1], system.Bodies[0].Position[2])

}
