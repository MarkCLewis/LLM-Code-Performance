/* The Computer Language Benchmarks Game
 * https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
 *
 * contributed by The Go Authors.
 * based on C program by Christoph Bauer
 * flag.Arg hack by Isaac Gouy
 * Modified by Antonio Petri
 */

package main

import (
	"flag"
	"math"
	"math/rand/v2"
	"strconv"
)

const (
	solarMass   = 4 * math.Pi * math.Pi
	daysPerYear = 365.24
)

type Position struct{ x, y, z float64 }
type Momentum struct{ x, y, z, m float64 }
type System struct {
	v []Momentum
	s []Position
}

func circular_orbits(n int) System {
	momentum := make([]Momentum, n)
	position := make([]Position, n)
	system := System{v: momentum, s: position}
	system.s[0] = Position{0.0, 0.0, 0.0}
	system.v[0] = Momentum{0.0, 0.0, 0.0, 1.0}

	for i := 1; i < n; i++ {
		d := float64(0.1) + (float64(i) * 5.0 / float64(n))
		v := math.Sqrt(1.0 / d)
		theta := rand.Float64() * 6.28
		x := d * math.Cos(theta)
		y := d * math.Sin(theta)
		vx := -v * math.Sin(theta)
		vy := v * math.Cos(theta)
		system.s[i] = Position{x, y, 0.0}
		system.v[i] = Momentum{vx, vy, 0.0, 1.0e-7}
	}
	return system
}

func energy(sys *System) float64 {
	N := len(sys.s)
	var e float64
	for i := 0; i < N; i++ {
		e += 0.5 * sys.v[i].m *
			(sys.v[i].x*sys.v[i].x + sys.v[i].y*sys.v[i].y + sys.v[i].z*sys.v[i].z)
		for j := i + 1; j < N; j++ {
			dx := sys.s[i].x - sys.s[j].x
			dy := sys.s[i].y - sys.s[j].y
			dz := sys.s[i].z - sys.s[j].z
			distance := math.Sqrt(dx*dx + dy*dy + dz*dz)
			e -= (sys.v[i].m * sys.v[j].m) / distance
		}
	}
	return e
}

func advance(dt float64, sys *System) {
	N := len(sys.s)
	for i := 0; i < N-1; i++ {
		_vx, _vy, _vz := sys.v[i].x, sys.v[i].y, sys.v[i].z

		for j := i + 1; j < N; j++ {

			dx := sys.s[i].x - sys.s[j].x
			dy := sys.s[i].y - sys.s[j].y
			dz := sys.s[i].z - sys.s[j].z

			dSquared := dx*dx + dy*dy + dz*dz
			distance := math.Sqrt(dSquared)
			mag := (dt / (dSquared * distance))
			mi := sys.v[i].m
			_vx -= dx * sys.v[j].m * mag
			_vy -= dy * sys.v[j].m * mag
			_vz -= dz * sys.v[j].m * mag
			sys.v[j].x += dx * mi * mag
			sys.v[j].y += dy * mi * mag
			sys.v[j].z += dz * mi * mag
		}
		sys.v[i].x, sys.v[i].y, sys.v[i].z = _vx, _vy, _vz
	}

	for i := 0; i < N; i++ {
		sys.s[i].x += dt * sys.v[i].x
		sys.s[i].y += dt * sys.v[i].y
		sys.s[i].z += dt * sys.v[i].z
	}
}

func main() {
	var steps int
	var n int
	flag.Parse()
	if flag.NArg() > 0 {
		steps, _ = strconv.Atoi(flag.Arg(0))
		n, _ = strconv.Atoi(flag.Arg(1))
	}

	sys := circular_orbits(n)

	// fmt.Printf("%.9f\n", energy())
	for i := 0; i < steps; i++ {
		advance(0.01, &sys)
	}
	// fmt.Printf("%.9f\n", energy())

}
