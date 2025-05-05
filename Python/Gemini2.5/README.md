## Gemini 2.5 on 5/3/2025

### Notes for parallel

Install Numba: If you don't have it, install it: pip install numba
Import Numba: Add import numba at the beginning of the script.
Decorate calculate_acceleration: Apply the Numba JIT (Just-In-Time) compiler decorator @numba.njit with parallelization enabled.
Modify Loop Structure for Parallelism: The original nested loop for i... for j in range(i+1, ...) is slightly tricky for automatic parallelization because of the paired updates (accelerations[i] += ..., accelerations[j] -= ...). A more parallel-friendly approach calculates the total force on each particle i independently in the outer loop.
Use numba.prange: Replace range with numba.prange in the outer loop of the modified calculate_acceleration to signal Numba that this loop can be parallelized.
(Optional but Recommended) Parallelize Potential Energy: The calculate_energy function also has an O(N 
2
 ) loop for potential energy. We can apply a similar Numba parallelization strategy there for faster energy checks, especially if you were to call it more frequently.

 ### Notes for kD-Tree

 Sat in "thinking" forever. I was able to pull up the code by expanding "Show thinking", but it never finished.

 