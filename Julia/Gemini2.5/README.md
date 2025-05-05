Yes, there are several Julia compiler and runtime flags that can potentially improve the performance of computationally intensive code like this N-body simulation. Here are the most relevant ones to consider:

Optimization Level (-O or --optimize)

-O3: This flag enables more aggressive optimizations than the default level (-O2). For CPU-bound numerical code, -O3 often provides noticeable performance gains. It might increase compilation time slightly, but that's usually negligible for a long-running simulation.
Usage: julia -O3 ... script.jl
Recommendation: Highly recommended to try.
Math Mode (--math-mode)

--math-mode=fast: This allows the compiler (LLVM) to perform optimizations that might slightly alter floating-point results according to strict IEEE 754 standards but can significantly speed up calculations. This includes using fused multiply-add (FMA) instructions where applicable and reordering operations (assuming associativity, which isn't strictly true for floats). For many simulations, the slight numerical difference is acceptable for the speed gain.
Usage: julia --math-mode=fast ... script.jl
Recommendation: Very likely beneficial for this type of code. Test for impact on results (e.g., energy drift) but often a good trade-off. (The default is --math-mode=ieee).
CPU Target (--cpu-target)

--cpu-target=native: This tells the compiler to optimize the generated code specifically for the CPU architecture of the machine you are running on.1 This allows it to use the most advanced instruction sets available (like AVX2, AVX512, FMA, etc.). This is the default setting in modern Julia versions, so you usually don't need to specify it explicitly, but it's good to be aware of it and ensure nothing is overriding it.   
1.
julia/Make.inc at master · JuliaLang/julia - GitHub

github.com

Usage (usually not needed): julia --cpu-target=native ... script.jl
Recommendation: Ensure this is the effective setting (default is usually fine).
Threads (-t or --threads)

-t auto or -t N: Essential for enabling the multi-threading used in the code (@threads). -t auto typically uses the number of logical cores available.
Usage: julia -t auto ... script.jl
Recommendation: Required for the threaded version.
Bounds Checking (--check-bounds)

--check-bounds=no: This globally disables array bounds checking. While this can improve performance by removing checks, it's generally unsafe. If there's an out-of-bounds access error in your code, it will likely lead to a crash or incorrect results without a clear error message. It's much safer to use the @inbounds macro selectively in performance-critical sections where you are certain the indices are valid (as done in the optimized code).
Usage: julia --check-bounds=no ... script.jl
Recommendation: Avoid unless absolutely necessary and code is rigorously tested; prefer using @inbounds locally.
Combining Flags:

You can combine these flags when running your script. A typical command line for running the optimized simulation aiming for maximum performance would be:

Bash

julia -t auto -O3 --math-mode=fast your_nbody_script.jl
Important Considerations:

Benchmarking: Always benchmark your code with and without these flags on your specific hardware to measure the actual performance impact. Sometimes -O3 might not yield much gain over -O2, or --math-mode=fast might have negligible effect if the bottleneck is elsewhere.
Numerical Accuracy: When using --math-mode=fast, carefully check if the numerical results (e.g., final positions, energy conservation) are still acceptable for your simulation's requirements.
Julia Version: Ensure you are using a reasonably recent version of Julia, as compiler optimizations improve over time.
In summary, -t auto -O3 --math-mode=fast is the combination most likely to provide the best runtime performance for this optimized N-body simulation code.