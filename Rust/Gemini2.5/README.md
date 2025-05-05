Yes, there are several flags and configuration options you can use with the Rust compiler (rustc) and Cargo to potentially squeeze out more performance, especially for CPU-intensive scientific computing tasks like this N-body simulation.

Here are the most relevant ones, typically configured either via the RUSTFLAGS environment variable or, more conveniently, in your project's .cargo/config.toml file:

1. Target CPU Features (-C target-cpu=native)

What it does: This is often the most impactful flag after --release. It tells the compiler (specifically the LLVM backend) to optimize the generated machine code specifically for the instruction set capabilities of the CPU you are compiling on. This allows it to use modern CPU instructions like AVX, AVX2, FMA (fused multiply-add), AVX512, etc., if your CPU supports them. These instructions can significantly speed up floating-point arithmetic, which is heavily used in your simulation. Without this, the compiler targets a more generic baseline architecture for broader compatibility.
How to use:
Environment variable: RUSTFLAGS="-C target-cpu=native" cargo build --release
.cargo/config.toml (Recommended):
Ini, TOML

[build]
rustflags = ["-C", "target-cpu=native"]
Caveat: The resulting binary might not run (or run optimally) on older CPUs that lack the specific instructions used. If you need portability, you might target a specific feature level instead (e.g., -C target-cpu=x86-64-v3).
2. Link-Time Optimization (LTO)

What it does: LTO allows the optimizer to perform optimizations across the entire program codebase (including all dependent crates) during the final linking stage, rather than just optimizing each crate individually. This can lead to better inlining across crate boundaries, more effective dead code elimination, and other global optimizations.
How to use (in Cargo.toml): Add this to your main Cargo.toml file (not .cargo/config.toml):
Ini, TOML

[profile.release]
lto = true # Recommended default, often similar to "thin"
# or
# lto = "thin" # Faster LTO, good balance.
# or
# lto = "fat"  # Slowest compile, potentially best optimization.
Caveat: LTO significantly increases link times (especially "fat" LTO). "thin" LTO offers a good compromise between compile time and runtime performance.
3. Codegen Units (codegen-units = 1) 1    
 1. 
Build Configuration - The Rust Performance Book

nnethercote.github.io


What it does: By default, Rust compiles crates in parallel using multiple "code generation units" to speed up compilation. However, this can prevent LLVM from performing certain optimizations (like inlining) across function boundaries within the same crate if those functions end up in different units. Setting codegen-units = 1 forces the entire crate to be treated as a single unit for code generation, giving LLVM maximum optimization scope within that crate.
How to use (in Cargo.toml): Add this to the [profile.release] section:
Ini, TOML

[profile.release]
codegen-units = 1
Caveat: This typically increases compilation time but can improve runtime performance, often working well in conjunction with LTO.
Combining Flags in Cargo.toml:

You would typically define these in the release profile within your Cargo.toml:

Ini, TOML

# In your Cargo.toml

[profile.release]
opt-level = 3      # Default for release, ensure it's set
debug = false      # Default for release, smaller binary, potentially faster
lto = "thin"       # Enable ThinLTO (good default)
codegen-units = 1  # Maximize intra-crate optimization scope
panic = 'abort'    # Optional: Slightly smaller/faster code by aborting on panic instead of unwinding
# strip = true     # Optional: Remove debug symbols for smaller binary size (no runtime perf impact)
And separately, enable native CPU targeting via .cargo/config.toml in your project's root directory (create the .cargo directory if it doesn't exist):

Ini, TOML

# In .cargo/config.toml (create this file/directory if needed)

[build]
rustflags = ["-C", "target-cpu=native"]
4. Profile-Guided Optimization (PGO)

What it does: PGO is a more advanced technique. You first compile your program with instrumentation (-C profile-generate), then run it with representative workloads to gather data about frequently executed code paths, branches taken, etc. Finally, you recompile the program using this profile data (-C profile-use), allowing the compiler to make much more informed optimization decisions (e.g., better inlining, code layout, branch prediction).
How to use: Requires multiple build steps and a representative workload. See the Rust PGO documentation for details.
Caveat: More complex workflow, but can provide significant gains (5-15%+) for complex applications.
5. Alternative Allocators (Less Likely Needed Now)

While Rust's default system allocator is good, sometimes specialized allocators like jemalloc or mimalloc can offer better performance for specific allocation patterns.
How: Add the dependency (e.g., tikv-jemallocator) and register it globally in your main.rs or lib.rs.
Relevance Here: Since you optimized the main allocation hotspot (tree nodes) using typed-arena, switching the global allocator is less likely to yield major gains for this specific code, but it's a general performance technique to be aware of.
Summary and Recommendation:

Definitely use --release.
Strongly recommend using -C target-cpu=native via .cargo/config.toml. This is often the biggest win for scientific code.
Enable LTO (e.g., lto = "thin") in Cargo.toml's [profile.release].
Consider setting codegen-units = 1 in Cargo.toml's [profile.release].
Benchmark! The impact of these flags can vary. Always measure the performance before and after applying them to confirm they are actually helping for your specific code and hardware. Sometimes, flags can interact in unexpected ways or even slightly regress performance.
Start with target-cpu=native and LTO, benchmark, then try adding codegen-units=1 and benchmark again.