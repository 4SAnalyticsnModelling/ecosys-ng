# ecosys-ng

The next generation ecosys model, `ecosys-ng`, is a modern re-implementation of the original ecosystem model **ECOSYS** developed by **Professor Dr. Robert Grant (University of Alberta)**, originally written in **Fortran-77**. The goal of this project is to preserve the scientific foundations and process-based structure of **ECOSYS** while modernizing its software architecture for contemporary large-scale, high-performance computing workflows.

This modernization focuses on improving **scalability, numerical robustness, and spatial representation**, while preserving the underlying biophysical formulations of the original **ECOSYS** model.

# Key modernizations

**Explicit spatial representation:**
 Geographic structure is represented using explicit latitude–longitude grids with lateral connectivity between neighboring cells, enabling more realistic spatial interactions and landscape-scale simulations.

**Out-of-core tiling for large domains:**
 Designed for very large spatial extents, the model supports tile-based, out-of-core execution, allowing simulations to scale beyond available memory on HPC systems.

**Robust nonlinear solvers:**
 Sub-hourly solutions of heat, water, solute, and gas transport use damped Picard–Newton (modified Newton) methods, improving convergence stability under strongly coupled and highly nonlinear conditions.

**Parallel execution:**
 Computational kernels are parallelized using multithreading, enabling efficient utilization of modern multi-core CPUs without external runtime dependencies.

**Improved error handling and safety:**
 `ecosys-ng` emphasizes explicit error handling, and clearer control flow, making the model more maintainable and easier to debug than legacy Fortran implementations.

⚠️  `ecosys-ng` is currently under active development. The code is incomplete and not yet ready for production use.

# Compile ecosys-ng
**(compile with Zig version 0.15.2)**

=> `zig build` (ecosys-ng binary within `./zig-out/ecosys-ng-bin`)

=> `zig build -p . --verbose --summary all` (ecosys-ng binary within `./ecosys-ng-bin`)

=> `zig build -p . test --verbose --summary all` (run ecosys-ng code tests during compilation; find ecosys-ng binary within `./ecosys-ng-bin`)
