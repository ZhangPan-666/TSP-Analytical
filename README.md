# TSP Analytical

A C++17 library for computing gravity and magnetic fields of homogeneous triangular spherical prisms. Each prism has three minor great-circle edges and two constant-radius faces. The library evaluates gravity potential and its spatial derivatives through third order, and the exterior magnetic field and its gradient for uniform magnetization. File-based drivers sum contributions from any number of prisms at each observation.

The [regional MMSI solver](MMSI/) evaluates nearby prisms analytically and integrates distant sources together. It computes gravity and magnetic fields for laterally varying source layers and provides serial, OpenMP, and MPI drivers.

[Read the PDF user guide](docs/TSP_CPP_User_Guide.pdf) for the binary format, component order, units, API, and verification workflow.

## Build and run

Use CMake 3.20 or later and a C++17 compiler. Windows x64 builds use AVX2. The serial targets require only the compiler and standard library. CMake enables separate OpenMP and MPI targets when it detects the corresponding runtime and development files. Keep floating-point optimizations that reassociate operations disabled.

**Windows PowerShell** with the Visual Studio C++ tools and CMake on `PATH`

```powershell
cmake -S . -B build -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
.\build\Release\Work1_GravityEstimate.exe
.\build\Release\Work1_MagneticEstimate.exe
```

**Linux** with a C++17 compiler and CMake on `PATH`

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/Work1_GravityEstimate
./build/Work1_MagneticEstimate
```

Run the drivers from the repository root. `TFG.ForPar` and `TFM.ForPar` each contain one synthetic prism and two observations. The examples produce `Result_TFG.dat` and `Result_TFM.dat`, respectively. Use a separate working directory for each experiment, containing its own parameter file, to keep inputs and outputs paired. The MATLAB writers and readers use the same file layout.

## Execution modes

| Mode | Full-field targets | Selected components |
| --- | --- | --- |
| Serial | `Work1_GravityEstimate`, `Work1_MagneticEstimate` | C++ array interfaces |
| OpenMP | `Work1_GravityEstimateOpenMP`, `Work1_MagneticEstimateOpenMP` | C++ array interfaces |
| MPI | `Work1_GravityEstimateMPI`, `Work1_MagneticEstimateMPI` | `Work2_GravityEstimateVzMPI`, `Work2_MagneticEstimateVectorMPI` |

Use `OMP_NUM_THREADS` to control an OpenMP run. Initialize the matching MPI environment before configuring CMake and launch an MPI target with that environment's `mpiexec`. The MPI drivers use one computational thread per rank. Selected MPI outputs retain the full file layout and fill skipped components with `NaN`.

## Call the C++ API

Include `_TSP.h`. This example uses fixed Cartesian coordinates for the observation and magnetization, with latitudes and longitudes in degrees provided separately to `MakePrism`.

```cpp
#include "_TSP.h"

int main() {
    const auto prism = TSP::MakePrism(60., 60., 30., 0., 45., 20., 1., 2., 1.);
    const TSP::Vec point{0.29, -0.27, 3.1};
    const TSP::Vec magnetization{0.3, -0.2, 0.7};
    const auto gravity = TSP::GravityCartesian(prism, point);
    const auto magnetic = TSP::MagneticCartesian(prism, point, magnetization);
    return (std::isfinite(gravity[0]) && std::isfinite(magnetic[0])) ? 0 : 1;
}
```

File-based outputs use a fixed local north, east, and radial-up basis at each observation. For gravity, the kernel uses `G = 1` and includes the source density. Magnetic quantities use the raw potential-derivative contractions. Apply the conversion factors in the guide before interpreting the files as physical fields. Use consistent lengths throughout an input. Complete third-order results are supported away from physical source interfaces.

## Verification

CTest runs the three C++ reference regressions. The Python workflow also checks input validation, binary layout, and superposition. It compares OpenMP and serial outputs byte for byte when the OpenMP targets are present. The MPI workflow exercises 1, 2, and 4 ranks and selected-output files.

```powershell
python tests/VerifyWorkflow.py --bin-dir build/Release
python tests/VerifyMPI.py --bin-dir build/Release --mpiexec mpiexec
```

For a Linux single-configuration build, pass `--bin-dir build`. Use the MPI command after building the MPI targets. The regressions retain their original thresholds. The analytic solver requests absolute and relative tolerances of `1e-16` by default. Use reference discrepancies to report achieved accuracy for a given calculation.

## Source and notices

`_TSP.h` is the ordinary TSP kernel from the `TSP_Analytical` project. See [third-party notices](THIRD_PARTY.md) for the embedded QD arithmetic and its license. Record the commit identifier, compiler settings, run mode, and input checksums in reproducibility records.

## 中文速查

在仓库根目录编译并运行。重力示例读取 `TFG.ForPar`，磁场示例读取 `TFM.ForPar`。结果依次写入 `Result_TFG.dat` 和 `Result_TFM.dat`。实际单位换算、坐标方向及完整分量顺序见 [PDF 使用手册](docs/TSP_CPP_User_Guide.pdf)。
