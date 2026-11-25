# Regional gravity and magnetic modelling with MMSI

The solver evaluates analytic triangular spherical prisms near each observation and uses adaptive quadrature over the distant part of a source layer. The input mesh supplies interface radii and density or magnetization at its vertices. The full gravity and magnetic drivers return 20 and 10 quantities, respectively. Two additional MPI drivers compute vertical gravity alone or the three magnetic induction components.

## Build

From this directory, use GNU Make with a C++17 compiler. Select the target for the execution mode you need.

```sh
make PROFILE=gnu serial
make PROFILE=gnu openmp
make PROFILE=gnu mpi
```

The MPI target uses `mpicxx`. The Makefile also supports the Intel Classic (`PROFILE=intel`) and Intel LLVM (`PROFILE=intel-llvm`) compiler families.

## Run the shell example

Run `Work1_SetMSIGParameters.m` or `Work1_SetMSIMParameters.m` in MATLAB from this directory to write a gravity or magnetic input file. Source radii are in metres, density is in kg/m³, and magnetization is in A/m in global Cartesian coordinates.

| Field | Serial driver | Input | Output |
| --- | --- | --- | --- |
| Gravity | `./Work1_MassShellIntegration` | `MSIG.ForPar` | `ResultG.dat` |
| Magnetic | `./Work1_MagnetizationShellIntegration` | `MSIM.ForPar` | `ResultM.dat` |

Read the binary outputs with `Work2_ReadGravityForwardResults.m` or `Work2_ReadMagneticForwardResults.m`. The readers give the component names in file order. Gravity outputs use unit gravitational constant. For an MPI run, launch the corresponding executable with `mpiexec`. The supplied `run.pbs` is a Slurm script and sets one thread per MPI rank.
