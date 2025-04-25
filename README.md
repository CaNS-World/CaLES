<!--- the logo -->
<img src="assets/img/CaLES-logo.png" height=*>

<<<<<<< HEAD
CaLES (Canonical Large-Eddy Simulation) is a GPU-accelerated finite-difference solver for large-eddy simulations (LES) of incompressible wall-bounded flows in massively parallel environments. Built upon the DNS solver CaNS, it employs low-storage third-order Runge-Kutta schemes for time integration and supports implicit Crank-Nicolson treatment of viscous terms. A fast direct solver based on eigenfunction expansions is used for Poisson/Helmholtz equations. CaLES implements the Smagorinsky model (with van Driest damping), the dynamic Smagorinsky model, and a logarithmic law wall model. It supports various flow configurations, including channel flows, square ducts, tri-periodic domains, and lid-driven cavities. More details can be found in the official documentation: [CaLES Documentation](https://cales-test.readthedocs.io/en/latest/).
=======
## Overview

**CaLES (Canonical Large-Eddy Simulation)** is a GPU-accelerated finite-difference solver designed for large-eddy simulations (LES) of incompressible wall-bounded flows in massively parallel environments. Built upon the existing direct numerical simulation (DNS) solver [CaNS](https://github.com/CaNS-World/CaNS), CaLES relies on low-storage, third-order Runge-Kutta schemes for temporal discretization, with the option to treat viscous terms via an implicit Crank-Nicolson scheme in one or three directions. A fast direct solver, based on eigenfunction expansions, is used to solve the discretized Poisson/Helmholtz equations. For turbulence modeling, the classical Smagorinsky model with van Driest near-wall damping and the dynamic Smagorinsky model are implemented, along with a logarithmic law wall model. GPU acceleration is achieved through OpenACC directives. Some examples of flows that CaLES can solve include:

 * periodic or developing channel
 * periodic or developing square duct
 * tri-periodic domain
 * lid-driven cavity

In the `x` and `y` directions, the grid is uniform and the solver supports the following homogeneous boundary conditions for the Poisson equation:

 * Neumann-Neumann
 * Dirichlet-Dirichlet
 * Neumann-Dirichlet
 * Periodic

In the `z` direction, the solver offers more flexibility by using Gauss elimination, allowing for a non-uniform grid.

CaLES allows for choosing an implicit temporal discretization of the momentum diffusion terms, either fully implicit or only along the `z` direction. This results in solving a 3D/1D Helmholtz equation per velocity component. In the fully implicit case, FFT-based solvers are also used, and the same options described above for pressure boundary conditions apply to the velocity. This implies that if the wall model is used, it can only be applied in the `z` direction, as homogeneous boundary conditions are required in the other two directions where FFTs are applied. However, this limitation may be irrelevant, as explicit time integration is commonly used for wall-modeled LES due to the large thickness of the first off-wall layer of cells. Additionally, this restriction does not apply when the viscous terms are handled implicitly only in the `z`-direction.

Please refer to the theory manual [docs/CaLES_theory.pdf](docs/CaLES_theory.pdf) for more information on the methodology of CaLES. The theory manual is also available on Overleaf at the following link: [https://www.overleaf.com/read/ggfmdsymppxm#9b973d](https://www.overleaf.com/read/ggfmdsymppxm#9b973d).

## Implementation

Some of the key implementation features of CaLES are:

 * MPI parallelization
 * FFTW guru interface/cuFFT used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, cosine, sine, etc) is automatically determined form the input file
 * [cuDecomp](https://github.com/NVIDIA/cuDecomp) pencil decomposition library for hardware-adaptive distributed memory calculations on many GPUs
 * [2DECOMP&FFT](https://github.com/xcompact3d/2decomp-fft) library used for performing global data transpositions on CPUs and some of the data I/O
 * GPU acceleration using OpenACC directives
 * A different canonical flow can be simulated by simply changing the input files

Please refer to the following references for further details:

 * [Xiao, Maochao, Alessandro Ceci, Pedro Costa, Johan Larsson, and Sergio Pirozzoli. "CaLES: A GPU-accelerated solver for large-eddy simulation of wall-bounded flows." arXiv preprint arXiv:2411.09364 (2024).](https://arxiv.org/abs/2411.09364)

 * [Costa, Pedro, Everett Phillips, Luca Brandt, and Massimiliano Fatica. "GPU acceleration of CaNS for massively-parallel direct numerical simulations of canonical fluid flows." Computers & Mathematics with Applications 81 (2021): 502-511.](https://doi.org/10.1016/j.camwa.2020.01.002)

 * [Costa, Pedro. "A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows." Computers & Mathematics with Applications 76 (2018): 1853-1862.](https://doi.org/10.1016/j.camwa.2018.07.034)

## Getting Started

### Installation

Since CaLES loads the external pencil decomposition libraries as Git Submodules, the repository should be cloned as follows:
```bash
git clone --recursive https://github.com/soaringxmc/CaLES.git
```
so the libraries are downloaded too. Alternatively, in case the repository has already been cloned without the Submodules (i.e., folders `cuDecomp` and `2decomp-fft` under `dependencies/` are empty), the following command can be used to update them:
```bash
git submodule update --init --recursive
```

The prerequisites for compiling CaLES are the following:

 * MPI
 * FFTW3/cuFFT library for CPU/GPU runs
 * The `nvfortran` compiler (for GPU runs)
 * CMake for compiling the cuDecomp library (for GPU runs)
 * NCCL and NVSHMEM (optional, may be exploited by the cuDecomp library)

For most systems, CaLES can be compiled from the root directory with the following commands `make libs && make`, which will compile the 2DECOMP&FFT/cuDecomp libraries and CaLES.

The `Makefile` in root directory is used to compile the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), and a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options, see [`INFO_COMPILING.md`](docs/INFO_COMPILING.md) for more details. Concerning the pre-processing options, the following are available:

 * `DEBUG`                    : performs some basic checks for debugging purposes
 * `TIMING`                   : wall-clock time per time step is computed
 * `IMPDIFF`                  : diffusion terms are integrated implicitly in time (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `IMPDIFF_1D`               : same as above, but with implicit diffusion only along z; for optimal parallel performance this option should be combined with `PENCIL_AXIS=3`
 * `PENCIL_AXIS`              : sets the default pencil direction, one of [1,2,3] for [X,Y,Z]-aligned pencils; X-aligned is the default and should be optimal for all cases except for Z implicit diffusion, where using Z-pencils is recommended
 * `SINGLE_PRECISION`         : calculation will be carried out in single precision (the default precision is double)
 * `GPU`                      : enable GPU-accelerated runs

The executable `cales` will be created in the `build/` directory.

### Running Simulations

The input file `input.nml` sets the physical and computational parameters. In the `examples/` folder are examples of input files for several canonical flows. See [`INFO_INPUT.md`](docs/INFO_INPUT.md) for a detailed description of the input file.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2- and 3-dimensional output files, respectively. The code should be recompiled after editing out?d.h90 files.

Run the executable with `mpirun` with a number of tasks. For example, to run the code with 2 tasks:
```bash
mpirun -n 2 cales
```

See [`INFO_VISU.md`](docs/INFO_VISU.md) for how to visualize the output files.

## Contributing

We appreciate any contributions and feedback that can improve CaLES. If you wish to contribute to the tool, please get in touch with the maintainers or open an Issue in the repository / a thread in Discussions. Pull Requests are welcome, but please propose/discuss the changes in an linked Issue first.

To contribute to CaLES, please follow the steps listed below:

1. Fork this repository to your GitHub account
2. Create a new branch in your forked repository
3. Make changes in your new branch
4. Submit your changes by creating a Pull Request to the `develop` branch of the original PyFR repository

Modifications to the `develop` branch are eventually merged into the master
branch for a new release.
>>>>>>> 2cf7d482595dc9a2b6ba78103c6ade7697d86c06
