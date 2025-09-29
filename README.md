# lbcpp
A C++ template library for constructing models of fluid flow with active tracers using the lattice Boltzmann method.

## Installation instructions

Clone the repository. See the examples in the `simulations` directory for setting up models.

## Example simulations

These are contained in directories within the `simulations` directory. 
Compile them by running 
```
make
```
at the commandline, and then run the executables.
NB you may have to create an `output` directory in each simulation's directory.

A summary of the examples are below;

* Taylor-Green `simulations/taylorgreen'
    * Taylor-Green decaying vortex flow.
    * The simulation `TaylorGreen.cpp` tests four different initialisation schemes, at both 32-bit and 64-bit floating point precision.
    * See `simulations/taylorgreen/TaylorGreen.pdf` for more information on the initialisation schemes and some analysis of the errors produced by the different methods.
    * The design of this benchmark test is based on the description in section 5.5 of Kruger et al., 2017: The Lattice Boltzmann method: principles and practice.

* Gaussian hill `simulations/gaussianhill`
    * Advection-diffusion of a patch of scalar field concentration with a Gaussian hill distribution.
    * Based on benchmark test described in section 8.6.1 of Kruger et al., 2017.

* Diffusion plate `simulations/diffusionplate`
    * Diffusion from a plate in a uniform flow
    * Based on the benchmark test described in section 8.6.3 of Kruger et al., 2017.

* Rayleigh-Benard convection
    * Convection of a fluid heated below by a constant temperature hot plate 
    and cooled from above by a constant temperature cold plate.
    * Three versions:
        * `simulations/rayleighbenard2dsingle`: two-dimensional, single processor
        * `simulations/rayleighbenard3dsingle`: three-dimensional, single processor
        * `simulations/rayleighbenard2dmpi`: two-dimensional, mutliple processors
    * The latter case is run as follows:
    
    Generate the executable by running the Makefile
    ```
    make
    ```
    Then run the executable using `mpirun` (or `mpiexec`),
    ```
    mpirun -n 4 ./RayleighBenard2DMPI
    ```
    The `-n` flag specifies the number of processors. 
    For the code as written, the number of processors must be a multiple of four.