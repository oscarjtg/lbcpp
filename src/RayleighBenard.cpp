#include <iostream> // For cout.
#include <cstdlib>  // For atof and exit.
#include <cmath>    // For pow and sqrt.
#include <mpi.h>
#include <sstream>
#include <iomanip>

#include "equilibria/Equilibria.h"
#include "lattice/LatticeD2Q5.h"
#include "lattice/LatticeD2Q9.h"
#include "macroscopic/Macroscopic2D.h"
#include "RayleighBenardStructs.h"

// Functions.
CharacteristicScales computeCharScales(DimensionlessNumbers dn, SimulationParameters sp, GridParameters grid);
void displayParameters(GridParameters grid, DimensionlessNumbers dn, SimulationParameters sp, CharacteristicScales scales);
DimensionlessNumbers getInputParameters(int argc, char* argv[]);
SimulationParameters goodParameters(DimensionlessNumbers dn, GridParameters grid);

// Templated functions.
template<class T>
inline T bodyforce(T temperature, SimulationParameters sp) {
    return static_cast<T>(sp.ag) * temperature;
}

template <class T>
void initialiseDistributions(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, GridParameters grid);

template<class T>
void initialiseMacroscopic(Macroscopic2D<T>& macros, GridParameters grid);

template<class T>
void runTimestep(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, GridParameters grid);

/*****************************************************************************
 *                                                                           *
 *                                   Main                                    *
 *                                                                           *
 ****************************************************************************/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int process_number, number_of_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_number);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);

    std::ostringstream oss;
    oss << "00000_" << process_number << "_" << number_of_ranks;
    std::string runId = oss.str();

    // Define constants.
    const int NX = 100;
    const int NY = 100;
    const int NT = 1000;
    // nx, ny, size
    GridParameters grid(NX, NY);

    // rayleigh and prandtl numbers, from command line or default values.
    DimensionlessNumbers dn = getInputParameters(argc, argv);

    // ag, omega_f, and omega_g in lattice units.
    SimulationParameters sp = goodParameters(dn, grid);

    // Calculate the characteristic scales of the simulated flow.
    CharacteristicScales scales = computeCharScales(dn, sp, grid);

    // Display all the parameters.
    displayParameters(grid, dn, sp, scales);
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate the arrays for the simulation.
    Macroscopic2D<double> macros(grid.nx, grid.ny);
    LatticeD2Q9<double> f(grid.nx, grid.ny);
    LatticeD2Q5<double> g(grid.nx, grid.ny);

    // Fill macroscopic arrays with initial conditions.
    initialiseMacroscopic(macros, grid);

    // Fill distribution function arrays with initial conditions.
    initialiseDistributions(f, g, macros, grid);

    // Save data to CSV.
    macros.WriteToCSV("output", runId);
    f.WriteToCSV("output", 'f', runId);
    g.WriteToCSV("output", 'g', runId);

    MPI_Barrier(MPI_COMM_WORLD);

    /*
    // The main loop.
    for (int t=0; t<NT; ++t)
    {
        runTimestep();
    }
    */

    // Close the MPI processes nicely (free the allocated memory, sever network connections, etc...)
    MPI_Finalize();

    return 0;
}


/*****************************************************************************
 *                                                                           *
 *                          Function definitions                             *
 *                                                                           *
 ****************************************************************************/

/**
 * @brief Computes the characteristic scales for the Rayleigh Benard problem.
 * 
 * @param dn A DimensionlessNumbers instance holding Ra and Pr for this flow.
 * @param sp A SimulationParameters instance for this simulation.
 * @param grid A GridParameters instance holding the grid characteristics.
 * @return The velocity, time, and length scales of the flow, 
 *         as an instance of the CharacteristicScales struct.
 */
CharacteristicScales computeCharScales(DimensionlessNumbers dn, SimulationParameters sp, GridParameters grid)
{
    // Temporary variables.
    double ny2 = pow(static_cast<double>(grid.ny), 2.0); // NY^2.
    double nu = ((1.0 / sp.omega_f) - 0.5) / 3.0; // Kinematic viscosity.
    double u = sqrt(sp.ag * grid.ny * 1.0); // delta T = 1.0.
    double t = (ny2 / nu) * sqrt(dn.prandtl / dn.rayleigh);
    double h = u * t;
    CharacteristicScales scales(u, t, h);
    return scales;
}

/**
 * @brief Display grid, dimensionless, and simulation parameters to screen.
 * 
 * @param dn A DimensionlessNumbers instance holding Ra and Pr for this flow.
 * @param sp A SimulationParameters instance for this simulation.
 * @param grid A GridParameters instance holding the grid characteristics.
 * @param scales A CharacteristicScales instance.
 */
void displayParameters(GridParameters grid, DimensionlessNumbers dn, SimulationParameters sp, CharacteristicScales scales)
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Grid parameters:\n";
    std::cout << "    nx               = " << grid.nx << "\n";
    std::cout << "    ny               = " << grid.ny << "\n";
    std::cout << "    size             = " << grid.size << "\n";
    std::cout << "Dimensionless numbers:\n";
    std::cout << "    rayleigh         = " << dn.rayleigh << "\n";
    std::cout << "    prandtl          = " << dn.prandtl << "\n";
    std::cout << "Simulation parameters:\n";
    std::cout << "    ag               = " << sp.ag << "\n";
    std::cout << "    omega_f          = " << sp.omega_f << "\n";
    std::cout << "    omega_g          = " << sp.omega_g << "\n";
    std::cout << "Characteristic scales:\n";
    std::cout << "    length scale     = " << scales.length << "\n";
    std::cout << "    time scale       = " << scales.time << "\n";
    std::cout << "    velocity scale   = " << scales.velocity << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

/** @brief Gets non-dimensional parameters from the command line, 
 * or sets default values.
 * 
 * @param argc The number of command line arguments.
 * @param argv Array of pointers to the command line arguments.
 * @return The Rayleigh and Prandtl numbers for the simulation 
 * as a DimensionlessNumbers instance.
 */
DimensionlessNumbers getInputParameters(int argc, char* argv[])
{
    if (argc == 3)
    {
        // Values from command line.
        DimensionlessNumbers dn(atof(argv[1]), atof(argv[2]));
        return dn;
    }
    else
    {
        // Default values (hard-coded).
        DimensionlessNumbers dn(1.0e3, 1.0);
        return dn;
    }
}

/** @brief Finds accurate and stable parameters for the simulation.
 * @details Iteratively tries combinations of parameter values that are
 * consistent with the dimensionless numbers governing the system
 * behaviour until the resulting parameters combinations lie within
 * a range that should meet accuracy and stability requirements.
 * Terminates program execution if no suitable parameters are found.
 * 
 * @param dn A DimensionlessNumbers instance containing the 
 *        dimensionless numbers.
 * @param grid A GridParameters instance containing the grid parameters.
 * @return A set of adequate parameter values 
 *         as a SimulationParameters instance.
 */
SimulationParameters goodParameters(DimensionlessNumbers dn, GridParameters grid)
{
    // Set maximum and minimum allowed relaxation time values.
    const double tau_max = 1.1;    // Much higher than this -> lose accuracy.
    const double tau_min = 0.5125; // Lower than this -> unstable.
    double nu, kappa, tau_f, tau_g;
    double alpha_g = 0.1 * 0.1 / static_cast<double>(grid.ny);
    double ny3 = pow(static_cast<double>(grid.ny), 3.0);
    int max_attempts = 1000; // Maximum number of iterations to find good parameter values.
    for (int it = 0; it < max_attempts; it++)
    {
        nu = sqrt(alpha_g * ny3 * dn.prandtl / dn.rayleigh); // Kinematic viscosity in lattice units.
        kappa = nu / dn.prandtl; // Thermal diffusivity in lattice units.
        tau_f = (0.5 + 3*nu); // Fluid relaxation time in lattice units.
        tau_g = (0.5 + 3*kappa); // Thermal relaxation time in lattice units.
        if (tau_f > tau_max || tau_g > tau_max)
        {
            // If relaxation time too large (low accuracy solution),
            // reduce gravity to reduce nu, kappa, and hence tau_f, tau_g.
            // (Disadvantage: more time steps required).
            alpha_g *= 0.999;
        }
        else if (tau_f < tau_min || tau_g < tau_min)
        {
            // If relaxation time too small (unstable),
            // increase gravity to increase nu, kappa, and hence tau_f, tau_g.
            // (Disadvantage: higher Mach number, larger compressibility errors, may also lead to instability).
            alpha_g *= 1.0001;
        }
        else
        {
            // Convert to relaxation rates, then return as SP instance.
            double omega_f = 1 / tau_f;
            double omega_g = 1 / tau_g;
            SimulationParameters sp(alpha_g, omega_f, omega_g);
            return sp;
        }
    }
    std::cerr << " !! Unable to find adequate parameter values.\n";
    std::cerr << " !! Terminating program :(\n\n";
    exit(EXIT_FAILURE); // Terminate the program.
}

/*****************************************************************************
 *                                                                           *
 *                      Templated function definitions                       *
 *                                                                           *
 ****************************************************************************/

/**
 * @brief Initialise distribution function arrays for Rayleigh Benard test.
 * @details Calculate the equilibrium deviaiton distribution functions
 *          corresponding to the macroscopic arrays.
 *
 * @param f The D2Q9-based fluid velocity deviation distribution functions.
 * @param g The D2Q5-based temperature deviation distribution functions.
 * @param macros An instance of the Macroscopic2D class.
 * @param grid An instance of the GridParameters class. 
 */
template <class T>
void initialiseDistributions(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, GridParameters grid)
{
    // Declare temporary variables.
    T r_, u_, v_, t_, dfeq, eu, u2, w;
    T const densAnom = static_cast<T>(0.0);

    // Loop over domain.
    for (int j = 0; j < grid.ny; j++)
    {
        for (int i = 0; i < grid.nx; i++)
        {
            // Grab data from macroscopic arrays.
            r_ = macros.GetDensity(i, j);
            u_ = macros.GetVelocityX(i, j);
            v_ = macros.GetVelocityY(i, j);
            t_ = macros.GetTemperature(i, j);
            u2 = u_ * u_ + v_ * v_;

            // D2Q9 velocity distribution functions.
            // (0, 0)
            eu = static_cast<T>(0.0);
            w = static_cast<T>(4.0/9.0);
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, eu, u2, w);
            f.SetF0(dfeq, i, j);

            // (1, 0)
            w = static_cast<T>(1.0/9.0);
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, u_, u2, w);
            f.SetF1(dfeq, i, j);

            // (-1, 0)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, -u_, u2, w);
            f.SetF2(dfeq, i, j);

            // (0, 1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, v_, u2, w);
            f.SetF3(dfeq, i, j);

            // (0, -1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, -v_, u2, w);
            f.SetF4(dfeq, i, j);

            // (1, 1)
            w = static_cast<T>(1.0/36.0);
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, u_ + v_, u2, w);
            f.SetF5(dfeq, i, j);

            // (-1, -1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, -u_ - v_, u2, w);
            f.SetF6(dfeq, i, j);

            // (1, -1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, u_ - v_, u2, w);
            f.SetF7(dfeq, i, j);

            // (-1, 1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, -u_ + v_, u2, w);
            f.SetF8(dfeq, i, j);

            // D2Q5 temperature distributions.
            // (0, 0)
            w = static_cast<T>(1.0/3.0);
            dfeq = computeSecondOrderDeviationEquilibrium(t_, densAnom, eu, u2, w);
            g.SetF0(dfeq, i, j);

            // (1, 0)
            w = static_cast<T>(1.0/6.0);
            dfeq = computeSecondOrderDeviationEquilibrium(t_, densAnom, u_, u2, w);
            g.SetF1(dfeq, i, j);

            // (-1, 0)
            dfeq = computeSecondOrderDeviationEquilibrium(t_, densAnom, -u_, u2, w);
            g.SetF2(dfeq, i, j);

            // (0, 1)
            dfeq = computeSecondOrderDeviationEquilibrium(t_, densAnom, v_, u2, w);
            g.SetF3(dfeq, i, j);

            // (0, -1)
            dfeq = computeSecondOrderDeviationEquilibrium(t_, densAnom, -v_, u2, w);
            g.SetF4(dfeq, i, j);
        }
    }
}

/**
 * @brief Initialise macroscopic arrays for Rayleigh Benard test.
 * @details Constant density, quiescent (zero velocity), linear temperature gradient
 *          from 1 at bottom plate to 0 at top plate.
 *          Apply a small sinusoidal temperature perturbation at the middle.
 * 
 * @param macros An instance of Macroscopic2D class.
 * @param grid An instance of the GridParameters class. 
 */
template<class T>
void initialiseMacroscopic(Macroscopic2D<T>& macros, GridParameters grid)
{
    // Declare temporary variables.
    T r_ = 1, u_ = 0, v_ = 0, t_;

    // Loop over domain.
    for (int j = 0; j < grid.ny; j++)
    {
        for (int i = 0; i < grid.nx; i++)
        {
            // Linear temperature profile. 1 at bottom, 0 at top.
            t_ = 1.0 - (static_cast<T>(j) + 0.5) / static_cast<T>(grid.ny);

            // Perturb the initial temperature profile at the middle of domain.
            T perturbation_size = 1e-2;
            if (j == grid.ny / 2)
            {
                T argument = static_cast<T>(i) * 2 * M_PI / static_cast<T>(grid.nx);
                t_ += perturbation_size * sin(argument);
            }

            // Store values in Macroscopic2D arrays.
            macros.SetDensity(r_, i, j);
            macros.SetVelocityX(u_, i, j);
            macros.SetVelocityY(v_, i, j);
            macros.SetTemperature(t_, i, j);
        }
    }
    std::cout << "Initialised macroscopic arrays.\n";
}

/**
 * @brief Runs one time step of the Lattice Boltzmann algorithm.
 *
 * @param f The D2Q9-based fluid velocity deviation distribution functions.
 * @param g The D2Q5-based temperature deviation distribution functions.
 * @param macros An instance of the Macroscopic2D class.
 * @param grid An instance of the GridParameters class. 
 */
template<class T>
void runTimestep(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, GridParameters grid)
{
    // "Stream" by swapping pointers.
    f.SwapPointers();
    g.SwapPointers();

    // Loop over domain.
    for (int j = 0; j < grid.ny; ++j)
    {
        for (int i = 0; i < grid.nx; ++i)
        {
            // Grab data from arrays to registers.
            T f0 = f.GetF0(i, j);
            T f1 = f.GetF1(i, j);
            T f2 = f.GetF2(i, j);
            T f3 = f.GetF3(i, j);
            T f4 = f.GetF4(i, j);
            T f5 = f.GetF5(i, j);
            T f6 = f.GetF6(i, j);
            T f7 = f.GetF7(i, j);
            T f8 = f.GetF8(i, j);

            T g0 = g.GetF0(i, j);
            T g1 = g.GetF1(i, j);
            T g2 = g.GetF2(i, j);
            T g3 = g.GetF3(i, j);
            T g4 = g.GetF4(i, j);

            // Boundary conditions.
            // Bottom boundary. u=0, T=1.
            if (j == 0)
            {
                f3 = f.GetBouncedF3(i, j);
                f5 = f.GetBouncedF5(i, j);
                f7 = f.GetBouncedF7(i, j);
            }
        }
    }
}