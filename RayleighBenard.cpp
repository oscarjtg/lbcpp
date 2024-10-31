#include <iostream> // For cout.
#include <cstdlib>  // For atof and exit.
#include <cmath>    // For pow and sqrt.
#include <mpi.h>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>

#include "equilibria/Equilibria.h"
#include "grid/GridParameters.h"
#include "lattice/LatticeD2Q5.h"
#include "lattice/LatticeD2Q9.h"
#include "macroscopic/Macroscopic2D.h"
#include "RayleighBenardStructs.h"

// Functions.
CharacteristicScales computeCharScales(DimensionlessNumbers dn, SimulationParameters sp, GridParameters grid);
void displayParameters(GridParameters grid, DimensionlessNumbers dn, SimulationParameters sp, CharacteristicScales scales);
DimensionlessNumbers getInputParameters(int argc, char* argv[]);
SimulationParameters goodParameters(DimensionlessNumbers dn, GridParameters grid);
//GridParameters setGridParameters(int process_number, int number_of_processes, int nx_full, int ny_full, int nx_sub_max, int ny_sub_max);

// Templated non-inline functions.
template <class T>
void initialiseDistributions(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, GridParameters grid);

template<class T>
void initialiseMacroscopic(Macroscopic2D<T>& macros, GridParameters grid);

template<class T>
void runTimestep(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, T* buffers[], GridParameters grid, SimulationParameters sp, bool saveMacros);

// Templated inline functions.
template<class T>
inline T bodyforce(T temperature, SimulationParameters sp) 
{
    return static_cast<T>(sp.ag) * temperature;
}

template<class T>
inline T calculateZerothMoment_D2Q9(T f0, T f1, T f2, T f3, T f4, T f5, T f6, T f7, T f8)
{
    return f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
}

template<class T>
inline T calculateFirstMoment_D2Q9(T a, T b, T c, T d, T e, T f)
{
    return a - b + c - d + e - f;
}

template<class T>
inline T calculateZerothMoment_D2Q5(T a, T b, T c, T d, T e)
{
    return a + b + c + d + e;
}

template<class T>
inline T addOne(T x)
{
    return static_cast<T>(1) + x;
}

template<class T>
inline T calculateVelocity(const T first_moment, const T density)
{
    return first_moment / density;
}

template<class T>
inline T calculateVelocity(const T first_moment, const T density, const T gravity)
{
    return (first_moment / density) + static_cast<T>(0.5) * gravity;
}

template<class T>
inline T computeForceTerm(const T forceProjection, const T velocityProjection, const T forceDotVelocity, const T latticeWeight)
{
    return latticeWeight * (static_cast<T>(3) * forceProjection + static_cast<T>(9) * forceProjection * velocityProjection - static_cast<T>(3) * forceDotVelocity);
}

template<class T>
inline T computeSRTCollision(T df, const T dfeq, const T omeg)
{
    return omeg * dfeq + (static_cast<T>(1) - omeg) * df;
}

template<class T>
inline T computeSRTCollision(T df, const T dfeq, const T F, const T omeg)
{
    return omeg * dfeq + (static_cast<T>(1) - omeg) * df + (static_cast<T>(1) - static_cast<T>(0.5)*omeg) * F;
}

/*****************************************************************************
 *                                                                           *
 *                                   Main                                    *
 *                                                                           *
 ****************************************************************************/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int process_number, number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_number);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    //std::ostringstream oss;
    //oss << "00000_" << process_number << "_" << number_of_processes;
    //std::string runId = oss.str();
    std::string runId = "00000";

    // Define constants.
    const int NX = 90;
    const int NY = 90;
    const int NT = 10;
    const int NROWS = 3;
    int NCOLUMNS = number_of_processes / NROWS;

    //assert(number_of_processes % NROWS == 0);
    //assert(NX % NCOLUMNS == 0);
    //assert(NY % NROWS == 0);

    if (number_of_processes % NROWS != 0 || NX % NCOLUMNS != 0 || NY % NROWS != 0)
    {
        std::cout << " !! Error, incompatible number of MPI processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    
    // nx, ny, size
    GridParameters grid(NX/NCOLUMNS, NY/NROWS, NCOLUMNS, NROWS, process_number%NCOLUMNS, process_number/NCOLUMNS, process_number/NCOLUMNS==0, process_number/NCOLUMNS==NROWS-1, process_number%NCOLUMNS==0, process_number%NCOLUMNS==NCOLUMNS-1);
    
    int rigt_proc_num = (process_number + 1 + number_of_processes) % number_of_processes;
    int left_proc_num = (process_number - 1 + number_of_processes) % number_of_processes;
    int upup_proc_num = (process_number + grid.nSubDomsX + number_of_processes) % number_of_processes;
    int down_proc_num = (process_number - grid.nSubDomsX + number_of_processes) % number_of_processes;

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
    Macroscopic2D<float> macros(grid.nx, grid.ny);
    LatticeD2Q9<float> f(grid.nx, grid.ny);
    LatticeD2Q5<float> g(grid.nx, grid.ny);

    // Allocate send/receive buffers for MPI.
    float* bsend_rigt = new float[4 * grid.ny]; // send to the process to the right of this one.
    float* brecv_left = new float[4 * grid.ny]; // receive from the process to the left of this one.
    float* bsend_left = new float[4 * grid.ny]; // etc...
    float* brecv_rigt = new float[4 * grid.ny];
    float* bsend_upup = new float[4 * grid.nx]; // note change from grid.ny to grid.nx.
    float* brecv_down = new float[4 * grid.nx];
    float* bsend_down = new float[4 * grid.nx];
    float* brecv_upup = new float[4 * grid.nx];

    // Define an array of pointers to the first elements of each array
    float* buffers[] = {
        bsend_rigt,
        brecv_left,
        bsend_left,
        brecv_rigt,
        bsend_upup,
        brecv_down,
        bsend_down,
        brecv_upup
    };

    // Fill macroscopic arrays with initial conditions.
    initialiseMacroscopic(macros, grid);

    // Fill distribution function arrays with initial conditions.
    initialiseDistributions(f, g, macros, grid);

    // Save data to CSV.
    macros.WriteToCSV("output", runId, 0, process_number);
    f.WriteToCSV("output", 'f', runId, 0, process_number);
    g.WriteToCSV("output", 'g', runId, 0, process_number);

    MPI_Barrier(MPI_COMM_WORLD);
    
    
    // The main loop.
    for (int t=0; t<NT; ++t)
    {
        bool save_data = t % 1 == 0;

        // Time step.
        runTimestep(f, g, macros, buffers, grid, sp, save_data);

        // Send/recv data.
        int tag = 0;
        MPI_Sendrecv(bsend_rigt, 4*grid.ny, MPI_FLOAT, rigt_proc_num, tag, brecv_left, 4*grid.ny, MPI_FLOAT, left_proc_num, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        tag = 1;
        MPI_Sendrecv(bsend_left, 4*grid.ny, MPI_FLOAT, left_proc_num, tag, brecv_rigt, 4*grid.ny, MPI_FLOAT, rigt_proc_num, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        tag = 2;
        MPI_Sendrecv(bsend_upup, 4*grid.nx, MPI_FLOAT, upup_proc_num, tag, brecv_down, 4*grid.nx, MPI_FLOAT, down_proc_num, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        tag = 3;
        MPI_Sendrecv(bsend_down, 4*grid.nx, MPI_FLOAT, down_proc_num, tag, brecv_upup, 4*grid.nx, MPI_FLOAT, upup_proc_num, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Save data periodically.
        if (save_data)
        {
            //f.WriteToCSV("output", 'f', runId, t, process_number);
            //g.WriteToCSV("output", 'g', runId, t, process_number);
            macros.WriteToCSV("output", runId, t, process_number);
        }
    }
    
    // Freeing the allocated memory.
    delete[] bsend_rigt;
    delete[] brecv_left;
    delete[] bsend_left;
    delete[] brecv_rigt;
    delete[] bsend_upup;
    delete[] brecv_down;
    delete[] bsend_down;
    delete[] brecv_upup;

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
    float ny2 = pow(static_cast<float>(grid.ny*grid.nSubDomsY), 2.0); // NY^2.
    float nu = ((1.0 / sp.omega_f) - 0.5) / 3.0; // Kinematic viscosity.
    float u = sqrt(sp.ag * grid.ny * grid.nSubDomsY * 1.0); // delta T = 1.0.
    float t = (ny2 / nu) * sqrt(dn.prandtl / dn.rayleigh);
    float h = u * t;
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
    std::cout << "    subDomCoordX     = " << grid.subDomCoordX << " / " << grid.nSubDomsX << "\n";
    std::cout << "    subDomCoordY     = " << grid.subDomCoordY << " / " << grid.nSubDomsY << "\n";
    std::cout << "    isBottom         = " << grid.isBottom << "\n";
    std::cout << "    isTop            = " << grid.isTop << "\n";
    std::cout << "    isLeft           = " << grid.isLeft << "\n";
    std::cout << "    isRight          = " << grid.isRight << "\n";
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
    const float tau_max = 1.1;    // Much higher than this -> lose accuracy.
    const float tau_min = 0.5125; // Lower than this -> unstable.
    float nu, kappa, tau_f, tau_g;
    float alpha_g = 0.1 * 0.1 / static_cast<float>(grid.ny);
    float ny3 = pow(static_cast<float>(grid.ny), 3.0);
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
            float omega_f = 1 / tau_f;
            float omega_g = 1 / tau_g;
            SimulationParameters sp(alpha_g, omega_f, omega_g);
            return sp;
        }
    }
    std::cerr << " !! Unable to find adequate parameter values.\n";
    std::cerr << " !! Terminating program :(\n\n";
    exit(EXIT_FAILURE); // Terminate the program.
}

/**
 * @brief Sets appropriate grid sizes and positions for subdomain.
 * 
 * @param process_number The id of this process. Should be less than 
 *                       number_of_processes!
 * @param number_of_processes The total number of processes.
 * @param nx_full The number of grid points in the x direction for the 
 *                whole domain.
 * @param ny_full The number of grid points in the y directions for the
 *                whole domain.
 * @param nx_sub_max The maximum number of x grid points of a subdomain.
 * @param ny_sub_max The maximum number of y grid points of a subdomain.
 * @return A GridParameters instance.
 */
/*
GridParameters setGridParameters(int process_number, int number_of_processes, int nx_full, int ny_full, int nx_sub_max, int ny_sub_max)
{
    // If we had unlimited processors:
    int nSubDomsX = nx_full / nx_sub_max;
    int nSubDomsY = ny_full / ny_sub_max;
    int number_of_subdomains = nSubDomsX * nSubDomsY;
    if (number_of_processes >= number_of_subdomains)
    {
        
    }
}
*/

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

            // (-1, 1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, -u_ + v_, u2, w);
            f.SetF7(dfeq, i, j);

            // (1, -1)
            dfeq = computeSecondOrderDeviationEquilibrium(r_, densAnom, u_ - v_, u2, w);
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
            t_ = 1.0 - (static_cast<T>(j) + grid.ny*grid.subDomCoordY + 0.5) / static_cast<T>(grid.ny*grid.nSubDomsY);

            // Perturb the initial temperature profile at the middle of domain.
            T perturbation_size = 1e-2;
            if (grid.subDomCoordY == 1 && j == grid.ny / 2)
            {
                T argument = static_cast<T>(i) * 2 * M_PI / static_cast<T>(grid.nx*grid.nSubDomsX);
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
void runTimestep(LatticeD2Q9<T>& f, LatticeD2Q5<T>& g, Macroscopic2D<T>& macros, T* buffers[], GridParameters grid, SimulationParameters sp, bool saveMacros)
{
    // "Stream" by swapping pointers.
    f.SwapPointers();
    g.SwapPointers();

    // Write in data streamed from neighbouring processes.
    // Left hand side. Receive f1, f5, f8, g1 from buffer1
    for (int j = 0; j < grid.ny; ++j)
    {
        T f1 = buffers[1][4*j];
        f.SetInitF1(f1, 0, j);
        T f5 = buffers[1][4*j+1];
        f.SetInitF5(f5, 0, j);
        T f8 = buffers[1][4*j+2];
        f.SetInitF8(f8, 0, j);
        T g1 = buffers[1][4*j+3];
        g.SetInitF1(g1, 0, j);
    }

    // Right hand side. Receive f2, f6, f7, g2 from buffer3
    for (int j = 0; j < grid.ny; ++j)
    {
        const int iright = grid.nx - 1;
        T f2 = buffers[3][4*j];
        f.SetInitF2(f2, iright, j);
        T f6 = buffers[3][4*j+1];
        f.SetInitF6(f6, iright, j);
        T f7 = buffers[3][4*j+2];
        f.SetInitF7(f7, iright, j);
        T g2 = buffers[3][4*j+3];
        g.SetInitF2(g2, iright, j);
    }

    // Bottom side. Receive f3, f5, f7, g3 from buffer5
    for (int i = 0; i < grid.nx; ++i)
    {
        T f3 = buffers[5][4*i];
        f.SetInitF3(f3, i, 0);
        T f5 = buffers[5][4*i+1];
        f.SetInitF5(f5, i, 0);
        T f7 = buffers[5][4*i+2];
        f.SetInitF7(f7, i, 0);
        T g3 = buffers[5][4*i+3];
        g.SetInitF3(g3, i, 0);
    }

    // Top side. Receive f4, f6, f8, g4 from buffer7
    for (int i = 0; i < grid.nx; ++i)
    {   
        T jtop = grid.ny - 1;
        T f4 = buffers[7][4*i];
        f.SetInitF4(f4, i, jtop);
        T f6 = buffers[7][4*i+1];
        f.SetInitF6(f6, i, jtop);
        T f8 = buffers[7][4*i+2];
        f.SetInitF8(f8, i, jtop);
        T g4 = buffers[7][4*i+3];
        g.SetInitF4(g4, i, 0);
    }

    // Loop over domain.
    for (int j = 0; j < grid.ny; ++j)
    {
        for (int i = 0; i < grid.nx; ++i)
        {
            // Grab data from arrays to registers.
            T df0 = f.GetF0(i, j);
            T df1 = f.GetF1(i, j);
            T df2 = f.GetF2(i, j);
            T df3 = f.GetF3(i, j);
            T df4 = f.GetF4(i, j);
            T df5 = f.GetF5(i, j);
            T df6 = f.GetF6(i, j);
            T df7 = f.GetF7(i, j);
            T df8 = f.GetF8(i, j);

            T dg0 = g.GetF0(i, j);
            T dg1 = g.GetF1(i, j);
            T dg2 = g.GetF2(i, j);
            T dg3 = g.GetF3(i, j);
            T dg4 = g.GetF4(i, j);

            // Boundary conditions.
            // Bottom boundary. u=0, T=1.
            if (grid.isBottom && j == 0)
            {
                df3 = f.GetBouncedF3(i, j);
                df5 = f.GetBouncedF5(i, j);
                df7 = f.GetBouncedF7(i, j);
                dg3 = -g.GetBouncedF3(i, j) + static_cast<T>(1.0/3.0); // 2*w_i
            }

            // Top boundary. u=0, T=0;
            if (grid.isTop && j == grid.ny-1)
            {
                df3 = f.GetBouncedF3(i, j);
                df5 = f.GetBouncedF5(i, j);
                df7 = f.GetBouncedF7(i, j);
                dg3 = -g.GetBouncedF3(i, j);   
            }

            // Do collision.
            // Assign to temp variables f0p, etc...
            // Calculate moments.
            T dfM00 = calculateZerothMoment_D2Q9(df0, df1, df2, df3, df4, df5, df6, df7, df8);
            T dfM10 = calculateFirstMoment_D2Q9(df1, df2, df5, df6, df8, df7);
            T dfM01 = calculateFirstMoment_D2Q9(df3, df4, df5, df6, df7, df8);
            T dgM00 = calculateZerothMoment_D2Q5(dg0, dg1, dg2, dg3, dg4);

            // Calculate macroscopic properties.
            T t_ = addOne(dgM00); // temperature.
            T gy_ = bodyforce(t_, sp); // gravitational acceleration.
            T r_ = addOne(dfM00); // density.
            T u_ = calculateVelocity(dfM10, r_); // velocity x-component.
            T v_ = calculateVelocity(dfM01, r_, gy_); // velocity y-component.
            T u2 = u_*u_ + v_*v_;

            if (saveMacros)
            {
                macros.SetDensity(r_, i, j);
                macros.SetVelocityX(u_, i, j);
                macros.SetVelocityY(v_, i, j);
                macros.SetTemperature(t_, i, j);
            }

            // Calculate equilibria.
            // D2Q9 velocity deviation DFs
            // (0, 0)
            T w = static_cast<T>(4.0 / 9.0);
            T zero_velocity = static_cast<T>(0.0);
            T dfeq0 = computeSecondOrderDeviationEquilibrium(r_, dfM00, zero_velocity, u2, w);

            // (1, 0)
            w = static_cast<T>(1.0 / 9.0);
            T dfeq1 = computeSecondOrderDeviationEquilibrium(r_, dfM00, u_, u2, w);

            // (-1, 0)
            T dfeq2 = computeSecondOrderDeviationEquilibrium(r_, dfM00, -u_, u2, w);

            // (0, 1)
            T dfeq3 = computeSecondOrderDeviationEquilibrium(r_, dfM00, v_, u2, w);

            // (0, -1)
            T dfeq4 = computeSecondOrderDeviationEquilibrium(r_, dfM00, -v_, u2, w);

            // (1, 1)
            w = static_cast<T>(1.0/36.0);
            T dfeq5 = computeSecondOrderDeviationEquilibrium(r_, dfM00, u_ + v_, u2, w);

            // (-1, -1)
            T dfeq6 = computeSecondOrderDeviationEquilibrium(r_, dfM00, -u_ - v_, u2, w);

            // (-1, 1)
            T dfeq7 = computeSecondOrderDeviationEquilibrium(r_, dfM00, -u_ + v_, u2, w);

            // (1, -1)
            T dfeq8 = computeSecondOrderDeviationEquilibrium(r_, dfM00, u_ - v_, u2, w);

            // D2Q5 temperature distributions.
            // (0, 0)
            w = static_cast<T>(1.0/3.0);
            T dgeq0 = computeSecondOrderDeviationEquilibrium(t_, dgM00, zero_velocity, u2, w);

            // (1, 0)
            w = static_cast<T>(1.0/6.0);
            T dgeq1 = computeSecondOrderDeviationEquilibrium(t_, dgM00, u_, u2, w);

            // (-1, 0)
            T dgeq2 = computeSecondOrderDeviationEquilibrium(t_, dgM00, -u_, u2, w);

            // (0, 1)
            T dgeq3 = computeSecondOrderDeviationEquilibrium(t_, dgM00, v_, u2, w);

            // (0, -1)
            T dgeq4 = computeSecondOrderDeviationEquilibrium(t_, dgM00, -v_, u2, w);

            // Calculate forcing terms (velocity DFs only).
            T F_dot_u = gy_ * v_;
            T zero_force = static_cast<T>(0);

            // (0, 0)
            w = static_cast<T>(4.0 / 9.0);
            T F0 = computeForceTerm(zero_force, zero_velocity, F_dot_u, w);

            // (1, 0)
            w = static_cast<T>(1.0 / 9.0);
            T F1 = computeForceTerm(zero_force, u_, F_dot_u, w);

            // (-1, 0)
            T F2 = computeForceTerm(zero_force, -u_, F_dot_u, w);

            // (0, 1)
            T F3 = computeForceTerm(gy_, v_, F_dot_u, w);

            // (0, -1)
            T F4 = computeForceTerm(-gy_, -v_, F_dot_u, w);

            // (1, 1)
            w = static_cast<T>(1.0 / 36.0);
            T F5 = computeForceTerm(gy_, u_ + v_, F_dot_u, w);

            // (-1, -1)
            T F6 = computeForceTerm(-gy_, -u_ - v_, F_dot_u, w);

            // (-1, 1)
            T F7 = computeForceTerm(gy_, -u_ + v_, F_dot_u, w);

            // (1, -1)
            T F8 = computeForceTerm(-gy_, u_ - v_, F_dot_u, w);

            // Collide (SRT) and store.
            df0 = computeSRTCollision(df0, dfeq0, F0, sp.omega_f);
            f.SetF0(df0, i, j);
            df1 = computeSRTCollision(df1, dfeq1, F1, sp.omega_f);
            f.SetF1(df1, i, j);
            df2 = computeSRTCollision(df2, dfeq2, F2, sp.omega_f);
            f.SetF2(df2, i, j);
            df3 = computeSRTCollision(df3, dfeq3, F3, sp.omega_f);
            f.SetF3(df3, i, j);
            df4 = computeSRTCollision(df4, dfeq4, F4, sp.omega_f);
            f.SetF4(df4, i, j);
            df5 = computeSRTCollision(df5, dfeq5, F5, sp.omega_f);
            f.SetF5(df5, i, j);
            df6 = computeSRTCollision(df6, dfeq6, F6, sp.omega_f);
            f.SetF6(df6, i, j);
            df7 = computeSRTCollision(df7, dfeq7, F7, sp.omega_f);
            f.SetF7(df7, i, j);
            df8 = computeSRTCollision(df8, dfeq8, F8, sp.omega_f);
            f.SetF8(df8, i, j);
            dg0 = computeSRTCollision(dg0, dgeq0, sp.omega_g);
            g.SetF0(dg0, i, j);
            dg1 = computeSRTCollision(dg1, dgeq1, sp.omega_g);
            g.SetF1(dg1, i, j);
            dg2 = computeSRTCollision(dg2, dgeq2, sp.omega_g);
            g.SetF2(dg2, i, j);
            dg3 = computeSRTCollision(dg3, dgeq3, sp.omega_g);
            g.SetF3(dg3, i, j);
            dg4 = computeSRTCollision(dg4, dgeq4, sp.omega_g);
            g.SetF4(dg4, i, j);

            // Write to buffers.
            // Streaming out of right side (f1, f5, f8, g1, buffer0):
            if (i == grid.nx-1)
            {
                buffers[0][4*j] = df1;
                buffers[0][4*j+1] = df5;
                buffers[0][4*j+2] = df8;
                buffers[0][4*j+3] = dg1;
            }
            // Streaming out of left side (f2, f6, f7, g2, buffer2):
            else if (i == 0)
            {
                buffers[2][4*j] = df2;
                buffers[2][4*j+1] = df6;
                buffers[2][4*j+2] = df7;
                buffers[2][4*j+3] = dg2;
            }
            // Streaming out of top side (f3, f5, f7, g3, buffer4):
            if (j == grid.ny-1)
            {
                buffers[4][4*i] = df3;
                buffers[4][4*i+1] = df5;
                buffers[4][4*i+2] = df7;
                buffers[4][4*i+3] = dg3;
            }
            // Streaming out of bottom side (f4, f6, f8, g4, buffer6):
            if (j == 0)
            {
                buffers[6][4*i] = df4;
                buffers[6][4*i+1] = df6;
                buffers[6][4*i+2] = df8;
                buffers[6][4*i+3] = dg4;
            }
        }
    }
}