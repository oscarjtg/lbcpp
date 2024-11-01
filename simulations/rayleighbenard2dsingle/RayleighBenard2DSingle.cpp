#include <cmath> // for pow and sqrt
#include <cstdlib> // for atof and exit
#include <iostream> // for cout

#include "domain/BoundaryInfo.h"
#include "domain/NodeInfo.h"
#include "evolver/ScalarEvolverSRT.h"
#include "lattice/LatticeSoAPull.h"
#include "macroscopic/MacroscopicVariable.h"

void getInputParameters(double& rayleigh_number, double& prandtl_number, int& nx, int& ny, int& nt, int argc, char* argv[]);

void findGoodModelParameters(double& ag, double& tau_f, double& tau_g, const double Pr, const double Ra, const int ny);

/***************************************************
 *                                                 *
 *                     Main                        *
 *                                                 *
 ***************************************************/

int main(int argc, char* argv[])
{
    double Ra, Pr; // Rayleigh and Prandtl numbers.
    int nx, ny, nz = 1, nt; // Number of grid points and time steps.
    getInputParameters(Ra, Pr, nx, ny, nt, argc, argv);

    double ag, tau_f, tau_g; // Model parameters.
    findGoodModelParameters(ag, tau_f, tau_g, Pr, Ra, ny);

    double kinematic_viscosity, thermal_diffusivity;
    kinematic_viscosity = (tau_f - 0.5) / 3.0;
    thermal_diffusivity = (tau_g - 0.5) / 3.0;

    double length_scale, velocity_scale, time_scale; //Characteristic scales.
    velocity_scale = sqrt(ag * ny * 1.0); // Delta T = 1.0
    time_scale = (ny*ny / kinematic_viscosity) * sqrt(Pr / Ra);
    length_scale = velocity_scale * time_scale;

    // Information about this run.
    std::string run_id = "testrun", save_path = "output";

    // Initialise arrays.
    const int nd = 2, nq_f = 9, nq_g = 5; // Set LB lattice model type.
    LatticeSoAPull<double> f(nx, ny, nz, nd, nq_f, "f", run_id, save_path);
    LatticeSoAPull<double> g(nx, ny, nz, nd, nq_g, "g", run_id, save_path);
    MacroscopicVariable<double> dens(nx, ny, nz, "r", run_id, save_path);
    MacroscopicVariable<double> velx(nx, ny, nz, "u", run_id, save_path);
    MacroscopicVariable<double> vely(nx, ny, nz, "v", run_id, save_path);
    MacroscopicVariable<double> velz(nx, ny, nz, "w", run_id, save_path);
    MacroscopicVariable<double> temp(nx, ny, nz, "t", run_id, save_path);
    MacroscopicVariable<double> Fx(nx, ny, nz, "Fx", run_id, save_path);
    MacroscopicVariable<double> Fy(nx, ny, nz, "Fy", run_id, save_path);
    MacroscopicVariable<double> Fz(nx, ny, nz, "Fz", run_id, save_path);

    // Print info to terminal.
    f.DisplayLatticeParameters();
    g.DisplayLatticeParameters();
    dens.DisplayInfo();
    velx.DisplayInfo();
    vely.DisplayInfo();
    temp.DisplayInfo();

    // Try swapping pointers.
    f.StreamDistributions();
    g.StreamDistributions();

    // Save data.
    std::cout << "Trying to save data\n";
    f.WriteToTextFile();
    dens.WriteToTextFile();

    std::cout << length_scale << thermal_diffusivity << "\n";

    FluidEvolverSRT<double> fluid_evolver;
    ScalarEvolverSRT<double> scalar_evolver;

    BoundaryInfo bdry(2);
    NodeInfo node(nx, ny, nz);

    // Run algorithm.
    for (int t = 0; t < nt; ++t)
    {
        // Perform one timestep of ADE LB algorithm for g using velocity u.
        // This updates g and temp.
        scalar_evolver.DoTimestep(g, temp, velx, vely, velz, node, bdry);

        // Compute force density F from updated temp.
        // TODO

        // Perform one timestep of the standard LB algorithm for f using F.
        // This updates f, dens, velx, vely, and velz.
        fluid_evolver.DoTimestep(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry);

        // Compute diagnostic.
    }
}

/***************************************************
 *                                                 *
 *           Helper function definitions           *
 *                                                 *
 ***************************************************/

void getInputParameters(double& rayleigh_number, double& prandtl_number, int& nx, int& ny, int& nt, int argc, char* argv[])
{
    // Default values.
    rayleigh_number = 1.0e4;
    prandtl_number = 1.0;
    nx = 100;
    ny = 100;
    nt = 100;

    if (argc >= 2)
    {
        rayleigh_number = atof(argv[1]);
    }
    if (argc >= 3)
    {
        prandtl_number = atof(argv[2]);
    }
    if (argc >= 4)
    {
        nx = atoi(argv[3]);
    }
    if (argc >= 5)
    {
        ny = atoi(argv[4]);
    }
    if (argc >= 6)
    {
        ny = atoi(argv[5]);
    }
    return;
}

void findGoodModelParameters(double& alpha_g, double& tau_f, double& tau_g, const double Pr, const double Ra, const int ny)
{
    // Set maximum and minimum allowed relaxation time values.
    const double tau_max = 1.1;    // Much higher than this -> lose accuracy.
    const double tau_min = 0.5125; // Lower than this -> unstable.
    double nu, kappa;
    alpha_g = 0.1 * 0.1 / static_cast<double>(ny);
    double ny3 = pow(static_cast<double>(ny), 3.0);
    int max_attempts = 1000; // Maximum number of iterations to find good parameter values.
    for (int it = 0; it < max_attempts; it++)
    {
        nu = sqrt(alpha_g * ny3 * Pr / Ra); // Kinematic viscosity in lattice units.
        kappa = nu / Pr; // Thermal diffusivity in lattice units.
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
            return;
        }
    }
    std::cerr << " !! Unable to find adequate parameter values.\n";
    std::cerr << " !! Terminating program :(\n\n";
    exit(EXIT_FAILURE); // Terminate the program.
}