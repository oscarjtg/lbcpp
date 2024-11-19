#include <algorithm> // for max
#include <cmath> // for pow, sqrt, and fabs
#include <cstdlib> // for atof and exit
#include <iostream> // for cout
#include <string>

#include "boundary/AllBoundaryRules.h"
#include "convergence/ConvergenceTester.h"
#include "diagnostics/DiagnosticNusselt.h"
#include "domain/BoundaryInfo.h"
#include "domain/NodeInfo.h"
#include "evolver/AllEvolvers.h"
#include "force/AllForces.h"
#include "lattice/LatticeSoAPull.h"
#include "macroscopic/MacroscopicVariable.h"

void getInputParameters(double& rayleigh_number, double& prandtl_number, int& nx, int& nz, int& nt, std::string& run_id, int argc, char* argv[]);

void findGoodModelParameters(double& ag, double& tau_f, double& tau_g, const double Pr, const double Ra, const int ny, const double csi_f, const double csi_g);

template <class T>
void printAverages(std::string& message, MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
                   MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz, MacroscopicVariable<T>& temp);

template <class T, int ND, int NQ_F, int NQ_G>
void saveData(int timstep, std::string& message, AbstractLattice<T, ND, NQ_F>& f, AbstractLattice<T, ND, NQ_G>& g, 
             MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
             MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz, MacroscopicVariable<T>& temp);
                   
template <class T>
void saveMacroscopic(int timstep, std::string& message, 
             MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
             MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz, MacroscopicVariable<T>& temp);
/***************************************************
 *                                                 *
 *                     Main                        *
 *                                                 *
 ***************************************************/

int main(int argc, char* argv[])
{
    using dist_type = double; // type alias for distribution functions

    double Ra, Pr; // Rayleigh and Prandtl numbers.
    int nx, ny=1, nz, nt; // Number of grid points and time steps.
    std::string run_id;
    getInputParameters(Ra, Pr, nx, nz, nt, run_id, argc, argv);

    // Information about this run.
    std::string save_path = "output";

    // Initialise arrays.
    const int nd = 2, nq_f = 9, nq_g = 5; // Set LB lattice model type.
    LatticeSoAPull<dist_type, nd, nq_f> f(nx, ny, nz, "f", run_id, save_path);
    LatticeSoAPull<dist_type, nd, nq_g> g(nx, ny, nz, "g", run_id, save_path);

    MacroscopicVariable<dist_type> dens(nx, ny, nz, "r", run_id, save_path);
    MacroscopicVariable<dist_type> velx(nx, ny, nz, "u", run_id, save_path);
    MacroscopicVariable<dist_type> vely(nx, ny, nz, "v", run_id, save_path);
    MacroscopicVariable<dist_type> velz(nx, ny, nz, "w", run_id, save_path);
    MacroscopicVariable<dist_type> temp(nx, ny, nz, "t", run_id, save_path);
    MacroscopicVariable<dist_type> Fx(nx, ny, nz, "Fx", run_id, save_path);
    MacroscopicVariable<dist_type> Fy(nx, ny, nz, "Fy", run_id, save_path);
    MacroscopicVariable<dist_type> Fz(nx, ny, nz, "Fz", run_id, save_path);

    double ag, tau_f, tau_g; // Model parameters.
    findGoodModelParameters(ag, tau_f, tau_g, Pr, Ra, nz, f.CSI(), g.CSI());

    double kinematic_viscosity, thermal_diffusivity;
    kinematic_viscosity = (tau_f - 0.5) / 3.0;
    thermal_diffusivity = (tau_g - 0.5) / 3.0;

    double length_scale, velocity_scale, time_scale; //Characteristic scales.
    velocity_scale = sqrt(ag * nz * 1.0); // Delta T = 1.0
    time_scale = (nz*nz / kinematic_viscosity) * sqrt(Pr / Ra);
    length_scale = velocity_scale * time_scale;

    std::cout << "tau_f               = " << tau_f << "\n";
    std::cout << "tau_g               = " << tau_g << "\n";
    std::cout << "alpha*gravity       = " << ag << "\n";
    std::cout << "Length scale        = " << length_scale << "\n";
    std::cout << "Kinematic viscosity = " << kinematic_viscosity << "\n";
    std::cout << "Diffusivity         = " << thermal_diffusivity << "\n";

    // Print info to terminal.
    f.DisplayLatticeParameters();
    g.DisplayLatticeParameters();
    dens.DisplayInfo();
    velx.DisplayInfo();
    vely.DisplayInfo();
    temp.DisplayInfo();

    // Set initial density to 1 everywhere.
    double rho0 = 1.0;
    dens.SetToConstantValue(rho0);

    // Set initial temperature profile.
    double temp_top = 0.5 / static_cast<double>(nz), temp_bot = 1.0 - temp_top;
    temp.SetLinearGradientZ(temp_bot, temp_top);

    int k_mid = nz / 2;
    for (int i = 0; i < nx; ++ i)
    {
        double eps = 1.0e-2;
        double x = 2 * M_PI * i / nx;
        double perturbation = eps * cos(x);
        temp.AddToValue(perturbation, i, 0, k_mid);
    }


    // Set force information.
    double ag_x = 0, ag_y = 0, ag_z = -std::fabs(ag), ref_temp = 0.5;
    LinearForceAnomalyTemp<dist_type> force_updater(ag_x, ag_y, ag_z, ref_temp, nx, ny, nz);
    force_updater.UpdateForce(Fx, Fy, Fz, temp);

    // Initialise boundary information.
    int num_boundaries = 2;
    BoundaryInfo<dist_type, nd, nq_f> bdry_info_f(num_boundaries);
    BoundaryInfo<dist_type, nd, nq_g> bdry_info_g(num_boundaries);
    NodeInfo node(nx, ny, nz);

    // Construct boundary rules.
    double uwall_x = 0.0, uwall_y = 0.0, uwall_z = 0.0;
    double temp_wall_top = 0.0, temp_wall_bot = 1.0;
    BdryRuleBounceBackTop<dist_type, nd, nq_f> bdry_top_f(&f, &dens, uwall_x, uwall_y, uwall_z);
    BdryRuleBounceBackBottom<dist_type, nd, nq_f> bdry_bot_f(&f, &dens, uwall_x, uwall_y, uwall_z);
    BdryRuleScalarDirichletTop<dist_type, nd, nq_g> bdry_top_g(&g, temp_wall_top, uwall_x, uwall_y, uwall_z);
    BdryRuleScalarDirichletBottom<dist_type, nd, nq_g> bdry_bot_g(&g, temp_wall_bot, uwall_x, uwall_y, uwall_z);
    
    // Add boundary rules.
    // Bottom boundary.
    int bdry_id_bot = 0;
    bdry_info_f.AddBoundaryRule(&bdry_bot_f, bdry_id_bot);
    bdry_info_g.AddBoundaryRule(&bdry_bot_g, bdry_id_bot);
    
    // Top boundary.
    int bdry_id_top = 1;
    bdry_info_f.AddBoundaryRule(&bdry_top_f, bdry_id_top);
    bdry_info_g.AddBoundaryRule(&bdry_top_g, bdry_id_top);

    // Set node info.
    node.SetBoundaryOnBottom(bdry_id_bot);
    node.SetBoundaryOnTop(bdry_id_top);

    FluidEvolverTRT<dist_type, nd, nq_f> fluid_evolver;
    fluid_evolver.SetMagicParameter(1./12.);
    fluid_evolver.SetKinematicViscosity(f, kinematic_viscosity);
    // Note Fx, Fy, Fz must be initialised (i.e. updated) before the next line!
    fluid_evolver.Initialise(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry_info_f);

    ScalarEvolverSRT<dist_type, nd, nq_g> scalar_evolver;
    scalar_evolver.SetScalarDiffusivity(g, thermal_diffusivity);
    scalar_evolver.Initialise(g, temp, velx, vely, velz, node, bdry_info_g);

    // Save initial data.
    std::string save_msg_before = "Initial data saved successfully.\n";
    saveData(0, save_msg_before, f, g, dens, velx, vely, velz, temp);
    
    std::string message_before = "Before simulation:";
    printAverages(message_before, dens, velx, vely, velz, temp);

    // Diagnostics.
    DiagnosticNusselt<dist_type> diagnoser(nx, ny, nz, thermal_diffusivity);
    double nusselt_number;
    int nNuSamples = 20;
    ConvergenceTester<dist_type> convergence_tester(nNuSamples);

    // Run algorithm.
    for (int t = 1; t < nt; ++t)
    {   
        // Perform one timestep of ADE LB algorithm for g using velocity u.
        // This updates g and temp.
        scalar_evolver.DoTimestep(g, temp, velx, vely, velz, node, bdry_info_g);

        // Compute force density F from updated temp.
        force_updater.UpdateForce(Fx, Fy, Fz, temp);

        // Perform one timestep of the standard LB algorithm for f using F.
        // This updates f, dens, velx, vely, and velz.
        fluid_evolver.DoTimestep(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry_info_f);

        // Compute diagnostic.
        const int nNusseltFrequency = 1000;
        if (t % nNusseltFrequency == 0)
        {
            std::cout << "timestep: " << t << "\n";
            // Print averages after.
            std::string message_mid = "";
            printAverages(message_mid, dens, velx, vely, velz, temp);

            // Compute nusselt number.
            nusselt_number = diagnoser.ComputeNusseltNumber(velz, temp);
            // Store it.
            convergence_tester.AddValueToList(nusselt_number);
        }

        // Save data.
        const int min_save_time = 1000;  // Do not save any more often than every <n> timesteps.
        const int max_num_of_saves = 10; // Do not save any more often over a whole simulation than this number of times.
        const int save_time = std::max(nt / max_num_of_saves, min_save_time);
        if (t % save_time == 0)
        {
            // Save data to file.
            std::string save_msg_mid = "";
            saveData(t, save_msg_mid, f, g, dens, velx, vely, velz, temp);
        }
        
        // End loop or not
        if (convergence_tester.HasConverged())
        {   
            std::cout << "Simulation ended after " << t << " timesteps.\n";
            break;
        }
    }

    // Print averages after.
    std::string message_after = "After simulation:";
    printAverages(message_after, dens, velx, vely, velz, temp);

    // Save data after.
    std::string save_msg_after = "Final data saved successfully.\n";
    saveData(nt, save_msg_after, f, g, dens, velx, vely, velz, temp);
}


/***************************************************
 *                                                 *
 *           Helper function definitions           *
 *                                                 *
 ***************************************************/

void getInputParameters(double& rayleigh_number, double& prandtl_number, int& nx, int& nz, int& nt, std::string& run_id, int argc, char* argv[])
{
    // Default values.
    rayleigh_number = 1.0e4;
    prandtl_number = 1.0;
    nx = 100;
    nz = 100;
    nt = 100;
    run_id = "testrun";

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
        nz = atoi(argv[4]);
    }
    if (argc >= 6)
    {
        nt = atoi(argv[5]);
    }
    if (argc >= 7)
    {
        run_id = argv[6];
    }
    std::cout << "Ra     = " << rayleigh_number << "\n";
    std::cout << "Pr     = " << prandtl_number << "\n";
    std::cout << "nx     = " << nx << "\n";
    std::cout << "nz     = " << nz << "\n";
    std::cout << "nt     = " << nt << "\n";
    std::cout << "run_id = " << run_id << "\n";
    return;
}

void findGoodModelParameters(double& alpha_g, double& tau_f, double& tau_g, const double Pr, const double Ra, const int ny, const double csi_f, const double csi_g)
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
        tau_f = (0.5 + csi_f*nu); // Fluid relaxation time in lattice units.
        tau_g = (0.5 + csi_g*kappa); // Thermal relaxation time in lattice units.
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

template <class T>
void printAverages(std::string& message, MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
                   MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz, MacroscopicVariable<T>& temp)
{
    std::cout << message << "\n";
    std::cout << "Average density = " << dens.ComputeAverage() << "\n";
    std::cout << "Average X velocity = " << velx.ComputeAverage() << "\n";
    std::cout << "Average Y velocity = " << vely.ComputeAverage() << "\n";
    std::cout << "Average Z velocity = " << velz.ComputeAverage() << "\n";
    std::cout << "Average temperature = " << temp.ComputeAverage() << "\n";
}

template <class T, int ND, int NQ_F, int NQ_G>
void saveData(int timestep, std::string& message, AbstractLattice<T, ND, NQ_F>& f, AbstractLattice<T, ND, NQ_G>& g, 
             MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
             MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz, MacroscopicVariable<T>& temp)
{
    f.WriteToTextFile(timestep);
    g.WriteToTextFile(timestep);
    dens.WriteToTextFile(timestep);
    velx.WriteToTextFile(timestep);
    vely.WriteToTextFile(timestep);
    velz.WriteToTextFile(timestep);
    temp.WriteToTextFile(timestep);
    std::cout << message;
}

template <class T>
void saveMacroscopic(int timestep, std::string& message, 
             MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
             MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz, MacroscopicVariable<T>& temp)
{
    dens.WriteToTextFile(timestep);
    velx.WriteToTextFile(timestep);
    vely.WriteToTextFile(timestep);
    velz.WriteToTextFile(timestep);
    temp.WriteToTextFile(timestep);
    std::cout << message << ". Timestep = " << timestep << "\n";
}