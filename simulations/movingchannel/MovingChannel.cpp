/**
 * Moving 2d channel.
 * 
 * Initialise constant velocity u0 = (0.05, 0, 0), constant concentration C=1
 * 
 * Periodic in x; walls move with speed u0.
 * 
 * Evolve both velocity and concentration DFs.
 */

#include "lbcpp.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <sstream>

struct Parameters
{
    int nx, ny, nz;
    double nu, D, u0, C0, gz;
    std::string run_id, savepath;
};

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& rho, VectorField<T>& vel, VectorField<T>& force, ScalarField<T>& conc, const Parameters& params);

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& rho, VectorField<T>& vel, ScalarField<T>& conc);

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& rho, VectorField<T>& vel, ScalarField<T>& conc);

template <typename T, int nd, int nqf, int nqg>
void RunLatticeBoltzmannModel(AbstractFluidEvolver<T, nd, nqf>& fluid_evolver, AbstractLattice<T, nd, nqf>& f, AbstractScalarEvolver<T, nd, nqg>& scalar_evolver, AbstractLattice<T, nd, nqg>& g, ScalarField<T>& rho, VectorField<T>& vel, VectorField<T>& force, ScalarField<T>& conc, NodeInfo& node, BoundaryInfo<T, nd, nqf>& bdryf, BoundaryInfo<T, nd, nqg>& bdryg, const Parameters& params);

int main()
{
    Parameters params;

    // Set parameters. All in lattice units.
    const int nd = 2, nqf = 9, nqg = 9, n = 50;
    params.nx = n, params.ny = 1, params.nz = n;
    params.nu = 0.07407;
    params.D = 0.07407; //, tau = 0.72221; // D2Q5 or D2Q9 lattice, cs = 1/sqrt(3).
    params.u0 = 0.05;
    params.C0 = 1.0;
    params.gz = 0.0001; // Gravity in z direction.

    //std::string run_id = "SRTxx_MEI_F64";
    std::string run_id = "MovingChannel_D2Q9_gz0001_50x1x50";
    std::string savepath = "output";
    params.run_id = run_id;
    params.savepath = savepath;
    
    // Set precision.
    using dfType = double; std::cout << "Double precision (F64)\n";
    //using dfType = float;  std::cout << "Single precision (F32)\n";
    FluidEvolverSRT<dfType, nd, nqf> fluid_evolver;
    ScalarEvolverSRT<dfType, nd, nqg> scalar_evolver; // SRTxx

    // Declare lattice.
    LatticeSoAPull<dfType, nd, nqf> f(params.nx, params.ny, params.nz, "f", run_id, savepath);
    LatticeSoAPull<dfType, nd, nqg> g(params.nx, params.ny, params.nz, "g", run_id, savepath);

    // Declare macroscopic variables.
    ScalarField<dfType> rho(params.nx, params.ny, params.nz, "r", run_id, savepath);
    VectorField<dfType> vel(params.nx, params.ny, params.nz, "u", run_id, savepath);
    VectorField<dfType> force(params.nx, params.ny, params.nz, "F", run_id, savepath);
    ScalarField<dfType> conc(params.nx, params.ny, params.nz, "c", run_id, savepath);

    // Declare node and boundary info.
    NodeInfo node(params.nx, params.ny, params.nz);
    BoundaryInfo<dfType, nd, nqf> bdryf(2);
    BoundaryInfo<dfType, nd, nqg> bdryg(2);

    // Declare the boundary rules.
    BdryRuleBounceBackBottom<dfType, nd, nqf> bdryf_bot(&f, &rho, params.u0, 0.0, 0.0);
    BdryRuleBounceBackTop<dfType, nd, nqf> bdryf_top(&f, &rho, params.u0, 0.0, 0.0);
    BdryRuleScalarDirichletBottom<dfType, nd, nqg> bdryg_bot(&g, params.C0, params.u0, 0.0, 0.0);
    BdryRuleScalarDirichletTop<dfType, nd, nqg> bdryg_top(&g, params.C0, params.u0, 0.0, 0.0);

    // Set the boundary conditions.
    // Add boundary rules to boundary info.
    int bdry_id_bot = 0, bdry_id_top = 1;
    bdryf.AddBoundaryRule(&bdryf_bot, bdry_id_bot);
    bdryf.AddBoundaryRule(&bdryf_top, bdry_id_top);
    bdryg.AddBoundaryRule(&bdryg_bot, bdry_id_bot);
    bdryg.AddBoundaryRule(&bdryg_top, bdry_id_top);

    // Set node info.
    node.SetBoundaryOnBottom(bdry_id_bot);
    node.SetBoundaryOnTop(bdry_id_top);

    RunLatticeBoltzmannModel(fluid_evolver, f, scalar_evolver, g, rho, vel, force, conc, node, bdryf, bdryg, params);

    std::cout << "######################################################\n";
    std::cout << "#                                                    #\n";
    std::cout << "#                      Finished                      #\n";
    std::cout << "#                                                    #\n";
    std::cout << "######################################################\n";
    return 0;
}

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& rho, VectorField<T>& vel, VectorField<T>& force, ScalarField<T>& conc, const Parameters& params)
{
    // Constant pressure everywhere.
    rho.SetToConstantValue(1.0);
    
    // Constant velocity (u0, 0, 0)
    vel.SetToConstantValue(params.u0, 0);
    vel.SetToConstantValue(0.0, 1);
    vel.SetToConstantValue(0.0, 2);

    // Constant force (0, 0, gz)
    force.SetToConstantValue(0.0);
    force.SetToConstantValue(params.gz, 2);

    // Constant initial concentration to C0 everywhere.
    conc.SetToConstantValue(params.C0);
}

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& rho, VectorField<T>& vel, ScalarField<T>& conc)
{   
    rho.WriteToTextFile(timestep);
    vel.WriteToTextFile(timestep);
    conc.WriteToTextFile(timestep);
    std::cout << message << ". Timestep = " << timestep << "\n";
}

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& rho, VectorField<T>& vel, ScalarField<T>& conc)
{
    std::cout << message << "\n";
    std::cout << "Average density       = " << rho.ComputeAverage() << "\n";
    std::cout << "Average velx          = " << vel.ComputeAverage(0) << "\n";
    std::cout << "Average vely          = " << vel.ComputeAverage(1) << "\n";
    std::cout << "Average velz          = " << vel.ComputeAverage(2) << "\n";
    std::cout << "Average concentration = " << conc.ComputeAverage() << "\n";
}

template <typename T, int nd, int nqf, int nqg>
void RunLatticeBoltzmannModel(AbstractFluidEvolver<T, nd, nqf>& fluid_evolver, AbstractLattice<T, nd, nqf>& f, AbstractScalarEvolver<T, nd, nqg>& scalar_evolver, AbstractLattice<T, nd, nqg>& g, ScalarField<T>& rho, VectorField<T>& vel, VectorField<T>& force, ScalarField<T>& conc, NodeInfo& node, BoundaryInfo<T, nd, nqf>& bdryf, BoundaryInfo<T, nd, nqg>& bdryg, const Parameters& params)
{
    std::cout << "~~~~~~~~~~~~~~~~~\n";
    std::cout << params.run_id << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~\n";

    g.SetRunID(params.run_id);
    conc.SetRunID(params.run_id);

    // Set kinematic viscosity and diffusivity.
    fluid_evolver.SetKinematicViscosity(f, params.nu);
    scalar_evolver.SetScalarDiffusivity(g, params.D);

    // Initialise the macroscopic variables.
    InitialiseMacroscopic(rho, vel, force, conc, params);

    // Initialise the distribution functions.
    fluid_evolver.Initialise(f, rho, vel, force, node, bdryf, "NEQ");
    scalar_evolver.Initialise(g, conc, vel, node, bdryg);

    // Print averages before.
    std::string message_before = "Before simulation:";
    PrintAverages(message_before, rho, vel, conc);

    // Save initial conditions.
    //std::string msg = "Start";
    //SaveMacroscopic(0, msg, conc, velx, velz);

    double tolerance = 1.0e-6, L2errornorm=1.0;
    int count = 0, N_compute_norm = 1000;
    auto start = std::chrono::high_resolution_clock::now();
    while (L2errornorm > tolerance)
    {
        // Run Lattice Boltzmann time step.
        fluid_evolver.DoTimestep(f, rho, vel, force, node, bdryf);
        scalar_evolver.DoTimestep(g, conc, vel, node, bdryg);
        count++;
        if (count % N_compute_norm == 0)
        {
            L2errornorm = ComputeL2ErrorDistribution(g);
            std::cout << "Timestep = " << count << ", error = " << L2errornorm << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    int total_updates = params.nx * params.ny * params.nz * count;
    double updates_per_second = total_updates / (1.0e6*duration.count());
    
    // Print averages after.
    std::string message_after = "After error < 1.0e-6:";
    PrintAverages(message_after, rho, vel, conc);
    SaveMacroscopic(count, message_after, rho, vel, conc);

    std::cout << "Total timesteps = " << count << "\n";
    std::cout << "Execution time  = " << duration.count() << " seconds" << "\n";
    std::cout << "Update rate     = " << updates_per_second << " MLUPS" << std::endl;

    ///////////////////////////////////
    // Continue to higher accuracy.
    ///////////////////////////////////
    tolerance = 1.0e-8;
    while (L2errornorm > tolerance)
    {
        // Run Lattice Boltzmann time step.
        fluid_evolver.DoTimestep(f, rho, vel, force, node, bdryf);
        scalar_evolver.DoTimestep(g, conc, vel, node, bdryg);
        count++;
        if (count % N_compute_norm == 0)
        {
            L2errornorm = ComputeL2ErrorDistribution(g);
            std::cout << "Timestep = " << count << ", error = " << L2errornorm << std::endl;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    total_updates = params.nx * params.ny * params.nz * count;
    updates_per_second = total_updates / (1.0e6*duration.count());
    
    // Print averages after.
    message_after = "After error < 1.0e-8:";
    PrintAverages(message_after, rho, vel, conc);
    SaveMacroscopic(count, message_after, rho, vel, conc);

    std::cout << "Total timesteps = " << count << "\n";
    std::cout << "Execution time  = " << duration.count() << " seconds" << "\n";
    std::cout << "Update rate     = " << updates_per_second << " MLUPS" << std::endl;

    ///////////////////////////////////
    // Continue to higher accuracy.
    ///////////////////////////////////
    tolerance = 1.0e-10;
    while (L2errornorm > tolerance)
    {
        // Run Lattice Boltzmann time step.
        fluid_evolver.DoTimestep(f, rho, vel, force, node, bdryf);
        scalar_evolver.DoTimestep(g, conc, vel, node, bdryg);
        count++;
        if (count % N_compute_norm == 0)
        {
            L2errornorm = ComputeL2ErrorDistribution(g);
            std::cout << "Timestep = " << count << ", error = " << L2errornorm << std::endl;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    total_updates = params.nx * params.ny * params.nz * count;
    updates_per_second = total_updates / (1.0e6*duration.count());
    
    // Print averages after.
    message_after = "After error < 1.0e-10:";
    PrintAverages(message_after, rho, vel, conc);
    SaveMacroscopic(count, message_after, rho, vel, conc);

    std::cout << "Total timesteps = " << count << "\n";
    std::cout << "Execution time  = " << duration.count() << " seconds" << "\n";
    std::cout << "Update rate     = " << updates_per_second << " MLUPS" << std::endl;
}
