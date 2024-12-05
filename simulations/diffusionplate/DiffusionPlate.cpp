/**
 * Decaying Taylor-Green vortex flow.
 * Tests fluid initialisation and bulk fluid dynamics in periodic domain (no solid boundaries)
 * nx, ny = 96, 72
 * tau = 0.8, nu = 0.1
 * u0 = 0.03
 * rho0 = 1
 * p0 = 0.
 * 
 * Initialise.
 * Evolve for one charaacteristic timescale
 * t_d = 1 / (nu * (kx^2 + ky^2))
 * where wavenumbers
 * kx = 2*pi / nx
 * ky = 2*pi / ny;
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
    double D, u0, Cp;
    std::string run_id, savepath;
};

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& conc, VectorField<T>& vel, const Parameters& params);

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& conc);

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& conc);

template <typename dfType, int nd, int nq>
void RunLatticeBoltzmannModel(AbstractScalarEvolver<dfType, nd, nq>& scalar_evolver, AbstractLattice<dfType, nd, nq>& g, ScalarField<dfType>& conc, VectorField<dfType>& vel, NodeInfo& node, BoundaryInfo<dfType, nd, nq>& bdry, const Parameters& params);

int main()
{
    Parameters params;

    // Set parameters. All in lattice units.
    const int nd = 2, nq = 5;
    params.nx = 800, params.ny = 1, params.nz = 80;
    params.D = 0.07407; //, tau = 0.72221; // D2Q5 or D2Q9 lattice, cs = 1/sqrt(3).
    params.u0 = 0.05;
    params.Cp = 1.0;

    //std::string run_id = "SRTxx_MEI_F64";
    std::string run_id = "DiffusionPlate160x16";
    std::string savepath = "output";
    params.savepath = savepath;
    
    /****************************************************************************************
     * 
     * Set precision
     * 
     ***************************************************************************************/
    using dfType = double; std::cout << "Double precision (F64)\n";
    //using dfType = float;  std::cout << "Single precision (F32)\n";
    ScalarEvolverSRT<dfType, nd, nq> scalar_evolver_SRTxx; // SRTxx

    // Declare lattice.
    LatticeSoAPull<dfType, nd, nq> g(params.nx, params.ny, params.nz, "g", run_id, savepath);

    // Declare macroscopic variables.
    ScalarField<dfType> conc(params.nx, params.ny, params.nz, "c", run_id, savepath);
    VectorField<dfType> vel(params.nx, params.ny, params.nz, "u", run_id, savepath);

    // Declare node and boundary info.
    NodeInfo node(params.nx, params.ny, params.nz);
    BoundaryInfo<dfType, nd, nq> bdry(4);

    // Declare the boundary rules.
    BdryRuleScalarDirichletBottom<dfType, nd, nq> bdry_bot_g(&g, params.Cp, 0, 0, 0);
    BdryRuleScalarDirichletLeft<dfType, nd, nq> bdry_lef_g(&g, 0.0, params.u0, 0, 0);
    BdryRuleScalarNeumannRight<dfType, nd, nq> bdry_rig_g(&g, &conc, params.u0, 0.0, 0.0);
    BdryRuleScalarNeumannTop<dfType, nd, nq> bdry_top_g(&g, &conc, params.u0, 0.0, 0.0);

    // Set the boundary conditions.
    // Add boundary rules to boundary info.
    int bdry_id_bot = 0, bdry_id_top = 1, bdry_id_lef = 2, bdry_id_rig = 3;
    bdry.AddBoundaryRule(&bdry_bot_g, bdry_id_bot);
    bdry.AddBoundaryRule(&bdry_top_g, bdry_id_top);
    bdry.AddBoundaryRule(&bdry_lef_g, bdry_id_lef);
    bdry.AddBoundaryRule(&bdry_rig_g, bdry_id_rig);

    // Set node info.
    node.SetBoundaryOnLeft(bdry_id_lef);
    node.SetBoundaryOnRight(bdry_id_rig);
    node.SetBoundaryOnBottom(bdry_id_bot);
    node.SetBoundaryOnTop(bdry_id_top);

    // SRTxx
    std::string run_id1 = "D2Q5_800x1x80_SRTxx_F64";
    params.run_id = run_id1;

    RunLatticeBoltzmannModel(scalar_evolver_SRTxx, g, conc, vel, node, bdry, params);

    std::cout << "######################################################\n";
    std::cout << "#                                                    #\n";
    std::cout << "#                      Finished                      #\n";
    std::cout << "#                                                    #\n";
    std::cout << "######################################################\n";
    return 0;
}

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& conc, VectorField<T>& vel, const Parameters& params)
{
    // Constant initial concentration zero everywhere.
    conc.SetToConstantValue(0.0);

    // Constant velocity (u0, 0, 0)
    vel.SetToConstantValue(params.u0, 0);
    vel.SetToConstantValue(0.0, 1);
    vel.SetToConstantValue(0.0, 2);
}

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& conc)
{
    conc.WriteToTextFile(timestep);
    std::cout << message << ". Timestep = " << timestep << "\n";
}

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& conc)
{
    std::cout << message << "\n";
    std::cout << "Average concentration = " << conc.ComputeAverage() << "\n";
}

template <typename dfType, int nd, int nq>
void RunLatticeBoltzmannModel(AbstractScalarEvolver<dfType, nd, nq>& scalar_evolver, AbstractLattice<dfType, nd, nq>& g, ScalarField<dfType>& conc, VectorField<dfType>& vel, NodeInfo& node, BoundaryInfo<dfType, nd, nq>& bdry, const Parameters& params)
{
    std::cout << "~~~~~~~~~~~~~~~~~\n";
    std::cout << params.run_id << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~\n";

    g.SetRunID(params.run_id);
    conc.SetRunID(params.run_id);

    // Set kinematic viscosity.
    scalar_evolver.SetScalarDiffusivity(g, params.D);

    InitialiseMacroscopic(conc, vel, params);

    /****************************************************************************************
     * 
     * Initialise the distribution functions.
     * 
     ***************************************************************************************/
    scalar_evolver.Initialise(g, conc, vel, node, bdry);

    // Print averages before.
    std::string message_before = "Before simulation:";
    PrintAverages(message_before, conc);

    // Save initial conditions.
    //std::string msg = "Start";
    //SaveMacroscopic(0, msg, conc, velx, velz);

    double tolerance = 1.0e-6, L2errornorm=1.0;
    int count = 0, N_compute_norm = 1000;
    auto start = std::chrono::high_resolution_clock::now();
    while (L2errornorm > tolerance)
    {
        // Run Lattice Boltzmann time step.
        scalar_evolver.DoTimestep(g, conc, vel, node, bdry);
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
    PrintAverages(message_after, conc);
    SaveMacroscopic(count, message_after, conc);

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
        scalar_evolver.DoTimestep(g, conc, vel, node, bdry);
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
    PrintAverages(message_after, conc);
    SaveMacroscopic(count, message_after, conc);

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
        scalar_evolver.DoTimestep(g, conc, vel, node, bdry);
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
    PrintAverages(message_after, conc);
    SaveMacroscopic(count, message_after, conc);

    std::cout << "Total timesteps = " << count << "\n";
    std::cout << "Execution time  = " << duration.count() << " seconds" << "\n";
    std::cout << "Update rate     = " << updates_per_second << " MLUPS" << std::endl;
}
