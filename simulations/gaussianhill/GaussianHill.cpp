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

double ComputeGaussianHillConcentration(const int i, const int j [[maybe_unused]], const int k, const double t, const double C0, const double sigma0, const double x0, const double y0, const double u0, const double v0, const double D);

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& conc, VectorField<T>& vel, const int nx, const int ny, const int nz, const double C0, const double sigma0, const double x0, const double y0, const double u0, const double v0, const double D);

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& conc);

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& conc);

template <typename dfType, int nd, int nq>
void RunLatticeBoltzmannModel(const int nx, const int ny, const int nz, const double C0, const double sigma0, const double x0, const double y0, const double u0, const double v0, const double D, const double dt, const double tf, std::string run_id, std::string savepath, AbstractScalarEvolver<dfType, nd, nq>& scalar_evolver, AbstractLattice<dfType, nd, nq>& g, ScalarField<dfType>& conc, VectorField<dfType>& vel, NodeInfo& node, BoundaryInfo<dfType, nd, nq>& bdry);

int main()
{
    // Set parameters. All in lattice units.
    const int nd = 2, nq = 5;
    const int nx = 512, ny = 1, nz = 512;
    const double D = 0.0043; //, tau = 0.5129; // D2Q5 or D2Q9 lattice, cs = 1/sqrt(3).
    const double u0 = 0.1, v0 = 0.1;
    const double C0 = 1;
    const double sigma0 = 10;
    const double x0 = 200, y0 = 200;
    const double dt = 1.0;
    const double tf = 201;

    //std::string run_id = "SRTxx_MEI_F64";
    std::string run_id = "GaussianHill";
    std::string savepath = "output";
    
    /****************************************************************************************
     * 
     * Set precision
     * 
     ***************************************************************************************/
    using dfType = double; std::cout << "Double precision (F64)\n";
    //using dfType = float;  std::cout << "Single precision (F32)\n";
    ScalarEvolverSRT<dfType, nd, nq> scalar_evolver_SRTxx; // SRTxx
    ScalarEvolverTRT<dfType, nd, nq> scalar_evolver_TRT12; // TRT
    ScalarEvolverTRT<dfType, nd, nq> scalar_evolver_TRT06; // TRT
    ScalarEvolverTRT<dfType, nd, nq> scalar_evolver_TRT04; // TRT
    scalar_evolver_TRT12.SetMagicParameter(1./12.);        // TRT12
    scalar_evolver_TRT06.SetMagicParameter(1./6.);       // TRT06
    scalar_evolver_TRT04.SetMagicParameter(1./4.);       // TRT04

    // Declare lattice.
    LatticeSoAPull<dfType, nd, nq> g(nx, ny, nz, "g", run_id, savepath);

    // Declare macroscopic variables.
    ScalarField<dfType> conc(nx, ny, nz, "c", run_id, savepath);
    VectorField<dfType> vel(nx, ny, nz, "u", run_id, savepath);

    // Declare node and boundary info.
    NodeInfo node(nx, ny, nz);
    BoundaryInfo<dfType, nd, nq> bdry(0);

    // SRTxx
    RunLatticeBoltzmannModel(nx, ny, nz, C0, sigma0, x0, y0, u0, v0, D, dt, tf, "D2Q5_SRTxx_F64", savepath, scalar_evolver_SRTxx, g, conc, vel, node, bdry);

    // TRT12
    RunLatticeBoltzmannModel(nx, ny, nz, C0, sigma0, x0, y0, u0, v0, D, dt, tf, "D2Q5_TRT12_F64", savepath, scalar_evolver_TRT12, g, conc, vel, node, bdry);

    // TRT06
    RunLatticeBoltzmannModel(nx, ny, nz, C0, sigma0, x0, y0, u0, v0, D, dt, tf, "D2Q5_TRT06_F64", savepath, scalar_evolver_TRT06, g, conc, vel, node, bdry);

    // TRT04
    RunLatticeBoltzmannModel(nx, ny, nz, C0, sigma0, x0, y0, u0, v0, D, dt, tf, "D2Q5_TRT04_F64", savepath, scalar_evolver_TRT04, g, conc, vel, node, bdry);

    std::cout << "######################################################\n";
    std::cout << "#                                                    #\n";
    std::cout << "#                      Finished                      #\n";
    std::cout << "#                                                    #\n";
    std::cout << "######################################################\n";
    return 0;
}

double ComputeGaussianHillConcentration(const int i, const int j [[maybe_unused]], const int k, const double t, const double C0, const double sigma0, const double x0, const double y0, const double u0, const double v0, const double D)
{
    const double x = static_cast<double>(i), y = static_cast<double>(k);
    const double sigmaD = std::sqrt(2.0 * D * t);
    const double arg_numerator = (x - x0 - u0*t)*(x - x0 - u0*t) + (y - y0 - v0*t)*(y - y0 - v0*t);
    const double arg_denominator = 2.0 * (sigma0 * sigma0 + sigmaD * sigmaD);
    const double factor = sigma0*sigma0 / (sigma0*sigma0 + sigmaD*sigmaD) * C0;
    return factor * std::exp(-arg_numerator / arg_denominator);
}

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& conc, VectorField<T>& vel, const int nx, const int ny, const int nz, const double C0, const double sigma0, const double x0, const double y0, const double u0, const double v0, const double D)
{
    vel.SetToConstantValue(u0, 0);
    vel.SetToConstantValue(0.0, 1);
    vel.SetToConstantValue(v0, 2);

    auto compute_concentration = std::bind(ComputeGaussianHillConcentration, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, C0, sigma0, x0, y0, u0, v0, D);

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                double C_ = compute_concentration(i, j, k);
                conc.SetValue(C_, i, j, k);
            }
        }
    }
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
void RunLatticeBoltzmannModel(const int nx, const int ny, const int nz, const double C0, const double sigma0, const double x0, const double y0, const double u0, const double v0, const double D, const double dt, const double tf, std::string run_id, std::string savepath, AbstractScalarEvolver<dfType, nd, nq>& scalar_evolver, AbstractLattice<dfType, nd, nq>& g, ScalarField<dfType>& conc, VectorField<dfType>& vel, NodeInfo& node, BoundaryInfo<dfType, nd, nq>& bdry)
{
    std::cout << "~~~~~~~~~~~~~~~~~\n";
    std::cout << run_id << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~\n";

    g.SetRunID(run_id);
    conc.SetRunID(run_id);

    // Set kinematic viscosity.
    scalar_evolver.SetScalarDiffusivity(g, D);

    InitialiseMacroscopic(conc, vel, nx, ny, nz, C0, sigma0, x0, y0, u0, v0, D);

    /****************************************************************************************
     * 
     * Initialise the distribution functions.
     * 
     ***************************************************************************************/
    scalar_evolver.Initialise(g, conc, vel, node, bdry);

    // Open file for writing.
    std::ofstream write_l2error(savepath+"/"+run_id+"_l2error.csv");
    write_l2error << "time,c" << std::endl;
    write_l2error.precision(6);
    write_l2error.setf(std::ios::scientific);
    write_l2error.setf(std::ios::showpos);

    // Compute error at 0th timestep.
    auto compute_concentration = std::bind(ComputeGaussianHillConcentration, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, C0, sigma0, x0, y0, u0, v0, D);
    double l2error_c = ComputeL2ErrorScalar(conc, compute_concentration);

    // Save to file.
    write_l2error << -0.5 << "," << l2error_c << std::endl;

    // Print averages before.
    std::string message_before = "Before simulation:";
    PrintAverages(message_before, conc);

    //std::string msg = "Start";
    //SaveMacroscopic(0, msg, conc, velx, velz);

    auto start = std::chrono::high_resolution_clock::now();

    for (int t = 0; t <= static_cast<int>(tf); t+=static_cast<int>(dt))
    {
        //std::cout << "t = " << t << std::endl;
        // Run Lattice Boltzmann time step.
        scalar_evolver.DoTimestep(g, conc, vel, node, bdry);

        //std::cout << "Lattice Boltzmann time step complete." << std::endl;
        if (t == 200)
        {
            // Print averages.
            std::string message = "t=200";
            //PrintAverages(message, conc);
            SaveMacroscopic(t, message, conc);
        }

        // Calculate L2 error norm for p, ux, uy.
        auto compute_concentration = std::bind(ComputeGaussianHillConcentration, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, C0, sigma0, x0, y0, u0, v0, D);

        l2error_c = ComputeL2ErrorScalar(conc, compute_concentration);

        //std::cout << "L2 error norm calculated" << std::endl;

        // Save to file.
        write_l2error << t << "," << l2error_c << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    int total_updates = nx*ny*nz*tf;
    double updates_per_second = total_updates / (1.0e6*duration.count());
    

    // Close files.
    write_l2error.close();

    // Print averages after.
    std::string message_after = "After simulation:";
    PrintAverages(message_after, conc);
    SaveMacroscopic(tf, message_after, conc);

    //msg = "End";
    //SaveMacroscopic(td, msg, conc, velx, velz);
    std::cout << "---------------------------------------------------\n";
    std::cout << "Final L2 error norms:\n";
    std::cout << tf << "," << l2error_c << "\n";
    std::cout << "---------------------------------------------------\n";

    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "Lattice updates per second = " << updates_per_second << " MLUPS" << std::endl;
}
