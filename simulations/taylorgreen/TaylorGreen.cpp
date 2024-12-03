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

#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <sstream>

double ComputeTaylorGreenPressure(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double p0, const double rho0, const double u0, const double kx, const double ky);

double ComputeTaylorGreenLBDensity(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double p0, const double rho0, const double u0, const double kx, const double ky);

double ComputeTaylorGreenVelocityX(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double u0, const double kx, const double ky);

double ComputeTaylorGreenVelocityY(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double u0, const double kx, const double ky);

double ComputeTaylorGreenSigmaXX(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double rho0, const double u0, const double kx, const double ky, const double nu);

double ComputeTaylorGreenSigmaXY(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double rho0, const double u0, const double kx, const double ky, const double nu);

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& dens, VectorField<T>& vel, VectorField<T>& force, TensorField<T>& sigma_lb, TensorField<T>& sigma_fd, double td, double p0, double rho0, double u0, double kx, double ky, double nu, int nx, int ny, int nz, bool ConstDensity=false);

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& dens, VectorField<T>& vel);

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& dens, VectorField<T>& vel);

template <typename dfType, int nd, int nq>
void RunLatticeBoltzmannModel(const int nx, const int ny, const int nz, const double nu, const double u0, const double rho0, const double p0, const double kx, const double ky, const double td, const double dt, std::string run_id, std::string savepath, std::string init_string, AbstractFluidEvolver<dfType, nd, nq>& fluid_evolver, AbstractLattice<dfType, nd, nq>& f, ScalarField<dfType>& dens, VectorField<dfType>& vel, VectorField<dfType>& force, TensorField<dfType>& sigma_lb, TensorField<dfType>& sigma_fd, NodeInfo& node, BoundaryInfo<dfType, nd, nq>& bdry);

int main()
{
    // Set parameters. All in lattice units.
    const int nd = 2, nq = 9;
    const int nx = 72, ny = 1, nz = 96;
    const double nu = 0.1; //, tau = 0.8; // D2Q9 lattice, cs = 1/sqrt(3).
    const double u0 = 0.03;
    const double rho0 = 1.0;
    const double p0 = 0.0;
    const double kx = 2*M_PI / static_cast<double>(nx);
    const double ky = 2*M_PI / static_cast<double>(nz);
    const double td = 1.0 / (nu * (kx*kx + ky*ky));    // evaluates to 840 time steps.
    const double dt = 1.0; // dx = 1.0, 

    //std::string run_id = "SRTxx_MEI_F64";
    std::string run_id = "TaylorGreen";
    std::string savepath = "output";
    //std::string init_string = "MEI";
    
    /****************************************************************************************
     * 
     * Set precision
     * 
     ***************************************************************************************/
    //using dfType = double; std::cout << "Double precision (F64)\n";
    using dfType = float;  std::cout << "Single precision (F32)\n";
    FluidEvolverSRT<dfType, nd, nq> fluid_evolver_SRTxx; // SRTxx
    FluidEvolverTRT<dfType, nd, nq> fluid_evolver_TRT12; // TRT
    FluidEvolverTRT<dfType, nd, nq> fluid_evolver_TRT06; // TRT
    FluidEvolverTRT<dfType, nd, nq> fluid_evolver_TRT04; // TRT
    fluid_evolver_TRT12.SetMagicParameter(1./12.);        // TRT12
    fluid_evolver_TRT06.SetMagicParameter(1./6.);       // TRT06
    fluid_evolver_TRT04.SetMagicParameter(1./4.);       // TRT04

    // Declare lattice.
    LatticeSoAPull<dfType, nd, nq> f(nx, ny, nz, "f", run_id, savepath);

    // Declare macroscopic variables.
    ScalarField<dfType> dens(nx, ny, nz, "r", run_id, savepath);
    VectorField<dfType> vel(nx, ny, nz, "u", run_id, savepath);
    VectorField<dfType> force(nx, ny, nz, "F", run_id, savepath);
    TensorField<dfType> sigma_lb(nx, ny, nz, "SL", run_id, savepath);
    TensorField<dfType> sigma_fd(nx, ny, nz, "SF", run_id, savepath);

    // Declare node and boundary info.
    NodeInfo node(nx, ny, nz);
    BoundaryInfo<dfType, nd, nq> bdry(0);

    // SRTxx
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "SRTxx_CEQ_F32", savepath, "CEQ", fluid_evolver_SRTxx, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "SRTxx_FEQ_F32", savepath, "FEQ", fluid_evolver_SRTxx, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "SRTxx_NEQ_F32", savepath, "NEQ", fluid_evolver_SRTxx, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "SRTxx_MEI_F32", savepath, "MEI", fluid_evolver_SRTxx, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);

    // TRT12
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT12_CEQ_F32", savepath, "CEQ", fluid_evolver_TRT12, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT12_FEQ_F32", savepath, "FEQ", fluid_evolver_TRT12, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT12_NEQ_F32", savepath, "NEQ", fluid_evolver_TRT12, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT12_MEI_F32", savepath, "MEI", fluid_evolver_TRT12, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);

    // TRT06
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT06_CEQ_F32", savepath, "CEQ", fluid_evolver_TRT06, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT06_FEQ_F32", savepath, "FEQ", fluid_evolver_TRT06, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT06_NEQ_F32", savepath, "NEQ", fluid_evolver_TRT06, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT06_MEI_F32", savepath, "MEI", fluid_evolver_TRT06, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);

    // TRT04
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT04_CEQ_F32", savepath, "CEQ", fluid_evolver_TRT04, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT04_FEQ_F32", savepath, "FEQ", fluid_evolver_TRT04, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT04_NEQ_F32", savepath, "NEQ", fluid_evolver_TRT04, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);
    RunLatticeBoltzmannModel(nx, ny, nz, nu, u0, rho0, p0, kx, ky, td, dt, "TRT04_MEI_F32", savepath, "MEI", fluid_evolver_TRT04, f, dens, vel, force, sigma_lb, sigma_fd, node, bdry);

    std::cout << "######################################################\n";
    std::cout << "#                                                    #\n";
    std::cout << "#                      Finished                      #\n";
    std::cout << "#                                                    #\n";
    std::cout << "######################################################\n";
    return 0;
}

double ComputeTaylorGreenPressure(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double p0, const double rho0, const double u0, const double kx, const double ky)
{
    const double x = static_cast<double>(i);
    const double y = static_cast<double>(k);

    const double exponential = exp(-2*t / td);
    const double cosines = (ky/kx)*cos(2*kx*x) + (kx/ky)*cos(2*ky*y);
    return p0 - 0.25*rho0*u0*u0 * cosines * exponential;
}

double ComputeTaylorGreenLBDensity(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double p0, const double rho0, const double u0, const double kx, const double ky)
{
    const double p = ComputeTaylorGreenPressure(i, j, k, t, td, p0, rho0, u0, kx, ky);
    return rho0 + 3.0 * (p - p0);
}

double ComputeTaylorGreenVelocityX(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double u0, const double kx, const double ky)
{
    const double x = static_cast<double>(i);
    const double y = static_cast<double>(k);
    return -u0 * sqrt(ky/kx) * cos(kx*x) * sin(ky*y) * exp(-t/td);
}

double ComputeTaylorGreenVelocityY(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double u0, const double kx, const double ky)
{
    const double x = static_cast<double>(i);
    const double y = static_cast<double>(k);
    return u0 * sqrt(kx/ky) * sin(kx*x) * cos(ky*y) * exp(-t/td);
}

double ComputeTaylorGreenSigmaXX(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double rho0, const double u0, const double kx, const double ky, const double nu)
{
    const double x = static_cast<double>(i);
    const double y = static_cast<double>(k);

    return 2.0 * rho0 * nu * u0 * std::sqrt(kx*ky) * std::sin(kx*x) * std::sin(ky*y) * std::exp(-t/td);
}

double ComputeTaylorGreenSigmaXY(const int i, const int j [[maybe_unused]], const int k, const double t, const double td, const double rho0, const double u0, const double kx, const double ky, const double nu)
{
    const double x = static_cast<double>(i);
    const double y = static_cast<double>(k);

    const double sqrts = std::sqrt(kx*kx*kx/ky) - std::sqrt(ky*ky*ky/kx);
    const double coses = std::cos(kx*x) * std::cos(ky*y);
    const double expon = std::exp(-t / td);
    
    return rho0 * nu * u0 * sqrts * coses * expon;
}

template <typename T>
void InitialiseMacroscopic(ScalarField<T>& dens, VectorField<T>& vel, VectorField<T>& force, TensorField<T>& sigma_lb, TensorField<T>& sigma_fd, double td, double p0, double rho0, double u0, double kx, double ky, double nu, int nx, int ny, int nz, bool ConstDensity)
{
    std::cout << "Initially constant density (1 = true, 0 = false): " << ConstDensity << "\n";

    vel.SetToConstantValue(0.0, 1);
    force.SetToConstantValue(0.0);
    sigma_lb.SetToConstantValue(0.0);
    sigma_fd.SetToConstantValue(0.0);

    auto compute_density = std::bind(ComputeTaylorGreenLBDensity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, p0, rho0, u0, kx, ky);
    auto compute_velx = std::bind(ComputeTaylorGreenVelocityX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, u0, kx, ky);
    auto compute_vely = std::bind(ComputeTaylorGreenVelocityY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, u0, kx, ky);
    auto compute_sxx = std::bind(ComputeTaylorGreenSigmaXX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, rho0, u0, kx, ky, nu);
    auto compute_sxy = std::bind(ComputeTaylorGreenSigmaXY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, rho0, u0, kx, ky, nu);

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                if (ConstDensity)
                {
                    dens.SetValue(rho0, i, j, k);
                }
                else
                {
                    double rho = compute_density(i, j, k);
                    dens.SetValue(rho, i, j, k);
                }

                double u = compute_velx(i, j, k);
                vel.SetValue(u, 0, i, j, k);

                double v = compute_vely(i, j, k);
                vel.SetValue(v, 2, i, j, k);

                double sxx = compute_sxx(i, j, k);
                sigma_lb.SetValue(sxx, 0, 0, i, j, k);
                sigma_lb.SetValue(-sxx, 2, 2, i, j, k); // sigma_yy = -sigma_xx
                sigma_fd.SetValue(sxx, 0, 0, i, j, k);
                sigma_fd.SetValue(-sxx, 2, 2, i, j, k); // sigma_yy = -sigma_xx

                double sxy = compute_sxy(i, j, k);
                sigma_lb.SetValue(sxy, 0, 2, i, j, k);
                sigma_lb.SetValue(sxy, 2, 0, i, j, k); // sigma_yx = sigma_xy
                sigma_fd.SetValue(sxy, 0, 2, i, j, k);
                sigma_fd.SetValue(sxy, 2, 0, i, j, k); // because stress tensor is symmetric.
            }
        }
    }
}

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, ScalarField<T>& dens, VectorField<T>& vel)
{
    dens.WriteToTextFile(timestep);
    vel.WriteToTextFile(timestep);
    std::cout << message << ". Timestep = " << timestep << "\n";
}

template <typename T>
void PrintAverages(std::string& message, ScalarField<T>& dens, VectorField<T>& vel)
{
    std::cout << message << "\n";
    std::cout << "Average density = " << dens.ComputeAverage() << "\n";
    std::cout << "Average X velocity = " << vel.ComputeAverage(0) << "\n";
    std::cout << "Average Y velocity = " << vel.ComputeAverage(1) << "\n";
    std::cout << "Average Z velocity = " << vel.ComputeAverage(2) << "\n";
}

template <typename dfType, int nd, int nq>
void RunLatticeBoltzmannModel(const int nx, const int ny, const int nz, const double nu, const double u0, const double rho0, const double p0, const double kx, const double ky, const double td, const double dt, std::string run_id, std::string savepath, std::string init_string, AbstractFluidEvolver<dfType, nd, nq>& fluid_evolver, AbstractLattice<dfType, nd, nq>& f, ScalarField<dfType>& dens, VectorField<dfType>& vel, VectorField<dfType>& force, TensorField<dfType>& sigma_lb, TensorField<dfType>& sigma_fd, NodeInfo& node, BoundaryInfo<dfType, nd, nq>& bdry)
{
    std::cout << "~~~~~~~~~~~~~~~~~\n";
    std::cout << run_id << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~\n";

    // Set kinematic viscosity.
    fluid_evolver.SetKinematicViscosity(f, nu);

    // Initialise macroscopic variables.
    bool set_uniform_density;
    if (init_string == "MEI" || init_string == "CEQ")
    {
        set_uniform_density = true;
    }
    else
    {
        set_uniform_density = false;
    }
    InitialiseMacroscopic(dens, vel, force, sigma_lb, sigma_fd, td, p0, rho0, u0, kx, ky, nu, nx, ny, nz, set_uniform_density);

    /****************************************************************************************
     * 
     * Initialise the distribution functions.
     * 
     ***************************************************************************************/
    fluid_evolver.Initialise(f, dens, vel, force, node, bdry, init_string);

    // Open file for writing.
    std::ofstream write_l2error(savepath+"/"+run_id+"_l2error.csv");
    write_l2error << "time,r,p,u,v,sxx_lb,sxy_lb,sxx_fd,sxy_fd" << std::endl;
    write_l2error.precision(6);
    write_l2error.setf(std::ios::scientific);
    write_l2error.setf(std::ios::showpos);

    // Compute error at 0th timestep.
    auto compute_density = std::bind(ComputeTaylorGreenLBDensity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, p0, rho0, u0, kx, ky);
    auto compute_velx = std::bind(ComputeTaylorGreenVelocityX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, u0, kx, ky);
    auto compute_vely = std::bind(ComputeTaylorGreenVelocityY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, u0, kx, ky);
    auto compute_sxx = std::bind(ComputeTaylorGreenSigmaXX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, rho0, u0, kx, ky, nu);
    auto compute_sxy = std::bind(ComputeTaylorGreenSigmaXY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, rho0, u0, kx, ky, nu);
    double l2error_r = ComputeL2ErrorScalar(dens, compute_density);
    double l2error_p = ComputeL2ErrorScalar(dens, compute_density, rho0);
    double l2error_u = ComputeL2ErrorVectorComponent(vel, 0, compute_velx);
    double l2error_v = ComputeL2ErrorVectorComponent(vel, 2, compute_vely);
    double l2error_sxx_lb = ComputeL2ErrorTensorComponent(sigma_lb, 0, 0, compute_sxx);
    double l2error_sxy_lb = ComputeL2ErrorTensorComponent(sigma_lb, 0, 2, compute_sxy);
    double l2error_sxx_fd = ComputeL2ErrorTensorComponent(sigma_fd, 0, 0, compute_sxx);
    double l2error_sxy_fd = ComputeL2ErrorTensorComponent(sigma_fd, 0, 2, compute_sxy);

    // Save to file.
    write_l2error << -0.5 << "," << l2error_r << "," << l2error_p << "," << l2error_u << "," << l2error_v << "," << l2error_sxx_lb << "," << l2error_sxy_lb << "," << l2error_sxx_fd << "," << l2error_sxy_fd << std::endl;

    // Print averages before.
    std::string message_before = "Before simulation:";
    PrintAverages(message_before, dens, vel);

    //std::string msg = "Start";
    //SaveMacroscopic(0, msg, dens, velx, velz);

    for (double t = 0.0; t <= td; t+=dt)
    {
        // Run Lattice Boltzmann time step.
        fluid_evolver.DoTimestep(f, dens, vel, sigma_lb, force, node, bdry);

        // Compute stress from finite differences.
        ComputeStressFiniteDifferencePeriodic(sigma_fd, vel, nu, rho0);

        if (t == 0.0)
        {
            // Print averages.
            std::string message = "t=0";
            PrintAverages(message, dens, vel);
        }

        // Calculate L2 error norm for p, ux, uy.
        auto compute_density = std::bind(ComputeTaylorGreenLBDensity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, p0, rho0, u0, kx, ky);
        auto compute_velx = std::bind(ComputeTaylorGreenVelocityX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, u0, kx, ky);
        auto compute_vely = std::bind(ComputeTaylorGreenVelocityY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, u0, kx, ky);
        auto compute_sxx = std::bind(ComputeTaylorGreenSigmaXX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, rho0, u0, kx, ky, nu);
        auto compute_sxy = std::bind(ComputeTaylorGreenSigmaXY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, rho0, u0, kx, ky, nu);

        l2error_r = ComputeL2ErrorScalar(dens, compute_density);
        l2error_p = ComputeL2ErrorScalar(dens, compute_density, rho0);
        l2error_u = ComputeL2ErrorVectorComponent(vel, 0, compute_velx);
        l2error_v = ComputeL2ErrorVectorComponent(vel, 2, compute_vely);
        l2error_sxx_lb = ComputeL2ErrorTensorComponent(sigma_lb, 0, 0, compute_sxx);
        l2error_sxy_lb = ComputeL2ErrorTensorComponent(sigma_lb, 0, 2, compute_sxy);
        l2error_sxx_fd = ComputeL2ErrorTensorComponent(sigma_fd, 0, 0, compute_sxx);
        l2error_sxy_fd = ComputeL2ErrorTensorComponent(sigma_fd, 0, 2, compute_sxy);

        // Save to file.
        write_l2error << t << "," << l2error_r << "," << l2error_p << "," << l2error_u << "," << l2error_v << "," << l2error_sxx_lb << "," << l2error_sxy_lb << "," << l2error_sxx_fd << "," << l2error_sxy_fd << std::endl;
    }

    // Close files.
    write_l2error.close();

    // Print averages after.
    std::string message_after = "After simulation:";
    PrintAverages(message_after, dens, vel);

    //msg = "End";
    //SaveMacroscopic(td, msg, dens, velx, velz);
    std::cout << "---------------------------------------------------\n";
    std::cout << "Final L2 error norms:\n";
    std::cout << td << "," << l2error_r << "," << l2error_p << "," << l2error_u << "," << l2error_v << "," << l2error_sxx_lb << "," << l2error_sxy_lb << "," << l2error_sxx_fd << "," << l2error_sxy_fd << "\n";
    std::cout << "---------------------------------------------------\n";
}
