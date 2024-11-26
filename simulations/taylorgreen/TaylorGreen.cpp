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

template <typename T>
void SaveMacroscopic(int timstep, std::string& message, 
             MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
             MacroscopicVariable<T>& velz);

template <typename T>
void PrintAverages(std::string& message, MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
                   MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz);

int main()
{
    // Set parameters. All in lattice units.
    using dfType = double; std::cout << "Double precision (F64)\n";
    //using dfType = float; std::cout << "Single precision (F32)\n";
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

    std::string run_id = "SRTxx_WEI_F64";
    std::string savepath = "output";

    std::cout << "~~~~~~~~~~~~~~~~~\n";
    std::cout << run_id << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~\n";

    // Declare arrays.
    LatticeSoAPull<dfType, nd, nq> f(nx, ny, nz, "f", run_id, savepath);
    MacroscopicVariable<dfType> dens(nx, ny, nz, "r", run_id, savepath);
    MacroscopicVariable<dfType> velx(nx, ny, nz, "u", run_id, savepath);
    MacroscopicVariable<dfType> vely(nx, ny, nz, "w", run_id, savepath);
    MacroscopicVariable<dfType> velz(nx, ny, nz, "v", run_id, savepath);
    MacroscopicVariable<dfType> Fx(nx, ny, nz, "Fx", run_id, savepath);
    MacroscopicVariable<dfType> Fy(nx, ny, nz, "Fx", run_id, savepath);
    MacroscopicVariable<dfType> Fz(nx, ny, nz, "Fx", run_id, savepath);

    /****************************************************************************************
     * 
     * Declare collision operator.
     * 
     ***************************************************************************************/
    FluidEvolverSRT<dfType, nd, nq> fluid_evolver; // SRTxx
    //FluidEvolverTRT<dfType, nd, nq> fluid_evolver; // TRT
    //fluid_evolver.SetMagicParameter(1./4.);        // TRT04
    //fluid_evolver.SetMagicParameter(1./12.);       // TRT12
    
    fluid_evolver.SetKinematicViscosity(f, nu);

    // Declare node and boundary info.
    NodeInfo node(nx, ny, nz);
    BoundaryInfo<dfType, nd, nq> bdry(0);

    // Initialise macroscopic variables.
    auto compute_density = std::bind(ComputeTaylorGreenLBDensity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, p0, rho0, u0, kx, ky);
    auto compute_velx = std::bind(ComputeTaylorGreenVelocityX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, u0, kx, ky);
    auto compute_vely = std::bind(ComputeTaylorGreenVelocityY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 0.0, td, u0, kx, ky);

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                double rho = compute_density(i, j, k);
                dens.SetValue(rho, i, j, k);

                double u = compute_velx(i, j, k);
                velx.SetValue(u, i, j, k);

                double v = compute_vely(i, j, k);
                velz.SetValue(v, i, j, k);
            }
        }
    }
    vely.SetToConstantValue(0.0);
    Fx.SetToConstantValue(0.0);
    Fy.SetToConstantValue(0.0);
    Fz.SetToConstantValue(0.0);

    /****************************************************************************************
     * 
     * Initialise the distribution functions.
     * 
     ***************************************************************************************/
    //fluid_evolver.InitialiseEquilibrium(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry); // FEQ
    //fluid_evolver.Initialise(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry);          // NEQ
    fluid_evolver.InitialiseWei(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry);       // WEI

    // Open file for writing.
    std::ofstream write_l2error(savepath+"/"+run_id+"_l2error.csv");
    write_l2error << "time,r,u,v" << std::endl;
    write_l2error.precision(6);
    write_l2error.setf(std::ios::scientific);
    write_l2error.setf(std::ios::showpos);

    // Compute error at 0th timestep.
    double l2error_r = ComputeL2ErrorScalar(dens, compute_density);
    double l2error_u = ComputeL2ErrorScalar(velx, compute_velx);
    double l2error_v = ComputeL2ErrorScalar(velz, compute_vely);

    // Save to file.
    write_l2error << -0.5 << "," << l2error_r << "," << l2error_u << "," << l2error_v << std::endl;

    // Print averages before.
    std::string message_before = "Before simulation:";
    PrintAverages(message_before, dens, velx, vely, velz);

    std::string msg = "Start";
    SaveMacroscopic(0, msg, dens, velx, velz);

    for (double t = 0.0; t <= td; t+=dt)
    {
        // Run Lattice Boltzmann time step.
        fluid_evolver.DoTimestep(f, dens, velx, vely, velz, Fx, Fy, Fz, node, bdry);

        if (t == 0.0)
        {
            // Print averages.
            std::string message = "t=0";
            PrintAverages(message, dens, velx, vely, velz);
        }

        // Calculate L2 error norm for p, ux, uy.
        auto compute_density = std::bind(ComputeTaylorGreenLBDensity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, p0, rho0, u0, kx, ky);
        auto compute_velx = std::bind(ComputeTaylorGreenVelocityX, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, u0, kx, ky);
        auto compute_vely = std::bind(ComputeTaylorGreenVelocityY, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, t, td, u0, kx, ky);

        l2error_r = ComputeL2ErrorScalar(dens, compute_density);
        l2error_u = ComputeL2ErrorScalar(velx, compute_velx);
        l2error_v = ComputeL2ErrorScalar(velz, compute_vely);

        // Save to file.
        write_l2error << t << "," << l2error_r << "," << l2error_u << "," << l2error_v << std::endl;
    }

    // Close files.
    write_l2error.close();

    // Print averages after.
    std::string message_after = "After simulation:";
    PrintAverages(message_after, dens, velx, vely, velz);

    msg = "End";
    SaveMacroscopic(td, msg, dens, velx, velz);

    std::cout << td << "," << l2error_r << "," << l2error_u << "," << l2error_v << std::endl;

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

template <typename T>
void SaveMacroscopic(int timestep, std::string& message, 
             MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
             MacroscopicVariable<T>& velz)
{
    dens.WriteToTextFile(timestep);
    velx.WriteToTextFile(timestep);
    velz.WriteToTextFile(timestep);
    std::cout << message << ". Timestep = " << timestep << "\n";
}

template <typename T>
void PrintAverages(std::string& message, MacroscopicVariable<T>& dens, MacroscopicVariable<T>& velx,
                   MacroscopicVariable<T>& vely, MacroscopicVariable<T>& velz)
{
    std::cout << message << "\n";
    std::cout << "Average density = " << dens.ComputeAverage() << "\n";
    std::cout << "Average X velocity = " << velx.ComputeAverage() << "\n";
    std::cout << "Average Y velocity = " << vely.ComputeAverage() << "\n";
    std::cout << "Average Z velocity = " << velz.ComputeAverage() << "\n";
}
