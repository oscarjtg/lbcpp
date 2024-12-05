#include "AbstractFluidEvolver.h"
#include "finitedifference/FiniteDifference.h"
#include "diagnostics/L2Error.h"
#include "equilibria/Equilibria.h"

#include <type_traits>
#include <stdexcept>

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::SetKinematicViscosity(AbstractLattice<T, ND, NQ>& f, T nu)
{
    T tau = f.CSI() * nu + 0.5;
    mOmega = 1/tau;
    std::cout << "Fluid evolver tau = " << tau << ", omega = " << mOmega << "\n";
}

inline int mat_idx(int a, int b) { return 3*a + b; }

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::Initialise(AbstractLattice<T, ND, NQ>& f,
                            ScalarField<T>& dens,
                            const VectorField<T>& vel,
                            const VectorField<T>& force,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry, 
                            const std::string method)
{
    if (method == "CEQ")
    {
        dens.SetToConstantValue(1.0);
        InitialiseFEQ(f, dens, vel, force, node, bdry);
    }
    else if (method == "FEQ")
    {
        InitialiseFEQ(f, dens, vel, force, node, bdry);
    }
    else if (method == "NEQ")
    {
        InitialiseNEQ(f, dens, vel, force, node, bdry);
    }
    else if (method == "MEI")
    {
        InitialiseMEI(f, dens, vel, force, node, bdry);
    }
    else
    {
        throw std::invalid_argument("Argument must be CEQ, FEQ, NEQ, or MEI.");
    }
}

/**
 * @brief f^eq + f^neq initialisation. 
 * Calculates velocity derivatives using centred finite differences on fluid nodes
 * and forward or backwards finiet differences on boundary nodes.
 * Note that initialised distributions correspond to pre-collision populations.
 * 
 * For greater accuracy, compute everything in doubles.
 */
template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::InitialiseNEQ(AbstractLattice<T, ND, NQ>& f,
                                    const ScalarField<T>& dens,
                                    const VectorField<T>& vel,
                                    const VectorField<T>& force,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T, ND, NQ>& bdry [[maybe_unused]])
{   
    std::cout << "Non-equilibrium initialisation (NEQ)\n";
    std::array<double, 9> DU;
    std::array<double, 9> Qi;
    std::array<double, 9> uFFu;
    std::array<double, 3> vec_u;
    std::array<double, 3> vec_F;
    static constexpr std::array<int, 3> DX = {1, 0, 0};
    static constexpr std::array<int, 3> DY = {0, 1, 0};
    static constexpr std::array<int, 3> DZ = {0, 0, 1};
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                if (!node.IsSolid(i, j, k) && !node.IsGas(i, j, k))
                {
                    // Grab values from arrays.
                    double r_  = dens.GetValue(i, j, k);
                    double u_  = vel.GetValue(0, i, j, k);
                    double v_  = vel.GetValue(1, i, j, k);
                    double w_  = vel.GetValue(2, i, j, k);
                    double Fx_ = force.GetValue(0, i, j, k);
                    double Fy_ = force.GetValue(1, i, j, k);
                    double Fz_ = force.GetValue(2, i, j, k);

                    // Values for equilibrium are different because of force term!
                    double u_eq = u_ - 0.5 * Fx_ / r_;
                    double v_eq = v_ - 0.5 * Fy_ / r_;
                    double w_eq = w_ - 0.5 * Fz_ / r_;
                    double usq  = (u_eq*u_eq + v_eq*v_eq + w_eq*w_eq) * f.CSI();

                    // Finite differences.
                    if (node.IsFluid(i, j, k))
                    {
                        for (int a = 0; a < 3; ++a)
                        {
                            for (int b = 0; b < 3; ++b)
                            {   // Centred difference for vel. b=0: x. b=1: y. b=2: z.
                                DU[mat_idx(a, b)] = ComputeCentredDiff(vel, b, i, j, k, DX[a], DY[a], DZ[a]);
                            }
                        }
                    }
                    else
                    {
                        for (int a = 0; a < 3; ++a)
                        {
                            bool is_fluid_forward  = node.IsFluid_NoWrapAllowed(i+DX[a], j+DY[a], k+DZ[a]);
                            bool is_fluid_backward = node.IsFluid_NoWrapAllowed(i-DX[a], j-DY[a], k-DZ[a]);
                            if (is_fluid_forward && is_fluid_backward)
                            {
                                for (int b = 0; b < 3; ++b)
                                {   // Centred difference.
                                    DU[mat_idx(a, b)] = ComputeCentredDiff(vel, b, i, j, k, DX[a], DY[a], DZ[a]);
                                }
                            }
                            else if (is_fluid_forward && !is_fluid_backward)
                            {
                                for (int b = 0; b < 3; ++b)
                                {// Forward difference.
                                    DU[mat_idx(a, 0)] = ComputeForwardDiff(vel, b, i, j, k, DX[a], DY[a], DZ[a]);
                                }
                            }
                            else if (!is_fluid_forward && is_fluid_backward)
                            {
                                for (int b = 0; b < 3; ++b)
                                {// Backward difference.
                                    DU[mat_idx(a, b)] = ComputeBackwardDiff(vel, b, i, j, k, DX[a], DY[a], DZ[a]);
                                }
                            }
                            else
                            {
                                DU[mat_idx(a, 0)] = 0.0;
                                DU[mat_idx(a, 1)] = 0.0;
                                DU[mat_idx(a, 2)] = 0.0;
                            }
                        }
                    }

                    // Compute uFFu matrix.
                    vec_u[0] = static_cast<double>(u_);
                    vec_u[1] = static_cast<double>(v_);
                    vec_u[2] = static_cast<double>(w_);
                    vec_F[0] = static_cast<double>(Fx_);
                    vec_F[1] = static_cast<double>(Fy_);
                    vec_F[2] = static_cast<double>(Fz_);
                    for (int a = 0; a < 3; ++a) // No need to unroll loops as initialisation only done once.
                    {
                        for (int b = 0; b < 3; ++b)
                        {
                            uFFu[mat_idx(a, b)] = vec_u[a] * vec_F[b] + vec_F[a] * vec_u[b]; 
                        }
                    }

                    // Check that moments sum correctly.
                    //T rstar_ = 0;
                    for (int q = 0; q < NQ; ++q)
                    {
                        // Compute Qi components, multiplied by 1/cs^2
                        Qi[mat_idx(0, 0)] = f.CX(q) * f.CX(q) * f.CSI() - 1.0;
                        Qi[mat_idx(1, 0)] = f.CY(q) * f.CX(q) * f.CSI();
                        Qi[mat_idx(2, 0)] = f.CZ(q) * f.CX(q) * f.CSI();
                        Qi[mat_idx(0, 1)] = f.CX(q) * f.CY(q) * f.CSI();
                        Qi[mat_idx(1, 1)] = f.CY(q) * f.CY(q) * f.CSI() - 1.0;
                        Qi[mat_idx(2, 1)] = f.CZ(q) * f.CY(q) * f.CSI();
                        Qi[mat_idx(0, 2)] = f.CX(q) * f.CZ(q) * f.CSI();
                        Qi[mat_idx(1, 2)] = f.CY(q) * f.CZ(q) * f.CSI();
                        Qi[mat_idx(2, 2)] = f.CZ(q) * f.CZ(q) * f.CSI() - 1.0;

                        // Compute tensor inner products 
                        // QDU = Qi(a,b)DU(a,b)
                        // QuFFU = Qi(a,b)uFFu(a,b).
                        // Recall Qi is already multiplied  by 1/cs^2.
                        double QDU = 0, QuFFu = 0;
                        for (int a = 0; a < 3; ++a)
                        {
                            for (int b = 0; b < 3; ++b)
                            {
                                QDU   += Qi[mat_idx(a, b)] * DU[mat_idx(a, b)];
                                QuFFu += Qi[mat_idx(a, b)] * uFFu[mat_idx(a, b)];
                            }
                        }

                        // Projection of force F on lattice velocity c, multiplied by 1/cs^2.
                        double cF = (Fx_ * f.CX(q) + Fy_ * f.CY(q) + Fz_ * f.CZ(q)) * f.CSI();

                        // Projection of equilibrium fluid velocity u_eq on lattice velocity c, multiplied by 1/cs^2.
                        double cu_eq = (u_eq * f.CX(q) + v_eq * f.CY(q) + w_eq * f.CZ(q)) * f.CSI();

                        // Compute equilibrium distribution function.
                        double feq = computeSecondOrderEquilibrium(r_, cu_eq, usq, static_cast<double>(f.W(q)));

                        // Compute non-equilibrium part.
                        double tau = 1.0 / mOmega;
                        double fneq = -f.W(q) * tau * r_ * QDU - 0.5 * f.W(q) * (cF + 0.5*QuFFu);

                        // Initialised distribution is sum of equilibrium and non-equlibrium parts.
                        double finit = feq + fneq;
                        f.SetCurrF(static_cast<T>(finit), q, i, j, k);
                        
                        //rstar_ += finit;
                    }
                    //std::cout << "before, after\n";
                    //std::cout << r_ << ", " << rstar_ << std::endl;
                }
                else
                {
                    for (int q = 0; q < NQ; ++q)
                    {
                        f.SetCurrF(0, q, i, j, k);
                    }
                }
            }
        }
    }
    f.StreamDistributions();
}

/**
 * @brief f^eq initialisation. Neglects components of deviatoric stress.
 * Note that initialised distributions correspond to pre-collision populations.
 */
template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::InitialiseFEQ(AbstractLattice<T, ND, NQ>& f,
                                    const ScalarField<T>& dens,
                                    const VectorField<T>& vel,
                                    const VectorField<T>& force,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T, ND, NQ>& bdry [[maybe_unused]])
{
    std::cout << "Equilibrium initialisation (FEQ)\n";
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                if (!node.IsSolid(i, j, k) && !node.IsGas(i, j, k))
                {
                    // Grab values from arrays.
                    T r_  = dens.GetValue(i, j, k);
                    T u_  = vel.GetValue(0, i, j, k);
                    T v_  = vel.GetValue(1, i, j, k);
                    T w_  = vel.GetValue(2, i, j, k);
                    T Fx_ = force.GetValue(0, i, j, k);
                    T Fy_ = force.GetValue(1, i, j, k);
                    T Fz_ = force.GetValue(2, i, j, k);

                    // Values for equilibrium are different because of force term!
                    T u_eq = u_ - 0.5 * Fx_ / r_;
                    T v_eq = v_ - 0.5 * Fy_ / r_;
                    T w_eq = w_ - 0.5 * Fz_ / r_;
                    T usq  = (u_eq*u_eq + v_eq*v_eq + w_eq*w_eq) * f.CSI();

                    // Calculate and set second order equilibrium distribution.
                    for (int q = 0; q < NQ; ++q)
                    {
                        T cu = (u_eq * f.CX(q) + v_eq * f.CY(q) + w_eq * f.CZ(q)) * f.CSI();
                        T fstar = computeSecondOrderEquilibrium(r_, cu, usq, f.W(q));
                        f.SetCurrF(fstar, q, i, j, k);
                    }
                }
                else // No distribution functions in solid and gas nodes.
                {
                    for (int q = 0; q < NQ; ++q)
                    {
                        f.SetCurrF(0, q, i, j, k);
                    }
                }
            }
        }
    }
    // Stream will be undone at beginning of DoTimestep(),
    // leaving the initialised values in the "pre-collision" positions.
    f.StreamDistributions();
}

/**
 * @brief Mei's consistent initialisation scheme.
 * Carries out lattice boltzmann timesteps with fixed velocity field
 * until the density (aka pressure) field converges. 
 * Note that initialised distributions correspond to pre-collision populations.
 */
template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::InitialiseMEI(AbstractLattice<T, ND, NQ>& f,
                                    ScalarField<T>& dens,
                                    const VectorField<T>& vel,
                                    const VectorField<T>& force,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T, ND, NQ>& bdry [[maybe_unused]])
{
    std::cout << "Mei's consistent initialisation scheme (Mei)\n";
    
    // Start from equilibrium initialisation.
    InitialiseFEQ(f, dens, vel, force, node, bdry);

    // Now employ Mei's iterative scheme until L2 error norm goes below tolerance.
    double tolerance;
    if (std::is_same<T, float>::value)
    {
        tolerance = 1.0e-7;
    }
    else if (std::is_same<T, double>::value)
    {
        tolerance = 1.0e-10;
    }
    else
    {
        tolerance = 1;
    }
    int count = 0;
    while (ComputeL2ErrorDistribution(f) > tolerance)
    {
        // Stream.
        f.StreamDistributions();

        for (int k = 0; k < node.GetNZ(); ++k)
        {
            for (int j = 0; j < node.GetNY(); ++j)
            {
                for (int i = 0; i < node.GetNX(); ++i)
                {
                    // Initialise local density (only density changes; velocity is fixed by initial condition).
                    T r_ = 0;

                    // Do timestep at this node (does not modify MacroscopicVariable arrays).
                    DoLocalTimestepIncompressible(f, vel, force, node, bdry, r_, i, j, k);

                    // Save only the density. Velocities remain unchanged during Mei's initialisation.
                    dens.SetValue(r_, i, j, k);
                }
            }
        }
        count++;
    }
    std::cout << "Completed Mei's initialization after " << count << " iterations.\n";
    f.StreamDistributions();
}

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::DoTimestep(AbstractLattice<T, ND, NQ>& f,
                                    ScalarField<T>& dens,
                                    VectorField<T>& vel,
                                    const VectorField<T>& force,
                                    NodeInfo& node,
                                    BoundaryInfo<T, ND, NQ>& bdry,
                                    const bool StoreMacros)
{
    // Stream.
    f.StreamDistributions();

    // Collide.
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                // Initialise local density and velocity.
                T r_ = 0, u_ = 0, v_ = 0, w_ = 0;

                // Do timestep at this node (does not modify MacroscopicVariable arrays).
                DoLocalTimestep(f, force, node, bdry, r_, u_, v_, w_, i, j, k);

                if (StoreMacros)
                {   // Save local macroscopic variables to arrays.
                    dens.SetValue(r_, i, j, k);
                    vel.SetValue(u_, 0, i, j, k);
                    vel.SetValue(v_, 1, i, j, k);
                    vel.SetValue(w_, 2, i, j, k);
                }
            }
        }
    }
}

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::DoTimestep(AbstractLattice<T, ND, NQ>& f,
                                    ScalarField<T>& dens,
                                    VectorField<T>& vel,
                                    TensorField<T>& sigma,
                                    const VectorField<T>& force,
                                    NodeInfo& node,
                                    BoundaryInfo<T, ND, NQ>& bdry,
                                    const bool StoreMacros)
{
    // Stream.
    f.StreamDistributions();

    // Initialise variables and arrays of local density, velocity, and stress
    T r_ = 0;
    static std::array<T, 3> u_vec = {0.0, 0.0, 0.0};
    static std::array<T, 9> s_mat = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Collide.
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                // Do timestep at this node (does not modify MacroscopicVariable arrays).
                DoLocalTimestep(f, force, node, bdry, r_, u_vec, s_mat, i, j, k);

                if (StoreMacros)
                {   // Save local macroscopic variables to arrays.
                    dens.SetValue(r_, i, j, k);
                    for (int a = 0; a < 3; ++a)
                    {
                        vel.SetValue(u_vec[a], a, i, j, k);
                    }
                    for (int b = 0; b < 3; ++b)
                    {
                        for (int a = 0; a < 3; ++a)
                        {   // note mat_idx and TensorField index in opposite ways.
                            sigma.SetValue(s_mat[mat_idx(a, b)], a, b, i, j, k); 
                        }
                    }
                }
            }
        }
    }
}

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::DoLocalTimestep(AbstractLattice<T, ND, NQ>& f,
                        const VectorField<T>& force,
                        NodeInfo& node,
                        BoundaryInfo<T, ND, NQ>& bdry, 
                        T& r_, T& u_, T& v_, T& w_, int i, int j, int k)
{
    static std::array<T, NQ> flocal;
    r_ = 0, u_ = 0, v_ = 0, w_ = 0;

    if (node.IsFluid(i, j, k))
    {
        // Fluid streaming step.

        // Grab local distribution function data
        // and add to local concentration.
        for (int q = 0; q < NQ; ++q)
        {
            T f_ = f.GetCurrF(q, i, j, k);
            r_ += f_;
            u_ += f_ * f.CX(q);
            v_ += f_ * f.CY(q);
            w_ += f_ * f.CZ(q);
            flocal[q] = f_;
        }
    }
    else if (node.IsBoundary(i, j, k))
    {
        // Boundary streaming step.
        //std::cout << "Boundary node (" << i << ", " << j << ", " << k << "):\n";

        // Work out which boundary condition rule you should be using.
        int bdry_id = node.GetBoundaryID(i, j, k);

        // Grab data according to boundary condition rule,
        // and add to local concentration.
        for (int q = 0; q < NQ; ++q)
        {
            T f_ = bdry.ComputeBdryRule(bdry_id, q, i, j, k);

            //std::cout << "f[" << q << "] = " << f_ << "\n";

            r_ += f_;
            u_ += f_ * f.CX(q);
            v_ += f_ * f.CY(q);
            w_ += f_ * f.CZ(q);
            flocal[q] = f_;
        }
    }
    else if (node.IsInterface(i, j, k))
    {
        // Interface LB step.
        // To be implemented.
        return;
    }
    else
    {
        // do nothing if it's a solid or gas node.
        return;
    }
    // Local collision step.

    // Grab local force.
    T Fx_ = force.GetValue(0, i, j, k);
    T Fy_ = force.GetValue(1, i, j, k);
    T Fz_ = force.GetValue(2, i, j, k);

    // Calculate local velocity (Guo's 2nd order scheme).
    u_ += 0.5 * Fx_;
    u_ /= r_;
    v_ += 0.5 * Fy_;
    v_ /= r_;
    w_ += 0.5 * Fz_;
    w_ /= r_;

    // Do collision.
    DoLocalCollision(computeSecondOrderEquilibrium<T>, f, flocal, r_, u_, v_, w_, Fx_, Fy_, Fz_, i, j, k);
}

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::DoLocalTimestep(AbstractLattice<T, ND, NQ>& f,
                        const VectorField<T>& force,
                        NodeInfo& node,
                        BoundaryInfo<T, ND, NQ>& bdry,
                        T& r_, std::array<T, 3>& u_vec, std::array<T, 9>& s_mat, int i, int j, int k)
{
    static std::array<T, NQ> flocal;
    static std::array<T, NQ> feq;
    static std::array<T, 3> F_vec;
    r_ = 0;
    T u_ = 0, v_ = 0, w_ = 0; // Setting these to zero is important!!

    if (node.IsFluid(i, j, k))
    {
        // Fluid streaming step.

        // Grab local distribution function data
        // and add to local concentration.
        for (int q = 0; q < NQ; ++q)
        {
            T f_ = f.GetCurrF(q, i, j, k);
            r_ += f_;
            u_ += f_ * f.CX(q);
            v_ += f_ * f.CY(q);
            w_ += f_ * f.CZ(q);
            flocal[q] = f_;
        }
    }
    else if (node.IsBoundary(i, j, k))
    {
        // Boundary streaming step.
        //std::cout << "Boundary node (" << i << ", " << j << ", " << k << "):\n";

        // Work out which boundary condition rule you should be using.
        int bdry_id = node.GetBoundaryID(i, j, k);

        // Grab data according to boundary condition rule,
        // and add to local concentration.
        for (int q = 0; q < NQ; ++q)
        {
            T f_ = bdry.ComputeBdryRule(bdry_id, q, i, j, k);

            //std::cout << "f[" << q << "] = " << f_ << "\n";

            r_ += f_;
            u_ += f_ * f.CX(q);
            v_ += f_ * f.CY(q);
            w_ += f_ * f.CZ(q);
            flocal[q] = f_;
        }
    }
    else if (node.IsInterface(i, j, k))
    {
        // Interface LB step.
        // To be implemented.
        return;
    }
    else
    {
        // do nothing if it's a solid or gas node.
        return;
    }
    // Local collision step.

    // Grab local force.
    T Fx_ = force.GetValue(0, i, j, k);
    T Fy_ = force.GetValue(1, i, j, k);
    T Fz_ = force.GetValue(2, i, j, k);

    // Calculate local velocity (Guo's 2nd order scheme).
    u_ += 0.5 * Fx_;
    u_ /= r_;
    v_ += 0.5 * Fy_;
    v_ /= r_;
    w_ += 0.5 * Fz_;
    w_ /= r_;

    // Calculate velocity magnitude squared, scaled by lattice speed of sound.
    T usq = (u_*u_ + v_*v_ + w_*w_) * f.CSI();

    // Compute feq.
    for (int q = 0; q < NQ; ++q)
    {
        // Compute c dot u, scaled by lattice speed of sound.
        T cu = (u_ * f.CX(q) + v_ * f.CY(q) + w_ * f.CZ(q)) * f.CSI();

        // Compute equilibrium distribution.
        feq[q] = computeSecondOrderEquilibrium(r_, cu, usq, f.W(q));
    }

    // Fill vectors.
    u_vec[0] = u_,  u_vec[1] = v_,  u_vec[2] = w_;
    F_vec[0] = Fx_, F_vec[1] = Fy_, F_vec[2] = Fz_;

    // Compute stress tensor. Must be computed using pre-collision distribution functions.
    for (int b = 0; b < 3; ++b)
    {
        for (int a = 0; a < 3; ++a)
        {
            double sum = 0;
            for (int q = 0; q < NQ; ++q)
            {
                T fneq = flocal[q] - feq[q];
                sum += static_cast<double>(fneq * f.GetCqa(q, a) * f.GetCqa(q, b));
            }
            double first_term = -(1 - 0.5*mOmega) * sum;
            double second_term = -0.5*(1 - 0.5*mOmega) * (F_vec[a]*u_vec[b] + u_vec[a]*F_vec[b]);
            s_mat[mat_idx(a, b)] = first_term + second_term;
        }
    }

    // Do collision.
    DoLocalCollision(computeSecondOrderEquilibrium<T>, f, flocal, r_, u_vec[0], u_vec[1], u_vec[2], F_vec[0], F_vec[1], F_vec[2], i, j, k);
}

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::DoLocalTimestepIncompressible(AbstractLattice<T, ND, NQ>& f,
                        const VectorField<T>& vel,
                        const VectorField<T>& force,
                        const NodeInfo& node,
                        const BoundaryInfo<T, ND, NQ>& bdry, 
                        T& r_, int i, int j, int k)
{
    static std::array<T, NQ> flocal;

    if (node.IsFluid(i, j, k))
    {
        // Fluid streaming step.

        // Grab local distribution function data
        // and add to local concentration.
        for (int q = 0; q < NQ; ++q)
        {
            T f_ = f.GetCurrF(q, i, j, k);
            r_ += f_;
            flocal[q] = f_;
        }
    }
    else if (node.IsBoundary(i, j, k))
    {
        // Boundary streaming step.
        //std::cout << "Boundary node (" << i << ", " << j << ", " << k << "):\n";

        // Work out which boundary condition rule you should be using.
        int bdry_id = node.GetBoundaryID(i, j, k);

        // Grab data according to boundary condition rule,
        // and add to local concentration.
        for (int q = 0; q < NQ; ++q)
        {
            T f_ = bdry.ComputeBdryRule(bdry_id, q, i, j, k);

            //std::cout << "f[" << q << "] = " << f_ << "\n";

            r_ += f_;
            flocal[q] = f_;
        }
    }
    else if (node.IsInterface(i, j, k))
    {
        // Interface LB step.
        // To be implemented.
        return;
    }
    else
    {
        // do nothing if it's a solid or gas node.
        return;
    }
    // Local collision step.

    // Grab velocity.
    T u_ = vel.GetValue(0, i, j, k);
    T v_ = vel.GetValue(1, i, j, k);
    T w_ = vel.GetValue(2, i, j, k);

    // Grab local force.
    T Fx_ = force.GetValue(0, i, j, k);
    T Fy_ = force.GetValue(1, i, j, k);
    T Fz_ = force.GetValue(2, i, j, k);

    // Calculate local velocity (Guo's 2nd order scheme).
    // Note we imposed the velocities, instead of computing from first moments of f, so this is not needed.
    // u_ += 0.5 * Fx_ / r_;
    // v_ += 0.5 * Fy_ / r_;
    // w_ += 0.5 * Fz_ / r_;

    // Do collision, with incompressible equilibrium.
    DoLocalCollision(computeSecondOrderIncompressibleEquilibrium<T>, f, flocal, r_, u_, v_, w_, Fx_, Fy_, Fz_, i, j, k);
}