#include "AbstractFluidEvolver.h"
#include "finitedifference/FiniteDifference.h"

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::SetKinematicViscosity(AbstractLattice<T, ND, NQ>& f, T nu)
{
    T tau = f.CSI() * nu + 0.5;
    mOmega = 1/tau;
    std::cout << "Fluid evolver tau = " << tau << ", omega = " << mOmega << "\n";
}

inline double mat_idx(int a, int b) { return 3*a + b; }

/**
 * @brief f^eq + f^neq initialisation. 
 * Calculates velocity derivatives using centred finite differences on fluid nodes
 * and forward or backwards finiet differences on boundary nodes.
 * Note that initialised distributions correspond to pre-collision populations.
 */
template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::Initialise(AbstractLattice<T, ND, NQ>& f,
                                    const MacroscopicVariable<T>& dens,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
                                    const MacroscopicVariable<T>& Fx,
                                    const MacroscopicVariable<T>& Fy,
                                    const MacroscopicVariable<T>& Fz,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T, ND, NQ>& bdry [[maybe_unused]])
{   
    std::array<double, 9> DU;
    std::array<double, 9> Qi;
    std::array<double, 9> uFFu;
    std::array<double, 3> vec_u;
    std::array<double, 3> vec_F;
    std::array<int, 3> DX = {1, 0, 0};
    std::array<int, 3> DY = {0, 1, 0};
    std::array<int, 3> DZ = {0, 0, 1};
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                if (!node.IsSolid(i, j, k) && !node.IsGas(i, j, k))
                {
                    // Grab values from arrays.
                    T r_ = dens.GetValue(i, j, k);
                    T u_ = velx.GetValue(i, j, k);
                    T v_ = vely.GetValue(i, j, k);
                    T w_ = velz.GetValue(i, j, k);
                    T Fx_ = Fx.GetValue(i, j, k);
                    T Fy_ = Fy.GetValue(i, j, k);
                    T Fz_ = Fz.GetValue(i, j, k);

                    // Values for equilibrium are different because of force term!
                    T u_eq = u_ - 0.5 * Fx_ / r_;
                    T v_eq = v_ - 0.5 * Fy_ / r_;
                    T w_eq = w_ - 0.5 * Fz_ / r_;
                    double usq = (u_eq*u_eq + v_eq*v_eq + w_eq*w_eq) * f.CSI();

                    // Finite differences.
                    if (node.IsFluid(i, j, k))
                    {
                        for (int a = 0; a < 3; ++a)
                        {
                            // Centred difference for velx.
                            DU[mat_idx(a, 0)] = ComputeCentredDiff(velx, i, j, k, DX[a], DY[a], DZ[a]);
                            // Centred difference for vely.
                            DU[mat_idx(a, 1)] = ComputeCentredDiff(vely, i, j, k, DX[a], DY[a], DZ[a]);
                            // Centred difference for velz.
                            DU[mat_idx(a, 2)] = ComputeCentredDiff(velz, i, j, k, DX[a], DY[a], DZ[a]);
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
                                // Centred difference for velx.
                                DU[mat_idx(a, 0)] = ComputeCentredDiff(velx, i, j, k, DX[a], DY[a], DZ[a]);
                                // Centred difference for vely.
                                DU[mat_idx(a, 1)] = ComputeCentredDiff(vely, i, j, k, DX[a], DY[a], DZ[a]);
                                // Centred difference for velz.
                                DU[mat_idx(a, 2)] = ComputeCentredDiff(velz, i, j, k, DX[a], DY[a], DZ[a]);
                            }
                            else if (is_fluid_forward && !is_fluid_backward)
                            {
                                // Forward difference for velx.
                                DU[mat_idx(a, 0)] = ComputeForwardDiff(velx, i, j, k, DX[a], DY[a], DZ[a]);
                                // Forward difference for vely.
                                DU[mat_idx(a, 1)] = ComputeForwardDiff(vely, i, j, k, DX[a], DY[a], DZ[a]);
                                // Forward difference for velz.
                                DU[mat_idx(a, 2)] = ComputeForwardDiff(velz, i, j, k, DX[a], DY[a], DZ[a]);
                            }
                            else if (!is_fluid_forward && is_fluid_backward)
                            {
                                // Backward difference for velx.
                                DU[mat_idx(a, 0)] = ComputeBackwardDiff(velx, i, j, k, DX[a], DY[a], DZ[a]);
                                // Backward difference for vely.
                                DU[mat_idx(a, 1)] = ComputeBackwardDiff(vely, i, j, k, DX[a], DY[a], DZ[a]);
                                // Backward difference for velz.
                                DU[mat_idx(a, 2)] = ComputeBackwardDiff(velz, i, j, k, DX[a], DY[a], DZ[a]);
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
                    T rstar_ = 0;
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
                        double feq = computeSecondOrderEquilibrium(r_, cu_eq, usq, f.W(q));

                        // Compute non-equilibrium part.
                        double tau = 1.0 / mOmega;
                        double fneq = -f.W(q) * tau * r_ * QDU - 0.5 * f.W(q) * (cF + 0.5*QuFFu);

                        // Initialised distribution is sum of equilibrium and non-equlibrium parts.
                        double finit = feq + fneq;
                        f.SetCurrF(static_cast<T>(finit), q, i, j, k);
                        
                        rstar_ += finit;
                    }
                    std::cout << "before, after\n";
                    std::cout << r_ << ", " << rstar_ << std::endl;
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
void AbstractFluidEvolver<T, ND, NQ>::InitialiseEquilibrium(AbstractLattice<T, ND, NQ>& f,
                                    const MacroscopicVariable<T>& dens,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
                                    const MacroscopicVariable<T>& Fx,
                                    const MacroscopicVariable<T>& Fy,
                                    const MacroscopicVariable<T>& Fz,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T, ND, NQ>& bdry [[maybe_unused]])
{
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                if (!node.IsSolid(i, j, k) && !node.IsGas(i, j, k))
                {
                    // Grab values from arrays.
                    T r_ = dens.GetValue(i, j, k);
                    T u_ = velx.GetValue(i, j, k);
                    T v_ = vely.GetValue(i, j, k);
                    T w_ = velz.GetValue(i, j, k);
                    T Fx_ = Fx.GetValue(i, j, k);
                    T Fy_ = Fy.GetValue(i, j, k);
                    T Fz_ = Fz.GetValue(i, j, k);

                    // Values for equilibrium are different because of force term!
                    T u_eq = u_ - 0.5 * Fx_ / r_;
                    T v_eq = v_ - 0.5 * Fy_ / r_;
                    T w_eq = w_ - 0.5 * Fz_ / r_;
                    T usq = (u_eq*u_eq + v_eq*v_eq + w_eq*w_eq) * f.CSI();

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
 * @brief Wei's consistent initialisation scheme.
 * Carries out lattice boltzmann timesteps with fixed velocity field
 * until the density (aka pressure) field converges. 
 * Note that initialised distributions correspond to pre-collision populations.
 */
template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::InitialiseWei(AbstractLattice<T, ND, NQ>& f,
                                    const MacroscopicVariable<T>& dens,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
                                    const MacroscopicVariable<T>& Fx,
                                    const MacroscopicVariable<T>& Fy,
                                    const MacroscopicVariable<T>& Fz,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T, ND, NQ>& bdry [[maybe_unused]])
{
    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                if (!node.IsSolid(i, j, k) && !node.IsGas(i, j, k))
                {
                    // Grab values from arrays.
                    T r_ = dens.GetValue(i, j, k);
                    T u_ = velx.GetValue(i, j, k);
                    T v_ = vely.GetValue(i, j, k);
                    T w_ = velz.GetValue(i, j, k);
                    T Fx_ = Fx.GetValue(i, j, k);
                    T Fy_ = Fy.GetValue(i, j, k);
                    T Fz_ = Fz.GetValue(i, j, k);
                    //std::cout << r_ << ", " << u_ << ", " << v_ << ", " << w_ << ", " << Fx_ << ", " << Fy_ << ", " << Fz_ << std::endl;

                    // Values for equilibrium are different because of force term!
                    T u_eq = u_ - 0.5 * Fx_ / r_;
                    T v_eq = v_ - 0.5 * Fy_ / r_;
                    T w_eq = w_ - 0.5 * Fz_ / r_;
                    T usq = (u_eq*u_eq + v_eq*v_eq + w_eq*w_eq) * f.CSI();

                    // Check that moments sum correctly.
                    //T rstar_ = 0;
                    for (int q = 0; q < NQ; ++q)
                    {
                        T cu = (u_eq * f.CX(q) + v_eq * f.CY(q) + w_eq * f.CZ(q)) * f.CSI();
                        T fstar = computeSecondOrderEquilibrium(r_, cu, usq, f.W(q));
                        f.SetCurrF(fstar, q, i, j, k);
                        //rstar_ += fstar;
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

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::DoTimestep(AbstractLattice<T, ND, NQ>& f,
                                    MacroscopicVariable<T>& dens,
                                    MacroscopicVariable<T>& velx,
                                    MacroscopicVariable<T>& vely,
                                    MacroscopicVariable<T>& velz,
                                    const MacroscopicVariable<T>& Fx,
                                    const MacroscopicVariable<T>& Fy,
                                    const MacroscopicVariable<T>& Fz,
                                    NodeInfo& node,
                                    BoundaryInfo<T, ND, NQ>& bdry)
{
    // Stream.
    f.StreamDistributions();

    // Collide.
    std::array<T, NQ> flocal;

    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                // Initialise local density and velocity.
                T r_ = 0, u_ = 0, v_ = 0, w_ = 0; 
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
                    continue;
                }
                else
                {
                    // do nothing if it's a solid or gas node.
                    continue;
                }
                // Local collision step.

                // Grab local force.
                T Fx_ = Fx.GetValue(i, j, k);
                T Fy_ = Fy.GetValue(i, j, k);
                T Fz_ = Fz.GetValue(i, j, k);

                // Calculate local velocity (Guo's 2nd order scheme).
                u_ += 0.5 * Fx_;
                u_ /= r_;
                v_ += 0.5 * Fy_;
                v_ /= r_;
                w_ += 0.5 * Fz_;
                w_ /= r_;

                // Save local moments to array.
                dens.SetValue(r_, i, j, k);
                velx.SetValue(u_, i, j, k);
                vely.SetValue(v_, i, j, k);
                velz.SetValue(w_, i, j, k);

                // Do collision.
                DoLocalCollision(f, flocal, r_, u_, v_, w_, Fx_, Fy_, Fz_, i, j, k);
            }
        }
    }
}