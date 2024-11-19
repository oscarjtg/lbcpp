#include "AbstractFluidEvolver.h"

template <typename T, int ND, int NQ>
void AbstractFluidEvolver<T, ND, NQ>::SetKinematicViscosity(AbstractLattice<T, ND, NQ>& f, T nu)
{
    T tau = f.CSI() * nu + 0.5;
    this->mOmega = 1/tau;
}

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
                    T vel_squared = u_eq*u_eq + v_eq*v_eq + w_eq*w_eq;

                    // Check that moments sum correctly.
                    T rstar_ = 0;
                    for (int q = 0; q < NQ; ++q)
                    {
                        T vel_projection = u_eq * f.CX(q) + v_eq * f.CY(q) + w_eq * f.CZ(q);
                        T fstar = computeSecondOrderEquilibrium(r_, vel_projection, vel_squared, f.W(q));
                        f.SetCurrFStar(fstar, q, i, j, k);
                        rstar_ += fstar;
                    }
                }
                else
                {
                    for (int q = 0; q < NQ; ++q)
                    {
                        f.SetCurrFStar(0, q, i, j, k);
                    }
                }
            }
        }
    }
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