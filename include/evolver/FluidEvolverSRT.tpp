#include <array>

#include "evolver/FluidEvolverSRT.h"
#include "equilibria/Equilibria.h"
#include "force/ForceSource.h"

template <class T>
void FluidEvolverSRT<T>::SetKinematicViscosity(T nu)
{
    T tau = 3 * nu + 0.5;
    this->mOmega = 1/tau;
}

template <class T>
void FluidEvolverSRT<T>::SetBulkViscosity(T eta)
{
    std::cout << "Bulk viscosity eta = ";
    std::cout << eta << "cannot be specified in SRT model.\n";
}

template <class T>
void FluidEvolverSRT<T>::Initialise(AbstractLattice<T>& f,
                                    const MacroscopicVariable<T>& dens,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
                                    const MacroscopicVariable<T>& Fx,
                                    const MacroscopicVariable<T>& Fy,
                                    const MacroscopicVariable<T>& Fz,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T>& bdry)
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
                    for (int q = 0; q < f.NQ(); ++q)
                    {
                        T vel_projection = u_eq * f.CX(q) + v_eq * f.CY(q) + w_eq * f.CZ(q);
                        T fstar = computeSecondOrderEquilibrium(r_, vel_projection, vel_squared, f.W(q));
                        f.SetCurrFStar(fstar, q, i, j, k);
                    }
                }
                else
                {
                    for (int q = 0; q < f.NQ(); ++ q)
                    {
                        f.SetCurrFStar(0, q, i, j, k);
                    }
                }
            }
        }
    }
}

template <class T>
void FluidEvolverSRT<T>::DoTimestep(AbstractLattice<T>& f,
                                    MacroscopicVariable<T>& dens,
                                    MacroscopicVariable<T>& velx,
                                    MacroscopicVariable<T>& vely,
                                    MacroscopicVariable<T>& velz,
                                    const MacroscopicVariable<T>& Fx,
                                    const MacroscopicVariable<T>& Fy,
                                    const MacroscopicVariable<T>& Fz,
                                    NodeInfo& node,
                                    BoundaryInfo<T>& bdry)
{
    // Stream.
    f.StreamDistributions();

    // Collide.
    const int nq = f.NQ();
    const int max_size = 27;
    std::array<T, max_size> flocal;

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
                    for (int q = 0; q < nq; ++q)
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

                    // Work out which boundary condition rule you should be using.
                    int bdry_id = node.GetBoundaryID(i, j, k);

                    // Grab data according to boundary condition rule,
                    // and add to local concentration.
                    for (int q = 0; q < nq; ++q)
                    {
                        T f_ = bdry.ComputeBdryRuleF(bdry_id, q, i, j, k);
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
                T Fx_ = velx.GetValue(i, j, k);
                T Fy_ = vely.GetValue(i, j, k);
                T Fz_ = velz.GetValue(i, j, k);

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

                // Calculate velocity magnitude squared.
                T usq = u_*u_ + v_*v_ + w_*w_;

                // Calculate u dot F
                T uF = u_*Fx_ + v_*Fy_ + w_*Fz_;

                // Do SRT collision.
                for (int q = 0; q < nq; ++q)
                {
                    // Compute c dot u.
                    T cu = u_ * f.CX(q) + v_ * f.CY(q) + w_ * f.CZ(q);
                    // Compute c dot F.
                    T cF = Fx_ * f.CX(q) + Fy_ * f.CY(q) + Fz_ * f.CZ(q);

                    // Compute equilibrium distribution.
                    T feq = computeSecondOrderEquilibrium(r_, cu, usq, f.W(q));
                    // Compute the force source term.
                    T force_source = computeSecondOrderForceSource(cu, cF, uF, f.W(q));

                    // Compute post-collision DF with SRT operator.
                    T fstar = mOmega * feq + (static_cast<T>(1.0) - mOmega) * flocal[q] + (static_cast<T>(1.0) - static_cast<T>(0.5) * mOmega) * force_source;

                    // Save the post-collision value.
                    f.SetCurrFStar(fstar, q, i, j, k);
                }
            }
        }
    }
}