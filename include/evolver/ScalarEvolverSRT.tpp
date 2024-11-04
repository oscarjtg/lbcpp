#include <array>

#include "evolver/ScalarEvolverSRT.h"
#include "equilibria/Equilibria.h"

template <class T>
void ScalarEvolverSRT<T>::SetScalarDiffusivity(T kappa)
{
    T tau = 3 * kappa + 0.5;
    this->mOmega = 1 / tau;
}

template <class T>
void ScalarEvolverSRT<T>::Initialise(AbstractLattice<T>& g,
                                    const MacroscopicVariable<T>& conc,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
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
                    T c_ = conc.GetValue(i, j, k);
                    T u_ = velx.GetValue(i, j, k);
                    T v_ = vely.GetValue(i, j, k);
                    T w_ = velz.GetValue(i, j, k);
                    T vel_squared = u_*u_ + v_*v_ + w_*w_;
                    for (int q = 0; q < g.NQ(); ++q)
                    {
                        T vel_projection = u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q);
                        T gstar = computeSecondOrderEquilibrium(c_, vel_projection, vel_squared, g.W(q));
                        g.SetCurrFStar(gstar, q, i, j, k);
                    }
                }
                else
                {
                    for (int q = 0; q < g.NQ(); ++ q)
                    {
                        g.SetCurrFStar(0, q, i, j, k);
                    }
                }
            }
        }
    }
}

template <class T>
void ScalarEvolverSRT<T>::DoTimestep(AbstractLattice<T>& g,
                                    MacroscopicVariable<T>& conc,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T>& bdry)
{
    // Stream.
    g.StreamDistributions();

    // Collide.
    const int nq = g.NQ();
    const int max_size = 27;
    std::array<T, max_size> glocal;

    for (int k = 0; k < node.GetNZ(); ++k)
    {
        for (int j = 0; j < node.GetNY(); ++j)
        {
            for (int i = 0; i < node.GetNX(); ++i)
            {
                // implementation here.
                if (node.IsFluid(i, j, k))
                {
                    // Fluid LB timestep.

                    // Grab local distribution function data
                    // and calculate local concentration.
                    T c_ = 0;
                    for (int q = 0; q < nq; ++q)
                    {
                        T g_ = g.GetCurrF(q, i, j, k);
                        c_ += g_;
                        glocal[q] = g_;
                    }
                    // Save local concentration to array.
                    conc.SetValue(c_, i, j, k);

                    // Grab local velocity.
                    T u_ = velx.GetValue(i, j, k);
                    T v_ = vely.GetValue(i, j, k);
                    T w_ = velz.GetValue(i, j, k);

                    // Calculate velocity magnitude squared.
                    T vel_squared = u_*u_ + v_*v_ + w_*w_;

                    // Do SRT collision.
                    for (int q = 0; q < nq; ++q)
                    {
                        T vel_projection = u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q);
                        T gstar = mOmega * computeSecondOrderEquilibrium(c_, vel_projection, vel_squared, g.W(q)) + (static_cast<T>(1.0) - mOmega) * glocal[q];
                        g.SetCurrFStar(gstar, q, i, j, k);
                    }

                }
                else if (node.IsBoundary(i, j, k))
                {
                    // Boundary LB timestep.

                    // Work out which boundary condition rule you should be using.
                    int bdry_id = node.GetBoundaryID(i, j, k);

                    // Grab data according to boundary condition rule,
                    // and calculate local concentration.
                    T c_ = 0;
                    for (int q = 0; q < nq; ++q)
                    {
                        T g_ = bdry.ComputeBdryRuleG(bdry_id, q, i, j, k);
                        c_ += g_;
                    }
                    // Save local concentration to array.
                    conc.SetValue(c_, i, j, k);

                    // Grab local velocity.
                    T u_ = velx.GetValue(i, j, k);
                    T v_ = vely.GetValue(i, j, k);
                    T w_ = velz.GetValue(i, j, k);

                    // Calculate velocity magnitude squared.
                    T vel_squared = u_*u_ + v_*v_ + w_*w_;

                    // Do SRT collision.
                    for (int q = 0; q < nq; ++q)
                    {
                        T vel_projection = u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q);
                        T gstar = mOmega * computeSecondOrderEquilibrium(c_, vel_projection, vel_squared, g.W(q)) + (static_cast<T>(1.0) - mOmega) * glocal[q];
                        g.SetCurrFStar(gstar, q, i, j, k);
                    }
                }
                else if (node.IsInterface(i, j, k))
                {
                    // Interface LB timestep.
                }
                else
                {
                    // do nothing if it's a solid node.
                }
            }
        }
    }
}