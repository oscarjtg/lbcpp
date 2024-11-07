#include <array>

#include "evolver/ScalarEvolverSRT1.h"
#include "equilibria/Equilibria.h"

template <class T>
void ScalarEvolverSRT1<T>::SetScalarDiffusivity(AbstractLattice<T>& g, T kappa)
{
    T tau = g.CSI() * kappa + 0.5;
    this->mOmega = 1 / tau;
}

template <class T>
void ScalarEvolverSRT1<T>::Initialise(AbstractLattice<T>& g,
                                    const MacroscopicVariable<T>& conc,
                                    const MacroscopicVariable<T>& velx,
                                    const MacroscopicVariable<T>& vely,
                                    const MacroscopicVariable<T>& velz,
                                    const NodeInfo& node,
                                    const BoundaryInfo<T>& bdry [[maybe_unused]])
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
                        T gstar = computeFirstOrderEquilibrium(c_, vel_projection, vel_squared, g.W(q));
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
void ScalarEvolverSRT1<T>::DoTimestep(AbstractLattice<T>& g,
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
                T c_ = 0; // Initialise local concentration.
                if (node.IsFluid(i, j, k))
                {
                    // Fluid streaming step.

                    // Grab local distribution function data
                    // and add to local concentration.
                    for (int q = 0; q < nq; ++q)
                    {
                        T g_ = g.GetCurrF(q, i, j, k);
                        c_ += g_;
                        glocal[q] = g_;
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
                        T g_ = bdry.ComputeBdryRule(bdry_id, q, i, j, k);
                        c_ += g_;
                        glocal[q] = g_;
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

                // Save local concentration to array.
                conc.SetValue(c_, i, j, k);

                // Grab local velocity.
                T u_ = velx.GetValue(i, j, k);
                T v_ = vely.GetValue(i, j, k);
                T w_ = velz.GetValue(i, j, k);

                // Calculate velocity magnitude squared, scaled by lattice speed of sound.
                T usq = (u_*u_ + v_*v_ + w_*w_) * g.CSI();

                // Do SRT collision.
                for (int q = 0; q < nq; ++q)
                {
                    // Compute velocity projection, scaled by lattice speed of sound.
                    T cu = (u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q)) * g.CSI();
                    // Compute equilibrium distribution.
                    T geq = computeFirstOrderEquilibrium(c_, cu, usq, g.W(q));
                    // Compute post-collision DF with SRT operator.
                    T gstar = mOmega * geq + (static_cast<T>(1.0) - mOmega) * glocal[q];

                    // Save the post-collision value.
                    g.SetCurrFStar(gstar, q, i, j, k);
                }
            }
        }
    }
}