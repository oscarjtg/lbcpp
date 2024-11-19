#include <array>

#include "evolver/ScalarEvolverSRT.h"
#include "equilibria/Equilibria.h"

template <typename T, int ND, int NQ>
void ScalarEvolverSRT<T, ND, NQ>::DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k)
{
    // Calculate velocity magnitude squared, scaled by lattice speed of sound.
    T usq = (u_*u_ + v_*v_ + w_*w_) * g.CSI();

    // Do SRT collision.
    for (int q = 0; q < NQ; ++q)
    {
        // Compute velocity projection, scaled by lattice speed of sound.
        T cu = (u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q)) * g.CSI();
        // Compute equilibrium distribution.
        T geq = computeSecondOrderEquilibrium(c_, cu, usq, g.W(q));
        // Compute post-collision DF with SRT operator.
        T gstar = this->mOmega * geq + (static_cast<T>(1.0) - this->mOmega) * glocal[q];

        // Save the post-collision value.
        g.SetCurrFStar(gstar, q, i, j, k);
    }
}