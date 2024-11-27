#include <array>

#include "evolver/FluidEvolverSRT.h"
#include "equilibria/Equilibria.h"
#include "force/ForceSource.h"

//#define DEBUG false

template <typename T, int ND, int NQ> 
void FluidEvolverSRT<T, ND, NQ>::DoLocalCollision(std::function<T(T, T, T, T)> ComputeEquilibrium, AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k)
{
    // Calculate velocity magnitude squared, scaled by lattice speed of sound.
    T usq = (u_*u_ + v_*v_ + w_*w_) * f.CSI();

    // Calculate u dot F, scaled by lattice speed of sound.
    T uF = (u_*Fx_ + v_*Fy_ + w_*Fz_) * f.CSI();

    /*
    if (DEBUG) {
        std::cout << "Before collision: ";
        std::cout << "r = " << r_ << ", u = " << u_ << ", v = " << v_ << ", w = " << w_ << std::endl;
    } std::array<T, NQ> feqlocal; */

    // Do SRT collision.
    for (int q = 0; q < NQ; ++q)
    {
        // Compute c dot u, scaled by lattice speed of sound.
        T cu = (u_ * f.CX(q) + v_ * f.CY(q) + w_ * f.CZ(q)) * f.CSI();
        // Compute c dot F, scaled by lattice speed of sound.
        T cF = (Fx_ * f.CX(q) + Fy_ * f.CY(q) + Fz_ * f.CZ(q)) * f.CSI();

        // Compute equilibrium distribution.
        T feq = ComputeEquilibrium(r_, cu, usq, f.W(q));
        // Compute the force source term.
        T force_source = computeSecondOrderForceSource(cu, cF, uF, f.W(q));

        // Compute post-collision DF with SRT operator.
        T fstar = this->mOmega * feq + (static_cast<T>(1.0) - this->mOmega) * flocal[q] + (static_cast<T>(1.0) - static_cast<T>(0.5) * this->mOmega) * force_source;

        // Save the post-collision value.
        f.SetCurrFStar(fstar, q, i, j, k);

        //if (DEBUG) { feqlocal[q] = feq; }
    }

    /*
    if (DEBUG)
    {
        T r_after=0.0, u_after=0.0, v_after=0.0, w_after=0.0;
        for (int q = 0; q < NQ; ++q)
        {
            r_after += feqlocal[q];
            u_after += feqlocal[q] * f.CX(q);
            v_after += feqlocal[q] * f.CY(q);
            w_after += feqlocal[q] * f.CZ(q);
        }

        std::cout << "After collision:  ";
        std::cout << "r = " << r_after << ", u = " << u_after;
        std::cout << ", v = " << v_after << ", w = " << w_after << std::endl;
    }
    */
}