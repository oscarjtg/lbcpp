#include "evolver/FluidEvolverTRT.h"
#include "equilibria/Equilibria.h"
#include "force/ForceSource.h"

template <typename T, int ND, int NQ>
void FluidEvolverTRT<T, ND, NQ>::SetKinematicViscosity(AbstractLattice<T, ND, NQ>& f, T nu)
{
    T tau_plus = f.CSI() * nu + 0.5;
    mOmegaPlus = static_cast<T>(1.0) / tau_plus;
    T tau_minus = mMagicParameter / (tau_plus - 0.5) + 0.5;
    mOmegaMinus = static_cast<T>(1.0) / tau_minus;
    // For initialisation, which uses mOmega sometimes.
    this->mOmega = mOmegaPlus;
    std::cout << "Fluid evolver parameters:\n";
    std::cout << "Magic parameter = " << mMagicParameter << "\n";
    std::cout << "tau+ = " << tau_plus << ", omega+ = " << mOmegaPlus << "\n";
    std::cout << "tau- = " << tau_minus << ", omega- = " << mOmegaMinus << "\n";
    std::cout << std::flush;
}

template <typename T, int ND, int NQ>
void FluidEvolverTRT<T, ND, NQ>::DoLocalCollision(std::function<T(T, T, T, T)> ComputeEquilibrium, AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k)
{
    // Calculate velocity magnitude squared, scaled by lattice speed of sound.
    T usq = (u_*u_ + v_*v_ + w_*w_) * f.CSI();

    // Calculate u dot F, scaled by lattice speed of sound.
    T uF = (u_*Fx_ + v_*Fy_ + w_*Fz_) * f.CSI();

    // Do TRT collision.
    for (int q = 0; q < NQ; ++q)
    {
        if (q > f.QRev(q))
        {
            continue; // skip loop.
        }
        // Compute c dot u, scaled by lattice speed of sound.
        T cu = (u_ * f.CX(q) + v_ * f.CY(q) + w_ * f.CZ(q)) * f.CSI();
        // Compute c dot F, scaled by lattice speed of sound.
        T cF = (Fx_ * f.CX(q) + Fy_ * f.CY(q) + Fz_ * f.CZ(q)) * f.CSI();

        // Compute equilibrium distribution.
        T feq = ComputeEquilibrium(r_, cu, usq, f.W(q));
        T feq_bar = ComputeEquilibrium(r_, -cu, usq, f.W(q));

        T feq_plus = 0.5 * (feq + feq_bar);
        T feq_minus = 0.5 * (feq - feq_bar);

        // Compute the force source term.
        T force_source = computeSecondOrderForceSource(cu, cF, uF, f.W(q));
        T force_source_bar = computeSecondOrderForceSource(-cu, -cF, uF, f.W(q));

        T force_source_plus = 0.5 * (force_source + force_source_bar);
        T force_source_minus = 0.5 * (force_source - force_source_bar);

        // Plus and minus f
        T f_plus = 0.5 * (flocal[q] + flocal[f.QRev(q)]);
        T f_minus = 0.5 * (flocal[q] - flocal[f.QRev(q)]);

        // Compute post-collision DF with TRT operator.
        T fstar = (
            flocal[q] 
            + mOmegaPlus * (feq_plus - f_plus) 
            + mOmegaMinus * (feq_minus - f_minus)
            + (static_cast<T>(1.0) - static_cast<T>(0.5) * mOmegaPlus) * force_source_plus
            + (static_cast<T>(1.0) - static_cast<T>(0.5) * mOmegaMinus) * force_source_minus
        );

        // Save the post-collision value.
        f.SetCurrFStar(fstar, q, i, j, k);

        if (q == f.QRev(q))
        {
            continue;
        }

        // Compute post-collision for fbar with TRT operator.
        T fstar_bar = (
            flocal[f.QRev(q)] 
            + mOmegaPlus * (feq_plus - f_plus) 
            - mOmegaMinus * (feq_minus - f_minus)
            + (static_cast<T>(1.0) - static_cast<T>(0.5) * mOmegaPlus) * force_source_plus
            - (static_cast<T>(1.0) - static_cast<T>(0.5) * mOmegaMinus) * force_source_minus
        );

        // Save the post-collision value.
        f.SetCurrFStar(fstar_bar, f.QRev(q), i, j, k);
    }
}