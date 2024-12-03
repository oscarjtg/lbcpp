#include "ScalarEvolverTRT.h"
#include "equilibria/Equilibria.h"

template <typename T, int ND, int NQ>
void ScalarEvolverTRT<T, ND, NQ>::SetScalarDiffusivity(AbstractLattice<T, ND, NQ>& g [[maybe_unused]], T kappa)
{
    T tau_plus = g.CSI() * kappa + 0.5;
    mOmegaPlus = static_cast<T>(1.0) / tau_plus;
    T tau_minus = mMagicParameter / (tau_plus - 0.5) + 0.5;
    mOmegaMinus = static_cast<T>(1.0) / tau_minus;
    // For initialisation, which uses mOmega sometimes.
    this->mOmega = mOmegaPlus;
    std::cout << "Scalar evolver parameters:\n";
    std::cout << "Magic parameter = " << mMagicParameter << "\n";
    std::cout << "tau+ = " << tau_plus << ", omega+ = " << mOmegaPlus << "\n";
    std::cout << "tau- = " << tau_minus << ", omega- = " << mOmegaMinus << "\n";
    std::cout << std::flush;
}

/*
template <typename T, int ND, int NQ>
void ScalarEvolverTRT<T, ND, NQ>::DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k)
{
    // Calculate velocity magnitude squared, scaled by lattice speed of sound.
    T usq = (u_*u_ + v_*v_ + w_*w_) * g.CSI();

    // Do TRT collision.
    for (int q = 0; q < NQ; ++q)
    {
        if (q > g.QRev(q))
        {
            continue; // skip loop.
        }

        // Compute c dot u, scaled by lattice speed of sound.
        T cu = (u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q)) * g.CSI();

        // Compute equilibrium distribution.
        T geq = computeSecondOrderEquilibrium(c_, cu, usq, g.W(q));
        T geq_bar = computeSecondOrderEquilibrium(c_, -cu, usq, g.W(q));

        T geq_plus = 0.5 * (geq + geq_bar);
        T geq_minus = 0.5 * (geq - geq_bar);

        // Plus and minus g
        T g_plus = 0.5 * (glocal[q] + glocal[g.QRev(q)]);
        T g_minus = 0.5 * (glocal[q] - glocal[g.QRev(q)]);

        // Compute post-collision DF with TRT operator.
        T gstar = (
            glocal[q] 
            + mOmegaPlus * (geq_plus - g_plus) 
            + mOmegaMinus * (geq_minus - g_minus)
        );

        // Save the post-collision value.
        g.SetCurrFStar(gstar, q, i, j, k);

        if (q == g.QRev(q))
        {
            continue;
        }

        // Compute post-collision DF with TRT operator.
        T gstar_bar = (
            glocal[g.QRev(q)] 
            + mOmegaPlus * (geq_plus - g_plus) 
            - mOmegaMinus * (geq_minus - g_minus)
        );

        // Save the post-collision value.
        g.SetCurrFStar(gstar_bar, g.QRev(q), i, j, k);
    }
}*/

template <typename T, int ND, int NQ>
void ScalarEvolverTRT<T, ND, NQ>::DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k)
{
    // Calculate velocity magnitude squared, scaled by lattice speed of sound.
    T usq = (u_*u_ + v_*v_ + w_*w_) * g.CSI();

    // Do TRT collision.
    for (int q = 0; q < NQ; ++q)
    {
        // Compute c dot u, scaled by lattice speed of sound.
        T cu = (u_ * g.CX(q) + v_ * g.CY(q) + w_ * g.CZ(q)) * g.CSI();

        // Compute equilibrium distribution.
        T geq = computeSecondOrderEquilibrium(c_, cu, usq, g.W(q));
        T geq_bar = computeSecondOrderEquilibrium(c_, -cu, usq, g.W(q));

        T geq_plus = 0.5 * (geq + geq_bar);
        T geq_minus = 0.5 * (geq - geq_bar);

        // Plus and minus g
        T g_plus = 0.5 * (glocal[q] + glocal[g.QRev(q)]);
        T g_minus = 0.5 * (glocal[q] - glocal[g.QRev(q)]);

        // Compute post-collision DF with TRT operator.
        T gstar = (
            glocal[q] 
            + mOmegaPlus * (geq_plus - g_plus) 
            + mOmegaMinus * (geq_minus - g_minus)
        );

        // Save the post-collision value.
        g.SetCurrFStar(gstar, q, i, j, k);
    }
}