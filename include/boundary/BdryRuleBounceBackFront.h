#ifndef BDRYRULEBOUNCEBACKFRONTDEF
#define BDRYRULEBOUNCEBACKFRONTDEF

#include "boundary/AbstractBdryRuleBounceBack.h"

template <typename T, int ND, int NQ>
class BdryRuleBounceBackFront : public AbstractBdryRuleBounceBack<T, ND, NQ>
{
public:
    BdryRuleBounceBackFront() = default;

    BdryRuleBounceBackFront(AbstractLattice<T, ND, NQ>* pDistribution)
        : AbstractBdryRuleBounceBack<T, ND, NQ>(pDistribution) {}

    BdryRuleBounceBackFront(AbstractLattice<T, ND, NQ>* pDistribution, ScalarField<T>* pDensity, T velx, T vely, T velz)
        : AbstractBdryRuleBounceBack<T, ND, NQ>(pDistribution, pDensity, velx, vely, velz) {}

    ~BdryRuleBounceBackFront() = default;

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EY(q) > 0)
        {
            // Bounce back rule for DFs coming from front wall.
            return this->compute_bounceback(q, i, j, k);
        }
        else
        {
            // Normal streaming.
            return (this->mpDistribution)->GetCurrF(q, i, j, k);
        }
    }
};

#endif // BDRYRULEBOUNCEBACKFRONTDEF