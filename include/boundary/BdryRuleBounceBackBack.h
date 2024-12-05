#ifndef BDRYRULEBOUNCEBACKBACKDEF
#define BDRYRULEBOUNCEBACKBACKDEF

#include "boundary/AbstractBdryRuleBounceBack.h"

template <typename T, int ND, int NQ>
class BdryRuleBounceBackBack : public AbstractBdryRuleBounceBack<T, ND, NQ>
{
public:
    BdryRuleBounceBackBack() = default;

    BdryRuleBounceBackBack(AbstractLattice<T, ND, NQ>* pDistribution)
        : AbstractBdryRuleBounceBack<T, ND, NQ>(pDistribution) {}

    BdryRuleBounceBackBack(AbstractLattice<T, ND, NQ>* pDistribution, ScalarField<T>* pDensity, T velx, T vely, T velz)
        : AbstractBdryRuleBounceBack<T, ND, NQ>(pDistribution, pDensity, velx, vely, velz) {}

    ~BdryRuleBounceBackBack() = default;

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EY(q) < 0)
        {
            // Bounce back rule for DFs coming from back wall.
            return this->compute_bounceback(q, i, j, k);
        }
        else
        {
            // Normal streaming.
            return (this->mpDistribution)->GetCurrF(q, i, j, k);
        }
    }
};

#endif // BDRYRULEBOUNCEBACKBACKDEF