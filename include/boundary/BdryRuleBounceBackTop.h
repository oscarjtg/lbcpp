#ifndef BDRYRULEBOUNCEBACKTOPDEF
#define BDRYRULEBOUNCEBACKTOPDEF

#include "boundary/AbstractBdryRuleBounceBack.h"

template <typename T, int ND, int NQ>
class BdryRuleBounceBackTop : public AbstractBdryRuleBounceBack<T, ND, NQ>
{
public:
    BdryRuleBounceBackTop() = default;

    BdryRuleBounceBackTop(AbstractLattice<T, ND, NQ>* pDistribution)
        : AbstractBdryRuleBounceBack<T, ND, NQ>(pDistribution) {}

    BdryRuleBounceBackTop(AbstractLattice<T, ND, NQ>* pDistribution, ScalarField<T>* pDensity, T velx, T vely, T velz)
        : AbstractBdryRuleBounceBack<T, ND, NQ>(pDistribution, pDensity, velx, vely, velz) {}

    ~BdryRuleBounceBackTop() = default;

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EZ(q) < 0)
        {
            // Bounce back rule for DFs coming from above.
            return this->compute_bounceback(q, i, j, k);
        }
        else
        {
            // Normal streaming.
            return (this->mpDistribution)->GetCurrF(q, i, j, k);
        }
    }
};

#endif // BDRYRULEBOUNCEBACKTOPDEF