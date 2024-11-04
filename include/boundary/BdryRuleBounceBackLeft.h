#ifndef BDRYRULEBOUNCEBACKLEFTDEF
#define BDRYRULEBOUNCEBACKLEFTDEF

#include "boundary/AbstractBdryRuleBounceBack.h"

template <class T>
class BdryRuleBounceBackLeft : public AbstractBdryRuleBounceBack<T>
{
public:
    BdryRuleBounceBackLeft() = default;

    BdryRuleBounceBackLeft(AbstractLattice<T>* pDistribution)
        : AbstractBdryRuleBounceBack<T>(pDistribution) {}

    BdryRuleBounceBackLeft(AbstractLattice<T>* pDistribution, MacroscopicVariable<T>* pDensity, T velx, T vely, T velz)
        : AbstractBdryRuleBounceBack<T>(pDistribution, pDensity, velx, vely, velz) {}

    ~BdryRuleBounceBackLeft() = default;

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EX(q) > 0)
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

#endif // BDRYRULEBOUNCEBACKLEFTDEF