#ifndef BDRYRULEBOUNCEBACKRIGHTDEF
#define BDRYRULEBOUNCEBACKRIGHTDEF

#include "boundary/AbstractBdryRuleBounceBack.h"

template <class T>
class BdryRuleBounceBackRight : public AbstractBdryRuleBounceBack<T>
{
public:
    BdryRuleBounceBackRight() = default;

    BdryRuleBounceBackRight(AbstractLattice<T>* pDistribution)
        : AbstractBdryRuleBounceBack<T>(pDistribution) {}

    BdryRuleBounceBackRight(AbstractLattice<T>* pDistribution, MacroscopicVariable<T>* pDensity, T velx, T vely, T velz)
        : AbstractBdryRuleBounceBack<T>(pDistribution, pDensity, velx, vely, velz) {}

    ~BdryRuleBounceBackRight() = default;

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EX(q) < 0)
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

#endif // BDRYRULEBOUNCEBACKRIGHTDEF