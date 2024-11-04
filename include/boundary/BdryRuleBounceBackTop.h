#ifndef BDRYRULEBOUNCEBACKTOPDEF
#define BDRYRULEBOUNCEBACKTOPDEF

#include "boundary/AbstractBdryRuleBounceBack.h"

template <class T>
class BdryRuleBounceBackTop : public AbstractBdryRuleBounceBack<T>
{
public:
    BdryRuleBounceBackTop() = default;

    BdryRuleBounceBackTop(AbstractLattice<T>* pDistribution)
        : AbstractBdryRuleBounceBack<T>(pDistribution) {}

    BdryRuleBounceBackTop(AbstractLattice<T>* pDistribution, MacroscopicVariable<T>* pDensity, T velx, T vely, T velz)
        : AbstractBdryRuleBounceBack<T>(pDistribution, pDensity, velx, vely, velz) {}

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