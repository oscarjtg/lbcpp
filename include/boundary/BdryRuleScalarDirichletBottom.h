#ifndef BDRYRULESCALARDIRICHLETBOTTOMDEF
#define BDRYRULESCALARDIRICHLETBOTTOMDEF

#include "boundary/AbstractBdryRuleAntiBounceBack.h"

template <class T>
class BdryRuleScalarDirichletBottom : public AbstractBdryRuleAntiBounceBack<T>
{
public:
    BdryRuleScalarDirichletBottom() = default;

    BdryRuleScalarDirichletBottom(AbstractLattice<T>* pDistribution, T wall_conc, T velx, T vely, T velz)
        : AbstractBdryRuleAntiBounceBack<T>(pDistribution, wall_conc, velx, vely, velz) {}

    ~BdryRuleScalarDirichletBottom() = default;

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EZ(q) > 0)
        {
            // Bounce back rule for DFs coming from above.
            return this->compute_antibounceback(q, i, j, k);
        }
        else
        {
            // Normal streaming.
            return (this->mpDistribution)->GetCurrF(q, i, j, k);
        }
    }
};

#endif // BDRYRULESCALARDIRICHLETBOTTOMDEF