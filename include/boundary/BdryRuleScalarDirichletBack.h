#ifndef BDRYRULESCALARDIRICHLETBACKDEF
#define BDRYRULESCALARDIRICHLETBACKDEF

#include "boundary/AbstractBdryRuleAntiBounceBack.h"

template <typename T, int ND, int NQ>
class BdryRuleScalarDirichletBack : public AbstractBdryRuleAntiBounceBack<T, ND, NQ>
{
public:
    BdryRuleScalarDirichletBack() = default;

    BdryRuleScalarDirichletBack(AbstractLattice<T, ND, NQ>* pDistribution, T wall_conc, T velx, T vely, T velz)
        : AbstractBdryRuleAntiBounceBack<T, ND, NQ>(pDistribution, wall_conc, velx, vely, velz) {}

    ~BdryRuleScalarDirichletBack() = default;

    T GetWallConc() const
    {
        return this->mWallConc;
    }

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EY(q) < 0)
        {
            // Bounce back rule for DFs coming from Back wall.
            T ans = this->compute_antibounceback(q, i, j, k);
            return ans;
        }
        else
        {
            // Normal streaming for DFs coming from fluid.
            T ans = (this->mpDistribution)->GetCurrF(q, i, j, k);
            return ans;
        }
    }
};

#endif // BDRYRULESCALARDIRICHLETBACKDEF