#ifndef BDRYRULESCALARDIRICHLETBOTTOMDEF
#define BDRYRULESCALARDIRICHLETBOTTOMDEF

#include "boundary/AbstractBdryRuleAntiBounceBack.h"

template <typename T, int ND, int NQ>
class BdryRuleScalarDirichletBottom : public AbstractBdryRuleAntiBounceBack<T, ND, NQ>
{
public:
    BdryRuleScalarDirichletBottom() = default;

    BdryRuleScalarDirichletBottom(AbstractLattice<T, ND, NQ>* pDistribution, T wall_conc, T velx, T vely, T velz)
        : AbstractBdryRuleAntiBounceBack<T, ND, NQ>(pDistribution, wall_conc, velx, vely, velz) {}

    ~BdryRuleScalarDirichletBottom() = default;

    T GetWallConc() const
    {
        return this->mWallConc;
    }

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        if ((this->mpDistribution)->EZ(q) > 0)
        {
            // Bounce back rule for DFs coming from above.
            //std::cout << "BdryRuleScalarDirichletBottom.GetDistributionValue(q, i, j, k):\n";
            //std::cout << "Computing bounce back for q = " << q;
            //std::cout << " at node (" << i << ", " << j << ", " << k << ")";
            //std::cout << " with wall conc " << this->mWallConc << "\n";
            T ans = this->compute_antibounceback(q, i, j, k);
            //std::cout << "Returns: " << ans << "\n";
            return ans;
        }
        else
        {
            // Normal streaming.
            T ans = (this->mpDistribution)->GetCurrF(q, i, j, k);
            ////std::cout << "BdryRuleScalarDirichletBottom.GetDistributionValue(q, i, j, k):\n";
            ////std::cout << "Grabbing streamed F for q = " << q;
            ////std::cout << " at node (" << i << ", " << j << ", " << k << ")";
            ////std::cout << "Returns: " << ans << "\n";
            return ans;
        }
    }
};

#endif // BDRYRULESCALARDIRICHLETBOTTOMDEF