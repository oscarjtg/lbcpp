#ifndef BDRYRULESCALARNEUMANNRIGHTDEF
#define BDRYRULESCALARNEUMANNRIGHTDEF

#include "boundary/AbstractBdryRuleAntiBounceBack.h"

template <typename T, int ND, int NQ>
class BdryRuleScalarNeumannRight : public AbstractBdryRuleAntiBounceBack<T, ND, NQ>
{
public:
    BdryRuleScalarNeumannRight() = default;

    BdryRuleScalarNeumannRight(AbstractLattice<T, ND, NQ>* pDistribution, ScalarField<T>* pConc, T velx, T vely, T velz)
        : AbstractBdryRuleAntiBounceBack<T, ND, NQ>(pDistribution, 0.0, velx, vely, velz), mpConc(pConc){}

    ~BdryRuleScalarNeumannRight() = default;

    void SetConcentrationPointer(ScalarField<T>* pConcentration)
    {
        mpConc = pConcentration;
    }

    T GetDistributionValue(int q, int i, int j, int k) const override
    {
        // Get concentration one grid point in direction of normal pointing out of wall 
        // For Right wall, the normal is (-1, 0, 0)
        this->mWallConc = mpConc->GetValue(i-1, j, k);

        //std::cout << "mWallConc = " << this->mWallConc << "\n";

        if ((this->mpDistribution)->EX(q) < 0)
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
private:
    ScalarField<T>* mpConc;
};

#endif // BDRYRULESCALARNEUMANNRIGHTDEF