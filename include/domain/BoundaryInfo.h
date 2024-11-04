#ifndef BOUNDARYINFODEF
#define BOUNDARYINFODEF

#include <vector>
#include <stdexcept>

#include "boundary/AbstractBoundaryRule.h"
#include "boundary/AllBoundaryRules.h"

template <class T>
class BoundaryInfo
{
public:
    BoundaryInfo(int num) : mNumberOfBdrys(num) {}

    ~BoundaryInfo() = default;

    void AddBoundaryRuleF(AbstractBoundaryRule<T>* bdry_rule, int bdry_id)
    {
        // Check that the bdry_id corresponds to the position of the new bdry_rule.
        if (bdry_id == static_cast<int>(mArrayF.size()))
        {
            mArrayF.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mArrayF.");
        }
    }

    void AddBoundaryRuleG(AbstractBoundaryRule<T>* bdry_rule, int bdry_id)
    {
        // Check that the bdry_id corresponds to the position of the new bdry_rule.
        if (bdry_id == static_cast<int>(mArrayG.size()))
        {
            mArrayG.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mArrayG.");
        }
    }

    void AddBoundaryRuleH(AbstractBoundaryRule<T>* bdry_rule, int bdry_id)
    {
        // Check that the bdry_id corresponds to the position of the new bdry_rule.
        if (bdry_id >= 0 && bdry_id == static_cast<int>(mArrayH.size()))
        {
            mArrayH.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mArrayH.");
        }
    }

    T ComputeBdryRuleF(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mArrayF.size()))
        {
            return mArrayF[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalid boundary ID for mArrayF.");
        }
    }

    T ComputeBdryRuleG(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mArrayG.size()))
        {
            return mArrayG[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalid boundary ID for mArrayG.");
        }
    }

    T ComputeBdryRuleH(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mArrayH.size()))
        {
            return mArrayH[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalud boundary ID for mArrayH.");
        }
    }

    T ComputeBdryRuleUnsafeF(int bdry_id, int q, int i, int j, int k) const
    {
        return mArrayF[bdry_id]->GetDistributionValue(q, i, j, k);
    }

    T ComputeBdryRuleUnsafeG(int bdry_id, int q, int i, int j, int k) const
    {
        return mArrayG[bdry_id]->GetDistributionValue(q, i, j, k);
    }

    T ComputeBdryRuleUnsafeH(int bdry_id, int q, int i, int j, int k) const
    {
        return mArrayH[bdry_id]->GetDistributionValue(q, i, j, k);
    }

private:
    int mNumberOfBdrys;

    std::vector<AbstractBoundaryRule<T>*> mArrayF;
    
    std::vector<AbstractBoundaryRule<T>*> mArrayG;

    std::vector<AbstractBoundaryRule<T>*> mArrayH;
};


#endif // BOUNDARYINFODEF