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
        if (bdry_id == static_cast<int>(mBdryRuleArrF.size()))
        {
            mBdryRuleArrF.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mBdryRuleArrF.");
        }
    }

    void AddBoundaryRuleG(AbstractBoundaryRule<T>* bdry_rule, int bdry_id)
    {
        // Check that the bdry_id corresponds to the position of the new bdry_rule.
        if (bdry_id == static_cast<int>(mBdryRuleArrG.size()))
        {
            mBdryRuleArrG.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mBdryRuleArrG.");
        }
    }

    void AddBoundaryRuleH(AbstractBoundaryRule<T>* bdry_rule, int bdry_id)
    {
        // Check that the bdry_id corresponds to the position of the new bdry_rule.
        if (bdry_id >= 0 && bdry_id == static_cast<int>(mBdryRuleArrH.size()))
        {
            mBdryRuleArrH.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mBdryRuleArrH.");
        }
    }

    T ComputeBdryRuleF(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mBdryRuleArrF.size()))
        {
            return mBdryRuleArrF[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalid boundary ID for mBdryRuleArrF.");
        }
    }

    T ComputeBdryRuleG(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mBdryRuleArrG.size()))
        {
            return mBdryRuleArrG[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalid boundary ID for mBdryRuleArrG.");
        }
    }

    T ComputeBdryRuleH(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mBdryRuleArrH.size()))
        {
            return mBdryRuleArrH[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalud boundary ID for mBdryRuleArrH.");
        }
    }

    T ComputeBdryRuleUnsafeF(int bdry_id, int q, int i, int j, int k) const
    {
        return mBdryRuleArrF[bdry_id]->GetDistributionValue(q, i, j, k);
    }

    T ComputeBdryRuleUnsafeG(int bdry_id, int q, int i, int j, int k) const
    {
        return mBdryRuleArrG[bdry_id]->GetDistributionValue(q, i, j, k);
    }

    T ComputeBdryRuleUnsafeH(int bdry_id, int q, int i, int j, int k) const
    {
        return mBdryRuleArrH[bdry_id]->GetDistributionValue(q, i, j, k);
    }

private:
    int mNumberOfBdrys;

    std::vector<AbstractBoundaryRule<T>*> mBdryRuleArrF;
    
    std::vector<AbstractBoundaryRule<T>*> mBdryRuleArrG;

    std::vector<AbstractBoundaryRule<T>*> mBdryRuleArrH;
};


#endif // BOUNDARYINFODEF