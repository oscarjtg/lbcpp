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

    void AddBoundaryRule(AbstractBoundaryRule<T>* bdry_rule, int bdry_id)
    {
        // Check that the bdry_id corresponds to the position of the new bdry_rule.
        if (bdry_id == static_cast<int>(mBdryRuleArr.size()))
        {
            mBdryRuleArr.push_back(bdry_rule);
        }
        else
        {
            throw std::invalid_argument("Invalid bdry_id for setting new bdry_rule for mBdryRuleArr.");
        }
    }

    T ComputeBdryRule(int bdry_id, int q, int i, int j, int k) const
    {
        if (bdry_id >= 0 && bdry_id < static_cast<int>(mBdryRuleArr.size()))
        {
            return mBdryRuleArr[bdry_id]->GetDistributionValue(q, i, j, k);
        }
        else
        {
            throw std::out_of_range("Invalid boundary ID for mBdryRuleArr.");
        }
    }

    T ComputeBdryRuleUnsafe(int bdry_id, int q, int i, int j, int k) const
    {
        return mBdryRuleArr[bdry_id]->GetDistributionValue(q, i, j, k);
    }

private:
    int mNumberOfBdrys;

    std::vector<AbstractBoundaryRule<T>*> mBdryRuleArr;
};


#endif // BOUNDARYINFODEF