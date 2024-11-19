#ifndef ABSTRACTBOUNDARYRULEDEF
#define ABSTRACTBOUNDARYRULEDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/MacroscopicVariable.h"

template <typename T, int ND, int NQ>
class AbstractBoundaryRule
{
public:
    AbstractBoundaryRule() = default;

    // Another constructor that sets mpDistribution
    AbstractBoundaryRule(AbstractLattice<T, ND, NQ>* pDistribution)
        : mpDistribution(pDistribution) {}

    virtual void SetDistributionPointer(AbstractLattice<T, ND, NQ>* pDistribution) {
        mpDistribution = pDistribution;
    }

    virtual T GetDistributionValue(int q, int i, int j, int k) const = 0;
protected:
    // Pointer to instance of AbstractLattice<T, ND, NQ> called mpDistribution
    AbstractLattice<T, ND, NQ>* mpDistribution = nullptr;
};

#endif // ABSTRACTBOUNDARYRULEDEF