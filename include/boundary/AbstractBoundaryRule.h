#ifndef ABSTRACTBOUNDARYRULEDEF
#define ABSTRACTBOUNDARYRULEDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/MacroscopicVariable.h"

template <class T>
class AbstractBoundaryRule
{
public:
    AbstractBoundaryRule() = default;

    // Another constructor that sets mpDistribution
    AbstractBoundaryRule(AbstractLattice<T>* pDistribution)
        : mpDistribution(pDistribution) {}

    virtual void SetDistributionPointer(AbstractLattice<T>* pDistribution) {
        mpDistribution = pDistribution;
    }

    virtual T GetDistributionValue(int q, int i, int j, int k) const = 0;
protected:
    // Pointer to instance of AbstractLattice<T> called mpDistribution
    AbstractLattice<T>* mpDistribution = nullptr;
};

#endif // ABSTRACTBOUNDARYRULEDEF