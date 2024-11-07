#ifndef ABSTRACTSCALAREVOLVERDEF
#define ABSTRACTSCALAREVOLVERDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/MacroscopicVariable.h"
#include "domain/BoundaryInfo.h"
#include "domain/NodeInfo.h"

template <class T>
class AbstractScalarEvolver
{
public:
    AbstractScalarEvolver() = default;

    virtual ~AbstractScalarEvolver() = default;

    virtual void SetScalarDiffusivity(AbstractLattice<T>& g, T kappa) = 0;

    virtual void Initialise(AbstractLattice<T>& g,
                            const MacroscopicVariable<T>& conc,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const NodeInfo& node,
                            const BoundaryInfo<T>& bdry) = 0;

    virtual void DoTimestep(AbstractLattice<T>& g,
                            MacroscopicVariable<T>& conc,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const NodeInfo& node,
                            const BoundaryInfo<T>& bdry) = 0;
};

#endif // ABSTRACTSCALAREVOLVERDEF