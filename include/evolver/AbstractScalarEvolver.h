#ifndef ABSTRACTSCALAREVOLVERDEF
#define ABSTRACTSCALAREVOLVERDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/AllFields.h"
#include "domain/BoundaryInfo.h"
#include "domain/NodeInfo.h"

template <typename T, int ND, int NQ>
class AbstractScalarEvolver
{
public:
    AbstractScalarEvolver() = default;

    virtual ~AbstractScalarEvolver() = default;

    virtual void SetScalarDiffusivity(AbstractLattice<T, ND, NQ>& g, T kappa);

    virtual void Initialise(AbstractLattice<T, ND, NQ>& g,
                            const MacroscopicVariable<T>& conc,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void DoTimestep(AbstractLattice<T, ND, NQ>& g,
                            MacroscopicVariable<T>& conc,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k) = 0;

protected:
    T mOmega;
};

#include "AbstractScalarEvolver.tpp"

#endif // ABSTRACTSCALAREVOLVERDEF