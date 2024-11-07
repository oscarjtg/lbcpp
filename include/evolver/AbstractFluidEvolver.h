#ifndef ABSTRACTFLUIDEVOLVERDEF
#define ABSTRACTFLUIDEVOLVERDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/MacroscopicVariable.h"
#include "domain/BoundaryInfo.h"
#include "domain/NodeInfo.h"

template <class T>
class AbstractFluidEvolver
{
public:
    AbstractFluidEvolver() = default;

    ~AbstractFluidEvolver() = default;

    virtual void SetKinematicViscosity(AbstractLattice<T>& f, T nu) = 0;

    virtual void SetBulkViscosity(AbstractLattice<T>& f, T eta) = 0;

    virtual void Initialise(AbstractLattice<T>& f,
                            const MacroscopicVariable<T>& dens,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            const NodeInfo& node,
                            const BoundaryInfo<T>& bdry) = 0;

    virtual void DoTimestep(AbstractLattice<T>& f,
                            MacroscopicVariable<T>& dens,
                            MacroscopicVariable<T>& velx,
                            MacroscopicVariable<T>& vely,
                            MacroscopicVariable<T>& velz,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            NodeInfo& node,
                            BoundaryInfo<T>& bdry) = 0;
};

#endif // ABSTRACTFLUIDEVOLVERDEF