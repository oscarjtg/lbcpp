#ifndef ABSTRACTFLUIDEVOLVERDEF
#define ABSTRACTFLUIDEVOLVERDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/MacroscopicVariable.h"
#include "domain/BoundaryInfo.h"
#include "domain/NodeInfo.h"

#include <functional>

template <typename T, int ND, int NQ>
class AbstractFluidEvolver
{
public:
    AbstractFluidEvolver() = default;

    ~AbstractFluidEvolver() = default;

    virtual void SetKinematicViscosity(AbstractLattice<T, ND, NQ>& f, T nu);

    virtual void Initialise(AbstractLattice<T, ND, NQ>& f,
                            const MacroscopicVariable<T>& dens,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void InitialiseEquilibrium(AbstractLattice<T, ND, NQ>& f,
                            const MacroscopicVariable<T>& dens,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void InitialiseWei(AbstractLattice<T, ND, NQ>& f,
                            MacroscopicVariable<T>& dens,
                            const MacroscopicVariable<T>& velx,
                            const MacroscopicVariable<T>& vely,
                            const MacroscopicVariable<T>& velz,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void DoTimestep(AbstractLattice<T, ND, NQ>& f,
                            MacroscopicVariable<T>& dens,
                            MacroscopicVariable<T>& velx,
                            MacroscopicVariable<T>& vely,
                            MacroscopicVariable<T>& velz,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            NodeInfo& node,
                            BoundaryInfo<T, ND, NQ>& bdry);

    virtual void DoLocalTimestep(AbstractLattice<T, ND, NQ>& f,
                            const MacroscopicVariable<T>& Fx,
                            const MacroscopicVariable<T>& Fy,
                            const MacroscopicVariable<T>& Fz,
                            NodeInfo& node,
                            BoundaryInfo<T, ND, NQ>& bdry,
                            T& r_, T& u_, T& v_, T& w_, int i, int j, int k);

    virtual void DoLocalTimestepIncompressible(AbstractLattice<T, ND, NQ>& f,
                        const MacroscopicVariable<T>& velx,
                        const MacroscopicVariable<T>& vely,
                        const MacroscopicVariable<T>& velz,
                        const MacroscopicVariable<T>& Fx,
                        const MacroscopicVariable<T>& Fy,
                        const MacroscopicVariable<T>& Fz,
                        const NodeInfo& node,
                        const BoundaryInfo<T, ND, NQ>& bdry, 
                        T& r_, int i, int j, int k);

    virtual void DoLocalCollision(std::function<T(T, T, T, T)> ComputeEquilibrium, AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k) = 0;

protected:
    T mOmega;
};

#include "AbstractFluidEvolver.tpp"

#endif // ABSTRACTFLUIDEVOLVERDEF