#ifndef ABSTRACTFLUIDEVOLVERDEF
#define ABSTRACTFLUIDEVOLVERDEF

#include "lattice/AbstractLattice.h"
#include "macroscopic/AllFields.h"
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
                            ScalarField<T>& dens,
                            const VectorField<T>& vel,
                            const VectorField<T>& force,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry, 
                            const std::string method="MEI");

    virtual void InitialiseNEQ(AbstractLattice<T, ND, NQ>& f,
                            const ScalarField<T>& dens,
                            const VectorField<T>& vel,
                            const VectorField<T>& force,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void InitialiseFEQ(AbstractLattice<T, ND, NQ>& f,
                            const ScalarField<T>& dens,
                            const VectorField<T>& vel,
                            const VectorField<T>& force,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void InitialiseMEI(AbstractLattice<T, ND, NQ>& f,
                            ScalarField<T>& dens,
                            const VectorField<T>& vel,
                            const VectorField<T>& force,
                            const NodeInfo& node,
                            const BoundaryInfo<T, ND, NQ>& bdry);

    virtual void DoTimestep(AbstractLattice<T, ND, NQ>& f,
                            ScalarField<T>& dens,
                            VectorField<T>& vel,
                            const VectorField<T>& force,
                            NodeInfo& node,
                            BoundaryInfo<T, ND, NQ>& bdry,
                            const bool StoreMacros=true);

    virtual void DoTimestep(AbstractLattice<T, ND, NQ>& f,
                            ScalarField<T>& dens,
                            VectorField<T>& vel,
                            TensorField<T>& sigma,
                            const VectorField<T>& force,
                            NodeInfo& node,
                            BoundaryInfo<T, ND, NQ>& bdry,
                            const bool StoreMacros=true);

    virtual void DoLocalTimestep(AbstractLattice<T, ND, NQ>& f,
                            const VectorField<T>& force,
                            NodeInfo& node,
                            BoundaryInfo<T, ND, NQ>& bdry,
                            T& r_, T& u_, T& v_, T& w_, int i, int j, int k);

    virtual void DoLocalTimestep(AbstractLattice<T, ND, NQ>& f,
                            const VectorField<T>& force,
                            NodeInfo& node,
                            BoundaryInfo<T, ND, NQ>& bdry,
                            T& r_, std::array<T, 3>& u_vec, std::array<T, 9>& s_mat, int i, int j, int k);

    virtual void DoLocalTimestepIncompressible(AbstractLattice<T, ND, NQ>& f,
                        const VectorField<T>& vel,
                        const VectorField<T>& force,
                        const NodeInfo& node,
                        const BoundaryInfo<T, ND, NQ>& bdry, 
                        T& r_, int i, int j, int k);

    virtual void DoLocalCollision(std::function<T(T, T, T, T)> ComputeEquilibrium, AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k) = 0;

protected:
    T mOmega;
};

#include "AbstractFluidEvolver.tpp"

#endif // ABSTRACTFLUIDEVOLVERDEF