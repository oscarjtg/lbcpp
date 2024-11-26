#ifndef FLUIDEVOLVERTRTDEF
#define FLUIDEVOLVERTRTDEF

#include "AbstractFluidEvolver.h"

template <typename T, int ND, int NQ>
class FluidEvolverTRT : public AbstractFluidEvolver<T, ND, NQ>
{
public:
    FluidEvolverTRT() {std::cout << "Two relaxation time collision operator (TRT)\n";};

    ~FluidEvolverTRT() = default;

    void SetKinematicViscosity(AbstractLattice<T, ND, NQ>& f, T nu) override;

    void SetMagicParameter(T magic_parameter) { mMagicParameter = magic_parameter; }

    void DoLocalCollision(std::function<T(T, T, T, T)> ComputeEquilibrium, AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k) override;

private:
    T mMagicParameter = 1./4., mOmegaPlus, mOmegaMinus;
};

#include "FluidEvolverTRT.tpp"

#endif // FLUIDEVOLVERTRTDEF