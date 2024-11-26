#ifndef FLUIDEVOLVERSRTDEF
#define FLUIDEVOLVERSRTDEF

#include "evolver/AbstractFluidEvolver.h"

template <typename T, int ND, int NQ>
class FluidEvolverSRT : public AbstractFluidEvolver<T, ND, NQ>
{
public:
    FluidEvolverSRT() {std::cout << "Single relaxation time collision operator (SRT)\n";};

    ~FluidEvolverSRT() = default;

    virtual void DoLocalCollision(std::function<T(T, T, T, T)> ComputeEquilibrium, AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k) override;
};

#include "evolver/FluidEvolverSRT.tpp"

#endif // FLUIDEVOLVERSRTDEF