#ifndef SCALAREVOLVERSRTDEF
#define SCALAREVOLVERSRTDEF

#include "evolver/AbstractScalarEvolver.h"

template <typename T, int ND, int NQ>
class ScalarEvolverSRT : public AbstractScalarEvolver<T, ND, NQ>
{
public:
    ScalarEvolverSRT() = default;

    ~ScalarEvolverSRT() = default;

    virtual void DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k) override;
};

#include "evolver/ScalarEvolverSRT.tpp"

#endif // SCALAREVOLVERSRTDEF