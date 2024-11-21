#ifndef SCALAREVOLVERTRTDEF
#define SCALAREVOLVERTRTDEF

#include "AbstractScalarEvolver.h"

template <typename T, int ND, int NQ>
class ScalarEvolverTRT : public AbstractScalarEvolver<T, ND, NQ>
{
public:
    ScalarEvolverTRT() = default;
    
    ~ScalarEvolverTRT() = default;

    void SetScalarDiffusivity(AbstractLattice<T, ND, NQ>& g, T kappa) override;

    void SetMagicParameter(T magic_parameter) { mMagicParameter = magic_parameter; }

    void DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k) override;

private:
    T mOmegaPlus, mOmegaMinus, mMagicParameter = 1./4.;
};

#include "ScalarEvolverTRT.tpp"

#endif // SCALAREVOLVERTRTDEF