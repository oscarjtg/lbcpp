#ifndef LOCALDISTRIBUTIOND2Q5DEF
#define LOCALDISTRIBUTIOND2Q5DEF

template<class T> class LocalDistributionD2Q5
{
public:
    // Distribution functions.
    T df0;
    T df1;
    T df2;
    T df3;
    T df4;
    // Equilibrium distribution functions.
    T ef0;
    T ef1;
    T ef2;
    T ef3;
    T ef4;
    // Moments.
    T m0;
    T m1;
    T m2;
    T m3;
    T m4;
    // Equilibrium moments.
    T em0;
    T em1;
    T em2;
    T em3;
    T em4;
    // Functions.
    void ComputeM0;
    void ComputeM1;
    void ComputeM2;
    void ComputeM3;
    void ComputeM4;
}

#endif // LOCALDISTRIBUTIOND2Q5DEF