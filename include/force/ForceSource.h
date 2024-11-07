#ifndef FORCESOURCEDEF
#define FORCESOURCEDEF

template <class T>
inline T computeSecondOrderForceSource(T cu, T cF, T uF, T w)
{
    return w * (cF + cF*cu - uF);
}

#endif // FORCESOURCEDEF