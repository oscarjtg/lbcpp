#ifndef FORCESOURCEDEF
#define FORCESOURCEDEF

template <class T>
inline T computeSecondOrderForceSource(T cu, T cF, T uF, T w)
{
    return w * (
        static_cast<T>(3.0)*cF
        + static_cast<T>(9.0)*cF*cu
        - static_cast<T>(3.0)*uF    
    );
}

#endif // FORCESOURCEDEF