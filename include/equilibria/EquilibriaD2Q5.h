#ifndef EQUILIBRIAD2Q5DEF
#define EQUILIBRIAD2Q5DEF

/**********************************************
 * 
 *  First order equilibria
 * 
 **********************************************/

template<class T>
inline T D2Q5_1stOrder_Feq0(T r_)
{
    return r_ / 3;
}

template<class T>
inline T D2Q5_1stOrder_Feq1(T r_, T u_)
{
    return r_ * (1 + 3*u_) / 6;
}

template<class T>
inline T D2Q5_1stOrder_Feq2(T r_, T u_)
{
    return r_ * (1 - 3*u_) / 6;
}

template<class T>
inline T D2Q5_1stOrder_Feq3(T r_, T v_)
{
    return r_ * (1 + 3*v_) / 6;
}

template<class T>
inline T D2Q5_1stOrder_Feq4(T r_, T v_)
{
    return r_ * (1 - 3*v_) / 6;
}

/**********************************************
 * 
 *  Second order equilibria
 * 
 **********************************************/

template<class T>
inline T D2Q5_2ndOrder_Feq0(T r_, T u2)
{
    return r_ * (2 - 3*u2) / 6;
}

/*
template<class T>
inline T D2Q5_2ndOrder_Feq1(T r_, T u_, T u2)
{
    return r_ * ()
}
*/

#endif // EQUILIBRIAD2Q5DEF