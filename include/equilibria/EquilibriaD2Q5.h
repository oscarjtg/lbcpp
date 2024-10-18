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

template<class T>
inline T D2Q5_2ndOrder_Feq1(T r_, T u_, T u2)
{
    return r_ * (2 + 6*u_ + 9*u_*u_ - 3*u2) / 12;
}

template<class T>
inline T D2Q5_2ndOrder_Feq2(T r_, T u_, T u2)
{
    return r_ * (2 - 6*u_ + 9*u_*u_ - 3*u2) / 12;
}

template<class T>
inline T D2Q5_2ndOrder_Feq3(T r_, T v_, T u2)
{
    return r_ * (2 + 6*v_ + 9*v_*v_ - 3*u2) / 12;
}

template<class T>
inline T D2Q5_2ndOrder_Feq4(T r_, T v_, T u2)
{
    return r_ * (2 - 6*v_ + 9*v_*v_ - 3*u2) / 12;
}

/**********************************************
 * 
 *  Second order equilibria perturbations
 * 
 **********************************************/



#endif // EQUILIBRIAD2Q5DEF