#include <iostream>
#include "lattice/LatticeD2Q5.h"

// Overridden constructor.
// Allocates memory. Prints size of data allocated.
template<class T>
LatticeD2Q5<T>::LatticeD2Q5(int size_x, int size_y)
    : mSizeX(size_x+1), mSizeY(size_y+1), mGridSize((size_x+1)*(size_y+1))
{
    mpF0 = new T [mGridSize];
    mpF1 = new T [mGridSize];
    mpF2 = new T [mGridSize];
    mpF3 = new T [mGridSize];
    mpF4 = new T [mGridSize];
    size_t alloc_size = 5 * mGridSize * sizeof(T);
    std::cout << "Allocated " << alloc_size << " bytes of memory (LatticeD2Q5)." << std::endl;
}

// Overriden destructor.
// Frees all dynamically-allocated memory.
template<class T>
LatticeD2Q5<T>::~LatticeD2Q5()
{
    delete[] mpF0;
    delete[] mpF1;
    delete[] mpF2;
    delete[] mpF3;
    delete[] mpF4;
    size_t delete_size = 5 * mGridSize * sizeof(T);
    std::cout << "Deleted " << delete_size << " bytes of memory (LatticeD2Q5)." << std::endl;
}

// Get data at (i, j), esoteric twist style!
template<class T>
T LatticeD2Q5<T>::GetF0(int i, int j) const
{
    return mpF0[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q5<T>::GetF1(int i, int j) const
{
    return mpF1[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q5<T>::GetF2(int i, int j) const
{
    return mpF2[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q5<T>::GetF3(int i, int j) const
{
    return mpF3[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q5<T>::GetF4(int i, int j) const
{
    return mpF4[scalar_index(i, j + 1)];
}

template<class T>
void LatticeD2Q5<T>::SetF0(T f_, int i, int j)
{
    mpF0[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetF1(T f_, int i, int j)
{
    mpF2[scalar_index(i + 1, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetF2(T f_, int i, int j)
{
    mpF1[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetF3(T f_, int i, int j)
{
    mpF4[scalar_index(i, j + 1)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetF4(T f_, int i, int j)
{
    mpF3[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetInitF0(T f_, int i, int j)
{
    mpF0[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetInitF1(T f_, int i, int j)
{
    mpF1[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetInitF2(T f_, int i, int j)
{
    mpF2[scalar_index(i + 1, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetInitF3(T f_, int i, int j)
{
    mpF3[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SetInitF4(T f_, int i, int j)
{  
    mpF4[scalar_index(i, j + 1)] = f_;
}

template<class T>
void LatticeD2Q5<T>::SwapPointers()
{
    std::swap(mpF1, mpF2);
    std::swap(mpF3, mpF4);
}

template<class T>
void LatticeD2Q5<T>::DisplayPointerInfo() const
{
    std::cout << "mpF0   = " << mpF0   << "\n";
    std::cout << "mpF1   = " << mpF1   << "\n";
    std::cout << "mpF2   = " << mpF2   << "\n";
    std::cout << "mpF3   = " << mpF3   << "\n";
    std::cout << "mpF4   = " << mpF4   << "\n";
}

template<class T>
void LatticeD2Q5<T>::DisplayInfo() const
{
    std::cout << "Lattice D2Q5\n";
    std::cout << "Size X = " << mSizeX << "\n";
    std::cout << "Size Y = " << mSizeY << "\n";
}

// Explicit template instantiation.
//template class LatticeD2Q5<_Float16>;
template class LatticeD2Q5<float>;
template class LatticeD2Q5<double>;