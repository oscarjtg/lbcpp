#include <iostream>

#include "macroscopic/Macroscopic2D.h"

/**
 * @brief Overridden constructor. Allocates memory for macroscopic arrays.
 * 
 * @param size_x The number of grid points along the x-axis.
 * @param size_y The number of grid points along the y-axis.
 */
template<class T>
Macroscopic2D<T>::Macroscopic2D(int size_x, int size_y)
: mSizeX(size_x), mSizeY(size_y), mGridSize(size_x * size_y)
{
    mpR = new T [mGridSize];
    mpU = new T [mGridSize];
    mpV = new T [mGridSize];
    mpT = new T [mGridSize];
    mpS = new T [mGridSize];
    size_t alloc_size = 5 * mGridSize * sizeof(T);
    std::cout << "Allocated " << alloc_size << " bytes of memory (Macroscopic2D)." << std::endl;
}

/**
 * @brief Overridden destructor. Frees the dynamically-allocated memory.
 */
template<class T>
Macroscopic2D<T>::~Macroscopic2D()
{
    delete[] mpR;
    delete[] mpU;
    delete[] mpV;
    delete[] mpT;
    delete[] mpS;
    size_t delete_size = 5 * mGridSize * sizeof(T);
    std::cout << "Deleted " << delete_size << " bytes of memory (Macroscopic2D)." << std::endl;
}

/**
 * @brief Get the fluid density at a given grid cell.
 * 
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 * @return The fluid density at grid cell (i, j).
 */
template<class T>
T Macroscopic2D<T>::GetDensity(int i, int j) const
{
    return mpR[scalar_index(i, j)];
}

/**
 * @brief Get the x-component of fluid velocity at a given grid cell.
 * 
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 * @return The x-component of fluid velocity at grid cell (i, j).
 */
template<class T>
T Macroscopic2D<T>::GetVelocityX(int i, int j) const
{
    return mpU[scalar_index(i, j)];
}

/**
 * @brief Get the y-component of fluid velocity at a given grid cell.
 * 
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 * @return The y-component of fluid velocity at grid cell (i, j).
 */
template<class T>
T Macroscopic2D<T>::GetVelocityY(int i, int j) const
{
    return mpV[scalar_index(i, j)];
}

/**
 * @brief Get the temperature at a given grid cell.
 * 
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 * @return The temperature at grid cell (i, j).
 */
template<class T>
T Macroscopic2D<T>::GetTemperature(int i, int j) const
{
    return mpT[scalar_index(i, j)];
}

/**
 * @brief Get the salt concentration at a given grid cell.
 * 
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 * @return The salt concentration at grid cell (i, j).
 */
template<class T>
T Macroscopic2D<T>::GetSalinity(int i, int j) const
{
    return mpS[scalar_index(i, j)];
}

/** 
 * @brief Store the density value at a given grid cell.
 * 
 * @param r_ The density value to be stored to an array.
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 */
template<class T>
void Macroscopic2D<T>::SetDensity(T r_, int i, int j)
{
    mpR[scalar_index(i, j)] = r_;
}

/** 
 * @brief Store the velocity x-component value at a given grid cell.
 * 
 * @param u_ The velocity x-component value to be stored to an array.
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 */
template<class T>
void Macroscopic2D<T>::SetVelocityX(T u_, int i, int j)
{
    mpU[scalar_index(i, j)] = u_;
}

/** 
 * @brief Store the velocity y-component value at a given grid cell.
 * 
 * @param v_ The velocity y-component value to be stored to an array.
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 */
template<class T>
void Macroscopic2D<T>::SetVelocityY(T v_, int i, int j)
{
    mpV[scalar_index(i, j)] = v_;
}

/** 
 * @brief Store the temperature value at a given grid cell.
 * 
 * @param t_ The density value to be stored to an array.
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 */
template<class T>
void Macroscopic2D<T>::SetTemperature(T t_, int i, int j)
{
    mpT[scalar_index(i, j)] = t_;
}

/** 
 * @brief Store the salt concentration value at a given grid cell.
 * 
 * @param s_ The salt concentration value to be stored to an array.
 * @param i The int index for the grid cell along the x-axis.
 * @param j The int index for the grid cell along the y-axis.
 */
template<class T>
void Macroscopic2D<T>::SetSalinity(T s_, int i, int j)
{
    mpT[scalar_index(i, j)] = s_;
}



// Explicit template instantiation.
//template class Macroscopic2D<_Float16>;
template class Macroscopic2D<float>;
template class Macroscopic2D<double>;