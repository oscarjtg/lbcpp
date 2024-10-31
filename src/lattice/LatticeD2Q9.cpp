#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // For std::setw and std::setfill
#include <sstream> // For std::ostringstream

#include "lattice/LatticeD2Q9.h"

/**********************************************
 * 
 *  Constructors and destructors
 * 
 **********************************************/

// Overridden constructor.
// Allocates memory.
template<class T>
LatticeD2Q9<T>::LatticeD2Q9(int size_x, int size_y)
    : mSizeX(size_x+1), mSizeY(size_y+1), mGridSize((size_x+1)*(size_y+1))
{
    mpF0 = new T [mGridSize];
    mpF1 = new T [mGridSize];
    mpF2 = new T [mGridSize];
    mpF3 = new T [mGridSize];
    mpF4 = new T [mGridSize];
    mpF5 = new T [mGridSize];
    mpF6 = new T [mGridSize];
    mpF7 = new T [mGridSize];
    mpF8 = new T [mGridSize];
    size_t alloc_size = 9 * mGridSize * sizeof(T);
    std::cout << "Allocated " << alloc_size << " bytes of memory (LatticeD2Q9)." << std::endl;
}

// Overridden destructor.
// Frees all dynamically-allocated memory.
template<class T>
LatticeD2Q9<T>::~LatticeD2Q9()
{
    delete[] mpF0;
    delete[] mpF1;
    delete[] mpF2;
    delete[] mpF3;
    delete[] mpF4;
    delete[] mpF5;
    delete[] mpF6;
    delete[] mpF7;
    delete[] mpF8;
    size_t delete_size = 9 * mGridSize * sizeof(T);
    std::cout << "Deleted " << delete_size << " bytes of memory (LatticeD2Q9)." << std::endl;
}

/**********************************************
 * 
 *  Standard getters
 * 
 **********************************************/

// Get data at (i, j), esoteric twist style!
template<class T>
T LatticeD2Q9<T>::GetF0(int i, int j) const
{
    return mpF0[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetF1(int i, int j) const
{
    return mpF1[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetF2(int i, int j) const
{
    return mpF2[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q9<T>::GetF3(int i, int j) const
{
    return mpF3[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetF4(int i, int j) const
{
    return mpF4[scalar_index(i, j + 1)];
}

template<class T>
T LatticeD2Q9<T>::GetF5(int i, int j) const
{
    return mpF5[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetF6(int i, int j) const
{
    return mpF6[scalar_index(i + 1, j + 1)];
}

template<class T>
T LatticeD2Q9<T>::GetF7(int i, int j) const
{
    return mpF7[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q9<T>::GetF8(int i, int j) const
{
    return mpF8[scalar_index(i, j + 1)];
}

/**********************************************
 * 
 *  Pre stream getters
 * 
 **********************************************/

// These should be used to get the postcollided values 
// that were at (i, j) before streaming,
// after having done the pointer swap.

template<class T>
T LatticeD2Q9<T>::GetPreStreamF0(int i, int j) const
{
    return mpF0[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF1(int i, int j) const
{
    return mpF2[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF2(int i, int j) const
{
    return mpF1[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF3(int i, int j) const
{
    return mpF4[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF4(int i, int j) const
{
    return mpF3[scalar_index(i, j + 1)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF5(int i, int j) const
{
    return mpF6[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF6(int i, int j) const
{
    return mpF5[scalar_index(i + 1, j + 1)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF7(int i, int j) const
{
    return mpF8[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q9<T>::GetPreStreamF8(int i, int j) const
{
    return mpF7[scalar_index(i, j + 1)];
}

/**********************************************
 * 
 *  Bounce back getters
 * 
 **********************************************/

// These should be used to get the DFs corresponding
// to half-way bounce back,
// after the pointer swap has been applied.

template<class T>
T LatticeD2Q9<T>::GetBouncedF0(int i, int j) const
{
    return mpF0[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF1(int i, int j) const
{
    return mpF2[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF2(int i, int j) const
{
    return mpF1[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF3(int i, int j) const
{
    return mpF4[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF4(int i, int j) const
{
    return mpF3[scalar_index(i, j + 1)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF5(int i, int j) const
{
    return mpF6[scalar_index(i, j)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF6(int i, int j) const
{
    return mpF5[scalar_index(i + 1, j + 1)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF7(int i, int j) const
{
    return mpF8[scalar_index(i + 1, j)];
}

template<class T>
T LatticeD2Q9<T>::GetBouncedF8(int i, int j) const
{
    return mpF7[scalar_index(i, j + 1)];
}

/**********************************************
 * 
 *  Standard setters
 * 
 **********************************************/

template<class T>
void LatticeD2Q9<T>::SetF0(T f_, int i, int j)
{
    mpF0[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF1(T f_, int i, int j)
{
    mpF2[scalar_index(i + 1, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF2(T f_, int i, int j)
{
    mpF1[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF3(T f_, int i, int j)
{
    mpF4[scalar_index(i, j + 1)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF4(T f_, int i, int j)
{
    mpF3[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF5(T f_, int i, int j)
{
    mpF6[scalar_index(i + 1, j + 1)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF6(T f_, int i, int j)
{
    mpF5[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF7(T f_, int i, int j)
{
    mpF8[scalar_index(i, j + 1)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetF8(T f_, int i, int j)
{
    mpF7[scalar_index(i + 1, j)] = f_;
}

/**********************************************
 * 
 *  Initial setters
 * 
 **********************************************/

template<class T>
void LatticeD2Q9<T>::SetInitF0(T f_, int i, int j)
{
    mpF0[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF1(T f_, int i, int j)
{
    mpF1[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF2(T f_, int i, int j)
{
    mpF2[scalar_index(i + 1, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF3(T f_, int i, int j)
{
    mpF3[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF4(T f_, int i, int j)
{  
    mpF4[scalar_index(i, j + 1)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF5(T f_, int i, int j)
{
    mpF5[scalar_index(i, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF6(T f_, int i, int j)
{
    mpF6[scalar_index(i + 1, j + 1)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF7(T f_, int i, int j)
{
    mpF7[scalar_index(i + 1, j)] = f_;
}

template<class T>
void LatticeD2Q9<T>::SetInitF8(T f_, int i, int j)
{
    mpF8[scalar_index(i, j + 1)] = f_;
}

/**********************************************
 * 
 *  Pointer swapping.
 * 
 **********************************************/

template<class T>
void LatticeD2Q9<T>::SwapPointers()
{
    std::swap(mpF1, mpF2);
    std::swap(mpF3, mpF4);
    std::swap(mpF5, mpF6);
    std::swap(mpF7, mpF8);
}

/**********************************************
 * 
 *  Display and saving utilities
 * 
 **********************************************/

template<class T>
void LatticeD2Q9<T>::DisplayPointerInfo() const
{
    std::cout << "mpF0   = " << mpF0   << "\n";
    std::cout << "mpF1   = " << mpF1   << "\n";
    std::cout << "mpF2   = " << mpF2   << "\n";
    std::cout << "mpF3   = " << mpF3   << "\n";
    std::cout << "mpF4   = " << mpF4   << "\n";
}

template<class T>
void LatticeD2Q9<T>::DisplayInfo() const
{
    std::cout << "Lattice D2Q5\n";
    std::cout << "Size X = " << mSizeX << "\n";
    std::cout << "Size Y = " << mSizeY << "\n";
}

/**
 * @brief Write the macroscopic data to CSV files with optional runId and timestep.
 * 
 * @param path The string path to the directory where the CSV files will be saved.
 * @param runId An optional string identifier for the run. Default is an empty string.
 * @param timestep An optional integer representing the timestep. Default is 0.
 */
template<class T>
void LatticeD2Q9<T>::WriteToCSV(const std::string& path, const char letter, const std::string& runId, const int timestep, const int process_number) const {

    // Construct base file paths with runId and timestep
    std::string basepath = construct_basepath(path, runId, timestep, letter, process_number);

    // Write each data array to its respective CSV file
    write_csv(basepath + "0.csv", &LatticeD2Q9<T>::GetF0);
    write_csv(basepath + "1.csv", &LatticeD2Q9<T>::GetF1);
    write_csv(basepath + "2.csv", &LatticeD2Q9<T>::GetF2);
    write_csv(basepath + "3.csv", &LatticeD2Q9<T>::GetF3);
    write_csv(basepath + "4.csv", &LatticeD2Q9<T>::GetF4);
    write_csv(basepath + "5.csv", &LatticeD2Q9<T>::GetF5);
    write_csv(basepath + "6.csv", &LatticeD2Q9<T>::GetF6);
    write_csv(basepath + "7.csv", &LatticeD2Q9<T>::GetF7);
    write_csv(basepath + "8.csv", &LatticeD2Q9<T>::GetF8);
}

/**********************************************
 * 
 *  Private member function definitions.
 * 
 **********************************************/

/**
 * @brief Construct the file name from path, runid, and timestep.
 * 
 * @param path The string path to the directory where the CSV files will be saved.
 * @param runId An optional string identifier for the run. Default is an empty string.
 * @param timestep An optional integer representing the timestep. Default is 0.
 * @return string basepath.
 */
 template<class T>
std::string LatticeD2Q9<T>::construct_basepath(const std::string& path, const std::string& runId, const int timestep, const char letter, const int process_number) const
{
    std::ostringstream oss;
    oss << path << "/" << runId << "_p" << std::setw(3) << std::setfill('0') << process_number << "_t" << std::setw(9) << std::setfill('0') << timestep << "_" << letter;
    return oss.str();
}

/**
 * @brief Write a given array to a file at a given path.
 *
 * @param path The path of the file (directory & filename).
 * @param mPtr Pointer to the array to be saved.
 */
 template<class T>
void LatticeD2Q9<T>::write_csv(const std::string& path, T (LatticeD2Q9<T>::*get_func)(int, int) const) const
{
    std::ofstream file(path);
    if (file.is_open()) 
    {
        for (int j = 0; j < mSizeY; ++j) 
        {
            for (int i = 0; i < mSizeX; ++i) 
            {
                file << (this->*get_func)(i, j);
                if (i < mSizeX - 1) 
                {
                    file << ",";
                }
            }
            file << "\n";
        }
        file.close();
    } 
    else 
    {
        throw std::runtime_error("Unable to open file: " + path);
    }
}

// Explicit template instantiation.
//template class LatticeD2Q9<_Float16>;
template class LatticeD2Q9<float>;
template class LatticeD2Q9<double>;