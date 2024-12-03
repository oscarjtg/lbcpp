#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include "macroscopic/VectorField.h"

template <typename T>
VectorField<T>::VectorField(int size_x, int size_y, int size_z, std::string id_string, std::string run_id, std::string savepath, int proc_num) 
    : mSizeX(size_x), mSizeY(size_y), mSizeZ(size_z), mGridSize(size_x*size_y*size_z), mIDString(id_string), mRunID(run_id), mSaveDirectory(savepath), mProcessNumber(proc_num), mArray(mDims*size_x*size_y*size_z) {}

template <typename T>
T VectorField<T>::GetValue(int a, int i, int j, int k) const
{
    return mArray[idx(a, i, j, k)];
}

template <typename T>
void VectorField<T>::SetValue(T value, int a, int i, int j, int k)
{
    mArray[idx(a, i, j, k)] = value;
}

template <typename T>
T VectorField<T>::GetValueWrap(int a, int i, int j, int k) const
{
    int ii = (i + mSizeX) % mSizeX;
    int jj = (j + mSizeY) % mSizeY;
    int kk = (k + mSizeZ) % mSizeZ;
    return mArray[idx(a, ii, jj, kk)];
}

template <typename T>
void VectorField<T>::SetValueWrap(T value, int a, int i, int j, int k)
{
    int ii = (i + mSizeX) % mSizeX;
    int jj = (j + mSizeY) % mSizeY;
    int kk = (k + mSizeZ) % mSizeZ;
    mArray[idx(a, ii, jj, kk)] = value;
}

template <typename T>
void VectorField<T>::AddToValue(T additive_value, int a, int i, int j, int k)
{
    mArray[idx(a, i, j, k)] += additive_value;
}

template <typename T>
void VectorField<T>::SetToConstantValue(T value, int a)
{
    for (int k = 0; k < mSizeZ; ++k)
    {
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                SetValue(value, a, i, j, k);
            }
        }
    }
}

template <typename T>
void VectorField<T>::SetToConstantValue(T value)
{
    for (int a = 0; a < mDims; ++a)
    {
        SetToConstantValue(value, a);
    }
}

template <typename T>
void VectorField<T>::SetLinearGradientZ(int a, T bottom, T top)
{
    for (int k = 0; k < mSizeZ; ++k)
    {
        T value = bottom + k * (top - bottom) / (mSizeZ - 1);
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                SetValue(value, a, i, j, k);
            }
        }
    }
}

template <typename T>
T VectorField<T>::ComputeAverage(int a) const
{
    double sum = 0.0;
    int total = 0;
    for (int k = 0; k < mSizeZ; ++k)
    {
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                sum += static_cast<double>(GetValue(a, i, j, k));
                total++;
            }
        }
    }
    double average = sum / static_cast<double>(total);
    return static_cast<T>(average);
}

template <typename T>
T VectorField<T>::ComputeAverage() const
{
    double sum = 0.0;
    int total = 0;

    for (int a = 0; a < mDims; ++a)
    {
        sum += ComputeAverage(a);
        total++;
    }

    double average = sum / static_cast<double>(total);
    return static_cast<T>(average);
}

template <typename T>
void VectorField<T>::DisplayInfo() const
{   
    std::cout << "--------------------------------------------\n";
    std::cout << "Type of variable:  " << mIDString << "\n";
    std::cout << "Run ID:            " << mRunID << "\n";
    std::cout << "Process number     " << mProcessNumber << "\n";
    std::cout << "Grid dimensions:   ";
    std::cout << mSizeX << " x " << mSizeY << " x " << mSizeZ << "\n";
    std::cout << "Vector dimensions: " << mDims << "\n";
    std::cout << "Saving to:         " << mSaveDirectory << "\n";
    std::cout << "--------------------------------------------\n";
}

template <typename T>
void VectorField<T>::WriteToTextFile(const int timestep, int a) const
{
    std::string basepath = construct_basepath(timestep, a);
    std::string fullpath = basepath + ".txt";
    std::ofstream file(fullpath);
    file << mSizeX << "," << mSizeY << "," << mSizeZ;
    if (file.is_open()) 
    {
        for (int k = 0; k < mSizeZ; ++k)
        {
            for (int j = 0; j < mSizeY; ++j) 
            {
                for (int i = 0; i < mSizeX; ++i) 
                {
                    if (i == 0 && j == 0 && k == 0)
                    {
                        file << "\n"; // End the first line.
                    }
                    else
                    {
                        file << ",";
                    }
                    file << GetValue(a, i, j, k);
                }
            }
        }
        file.close();
    } 
    else 
    {
        throw std::runtime_error("Unable to open file: " + fullpath);
    }
}

template <typename T>
void VectorField<T>::WriteToTextFile(const int timestep) const
{
    for (int a = 0; a < mDims; ++a)
    {
        WriteToTextFile(timestep, a);
    }
}

template <typename T>
void VectorField<T>::SetRunID(std::string run_id)
{
    mRunID = run_id;
}

/**
 * @brief Construct the file name from path, runid, and timestep.
 * 
 * @param timestep An optional integer representing the timestep. Default is 0.
 * @return string basepath.
 */
 template<typename T>
std::string VectorField<T>::construct_basepath(const int timestep, int a) const
{
    // Construct ID based on a.
    std::string id_string;
    if (a == 0)
    {
        id_string = mIDString + "x";
    }
    else if (a == 1)
    {
        id_string = mIDString + "y";
    }
    else if (a == 2)
    {
        id_string = mIDString + "z";
    }
    else 
    {
        throw std::invalid_argument("a must be 0, 1, or 2");
    }

    // Construct the path.
    std::ostringstream oss;
    oss << mSaveDirectory << "/" << mRunID << "_p" << std::setw(3) << std::setfill('0') << mProcessNumber << "_t" << std::setw(9) << std::setfill('0') << timestep << "_" << id_string;
    return oss.str();
}
