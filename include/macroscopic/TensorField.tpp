#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include "macroscopic/TensorField.h"

template <typename T>
TensorField<T>::TensorField(int size_x, int size_y, int size_z, std::string id_string, std::string run_id, std::string savepath, int proc_num) 
    : mSizeX(size_x), mSizeY(size_y), mSizeZ(size_z), mGridSize(size_x*size_y*size_z), mIDString(id_string), mRunID(run_id), mSaveDirectory(savepath), mProcessNumber(proc_num), mArray(mRows*mCols*size_x*size_y*size_z) {}

template <typename T>
T TensorField<T>::GetValue(int a, int b, int i, int j, int k) const
{
    return mArray[idx(a, b, i, j, k)];
}

template <typename T>
void TensorField<T>::SetValue(T value, int a, int b, int i, int j, int k)
{
    mArray[idx(a, b, i, j, k)] = value;
}

template <typename T>
T TensorField<T>::GetValueWrap(int a, int b, int i, int j, int k) const
{
    int ii = (i + mSizeX) % mSizeX;
    int jj = (j + mSizeY) % mSizeY;
    int kk = (k + mSizeZ) % mSizeZ;
    return mArray[idx(a, b, ii, jj, kk)];
}

template <typename T>
void TensorField<T>::SetValueWrap(T value, int a, int b, int i, int j, int k)
{
    int ii = (i + mSizeX) % mSizeX;
    int jj = (j + mSizeY) % mSizeY;
    int kk = (k + mSizeZ) % mSizeZ;
    mArray[idx(a, b, ii, jj, kk)] = value;
}

template <typename T>
void TensorField<T>::AddToValue(T additive_value, int a, int b, int i, int j, int k)
{
    mArray[idx(a, b, i, j, k)] += additive_value;
}

template <typename T>
void TensorField<T>::SetToConstantValue(T value, int a, int b)
{
    for (int k = 0; k < mSizeZ; ++k)
    {
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                SetValue(value, a, b, i, j, k);
            }
        }
    }
}

template <typename T>
void TensorField<T>::SetToConstantValue(T value)
{
    for (int a = 0; a < mRows; ++a)
    {
        for (int b = 0; b < mCols; ++b)
        {
            SetToConstantValue(value, a, b);
        }
    }
}

template <typename T>
void TensorField<T>::SetLinearGradientZ(int a, int b, T bottom, T top)
{
    for (int k = 0; k < mSizeZ; ++k)
    {
        T value = bottom + k * (top - bottom) / (mSizeZ - 1);
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                SetValue(value, a, b, i, j, k);
            }
        }
    }
}

template <typename T>
T TensorField<T>::ComputeAverage() const
{
    double sum = 0.0;
    int total = 0;

    for (int a = 0; a < mRows; ++a)
    {
        for (int b = 0; b < mCols; ++b)
        {
            sum += ComputeAverage(a, b);
            total++;
        }
    }

    double average = sum / static_cast<double>(total);
    return static_cast<T>(average);
}

template <typename T>
T TensorField<T>::ComputeAverage(int a, int b) const
{
    double sum = 0.0;
    long total = 0;
    for (int k = 0; k < mSizeZ; ++k)
    {
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                sum += static_cast<double>(GetValue(a, b, i, j, k));
                total++;
            }
        }
    }
    double average = sum / static_cast<double>(total);
    return static_cast<T>(average);
}

template <typename T>
void TensorField<T>::DisplayInfo() const
{   
    std::cout << "--------------------------------------------\n";
    std::cout << "Type of variable:  " << mIDString << "\n";
    std::cout << "Run ID:            " << mRunID << "\n";
    std::cout << "Process number     " << mProcessNumber << "\n";
    std::cout << "Grid dimensions:   ";
    std::cout << mSizeX << " x " << mSizeY << " x " << mSizeZ << "\n";
    std::cout << "Tensor dimensions: " << mRows << " x " << mCols << "\n";
    std::cout << "Saving to:         " << mSaveDirectory << "\n";
    std::cout << "--------------------------------------------\n";
}

template <typename T>
void TensorField<T>::WriteToTextFile(const int timestep, int a, int b) const
{
    std::string basepath = construct_basepath(timestep, a, b);
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
                    file << GetValue(a, b, i, j, k);
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
void TensorField<T>::SetRunID(std::string run_id)
{
    mRunID = run_id;
}

template <typename T>
void TensorField<T>::WriteToTextFile(const int timestep) const
{
    for (int a = 0; a < mRows; ++a)
    {
        WriteToTextFile(timestep, a);
    }
}

/**
 * @brief Construct the file name from path, runid, and timestep.
 * 
 * @param timestep An optional integer representing the timestep. Default is 0.
 * @return string basepath.
 */
 template<typename T>
std::string TensorField<T>::construct_basepath(const int timestep, int a, int b) const
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

    // Construct ID based on b.
    if (b == 0)
    {
        id_string += "x";
    }
    else if (b == 1)
    {
        id_string += "y";
    }
    else if (b == 2)
    {
        id_string += "z";
    }
    else 
    {
        throw std::invalid_argument("b must be 0, 1, or 2");
    }

    // Construct the path.
    std::ostringstream oss;
    oss << mSaveDirectory << "/" << mRunID << "_p" << std::setw(3) << std::setfill('0') << mProcessNumber << "_t" << std::setw(9) << std::setfill('0') << timestep << "_" << id_string;
    return oss.str();
}
