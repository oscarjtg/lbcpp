#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "macroscopic/MacroscopicVariable.h"

template <class T>
MacroscopicVariable<T>::MacroscopicVariable(int size_x, int size_y, int size_z, std::string id_string, std::string run_id, std::string savepath, int proc_num) 
    : mSizeX(size_x), mSizeY(size_y), mSizeZ(size_z), mGridSize(size_x*size_y*size_z), mIDString(id_string), mRunID(run_id), mSaveDirectory(savepath), mProcessNumber(proc_num), mArray(size_x*size_y*size_z) {}

template <class T>
T MacroscopicVariable<T>::GetValue(int i, int j, int k) const
{
    return mArray[idx(i, j, k)];
}

template <class T>
void MacroscopicVariable<T>::SetValue(T value, int i, int j, int k)
{
    mArray[idx(i, j, k)] = value;
}

template <class T>
void MacroscopicVariable<T>::AddToValue(T additive_value, int i, int j, int k)
{
    mArray[idx(i, j, k)] += additive_value;
}

template <class T>
void MacroscopicVariable<T>::SetToConstantValue(T value)
{
    for (int k = 0; k < mSizeZ; ++k)
    {
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                this->SetValue(value, i, j, k);
            }
        }
    }
}

template <class T>
void MacroscopicVariable<T>::SetLinearGradientZ(T bottom, T top)
{
    for (int k = 0; k < mSizeZ; ++k)
    {
        T value = bottom + k * (top - bottom) / (mSizeZ - 1);
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                this->SetValue(value, i, j, k);
            }
        }
    }
}

template <class T>
T MacroscopicVariable<T>::ComputeAverage() const
{
    double sum = 0.0;
    int total = 0;
    for (int k = 0; k < mSizeZ; ++k)
    {
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                sum += static_cast<double>(this->GetValue(i, j, k));
                total++;
            }
        }
    }
    double average = sum / static_cast<double>(total);
    return static_cast<T>(average);
}

template <class T>
void MacroscopicVariable<T>::DisplayInfo() const
{   
    std::cout << "--------------------------------------------\n";
    std::cout << "Type of variable: " << mIDString << "\n";
    std::cout << "Run ID:           " << mRunID << "\n";
    std::cout << "Process number    " << mProcessNumber << "\n";
    std::cout << "Grid dimensions:  ";
    std::cout << mSizeX << " x " << mSizeY << " x " << mSizeZ << "\n";
    std::cout << "Saving to:        " << mSaveDirectory << "\n";
    std::cout << "--------------------------------------------\n";
}

template <class T>
void MacroscopicVariable<T>::WriteToTextFile(const int timestep) const
{
    std::string basepath = construct_basepath(timestep);
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
                    file << GetValue(i, j, k);
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

/**
 * @brief Construct the file name from path, runid, and timestep.
 * 
 * @param timestep An optional integer representing the timestep. Default is 0.
 * @return string basepath.
 */
 template<class T>
std::string MacroscopicVariable<T>::construct_basepath(const int timestep) const
{
    std::ostringstream oss;
    oss << mSaveDirectory << "/" << mRunID << "_p" << std::setw(3) << std::setfill('0') << mProcessNumber << "_t" << std::setw(9) << std::setfill('0') << timestep << "_" << mIDString;
    return oss.str();
}
