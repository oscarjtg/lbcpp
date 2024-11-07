#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <array>

#include "lattice/AbstractLattice.h"

// Overriden constructor
template <class T>
AbstractLattice<T>::AbstractLattice(int size_x, int size_y, int size_z, int grid_size, int dims, int vels, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : mSizeX(size_x), mSizeY(size_y), mSizeZ(size_z), mGridSize(grid_size), mDims(dims), mVels(vels), mIDString(id_string), mRunID(run_id), mSaveDirectory(savepath), mProcessNumber(proc_num)
    {
        if (dims == 2 && size_y != 1)
        {
            throw std::invalid_argument("Invalid size_y (number of grid points along y direction) for dims = 2 grid. Should be 1.");
        }
        else if (dims == 2 && vels == 4)
        {
            mEX = {1, 0, -1, 0};
            mEY = {0, 0, 0, 0};
            mEZ = {0, 1, 0, -1};
            mW = {1./4., 1./4., 1./4., 1./4.};
            mQRev = {2, 3, 0, 1};
            mCSI = 2;
        }
        else if (dims == 2 && vels == 5)
        {
            mEX = {0, 1, 0, -1, 0};
            mEY = {0, 0, 0, 0, 0};
            mEZ = {0, 0, 1, 0, -1};
            mW = {1./3., 1./6., 1./6., 1./6., 1./6.};
            mQRev = {0, 3, 4, 1, 2};
            mCSI = 3;
        }
        else if (dims == 2 && vels == 9)
        {
            mEX = {0, 1, 0, -1, 0, 1, -1, -1, 1};
            mEY = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            mEZ = {0, 0, 1, 0, -1, 1, 1, -1, -1};
            mW = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
            mQRev = {0, 3, 4, 1, 2, 7, 8, 5, 6};
            mCSI = 3;
        }
        else if (dims == 3 && vels == 6)
        {
            mEX = {1, -1, 0, 0, 0, 0};
            mEY = {0, 0, 1, -1, 0, 0};
            mEZ = {0, 0, 0, 0, 1, -1};
            mW = {1./6., 1./6., 1./6., 1./6., 1./6., 1./6.};
            mQRev = {1, 0, 3, 2, 5, 4};
            mCSI = 2;
        }
        else if (dims == 3 && vels == 7)
        {
            mEX = {0, 1, -1, 0, 0, 0, 0};
            mEY = {0, 0, 0, 1, -1, 0, 0};
            mEZ = {0, 0, 0, 0, 0, 1, -1};
            mW = {1./4., 1./8., 1./8., 1./8., 1./8., 1./8., 1./8.};
            mQRev = {0, 2, 1, 4, 3, 6, 5};
            mCSI = 2;
        }
        else if (dims == 3 && vels == 15)
        {
            mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
            mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1};
            mEZ = {0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1};
            mW = {2./9., 1./9., 1./9., 1./9., 1./9., 1./9., 1./9., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
            mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};
            mCSI = 3;
        }
        else if (dims == 3 && vels == 19)
        {
            mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0};
            mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
            mEZ = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
            mW = {1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,  1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
            mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
            mCSI = 3;
        }
        else if (dims == 3 && vels == 27)
        {
            mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
            mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1};
            mEZ = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1};
            mW = {8./27., 2./27., 2./27., 2./27., 2./27., 2./27., 2./27., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.};
            mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
            mCSI = 3;
        }
        else
        {
            throw std::invalid_argument("Invalid combination of dims and vels");
        }
    }

template <class T>
void AbstractLattice<T>::DisplayLatticeParameters() const
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Type of variable:  " << mIDString << "\n";
    std::cout << "Run ID:            " << mRunID << "\n";
    std::cout << "Process number     " << mProcessNumber << "\n";
    std::cout << "Grid dimensions:   ";
    std::cout << mSizeX << " x " << mSizeY << " x " << mSizeZ << "\n";
    std::cout << "Lattice type:      " << "D" << mDims << "Q" << mVels << "\n";
    std::cout << "Lattice 1 / c_s^2: " << mCSI << "\n";
    std::cout << "i, ex, ey, ez, w, irev\n";
    double sum = 0.0;
    bool all_okay = true;
    for (int i = 0; i < mVels; ++i)
    {
        // Display values.
        std::cout << i << ", ";
        std::cout << mEX[i] << ", ";
        std::cout << mEY[i] << ", ";
        std::cout << mEZ[i] << ", ";
        std::cout << mW[i] << ", ";
        std::cout << mQRev[i] << "\n";

        // Increment sum of lattice weights.
        sum += mW[i];

        // Check that e_ibar = -e_i.
        if (all_okay && (
            (mEX[i] != -mEX[mQRev[i]]) 
            || (mEY[i] != -mEY[mQRev[i]]) 
            || (mEZ[i] != -mEZ[mQRev[i]])
            ))
        {
            all_okay = false;
        }
    }
    std::cout << "Sum of lattice weights = " << sum << "\n";
    std::cout << "All lattice values consistent = ";
    std::cout << std::boolalpha << all_okay << "\n";
    std::cout << "Saving to: " << mSaveDirectory << "\n";
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

/**
 * @brief Construct the file name from path, runid, and timestep.
 * 
 * @param timestep An optional integer representing the timestep. Default is 0.
 * @return string basepath.
 */
 template<class T>
std::string AbstractLattice<T>::construct_basepath(const int timestep) const
{
    std::ostringstream oss;
    oss << mSaveDirectory << "/" << mRunID << "_p" << std::setw(3) << std::setfill('0') << mProcessNumber << "_t" << std::setw(9) << std::setfill('0') << timestep << "_" << mIDString;
    return oss.str();
}

/**
 * @brief Write a given array to a file at a given path.
 *
 * @param path The path of the file (directory & filename).
 * @param mPtr Pointer to the array to be saved.
 */
template<class T>
void AbstractLattice<T>::write_txt(const std::string& path, T (AbstractLattice<T>::*get_func)(int, int, int, int) const, int q) const
{
    std::string fullpath = path + std::to_string(q) + ".txt";
    std::ofstream file(fullpath);
    if (file.is_open()) 
    {
        file << mSizeX << "," << mSizeY << "," << mSizeZ;
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
                    file << (this->*get_func)(q, i, j, k);
                }
            }
        }
        file.close();
    } 
    else 
    {
        throw std::runtime_error("Unable to open file: " + path);
    }
}
