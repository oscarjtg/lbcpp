#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <array>

#include "lattice/AbstractLattice.h"

// Overriden constructor
template <typename T, int ND, int NQ>
AbstractLattice<T, ND, NQ>::AbstractLattice(int size_x, int size_y, int size_z, int grid_size, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : mSizeX(size_x), mSizeY(size_y), mSizeZ(size_z), mGridSize(grid_size), mIDString(id_string), mRunID(run_id), mSaveDirectory(savepath), mProcessNumber(proc_num)
    {
        if (ND == 2 && size_y != 1)
        {
            throw std::invalid_argument("Invalid size_y (number of grid points along y direction) for ND = 2 grid. Should be 1.");
        }
        else if (ND == 2 && NQ == 4)
        {
            // do nothing
        }
        else if (ND == 2 && NQ == 5)
        {
            // do nothing
        }
        else if (ND == 2 && NQ == 9)
        {
            // do nothing
        }
        else if (ND == 3 && NQ == 6)
        {
            // do nothing
        }
        else if (ND == 3 && NQ == 7)
        {
            // do nothing
        }
        else if (ND == 3 && NQ == 15)
        {
            // do nothing
        }
        else if (ND == 3 && NQ == 19)
        {
            // do nothing
        }
        else if (ND == 3 && NQ == 27)
        {
            // do nothing
        }
        else
        {
            throw std::invalid_argument("Invalid combination of ND and NQ");
        }
    }

template <typename T, int ND, int NQ>
void AbstractLattice<T, ND, NQ>::DisplayLatticeParameters() const
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Type of variable:  " << mIDString << "\n";
    std::cout << "Run ID:            " << mRunID << "\n";
    std::cout << "Process number     " << mProcessNumber << "\n";
    std::cout << "Grid dimensions:   ";
    std::cout << mSizeX << " x " << mSizeY << " x " << mSizeZ << "\n";
    std::cout << "Lattice type:      " << "D" << ND << "Q" << NQ << "\n";
    std::cout << "Lattice 1 / c_s^2: " << mLP.mCSI << "\n";
    std::cout << "i, ex, ey, ez, w, irev\n";
    double sum = 0.0;
    bool all_okay = true;
    for (int i = 0; i < NQ; ++i)
    {
        // Display values.
        std::cout << i << ", ";
        std::cout << mLP.mEX[i] << ", ";
        std::cout << mLP.mEY[i] << ", ";
        std::cout << mLP.mEZ[i] << ", ";
        std::cout << mLP.mW[i] << ", ";
        std::cout << mLP.mQRev[i] << "\n";

        // Increment sum of lattice weights.
        sum += mLP.mW[i];

        // Check that e_ibar = -e_i.
        if (all_okay && (
            (mLP.mEX[i] != -mLP.mEX[mLP.mQRev[i]]) 
            || (mLP.mEY[i] != -mLP.mEY[mLP.mQRev[i]]) 
            || (mLP.mEZ[i] != -mLP.mEZ[mLP.mQRev[i]])
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
 template <typename T, int ND, int NQ>
std::string AbstractLattice<T, ND, NQ>::construct_basepath(const int timestep) const
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
template <typename T, int ND, int NQ>
void AbstractLattice<T, ND, NQ>::write_txt(const std::string& path, T (AbstractLattice<T, ND, NQ>::*get_func)(int, int, int, int) const, int q) const
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
