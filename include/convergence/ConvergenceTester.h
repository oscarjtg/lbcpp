#ifndef CONVERGENCETESTERDEF
#define CONVERGENCETESTERDEF

#include <cmath>

#include "util/FixedSizeContainer.h"

template <class T>
class ConvergenceTester
{
public:
    ConvergenceTester(int n) : mNumElements(n), mValues(n) {}

    ~ConvergenceTester() = default;

    void AddValueToList(T value)
    {
        mValues.AddElement(value);
    }

    T AverageValue()
    {
        return mValues.CalculateAverage();
    }

    bool HasConverged(double absolute_tolerance = 1.0e-3)
    {
        if (mValues.GetElement(mNumElements-1) == 0)
        {
            return false;
        }
        double difference, maximum_difference = 0.0;
        double current_value = mValues.GetElement(0);
        for (int i = 1; i < mNumElements; ++i)
        {
            difference = std::fabs(current_value);
            if (maximum_difference < difference)
            {
                maximum_difference = difference;
            }
        }
        return maximum_difference < absolute_tolerance;
    }

private:
    int mNumElements;
    FixedSizeContainer<double> mValues;
};

#endif // CONVERGENCETESTERDEF