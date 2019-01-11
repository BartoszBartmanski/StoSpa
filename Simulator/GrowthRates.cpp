//
// Created by bartmanski on 11/01/19.
//

#include "GrowthRates.hpp"

Exponential::Exponential(unsigned dims, double rate)
{
    assert(rate >= 0);
    assert(dims == 1 or dims == 2);
    mRateConstant = rate;
    mDims = dims;
}

double Exponential::GetFunction(const double& time)
{
    return exp(mRateConstant * time * mDims);
}