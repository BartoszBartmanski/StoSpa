//
// Created by bartmanski on 11/01/19.
//

#include "GrowthRates.hpp"

Exponential::Exponential(double rate)
{
    assert(rate >= 0);
    mRateConstant = rate;
}

double Exponential::GetFunction(const double& time)
{
    return exp(mRateConstant * time);
}