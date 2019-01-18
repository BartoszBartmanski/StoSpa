//
// Created by bartmanski on 11/01/19.
//

#ifndef STOSPA_GROWTHRATES_HPP
#define STOSPA_GROWTHRATES_HPP

#include <cassert>
#include <cmath>

class GrowthRate
{
protected:
    double mRateConstant=0;

public:
    GrowthRate() = default;

    virtual ~GrowthRate()= default;

    virtual double GetFunction(const double& time)=0;
};

class Exponential : public GrowthRate
{
public:
    explicit Exponential(double rate);

    double GetFunction(const double& time) override;
};


#endif //STOSPA_GROWTHRATES_HPP
