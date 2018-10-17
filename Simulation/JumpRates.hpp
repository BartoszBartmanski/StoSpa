//
// Created by bartmanski on 16/10/18.
//

#ifndef STOSPA_JUMPRATES_HPP
#define STOSPA_JUMPRATES_HPP


#include "cmath"
#include "algorithm"
#include "cassert"

using namespace std;

class JumpRate
{
protected:

    double mKappa = 1.0;

    double mH = 1.0;
public:
    JumpRate() = default;

    virtual ~JumpRate()= default;

    virtual double GetLambda0()=0;

    virtual double GetLambda1()=0;

    virtual double GetLambda2()=0;

    virtual double GetLambda3()=0;

};

class FDM : public JumpRate
{
private:
    double mAlpha = 0.0;

public:
    FDM(double kappa, double length, double alpha);

    double GetLambda0() override;

    double GetLambda1() override;

    double GetLambda2() override;

    double GetLambda3() override;

};

class FEM : public JumpRate
{
public:
    FEM(double kappa, double length);

    double GetLambda0() override;

    double GetLambda1() override;

    double GetLambda2() override;

    double GetLambda3() override;

};

class FVM : public JumpRate
{
public:
    FVM(double kappa, double length);

    double GetLambda0() override;

    double GetLambda1() override;

    double GetLambda2() override;

    double GetLambda3() override;

};

class FET : public JumpRate
{
protected:
    double mBetaX=0.0;

    double mBetaY=0.0;

    unsigned mTruncOrder=1000;

    double mLambda0=1.0;

    double mTheta1=0.25;

    double mTheta2=0.0;

    double mTheta3=0.25;

public:
    FET() = default;

    FET(double kappa, double length, double beta_x, double beta_y, unsigned truncation_order);

    virtual double GetTheta1();

    virtual double GetTheta2();

    virtual double GetTheta3();

    double GetLambda0() override;

    double GetLambda1() override;

    double GetLambda2() override;

    double GetLambda3() override;

};

class FETUniform : public FET
{
public:
    FETUniform(double kappa, double length, double beta_x, double beta_y, unsigned truncation_order);

    double GetTheta1() override;

    double GetTheta2() override;

    double GetTheta3() override;

    double GetLambda0() override;

};


#endif //STOSPA_JUMPRATES_HPP
