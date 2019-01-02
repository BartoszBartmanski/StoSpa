//
// Created by bartmanski on 16/10/18.
//

#include "JumpRates.hpp"

double JumpRate::GetLambda(vector<int> direction)
{
    double value;
    if (direction[1] == 0)  // Horizontal jumps
    {
        value = this->GetLambda1();
    }
    else if (direction[0] == 0)  // Vertical jumps
    {
        value = this->GetLambda3();
    }
    else  // Diagonal jumps
    {
        value = this->GetLambda2();
    }
    return value;
}

JumpRate1d::JumpRate1d(vector<double> voxel_dims)
{
    assert(!voxel_dims.empty());
    mH = voxel_dims[0];
}

double JumpRate1d::GetLambda0()
{
    double lambda_0 = 2.0 / pow(mH, 2);
    return lambda_0;
}

double JumpRate1d::GetLambda1()
{
    double lambda_1 = 1.0 / pow(mH, 2);
    return lambda_1;
}

double JumpRate1d::GetLambda2()
{
    return 0.0;
}

double JumpRate1d::GetLambda3()
{
    return 0.0;
}

string JumpRate::GetMethod()
{
    return mMethod;
}

FDM::FDM(vector<double> voxel_dims, double alpha)
{
    assert(voxel_dims.size() == 2);
    mMethod = "fdm";
    mKappa = voxel_dims[0]/voxel_dims[1];
    mH = voxel_dims[1];

    double upper_bound = min(mKappa, 1.0/mKappa);
    if (alpha <= 0.0 or alpha >= upper_bound)
    {
        throw runtime_error("Parameter alpha is not within a suitable range!");
    }
    mAlpha = alpha;
}

double FDM::GetLambda0()
{
    return 2.0 * (pow(mKappa, 2) - mAlpha * mKappa + 1) / pow(mKappa * mH, 2);
}

double FDM::GetLambda1()
{
    return (1.0 - mKappa * mAlpha) / pow(mKappa * mH, 2);
}

double FDM::GetLambda2()
{
    return mAlpha / (2.0 * mKappa * pow( mH, 2));
}

double FDM::GetLambda3()
{
    return (pow(mKappa, 2) - mKappa * mAlpha) / pow(mKappa * mH, 2);
}

FEM::FEM(vector<double> voxel_dims)
{
    // Check that the parameters satisfy the discrete maximum principle
    assert(voxel_dims.size() == 2);
    mMethod = "fem";
    mKappa = voxel_dims[0] / voxel_dims[1];
    mH = voxel_dims[1];

    assert(mKappa >= 1.0/sqrt(2) and mKappa <= sqrt(2));
}

double FEM::GetLambda0()
{
    return 4.0 * (pow(mKappa, 2) + 1) / (3.0 * pow(mKappa * mH, 2));
}

double FEM::GetLambda1()
{
    return (2.0 - pow(mKappa, 2)) / (3.0 * pow(mKappa * mH, 2));
}

double FEM::GetLambda2()
{
    return (pow(mKappa, 2) + 1) / (6.0 * pow(mKappa * mH, 2));
}

double FEM::GetLambda3()
{
    return (2.0 * pow(mKappa, 2) - 1) / (3.0 * pow(mKappa * mH, 2));
}

FVM::FVM(vector<double> voxel_dims)
{
    assert(voxel_dims.size() == 2);
    mMethod = "fvm";
    mKappa = voxel_dims[0] / voxel_dims[1];
    mH = voxel_dims[1];
}

double FVM::GetLambda0()
{
    return 2.0 * (pow(mKappa, 2) + 1) / (pow(mKappa * mH, 2));
}

double FVM::GetLambda1()
{
    return 1.0 / pow(mKappa * mH, 2);
}

double FVM::GetLambda2()
{
    return 0.0;
}

double FVM::GetLambda3()
{
    return 1.0 / pow(mH, 2);
}

FET::FET(vector<double> voxel_dims, double beta_x, double beta_y, unsigned truncation_order)
{
    assert(voxel_dims.size() == 2);
    assert(beta_x >= 0.0 and beta_x <= 1.0);
    assert(beta_y >= 0.0 and beta_y <= 1.0);
    mMethod = "fet";
    mTruncOrder = truncation_order;
    mKappa = voxel_dims[0] / voxel_dims[1];
    mH = voxel_dims[1];
    mBetaX = beta_x;
    mBetaY = beta_y;
    mLambda0 = GetLambda0();
    mTheta1 = GetTheta1();
    mTheta3 = GetTheta3();
    mTheta2 = 0.25*(1.0 - 2*mTheta1 - 2*mTheta3);
}

FET::FET(vector<double> voxel_dims, vector<double> beta, unsigned truncation_order)
{
    assert(voxel_dims.size() == 2);
    assert(beta.size() == 2);
    assert(beta[0] >= 0.0 and beta[0] <= 1.0);
    assert(beta[1] >= 0.0 and beta[1] <= 1.0);
    mMethod = "fet";
    mTruncOrder = truncation_order;
    mKappa = voxel_dims[0] / voxel_dims[1];
    mH = voxel_dims[1];
    mBetaX = beta[0];
    mBetaY = beta[1];
    mLambda0 = GetLambda0();
    mTheta1 = GetTheta1();
    mTheta3 = GetTheta3();
    mTheta2 = 0.25*(1.0 - 2*mTheta1 - 2*mTheta3);
}

double FET::GetLambda0()
{
    double e_t = 0;

    // First we calculate the expected time to leave the voxel, then we invert this to calculate the propensity
    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2.0*j-1.0), 2) + pow((2.0*k-1.0), 2));
            double top = pow((-1), (j + k));
            double bottom = ((2.0*k-1.0) * (2.0*j-1.0) * sum_squares);

            e_t += top/bottom;
        }
    }

    return pow(M_PI, 4) / (64.0 * pow(mKappa * mH, 2) * e_t);
}

double FET::GetTheta1()
{
    double theta_1 = 0;

    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2*j-1), 2) + pow((2*k-1), 2));

            double val = 8.0 * pow((-1), (k + 1)) * (2.0 * k - 1.0);
            val /= (pow(M_PI, 2) * sum_squares * (2.0 * j - 1.0));
            theta_1 += sin((j - 0.5) * M_PI * mBetaY) * val;
        }
    }

    return theta_1;
}

double FET::GetTheta2()
{
    double theta2 = 0.0;

    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2*j-1), 2) + pow((2*k-1), 2));
            double val1 = pow(mKappa, 2) * (2*j-1)*(sin(M_PI*(mBetaX*k-0.5*mBetaX+k))+1.0) / (2 * k - 1);
            double val2 = (2*k-1)*(sin(M_PI*(mBetaY*j-0.5*mBetaY+j))+1.0) / (2 * j - 1);
            theta2 += 4.0*pow((-1), (j+k)) * (val1 + val2) / (pow(M_PI, 2) * sum_squares);
        }
    }

    return theta2;
}

double FET::GetTheta3()
{
    double theta3 = 0.0;

    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2*j-1), 2) + pow((2*k-1), 2));

            double val = 8.0 * pow((-1), (j+1)) * (2 * j - 1) * pow(mKappa, 2) / (pow(M_PI, 2) * sum_squares * (2 * k - 1));
            theta3 += sin((k - 0.5) * M_PI * mBetaX) * val;
        }
    }

    return theta3;
}

double FET::GetLambda1()
{
    return mTheta1 * mLambda0;
}

double FET::GetLambda2()
{
    return mTheta2 * mLambda0;
}

double FET::GetLambda3()
{
    return mTheta3 * mLambda0;
}


FETUniform::FETUniform(vector<double> voxel_dims, double beta_x, double beta_y, unsigned truncation_order)
{
    assert(voxel_dims.size() == 2);
    assert(beta_x >= 0.0 and beta_x <= 1.0);
    assert(beta_y >= 0.0 and beta_y <= 1.0);

    mMethod = "fetU";
    mTruncOrder = truncation_order;
    mKappa = voxel_dims[0] / voxel_dims[1];
    mH = voxel_dims[1];
    mBetaX = beta_x;
    mBetaY = beta_y;
    mLambda0 = GetLambda0();
    mTheta1 = GetTheta1();
    mTheta3 = GetTheta3();
    mTheta2 = 0.25*(1.0 - 2*mTheta1 - 2*mTheta3);
}

FETUniform::FETUniform(vector<double> voxel_dims, vector<double> beta, unsigned truncation_order)
{
    assert(voxel_dims.size() == 2);
    assert(beta.size() == 2);
    assert(beta[0] >= 0.0 and beta[0] <= 1.0);
    assert(beta[1] >= 0.0 and beta[1] <= 1.0);

    mMethod = "fetU";
    mTruncOrder = truncation_order;
    mKappa = voxel_dims[0] / voxel_dims[1];
    mH = voxel_dims[1];
    mBetaX = beta[0];
    mBetaY = beta[1];
    mLambda0 = GetLambda0();
    mTheta1 = GetTheta1();
    mTheta3 = GetTheta3();
    mTheta2 = 0.25*(1.0 - 2*mTheta1 - 2*mTheta3);
}

double FETUniform::GetLambda0()
{
    double e_t = 0;

    // First we calculate the expected time to leave the voxel, then we invert this to calculate the propensity
    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2.0*j-1.0), 2) + pow((2.0*k-1.0), 2));
            e_t += 1.0/(pow((2.0*k-1.0) * (2.0*j-1.0), 2) * sum_squares);
        }
    }

    return pow(M_PI, 6) / (64.0 * pow(mKappa * mH, 2) * e_t);
}

double FETUniform::GetTheta1()
{

    double theta_1 = 0;

    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2*j-1), 2) + pow((2*k-1), 2));

            double val = 32.0 * pow((-1), (j + 1));
            val /= (pow(M_PI, 4) * sum_squares * pow((2.0 * j - 1.0), 2));
            theta_1 += sin((j - 0.5) * M_PI * mBetaY) * val;
        }
    }

    return theta_1;
}

double FETUniform::GetTheta2()
{
    double theta2 = 0.0;

    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2*j-1), 2) + pow((2*k-1), 2));
            double val1 = pow(mKappa, 2) * (sin(M_PI*(mBetaX*k-0.5*mBetaX+k))+1.0) / pow((2 * k - 1), 2);
            double val2 = (sin(M_PI*(mBetaY*j-0.5*mBetaY+j))+1.0) / pow((2 * j - 1), 2);
            theta2 += 16.0 * (val1 + val2) / (pow(M_PI, 4) * sum_squares);
        }
    }

    return theta2;
}

double FETUniform::GetTheta3()
{
    double theta3 = 0.0;

    for (unsigned j = 1; j < mTruncOrder + 1; j++)
    {
        for (unsigned k = 1; k < mTruncOrder + 1; k++)
        {
            double sum_squares = (pow(mKappa * (2*j-1), 2) + pow((2*k-1), 2));

            double val = 32.0 * pow((-1), (k+1)) * pow(mKappa, 2);
            val /= (pow(M_PI, 4) * sum_squares * pow((2 * k - 1), 2));
            theta3 += sin((k - 0.5) * M_PI * mBetaX) * val;
        }
    }

    return theta3;
}
