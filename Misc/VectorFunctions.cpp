//
// Created by bartmanski on 07/04/17.
//

#include "VectorFunctions.hpp"

vector<double> operator+(const double& alpha, const vector<double>& v)
{
    vector<double> output(v.size());
    for (unsigned i=0; i<v.size(); i++)
    {
        output[i] = v[i] + alpha;
    }

    return output;
}

vector<double> operator+(const vector<double>& v, const double& alpha)
{
    vector<double> output(v.size());
    for (unsigned i=0; i<v.size(); i++)
    {
        output[i] = v[i] + alpha;
    }

    return output;
}

vector<double> operator+(const vector<double>& v1, const vector<double>& v2)
{
    assert(v1.size() == v2.size());
    vector<double> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] + v2[i];
    }

    return output;
}

vector<unsigned> operator+(const vector<unsigned>& v1, const vector<unsigned>& v2)
{
    assert(v1.size() == v2.size());
    vector<unsigned> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] + v2[i];
    }

    return output;
}

vector<double> operator+(const vector<double>& v1, const vector<unsigned>& v2)
{
    assert(v1.size() == v2.size());
    vector<double> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] + v2[i];
    }

    return output;
}

vector<double> operator+(const vector<unsigned>& v1, const vector<double>& v2)
{
    assert(v1.size() == v2.size());
    vector<double> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] + v2[i];
    }

    return output;
}

vector<vector<unsigned> > operator*(const unsigned& alpha, const vector<vector<unsigned>>& m)
{
    vector<vector<unsigned> > output(m.size());
    for (unsigned i = 0; i < m.size(); i++)
    {
        output[i].resize(m[i].size());
    }
    for (unsigned i = 0; i < m.size(); i++)
    {
        for (unsigned j = 0; j < m[i].size(); j++)
        {
            output[i][j] = alpha * m[i][j];
        }
    }

    return output;
}

vector<double> operator/ (const vector<double>& v, const double& alpha)
{
    assert(alpha > 0);
    vector<double> output(v.size());
    for (unsigned i=0; i<v.size(); i++)
    {
        output[i] = v[i] / alpha;
    }
    return output;
}

vector<unsigned> floor_div(const vector<unsigned>& v, const unsigned& alpha)
{
    assert(alpha > 0);
    vector<unsigned> output(v.size());
    for (unsigned i=0; i<v.size(); i++)
    {
        output[i] = v[i] / alpha;
    }
    return output;
}

vector<double> linspace(double a, double b, unsigned n)
{
    vector<double> array;
    array.resize(n);
    double step = 0;
    if (n>1) { step = (b - a) / (n - 1); }
    for (unsigned i=0; i<n; i++)
    {
        array[i] = a + step * i;
    }
    return array;
}

vector<unsigned> ones(unsigned length)
{
    vector<unsigned> output = vector<unsigned>(length, 1);
    return output;
}

vector<vector<unsigned>> ones(unsigned num_cols, unsigned num_rows)
{
    vector<vector<unsigned> > output(num_rows);
    for (vector<unsigned> &row : output)
    {
        row = vector<unsigned>(num_cols, 1);
    }

    return output;
}

vector<vector<unsigned>> ones(vector<unsigned> shape)
{
    assert(!shape.empty());
    if (shape.size() == 1) { shape.push_back(1); }
    vector<vector<unsigned> > output(shape[1]);
    for (vector<unsigned> &row : output)
    {
        row = vector<unsigned>(shape[0], 1);
    }

    return output;
}
