//
// Created by bartmanski on 07/04/17.
//

#ifndef  STOSPA_VECTORFUNCTIONS_HPP
#define  STOSPA_VECTORFUNCTIONS_HPP

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <cassert>
#include <memory>

using namespace std;

/**
 * Produces a vector of linearly spaced points from the given lower bound to the given upper bound.
 * @param a - lower bound
 * @param b - upper bound
 * @param n - number of points in the vector (including the lower and upper bounds)
 * @return vector of linearly spaced points
 */
vector<double> linspace(double a, double b, unsigned n);

/**
 * Returns a vector or a matrix filled with ones
 * @param length - length of the vector
 * @return a vector or vector of vectors
 */
vector<unsigned> ones(unsigned length);
vector<vector<unsigned>> ones(unsigned num_cols, unsigned num_rows);
vector<vector<unsigned>> ones(vector<unsigned> shape);

/**
 * Prints a vector
 * @param vec - vector to be printed
 * @param name - name of the vector
 */
template<typename T>
void print_array(vector<T> vec, const string& name="vector")
{
    cout << endl << name << " = ";

    for (T elem : vec)
    {
        cout << elem << " ";
    }
    cout << endl;
}

/**
 * Prints a matrix
 * @param matrix - matrix to be printed
 * @param name - name of the matrix
 */
template<typename T>
void print_array(vector<vector<T>> matrix, const string& name="matrix")
{
    cout << endl << name << " = " << endl;
    for (const vector<double>& row : matrix)
    {
        for (double elem : row)
        {
            cout << elem << " ";
        }
        cout << endl;
    }
    cout << endl;
}

/**
 * Saves the vector to a given file. Adds to the file if the file exists already.
 * @param vec - vector to be saved
 * @param path_to_file - path to the file
 */
template<typename T>
void save_vector(const vector<T>& vec, const string& path_to_file)
{
    ofstream vec_file(path_to_file, ios_base::app);
    ostream_iterator<T> output_iterator(vec_file, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    vec_file << endl;
    vec_file.close();
}

template<typename T>
void save_vector(const vector<T>& vec, const unique_ptr<ofstream>& handle)
{
    ostream_iterator<T> output_iterator(*handle, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    *handle << endl;
}

/**
 * Sum of all the elements in a vector.
 * @param vec - vector for which to sum all the elements
 * @return sum of all the elements of the vector
 */
template<typename T>
T vec_sum(const vector<T>& vec)
{
    T sum = 0.0;
    for (T elem : vec)
    {
        sum += elem;
    }

    return sum;
}

/**
 * Finds the smallest element in a vector
 * @param vec - vector in which to find the minimum
 * @return a pair of the value and the index of the smallest element
 */
template<typename T>
pair<T, unsigned> find_min(const vector<T>& vec)
{
    auto result = min_element(vec.begin(), vec.end());
    unsigned index = unsigned(distance(vec.begin(), result));
    auto val_idx = make_pair(*result, index);
    return val_idx;
}

/**
 * Checks whether two vectors are equal
 * @param v1 - first vector in the comparison
 * @param v2 - second vector in the comparison
 * @return boolean whether the vectors are equal
 */
bool operator==(const vector<double>& v1, const vector<double>& v2);
bool operator==(const vector<unsigned>& v1, const vector<unsigned>& v2);
bool operator==(const vector<double>& v1, const vector<unsigned>& v2);
bool operator==(const vector<unsigned>& v1, const vector<double>& v2);

/**
 * Addition operator between a constant and a vector.
 * @param alpha - constant to be added to every element in a vector
 * @param v - vector to which to add the constant
 * @return vector to which a constant has been added to every element
 */
vector<double> operator+(const double& alpha, const vector<double>& v);
vector<double> operator+(const vector<double>& v, const double& alpha);

/**
 * Addition between two vectors
 * @param v1 - first vector in the sum
 * @param v2 - second vector in the sum
 * @return sum of two vectors
 */
vector<double> operator+(const vector<double>& v1, const vector<double>& v2);
vector<unsigned> operator+(const vector<unsigned>& v1, const vector<unsigned>& v2);
vector<double> operator+(const vector<double>& v1, const vector<unsigned>& v2);
vector<double> operator+(const vector<unsigned>& v1, const vector<double>& v2);

/**
 * Difference between two vectors
 * @param v1 - the first vector in the difference
 * @param v2 - the second vector in the difference
 * @return difference between vectors
 */
vector<double> operator-(const vector<double>& v1, const vector<double>& v2);
vector<unsigned> operator-(const vector<unsigned>& v1, const vector<unsigned>& v2);

/**
 * Multiplication between vectors
 * @param alpha - constant by which to multiply the vector
 * @param v - vector that is going to be multiplied by a constant
 * @return vector
 */
vector<double> operator*(const double& alpha, const vector<double>& v);
vector<unsigned> operator*(const unsigned& alpha, const vector<unsigned>& v);

/**
 * Multiplication between a constant and a matrix
 * @param alpha - constant by which to multiply the matrix
 * @param m - matrix that will be multiplied by a constant
 * @return matrix
 */
vector<vector<unsigned>> operator*(const unsigned& alpha, const vector<vector<unsigned>>& m);

/**
 * Element-by-element multiplication between two vectors
 * @param v1 - first vector
 * @param v2 - second vector
 * @return vector
 */
vector<double> operator*(const vector<double>& v1, const vector<double>& v2);

/**
 * Division of a vector by a constant
 * @param v - vector to be divided by a constant
 * @param alpha - constant by which to divide the vector
 * @return vector
 */
vector<double> operator/(const vector<double>& v, const double& alpha);

/**
 * Floor division of a vector by a constant
 * @param v - vector to be divided
 * @param alpha - constant by which to divide the vector
 * @return vector
 */
vector<unsigned> floor_div(const vector<unsigned>& v, const unsigned& alpha);

#endif //STOSPA_VECTORFUNCTIONS_HPP
