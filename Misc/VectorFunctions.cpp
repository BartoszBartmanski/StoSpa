//
// Created by bartmanski on 07/04/17.
//

#include "VectorFunctions.hpp"

void print_array(vector<double> vec, string name)
{
    cout << endl << name << " = ";

    for (double elem : vec)
    {
        cout << elem << " ";
    }
    cout << endl;
}

void print_array(vector<unsigned> vec, string name)
{
    cout << endl << name << " = ";

    for (unsigned elem : vec)
    {
        cout << elem << " ";
    }
    cout << endl;
}
void print_array(vector<int> vec, string name)
{
    cout << endl << name << " = ";

    for (int elem : vec)
    {
        cout << elem << " ";
    }
    cout << endl;
}

void print_array(vector<vector<double>> matrix, string name)
{
    cout << endl << name << " = " << endl;
    for (vector<double> row : matrix)
    {
        for (double elem : row)
        {
            cout << elem << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_array(vector<vector<unsigned>> matrix, string name)
{
    cout << endl << name << " = " << endl;
    for (vector<unsigned> row : matrix)
    {
        for (unsigned elem : row)
        {
            cout << elem << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_array(vector<vector<int>> matrix, string name)
{
    cout << endl << name << " = " << endl;
    for (vector<int> row : matrix)
    {
        for (int elem : row)
        {
            cout << elem << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void save_vector(vector<double> vec, string path_to_file)
{
    ofstream vec_file(path_to_file, ios_base::app);
    ostream_iterator<double> output_iterator(vec_file, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    vec_file << endl;
    vec_file.close();
}

void save_vector(vector<long double> vec, string path_to_file)
{
    ofstream vec_file(path_to_file, ios_base::app);
    ostream_iterator<long double> output_iterator(vec_file, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    vec_file << endl;
    vec_file.close();
}

void save_vector(vector<unsigned> vec, string path_to_file)
{
    ofstream vec_file(path_to_file, ios_base::app);
    ostream_iterator<unsigned> output_iterator(vec_file, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    vec_file << endl;
    vec_file.close();
}

void save_vector(vector<double> vec, const unique_ptr<ofstream>& handle)
{
    ostream_iterator<unsigned> output_iterator(*handle, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    *handle << endl;
}

void save_vector(vector<long double> vec, const unique_ptr<ofstream>& handle)
{
    ostream_iterator<unsigned> output_iterator(*handle, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    *handle << endl;
}

void save_vector(vector<unsigned> vec, const unique_ptr<ofstream>& handle)
{
    ostream_iterator<unsigned> output_iterator(*handle, " ");
    copy(vec.begin(), vec.end(), output_iterator);
    *handle << endl;
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

double vec_sum(vector<double> vec)
{
    double sum = 0.0;
    for (double elem : vec)
    {
        sum += elem;
    }

    return sum;
}

unsigned vec_sum(vector<unsigned> vec)
{
    unsigned sum = 0;
    for (unsigned elem : vec)
    {
        sum += elem;
    }

    return sum;
}

vector<double> add(const vector<double>& vec1, const vector<double>& vec2)
{
    // Check the input is sensible
    assert(vec1.size() == vec2.size());

    vector<double> output(vec1.size());
    for (unsigned i=0; i < vec1.size(); i++)
    {
        output[i] = vec1[i] + vec2[i];
    }
    return output;
}

vector<unsigned> add(const vector<unsigned>& vec1, const vector<unsigned>& vec2)
{
    // Check the input is sensible
    assert(vec1.size() == vec2.size());

    vector<unsigned> output(vec1.size());
    for (unsigned i=0; i < vec1.size(); i++)
    {
        output[i] = vec1[i] + vec2[i];
    }
    return output;
}

vector<int> add(const vector<int>& vec1, const vector<int>& vec2)
{
    // Check the input is sensible
    assert(vec1.size() == vec2.size());

    vector<int> output(vec1.size());
    for (unsigned i=0; i < vec1.size(); i++)
    {
        output[i] = vec1[i] + vec2[i];
    }
    return output;
}

double find_min(vector<double> vec)
{
    auto result = min_element(vec.begin(), vec.end());

    return *result;
}

unsigned find_min(vector<unsigned> vec)
{
    auto result = min_element(vec.begin(), vec.end());

    return *result;
}

double find_min(vector<vector<double>> matrix)
{
    vector<unsigned> index_of_min(2);
    vector<double> min_elements_for_each_row(matrix.size());
    vector<vector<double>::iterator> iterators_for_each_species(matrix.size());

    for (unsigned row=0; row < matrix.size(); row++)
    {
        auto result = min_element(matrix[row].begin(), matrix[row].end());
        min_elements_for_each_row[row] = *result;
        iterators_for_each_species[row] = result;
    }
    auto result = min_element(min_elements_for_each_row.begin(),
                                                  min_elements_for_each_row.end());
    return *result;
}

unsigned find_min(vector<vector<unsigned>> matrix)
{
    vector<unsigned> index_of_min(2);
    vector<unsigned> min_elements_for_each_row(matrix.size());
    vector<vector<unsigned>::iterator> iterators_for_each_species(matrix.size());

    for (unsigned row=0; row < matrix.size(); row++)
    {
        auto result = min_element(matrix[row].begin(), matrix[row].end());
        min_elements_for_each_row[row] = *result;
        iterators_for_each_species[row] = result;
    }
    auto result = min_element(min_elements_for_each_row.begin(), min_elements_for_each_row.end());
    return *result;
}

unsigned find_min_index(vector<double> vec)
{
    auto result = min_element(vec.begin(), vec.end());

    unsigned index = unsigned(distance(vec.begin(), result));

    return index;
}

unsigned find_min_index(vector<unsigned> vec)
{
    auto result = min_element(vec.begin(), vec.end());

    unsigned index = unsigned(distance(vec.begin(), result));

    return index;
}

vector<unsigned> find_min_index(vector<vector<double>> matrix)
{
    vector<unsigned> index_of_min(2);
    vector<double> min_elements_for_each_row(matrix.size());
    vector<vector<double>::iterator> iterators_for_each_species(matrix.size());

    for (unsigned row=0; row < matrix.size(); row++)
    {
        auto result = min_element(matrix[row].begin(), matrix[row].end());
        min_elements_for_each_row[row] = *result;
        iterators_for_each_species[row] = result;
    }
    auto result = min_element(min_elements_for_each_row.begin(),
                                                  min_elements_for_each_row.end());

    index_of_min[0] = unsigned(distance(min_elements_for_each_row.begin(), result));
    index_of_min[1] = unsigned(distance(matrix[index_of_min[0]].begin(), iterators_for_each_species[index_of_min[0]]));

    return index_of_min;
}

vector<unsigned> find_min_index(vector<vector<unsigned>> matrix)
{
    vector<unsigned> index_of_min(2);
    vector<unsigned> min_elements_for_each_row(matrix.size());
    vector<vector<unsigned>::iterator> iterators_for_each_species(matrix.size());

    for (unsigned row=0; row < matrix.size(); row++)
    {
        auto result = min_element(matrix[row].begin(), matrix[row].end());
        min_elements_for_each_row[row] = *result;
        iterators_for_each_species[row] = result;
    }
    auto result = min_element(min_elements_for_each_row.begin(), min_elements_for_each_row.end());

    index_of_min[0] = unsigned(distance(min_elements_for_each_row.begin(), result));
    index_of_min[1] = unsigned(distance(matrix[index_of_min[0]].begin(), iterators_for_each_species[index_of_min[0]]));

    return index_of_min;
}

bool operator==(const vector<double>& v1, const vector<double>& v2)
{
    assert(v1.size() == v2.size());
    for (unsigned i=0; i < v1.size(); i++)
    {
        if (v1[i] != v2[i])
            return false;

    }
    return true;
}

bool operator==(const vector<unsigned>& v1, const vector<unsigned>& v2)
{
    assert(v1.size() == v2.size());
    for (unsigned i=0; i < v1.size(); i++)
    {
        if (v1[i] != v2[i])
            return false;

    }
    return true;
}

bool operator==(const vector<double>& v1, const vector<unsigned>& v2)
{
    assert(v1.size() == v2.size());
    for (unsigned i=0; i < v1.size(); i++)
    {
        if (v1[i] != v2[i])
            return false;

    }
    return true;
}

bool operator==(const vector<unsigned>& v1, const vector<double>& v2)
{
    assert(v1.size() == v2.size());
    for (unsigned i=0; i < v1.size(); i++)
    {
        if (v1[i] != v2[i])
            return false;

    }
    return true;
}

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

vector<double> operator-(const vector<double>& v1, const vector<double>& v2)
{
    assert(v1.size() == v2.size());
    vector<double> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] - v2[i];
    }

    return output;
}

vector<unsigned> operator-(const vector<unsigned>& v1, const vector<unsigned>& v2)
{
    assert(v1.size() == v2.size());
    vector<unsigned> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] - v2[i];
    }

    return output;
}

vector<double> operator* (const double& alpha, const vector<double>& v)
{
    vector<double> output(v.size());
    for (unsigned i=0; i<v.size(); i++)
    {
        output[i] = alpha * v[i];
    }
    return output;
}

vector<unsigned> operator* (const unsigned& alpha, const vector<unsigned>& v)
{
    vector<unsigned> output(v.size());
    for (unsigned i=0; i<v.size(); i++)
    {
        output[i] = alpha * v[i];
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

vector<double> operator* (const vector<double>& v1, const vector<double>& v2)
{
    assert(v1.size() == v2.size());
    vector<double> output(v1.size());
    for (unsigned i=0; i<v1.size(); i++)
    {
        output[i] = v1[i] * v2[i];
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
