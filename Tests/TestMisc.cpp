//
// Created by bartosz on 25/12/17.
//

#include "catch.hpp"
#include "VectorFunctions.hpp"
#include "Grid.hpp"
#include "Utilities.hpp"

TEST_CASE("Test Grids.*pp")
{
    Grid grid = Grid(2, 0.1, 5);

    SECTION("Checkd voxels")
    {
        for (const vector<unsigned>& species : grid.voxels)
        {
            for (auto voxel : species)
            {
                REQUIRE(voxel == 0);
            }
        }
    }

    SECTION("Check time increments")
    {
        double inf = numeric_limits<double>::infinity();
        for (const double& voxel : grid.next_reaction_time)
        {
            REQUIRE(voxel == inf);
        }
    }

}

TEST_CASE("Test Useful.*pp")
{
    SECTION("Check pad_string function")
    {
        string test = "1";
        test = pad_string(test, '0', 3);
        REQUIRE(test == "001");

        test = "002";
        test = pad_string(test, '0', 2);
        REQUIRE(test == "002");
    }

    SECTION("Check get_substring function")
    {
        string test = "This is a test case\n";
        test = get_substring(test, "test ", "\n");
        REQUIRE(test == "case\n");
    }

    SECTION("Check check_dir function")
    {
        REQUIRE(check_dir("/tmp/"));
        REQUIRE(!check_dir("/tm/"));
    }
    SECTION("Check get_list_of_files function")
    {
        vector<string> files = get_list_of_files("/home/");
        REQUIRE(files.size() >= 3);
        REQUIRE(files[0] == ".");
        REQUIRE(files[1] == "..");
    }

    SECTION("Check get_dir function")
    {
        string test = get_dir();
        REQUIRE(check_dir(test));
    }

    SECTION("Check split function")
    {
        string a_string = "1000.0";
        vector<string> test = split(a_string, '.');

        REQUIRE(test[0] == "1000");
        REQUIRE(test[1] == "0");
    }

    SECTION("Check update_path function")
    {
        string tmp_name = "/tmp/test_file_000.dat";
        ofstream test(tmp_name);
        test << "This is a test file" << endl;
        test.close();
        string dir_name = "/tmp/";
        string file_name = "test_file";
        string path = update_path(dir_name, file_name, 0);
        REQUIRE(dir_name == "/tmp/");
        REQUIRE(file_name == "test_file");
        REQUIRE(path == "/tmp/test_file_001.dat");
        remove(tmp_name.c_str());
    }
}

TEST_CASE("Test VectorFunctions.*pp")
{
    string ans;

    SECTION("Check save_vector functions")
    {
        vector<double> test = {0, 1, 2, 3};
        save_vector(test, "./test.dat");

        ifstream test_file ("./test.dat");
        REQUIRE(test_file.is_open());
        if (test_file.is_open())
        {
            string line;
            getline(test_file, line);
            vector<string> vec = split(line, ' ');

            REQUIRE(vec.size() == test.size());
            for ( unsigned i=0; i < test.size(); i++ )
            {
                REQUIRE(stod(vec[i]) == test[i]);
            }

        }
        test_file.close();
        remove("test.dat");
    }

    SECTION("Check linspace function")
    {
        vector<double> test = linspace(0.0, 1.0, 10);
        for (unsigned i=0; i < 10; i++)
        {
            REQUIRE((test[i] - (0.0 + i * (1.0 - 0.0)/9)) < 0.00001);
        }
    }

    SECTION("Check ones function")
    {
        vector<unsigned> test_1d = ones(10);
        REQUIRE(test_1d.size() == 10);
        for (unsigned elem : test_1d)
        {
            REQUIRE(elem == 1);
        }

        vector<vector<unsigned>> test_2d = ones(10, 2);
        REQUIRE(test_2d.size() == 2);
        for (vector<unsigned> row : test_2d)
        {
            REQUIRE(row.size() == 10);
            for (unsigned elem : row)
            {
                REQUIRE(elem == 1);
            }
        }

    }

    SECTION("Check vec_sum function")
    {
        vector<double> test = {0, 1, 2, 3};
        double result = vec_sum(test);
        REQUIRE(result == 6.0);
    }

    SECTION("Check add function")
    {
        vector<double> test_double = {0.0, 1.0, 2.0, 3.0};
        vector<double> result_double = add(test_double, test_double);
        for (unsigned i=0; i < test_double.size(); i++)
        {
            REQUIRE(result_double[i] == 2.0 * test_double[i]);
        }

        vector<int> test_int = {1, 2, 3, 4};
        vector<int> result_int = add(test_int, test_int);
        for (unsigned i=0; i < test_int.size(); i++)
        {
            REQUIRE(result_int[i] == 2 * test_int[i]);
        }

        vector<unsigned> test_unsigned = {2, 3, 4, 5};
        vector<unsigned> result_unsigned = add(test_unsigned, test_unsigned);
        for (unsigned i=0; i < test_double.size(); i++)
        {
            REQUIRE(result_unsigned[i] == 2 * test_unsigned[i]);
        }

    }

    SECTION("Check find_min and find_min_index functions")
    {
        vector<double> test_double = {3, 0, 2, 1};
        double test_min_double = find_min(test_double);
        REQUIRE(test_min_double == 0.0);
        unsigned test_min_index_double = find_min_index(test_double);
        REQUIRE(test_min_index_double == 1);

        vector<unsigned> test_unsigned = {3, 0, 2, 1};
        unsigned test_min_unsigned = find_min(test_unsigned);
        REQUIRE(test_min_unsigned == 0);
        unsigned test_min_index_unsigned = find_min_index(test_unsigned);
        REQUIRE(test_min_index_unsigned == 1);

        vector<vector<double>> test_mat_double = {{4, 1}, {2, 6}};
        test_min_double = find_min(test_mat_double);
        REQUIRE(test_min_double == 1.0);
        vector<unsigned> test_min_index_mat_double = find_min_index(test_mat_double);
        REQUIRE(test_min_index_mat_double[0] == 0);
        REQUIRE(test_min_index_mat_double[1] == 1);

        vector<vector<unsigned>> test_mat_unsigned = {{4, 1}, {2, 6}};
        test_min_unsigned = find_min(test_mat_unsigned);
        REQUIRE(test_min_unsigned == 1);
        vector<unsigned> test_min_index_mat_unsigned = find_min_index(test_mat_unsigned);
        REQUIRE(test_min_index_mat_unsigned[0] == 0);
        REQUIRE(test_min_index_mat_unsigned[1] == 1);
    }

    SECTION("Check opertors")
    {
        vector<double> result_double;
        vector<unsigned> result_unsigned;

        vector<double> test_double_1 = {0, 1, 2};
        vector<double> test_double_2 = {5, 2, 3};
        vector<double> test_double_3 = {1, 2, 3};
        vector<double> test_double_4 = {5, 1, 1};
        vector<double> test_double_5 = {0, 2, 4};
        vector<double> test_double_6 = {0, 2, 6};

        vector<unsigned> test_unsigned_1 = {0, 1, 2};
        vector<unsigned> test_unsigned_2 = {5, 2, 3};
        vector<unsigned> test_unsigned_3 = {1, 2, 3};
        vector<unsigned> test_unsigned_4 = {5, 1, 1};
        vector<unsigned> test_unsigned_5 = {0, 2, 4};
        vector<unsigned> test_unsigned_6 = {0, 2, 6};

        REQUIRE(test_double_1 == test_double_1);
        REQUIRE(!(test_double_1 == test_double_2));
        REQUIRE(test_double_1 == test_unsigned_1);
        REQUIRE(test_unsigned_1 == test_double_1);
        REQUIRE(test_unsigned_1 == test_unsigned_1);

        result_double = 1.0 + test_double_1;
        REQUIRE(result_double == test_double_3);
        result_double = test_double_1 + 1.0;
        REQUIRE(result_double == test_double_3);
        
        result_double = test_double_1 + test_double_4;
        REQUIRE(result_double == test_double_2);
        result_unsigned = test_unsigned_1 + test_unsigned_4;
        REQUIRE(result_unsigned == test_unsigned_2);
        result_double = test_double_1 + test_unsigned_4;
        REQUIRE(result_double == test_double_2);
        result_double = test_unsigned_4 + test_double_1;
        REQUIRE(result_double == test_double_2);

        result_double = test_double_2 - test_double_4;
        REQUIRE(result_double == test_double_1);
        result_unsigned = test_unsigned_2 - test_unsigned_4;
        REQUIRE(result_unsigned == test_unsigned_1);

        result_double = 2 * test_double_1;
        REQUIRE(result_double == test_double_5);
        result_unsigned = 2 * test_unsigned_1;
        REQUIRE(result_unsigned == test_unsigned_5);

        vector<vector<unsigned>> matrix_1 = {{0, 1}, {2, 3}};
        vector<vector<unsigned>> matrix_2 = {{0, 2}, {4, 6}};
        vector<vector<unsigned>> matrix_result;
        matrix_result = 2 * matrix_1;
        REQUIRE(matrix_result[0][0] == matrix_2[0][0]);
        REQUIRE(matrix_result[0][1] == matrix_2[0][1]);
        REQUIRE(matrix_result[1][0] == matrix_2[1][0]);
        REQUIRE(matrix_result[1][1] == matrix_2[1][1]);

        result_double = test_double_1 * test_double_3;
        REQUIRE(result_double == test_double_6);
        result_double = test_double_5 / 2;
        REQUIRE(result_double == test_double_1);

        REQUIRE(mod(5, 3) == 2);
        REQUIRE(mod(-2, 5) == 3);

    }
}