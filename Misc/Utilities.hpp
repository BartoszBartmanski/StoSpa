#ifndef STOSPA_UTILITIES_HPP
#define STOSPA_UTILITIES_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <unistd.h>
#include <dirent.h>
#include "Version.hpp"

using namespace std;


/**
 * Proper modulus operator
 * @param a - number on which to apply the modulus operation
 * @param b - base number
 * @return an integer
 */
int mod (const int& a, const int& b);

/**
 * Calculate the greatest common divisor of two numbers
 * @param a - first number
 * @param b - second number
 * @return greatest common divisor of two numbers
 */
int gcd(int a, int b);

/**
 * Pads the given string with the given length of given character.
 * @param s - string to be padded
 * @param pad - character to pad the string with
 * @param length - total length of the string after padding
 * @return changed string
 */
string pad_string(string s, char pad='0', unsigned length=2);

/**
 * Get a substring from a given string from the string start to string end.
 * @param line - string from which to get the substring
 * @param start - string from which to start to find the substring
 * @param end - string to which stop searching for the substring
 * @return substring
 */
string get_substring(string line, string start="", string end="\n");

/**
 * Converts a character array to a string
 * @param num - number of pointers ot character arrays
 * @param array - an array of pointers to character arrays
 * @return string
 */
string arr_to_str(int num, const char** array);

/**
 * Splits a string based on the spacing character
 * @param input - the string to be split
 * @param separator - character that is used as a separation
 * @return vector of strings that results from the separating the original string
 */

template<typename T=string>
vector<T> split(const string &input, char separator=',');

/**
 * Checks whether the given path is valid.
 * @param path - string that is to be checked
 * @return boolean variable indicating whether the path is valid.
 */
bool check_dir(string path);

/**
 * Returns a list of files or directories for the given path
 * @param path - given path to the directory
 * @param pattern - only get the file names that contain the string pattern
 * @param cutoff - the number of characters to keep when adding file name to the return vector
 * @return vector<string>
 */
vector<string> get_list_of_files(const string& path, const string& pattern="", int cutoff=-1);
/**
 * Checks whether the directory exists and updates the file_name accordingly
 * @param dir_name - name of directory
 * @param file_name - name of file
 * @param start_index - the extension of the file to be added to file name
 */
string update_path(string dir_name, string file_name, unsigned start_index);

class Progress
{
private:
    /** Total number of steps that need to be taken to finish the simulation. */
    unsigned mNumSteps;

    unsigned mCurrentStep;

public:
    /**
     * Constructor.
     * @param num_steps - total number of steps
     */
    explicit Progress(unsigned num_steps);

    /**
     * Default destructor.
     */
    ~Progress()= default;

    /**
     * Displays the progress of the simulation.
     * @param step - the index of the current step.
     */
    void Show();

    /**
     * Resets the counters
     * @return
     */
    void Reset();

};

#endif //STOCHASTIC_DIFFUSION_UTILITIES_HPP
