//
// Created by bartosz on 13/07/17.
//

#include "Utilities.hpp"

int mod (const int& a, const int& b)
{
    if(b < 0) //you can check for b == 0 separately and do what you want
        return mod(a, -b);
    int ret = a % b;
    if(ret < 0)
        ret+=b;
    return ret;
}

int gcd(int a, int b)
{
    return b == 0 ? a : gcd(b, a % b);
}

string pad_string(string s, char pad, unsigned length)
{
    std::ostringstream o;
    o << setw(length) << setfill(pad) << s;
    string output = o.str();

    return output;
}

string get_substring(string line, string start, string end)
{
    return line.substr(line.find(start)+start.size(), line.find(end));
}

string arr_to_str(int num, const char** array)
{
    string output;
    for (int i=0; i < num; i++)
    {
        output += string(array[i]) + " ";
    }

    return output;
}

template<>
vector<string> split(const string &input, char separator)
{
    istringstream ss(input);
    string token;

    vector<string> output;
    while(getline(ss, token, separator))
    {
        output.push_back(token);
    }

    return output;
}

template<>
vector<double> split(const string &input, char separator)
{
    istringstream ss(input);
    string token;

    vector<double> output;
    while(getline(ss, token, separator))
    {
        output.push_back(stod(token));
    }

    return output;
}

template<>
vector<unsigned> split(const string &input, char separator)
{
    istringstream ss(input);
    string token;

    vector<unsigned> output;
    while(getline(ss, token, separator))
    {
        output.push_back(unsigned(stoi(token)));
    }

    return output;
}

bool check_dir(string path)
{
    DIR *dir = opendir (path.c_str());

    return dir != nullptr;
}

vector<string> get_list_of_files(const string& path, const string& pattern, int cutoff)
{
    DIR *dir;
    struct dirent *ent;
    vector<string> files;
    if ((dir = opendir(path.c_str())) != nullptr)
    {
        while ((ent = readdir(dir)) != nullptr)
        {
            string file = ent->d_name;
            if (file.find(pattern) == 0)
            {
                if (cutoff > -1) { file = file.substr(0, pattern.size()+cutoff); }
                files.emplace_back(file);
            }
        }
        closedir(dir);
    }

    sort(files.begin(), files.end());

    return files;
}

string update_path(string dir_name, string file_name, unsigned start_index)
{
    // Check that the path dir_name is valid
    if (!check_dir(dir_name)) {
        cout << "Invalid path specified!" << endl;
        while (!check_dir(dir_name)) {
            cout << "Specify the path of where to save the results of this simulation:" << endl;
            cin >> dir_name;
        }
    }

    vector<string> names = get_list_of_files(dir_name, file_name, 4);

    if (names.empty())
    {
        file_name += "_" + pad_string(to_string(start_index), '0', 3);
    }
    else
    {
        for (unsigned i = start_index; i < 1000; i++) {
            string tmp = file_name + "_" + pad_string(to_string(i), '0', 3);
            if (find(names.begin(), names.end(), tmp) == names.end()) {
                file_name = tmp;
                break;
            }
        }
    }

    string path;
    if (dir_name.back() != '/')
    {
         dir_name += "/";
    }
    path = dir_name + file_name + ".dat";

    return path;
}

Progress::Progress(unsigned num_steps)
{
    mNumSteps = num_steps;
    mCurrentStep = 0;
    cout << setprecision(2) << fixed << "\rProgress: " << 0.0 << "%" << flush;
}

void Progress::Show()
{
    mCurrentStep++;
    cout << setprecision(2) << fixed << "\rProgress: " << 100.0 * mCurrentStep / double(mNumSteps) << "%" << flush;
    if (mCurrentStep >= mNumSteps)
    {
        cout << endl;
    }
}

void Progress::Reset()
{
    mCurrentStep = 0;
}

void Progress::End(string path_to_file)
{
    cout << "Data saved in " << path_to_file << endl;
}
