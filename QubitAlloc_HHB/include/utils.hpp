#ifndef __UTILS__
    #define __UTILS__ 

    #include <iostream>
    #include <cassert>
    #include <limits>
    #include <vector>
    #include <queue>
    #include <stack>
    #include <algorithm>
    #include <chrono>
    #include <fstream>
    #include <string>
    #include <sstream>
    #include <omp.h>

    using namespace std;

    const int INF = std::numeric_limits<int>::max();

    /* overload operator << to display a 1D vector */
    template <typename T>
    std::ostream& operator<< (std::ostream& os, const vector<T>& vec)
    {
        os << "[ ";
        for (size_t i = 0; i < vec.size(); ++i)
        {
            os << vec[i] << " ";
        }
        os << "]";
        return os;
    }

    /* overload operator << to display a 2D matrix */
    template <typename T>
    std::ostream& operator<< (std::ostream& os, const vector<vector<T>>& matrix)
    {
        for (const auto &row : matrix)
            os << row << "\n";

        return os;
    }

    double ComputeTotalNumberOfNodes (int n, int m);

    vector<int> Arange (int n);

    int idx4D (int i, int j, int k, int l, int n);

#endif
