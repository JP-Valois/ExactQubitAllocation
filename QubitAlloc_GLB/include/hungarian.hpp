#ifndef __HUNGARIAN__
    #define __HUNGARIAN__ 

    #include "utils.hpp"

    int Hungarian (const vector<int>& L, int n, int m);
    int HungarianWithDuals(const std::vector<int>& L, int n, int m, std::vector<int>& rowDual, std::vector<int>& colDual);

#endif
