#ifndef __OBJECTIVE__
    #define __OBJECTIVE__ 

    #include "utils.hpp"

    int Objective (const vector<int>& map, const vector<vector<int>>& D, const vector<vector<int>>& F, int n);

    bool ValidSolutions (const vector<vector<int>>& op_sols, const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int min_cost);

#endif
