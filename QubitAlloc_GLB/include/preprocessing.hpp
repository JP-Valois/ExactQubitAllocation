#ifndef __PREPROCESSING__
    #define __PREPROCESSING__ 

    #include "objective.hpp"
    #include "utils.hpp"

    int Load_distanceMatrix (vector<vector<int>>& D, const std::string& file_path);

    int Load_interactionMatrix (vector<vector<int>>& M, const std::string& file_path, int m);

    vector<vector<int>> LinearCouplingDistanceMatrix (int n);

    vector<int> RowwiseSum (const vector<vector<int>>& F, const int m);
    
    vector<int> Prioritization (const vector<vector<int>>& F, const int n, const int m);

    int GreedyAllocation (const vector<vector<int>>& D, const vector<vector<int>>& F, const vector<int>& sort_arr, const int n, const int m, vector<int>& alloc);

#endif
