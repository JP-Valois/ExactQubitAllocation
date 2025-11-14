#ifndef __HEURISTICS__
    #define __HEURISTICS__ 

    #include "objective.hpp"

    vector<int> RowwiseSum (const vector<vector<int>>& F, const int m);
    
    vector<int> Prioritization (const vector<vector<int>>& F, const int n, const int m);

    int GreedyAllocation (const vector<vector<int>>& D, const vector<vector<int>>& F, const vector<int>& sort_arr, const int n, const int m, vector<int>& alloc);

#endif
