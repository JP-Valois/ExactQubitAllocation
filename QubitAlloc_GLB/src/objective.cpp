#include "../include/objective.hpp"


int Objective(const vector<int>& mapping, const vector<vector<int>>& D, const vector<vector<int>>& F, int n)
{
    int cost = 0;
    for (int i = 0; i < n; ++i)
    {
        if (mapping[i] == -1)
            continue;

        for (int j = 0; j < n; ++j)
        {
            if (mapping[j] == -1)
                continue;
                
            cost += F[i][j] * D[mapping[i]][mapping[j]];
        }
    }
    return cost;
}


bool ValidSolutions (const vector<vector<int>>& op_sols, const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int min_cost)
{
    for (vector<int> map : op_sols)
    {
        if (Objective(map, D, F, n) != min_cost)
            return false;
    }

    return true;
}
