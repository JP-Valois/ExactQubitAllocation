#include "../include/objective.hpp"


int ObjectiveFunction (const vector<int>& mapping, const vector<vector<int>>& D, const vector<vector<int>>& F, int n)
{
    int route_cost = 0;

    int i, j;

    for (i = 0; i < n; ++i)
    {
        for (j = i; j < n; ++j)
        {
            route_cost += F[i][j] * D[mapping[i]][mapping[j]];
        }
    }      

    return 2*route_cost;
}


bool ValidSolutions (const vector<vector<int>>& op_sols, const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int min_cost)
{
    for (vector<int> map : op_sols)
    {
        if (ObjectiveFunction(map, D, F, n) != min_cost)
            return false;
    }

    return true;
}
