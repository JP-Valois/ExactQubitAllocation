#include "../include/heuristics.hpp"


vector<int> RowwiseSum (const vector<vector<int>>& M, const int m)
{
    vector<int> sM(m, 0);

    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < m; ++j)
            sM[i] += M[i][j];

    return sM;
}


vector<int> RowwiseNumZeros (const vector<vector<int>>& D, const int m)
{
    vector<int> nzD(m, 0);

    int i, j;

    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            if (D[i][j] == 0)
                nzD[i] ++;
        }
    }

    return nzD;
}


/* rank logical qubits (facilities) based on their interaction count */
vector<int> Prioritization (const vector<vector<int>>& F, int n, int m)
{
    vector<int> priority;
    priority.reserve(n);

    vector<int> sF = RowwiseSum(F, n);

    int i, j, min_inter, min_inter_index;

    for (i = 0; i < n; ++i)
    {
        min_inter = sF[0];
        min_inter_index = 0;

        for (j = 1; j < n; ++j)
        {
            if (sF[j] < min_inter)
            {
                min_inter = sF[j];
                min_inter_index = j;
            }
        }

        priority.insert(priority.begin(), min_inter_index);

        sF[min_inter_index] = INF;

        for (j = 0; j < n; ++j)
        {
            if (sF[j] != INF)
                sF[j] -= F[j][min_inter_index];
        }
    }

    for (i = 0; i < n; ++i)
        assert(priority[i] < n && "Error: Logical prioritization failure.");

    return priority;
}

vector<int> Prioritization_loc_basic (int m)
{
    vector<int> priority(m, 0);
    for (int i = 0; i < m-1; ++i)
        priority[i] = m - 1 - i;

    return priority;
}

/* rank physical qubits (locations) based on their overall distances */
vector<int> Prioritization_loc_dist (const vector<vector<int>>& D, int m)
{
    vector<int> priority;
    priority.reserve(m);

    vector<int> sD = RowwiseSum(D, m);

    int i, j, min_dist, min_dist_index;

    for (i = 0; i < m; ++i)
    {
        min_dist = sD[0];
        min_dist_index = 0;

        for (j = 1; j < m; ++j)
        {
            if (sD[j] < min_dist)
            {
                min_dist = sD[j];
                min_dist_index = j;
            }
        }

        priority.insert(priority.begin(), min_dist_index);
        //priority.push_back(min_dist_index);

        sD[min_dist_index] = INF;
    }

    return priority;
}


/* rank physical qubits (locations) based on their connectivity degree */
vector<int> Prioritization_loc_connec (const vector<vector<int>>& D, int m)
{
    vector<int> priority;
    priority.reserve(m);

    vector<int> nzD = RowwiseNumZeros(D, m);

    int i, j, min_connec, min_connec_index;

    for (i = 0; i < m; ++i)
    {
        min_connec = nzD[0];
        min_connec_index = 0;

        for (j = 1; j < m; ++j)
        {
            if (nzD[j] < min_connec)
            {
                min_connec = nzD[j];
                min_connec_index = j;
            }
        }

        //priority.insert(priority.begin(), min_connec_index);
        priority.push_back(min_connec_index);

        nzD[min_connec_index] = INF;
    }

    return priority;
}


int GreedyAllocation (const vector<vector<int>>& D, const vector<vector<int>>& F, const vector<int>& priority, int n, int m, vector<int>& alloc)
{
    alloc.clear();
    int route_cost = INF;

    int i, j, k, l, p, q, l_min{0};
    int route_cost_temp, cost_incre, min_cost_incre;

    for (j = 0; j < m; ++j)
    {
        vector<int> alloc_temp(n, -1);
        vector<bool> available(m, true);

        alloc_temp[priority[0]] = j;
        available[j] = false;

        // for each logical qubit (after the first one)
        for (p = 1; p < n; ++p)
        {
            k = priority[p];

            min_cost_incre = INF;

            // find physical qubit with least increasing route cost
            for (l = 0; l < m; ++l)
            {
                if (available[l])
                {
                    cost_incre = 0;
                    for (q = 0; q < p; ++q)
                    {
                        i = priority[q];
                        cost_incre += F[i][k] * D[alloc_temp[i]][l];
                    }

                    if (cost_incre < min_cost_incre)
                    {
                        l_min = l;
                        min_cost_incre = cost_incre;
                    }
                }
            }

            alloc_temp[k] = l_min;
            available[l_min] = false;
        }

        route_cost_temp = ObjectiveFunction(alloc_temp, D, F, n);

        if (route_cost_temp < route_cost)
        {
            alloc = alloc_temp;
            route_cost = route_cost_temp;
        }
    }

    //std::cout << alloc << std::endl;
    //std::cout << ObjectiveFunction(alloc, D, F, n) << std::endl;

    return route_cost;
}
