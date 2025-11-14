#include "../include/preprocessing.hpp"


int Load_distanceMatrix (vector<vector<int>>& D, const std::string& file_path)
{
    std::ifstream file;
    file.open(file_path, std::ios::in);

    if (!file.is_open())
    {
        cerr << "Error opening file: " << file_path << endl;
        exit(1);
    }

    std::string line;
    int m = 0;

    // read the first line to get the input matrix size (physical qubits)
    if (getline(file, line))
    {
        istringstream iss(line);
        iss >> m;
    }

    D.assign(m, vector<int>(m, 0));

    // skip the second line
    getline(file, line);
    
    // read data
    int row{0}, col{0};
    while (getline(file, line))
    {
        std::istringstream iss(line);

        int value;

        while (iss >> value && col < m)
        {
            D[row][col] = value;
            ++col;
        }
        ++row;

        col = 0;
    }

    file.close();

    return m;
}


int Load_interactionMatrix (vector<vector<int>>& F, const std::string& file_path, int m)
{
    F.assign(m, vector<int>(m, 0));

    std::ifstream file;
    file.open(file_path, std::ios::in);

    if (!file.is_open())
    {
        cerr << "Error opening file: " << file_path << endl;
        exit(1);
    }

    std::string line;
    int n = 0;

    // read the first line to get the input matrix size (logical qubits)
    if (getline(file, line))
    {
        istringstream iss(line);
        iss >> n;
    }

    // assertion
    assert(n <= m && "Error: There must be at least as many physical qubits as logical qubits.");

    // skip the second line
    getline(file, line);
    
    // read data
    int row{0}, col{0};
    while (getline(file, line) && row < m)
    {
        std::istringstream iss(line);

        int value;

        while (iss >> value && col < m)
        {
            F[row][col] = value;
            ++col;
        }
        ++row;

        col = 0;
    }

    file.close();

    return n;
}


vector<vector<int>> LinearCouplingDistanceMatrix (int n)
{
    vector<vector<int>> D(n, vector<int>(n, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                D[i][j] = abs(i - j) - 1;
            }
        }
    }

    return D;
}


vector<int> RowwiseSum (const vector<vector<int>>& F, const int m)
{
    vector<int> sF(m, 0);

    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < m; ++j)
            sF[i] += F[i][j];

    return sF;
}


vector<int> Prioritization (const vector<vector<int>>& F, int n, int m)
{
    vector<int> priority;
    priority.reserve(m);

    vector<int> sF = RowwiseSum(F, m);

    int i, j, min_inter, min_inter_index;

    for (i = 0; i < m; ++i)
    {
        min_inter = sF[0];
        min_inter_index = 0;

        for (j = 1; j < m; ++j)
        {
            if (sF[j] < min_inter)
            {
                min_inter = sF[j];
                min_inter_index = j;
            }
        }

        priority.insert(priority.begin(), min_inter_index);
        //priority.push_back(min_inter_index)

        sF[min_inter_index] = INF;

        for (j = 0; j < m; ++j)
        {
            if (sF[j] != INF)
                sF[j] -= F[j][min_inter_index];
        }
    }

    //std::reverse(priority);

    for (i = 0; i < n; ++i)
        assert(priority[i] < n && "Error: Logical prioritization failure.");

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

        route_cost_temp = Objective(alloc_temp, D, F, n);

        if (route_cost_temp < route_cost)
        {
            alloc = alloc_temp;
            route_cost = route_cost_temp;
        }
    }

    return route_cost;
}
