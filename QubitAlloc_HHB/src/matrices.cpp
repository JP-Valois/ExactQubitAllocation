#include "../include/matrices.hpp"


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


CostMatrix CostMatrix::Assemble (const vector<vector<int>>& D, const vector<vector<int>>& F, int m)
{
    vector<int> costs(m*m*m*m, 0);
    vector<int> leader(m*m, 0);

    const int inf = INF / 2;

    int i, j, k, l;

    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            for (k = 0; k < m; ++k)
            {
                for (l = 0; l < m; ++l)
                {
                    if ((k == i) ^ (l == j))
                        costs[idx4D(i, j, k, l, m)] = inf;
                    else
                        costs[idx4D(i, j, k, l, m)] = F[i][k] * D[j][l];

                    //if (!((k == i) != (l == j)))
                    //    costs[idx4D(i, j, k, l, m)] = F[i][k] * D[j][l];
                }
            }
            leader[i*m + j] = costs[idx4D(i, j, i, j, m)];
        }
    }

    return CostMatrix{costs, leader, m};
}


void CostMatrix::distributeLeader ()
{
    vector<int>& C = this -> costs;
    vector<int>& L = this -> leader;
    const int n = this -> size;

    int i, j, k, l;
    int leader_cost, leader_cost_div, leader_cost_rem, val;

    if (n == 1)
    {
        C[0] = 0;
        L[0] = 0;

        return;
    }

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            leader_cost = L[i*n + j];

            C[idx4D(i, j, i, j, n)] = 0;
            L[i*n + j] = 0;

            if (leader_cost == 0)
            {
                continue;
            }

            leader_cost_div = leader_cost / (n - 1);
            leader_cost_rem = leader_cost % (n - 1);

            for (k = 0; k < n; ++k)
            {
                if (k == i)
                    continue;

                val = leader_cost_div + (k < leader_cost_rem || (k == leader_cost_rem && i < k));
                
                for (l = 0; l < n; ++l)
                {
                    if (l != j)
                        C[idx4D(i, j, k, l, n)] += val;
                }
            }
        }
    }
}


void CostMatrix::distributeLeader_OMP ()
{
    vector<int>& C = this -> costs;
    vector<int>& L = this -> leader;
    const int n = this -> size;

    int i, j, k, l;
    int leader_cost, leader_cost_div, leader_cost_rem, val;

    if (n == 1)
    {
        C[0] = 0;
        L[0] = 0;

        return;
    }

    const int nb_threads = std::min(n*n, omp_get_max_threads());

    #pragma omp parallel for collapse(2) default(none) private(i, j, k, l, leader_cost, leader_cost_div, leader_cost_rem, val) \
            shared(C, L, n) num_threads(nb_threads)
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            leader_cost = L[i*n + j];

            C[idx4D(i, j, i, j, n)] = 0;
            L[i*n + j] = 0;

            if (leader_cost == 0)
            {
                continue;
            }

            leader_cost_div = leader_cost / (n - 1);
            leader_cost_rem = leader_cost % (n - 1);

            for (k = 0; k < n; ++k)
            {
                if (k == i)
                    continue;

                val = leader_cost_div + (k < leader_cost_rem || (k == leader_cost_rem && i < k));
                
                for (l = 0; l < n; ++l)
                {
                    if (l != j)
                        C[idx4D(i, j, k, l, n)] += val;
                }
            }
        }
    }
}


void CostMatrix::halveComplementary ()
{
    vector<int>& C = this -> costs;
    const int n = this -> size;

    int i, j, k, l;
    int cost_sum;

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = i; k < n; ++k)
            {
                for (l = 0; l < n; ++l)
                {
                    if ((k != i) && (l != j))
                    {
                        cost_sum = C[idx4D(i, j, k, l, n)] + C[idx4D(k, l, i, j, n)];
                        C[idx4D(i, j, k, l, n)] = cost_sum / 2;
                        C[idx4D(k, l, i, j, n)] = cost_sum / 2;

                        if (cost_sum % 2 == 1)
                        {
                            if ((i + j + k + l) % 2 == 0) // total index parity for balance
                                C[idx4D(i, j, k, l, n)] += 1;
                            else
                                C[idx4D(k, l, i, j, n)] += 1;
                        }
                    }
                }
            }
        }
    }
}


void CostMatrix::halveComplementary_OMP ()
{
    vector<int>& C = this -> costs;
    const int n = this -> size;

    int i, j, k, l;
    int cost_sum;

    const int nb_threads = std::min(n*n, omp_get_max_threads());

    #pragma omp parallel for collapse(2) default(none) private(i, j, k, l, cost_sum) shared(C, n) num_threads(nb_threads)
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = i; k < n; ++k)
            {
                for (l = 0; l < n; ++l)
                {
                    if ((k != i) && (l != j))
                    {
                        cost_sum = C[idx4D(i, j, k, l, n)] + C[idx4D(k, l, i, j, n)];
                        C[idx4D(i, j, k, l, n)] = cost_sum / 2;
                        C[idx4D(k, l, i, j, n)] = cost_sum / 2;

                        if (cost_sum % 2 == 1)
                        {
                            if ((i + j + k + l) % 2 == 0) // total index parity for balance
                                C[idx4D(i, j, k, l, n)] += 1;
                            else
                                C[idx4D(k, l, i, j, n)] += 1;
                        }
                    }
                }
            }
        }
    }
}


CostMatrix CostMatrix::reduce (int i, int j) const
{
    const int n = this -> size;
    const int m = n - 1;

    assert(n > 0 && "Cannot reduce problem of size 0.");
    assert(std::min(i, j) >= 0 && std::max(i, j) < n && "Invalid reduction indices.");

    const vector<int>& C = this -> costs;
    vector<int> L_copy = this -> leader;
    
    vector<int> C_new;
    vector<int> L_new;

    int x, y, k, l;
    int x2, y2, k2, l2;

    // Updating the leader
    for (x = 0; x < n; ++x)
    {
        if (x == i)
            continue;

        for (y = 0; y < n; ++y)
        {
            if (y != j)
            {
                L_copy[x*n + y] += (C[idx4D(x, y, i, j, n)] + C[idx4D(i, j, x, y, n)]);
            }
        }
    }

    // clearing + reallocating
    C_new.assign(m*m*m*m, 0);
    L_new.assign(m*m, 0);

    // reducing the matrix
    x2 = 0;
    for (x = 0; x < n; ++x)
    {
        if (x == i)
            continue;

        y2 = 0;
        for (y = 0; y < n; ++y)
        {
            if (y == j)
                continue;

            // copy C_xy into C_x2y2
            k2 = 0;
            for (k = 0; k < n; ++k)
            {
                if (k == i)
                    continue;

                l2 = 0;
                for (l = 0; l < n; ++l)
                {
                    if (l == j)
                        continue;

                    C_new[idx4D(x2, y2, k2, l2, m)] = C[idx4D(x, y, k, l, n)];
                    ++l2;
                }
                ++k2;
            }

            L_new[x2*m + y2] = L_copy[x*n + y];
            ++y2;
        }
        ++x2;
    }

    return CostMatrix{C_new, L_new, m};
}
