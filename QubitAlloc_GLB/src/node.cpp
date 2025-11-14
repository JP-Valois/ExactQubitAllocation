#include "../include/node.hpp"
#include "../include/hungarian.hpp"


Node Node::Root (int n, int m)
{
    vector<int> mapping(n, -1);
    vector<bool> available(m, true);

    return Node{mapping, available, 0};
}


vector<Node> Node::decompose (const vector<int>& priority, int n, int m)
{
    vector<Node> children;

    vector<int> map = this->mapping;
    vector<bool>& av = this->available;
    const int dp = this->depth;

    // next logical qubit q_i to assign
    int i = priority[dp];

    // iterate over available physical qubits
    for (int j = m - 1; j >= 0; --j)
    {
        if (!av[j])
            continue; // skip if not available

        // assign q_i to P_j
        map[i] = j;
        av[j] = false;

        // insert in children vector
        children.push_back(Node{map, av, dp+1});

        // restore data
        map[i] = -1;
        av[j] = true;
    }

    return children;
}


// unoptimized
// partial_mapping[i] = -1 means facility i is unassigned (n = facilities, m = locations)
// complexity = O(u^2 r^2), if n ~ m, complexity = O(M^4) with M = m - <depth of the node>
/*
vector<int> Node::Assemble_LAP_v0 (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m)
{
    vector<int> partial_mapping = this->mapping;
    vector<bool> av = this->available;

    // Identify assigned and unassigned facilities/locations
    vector<int> assigned_fac, unassigned_fac;
    vector<int> assigned_loc, unassigned_loc;

    for (int i = 0; i < n; ++i)
    {
        if (partial_mapping[i] != -1)
        {
            assigned_fac.push_back(i);
            assigned_loc.push_back(partial_mapping[i]);
        }
        else
            unassigned_fac.push_back(i);
    }
    for (int k = 0; k < m; ++k)
    {
        if (av[k])
            unassigned_loc.push_back(k);
    }

    int u = (int) unassigned_fac.size(); // number of unassigned facilities (rows)
    int r = (int) unassigned_loc.size(); // number of available locations (columns)

    assert(u > 0 && r >= u && "Error: Bounding failure.");

    vector<int> L(u*r, 0); // reduced L-matrix (r x r)

    int val, min_val;

    // Build reduced L-matrix
    for (int i_idx = 0; i_idx < u; ++i_idx)
    {
        int i = unassigned_fac[i_idx];

        for (int k_idx = 0; k_idx < r; ++k_idx)
        {
            int k = unassigned_loc[k_idx];
            int cost = 0;

            // Interaction with other unassigned facilities
            for (int j_idx = 0; j_idx < u; ++j_idx)
            {
                int j = unassigned_fac[j_idx];

                if (i == j)
                    continue;

                min_val = INF;

                for (int l_idx = 0; l_idx < r; ++l_idx)
                {
                    if (k_idx == l_idx)
                        continue;

                    int l = unassigned_loc[l_idx];

                    val = F[i][j] * D[k][l];

                    if (val < min_val)
                        min_val = val; 
                }
                cost += min_val;
            }

            // Interaction with assigned facilities
            for (int a_idx = 0; a_idx < (int) assigned_fac.size(); ++a_idx)
            {
                int j = assigned_fac[a_idx];
                int l = partial_mapping[j];

                cost += F[i][j] * D[k][l];
            }

            L[i_idx*r + k_idx] = cost;
        }
    }

    return L;
}*/


// optimized
// complexity = O(u^2 r), if n ~ m, complexity = O(M^3) with M = m - <depth of the node>
vector<int> Node::Assemble_LAP (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m)
{
    vector<int> partial_mapping = this->mapping;
    vector<bool> av = this->available;
    int dp = this->depth;

    // Identify assigned and unassigned facilities/locations
    vector<int> assigned_fac, unassigned_fac;
    vector<int> assigned_loc, unassigned_loc;

    for (int i = 0; i < n; ++i)
    {
        if (partial_mapping[i] != -1)
        {
            assigned_fac.push_back(i);
            assigned_loc.push_back(partial_mapping[i]);
        }
        else
            unassigned_fac.push_back(i);
    }
    for (int k = 0; k < m; ++k)
    {
        if (av[k])
            unassigned_loc.push_back(k);
    }

    int u = (int) unassigned_fac.size(); // number of unassigned facilities (rows)
    int r = (int) unassigned_loc.size(); // number of available locations (columns)

    assert(u > 0 && r >= u && "Error: Bounding failure.");
    assert(u == n - dp && r == m - dp && "Error: Bounding failure.");

    vector<int> L(u*r, 0);

    // Precompute for each location the two smallest distances to other locations
    struct MinPair { int min1, idx1, min2; };
    vector<MinPair> best(r);

    for (int k_idx = 0; k_idx < r; ++k_idx)
    {
        int k = unassigned_loc[k_idx];
        int min1 = INF, idx1 = -1, min2 = INF;

        for (int l_idx = 0; l_idx < r; ++l_idx)
        {
            if (k_idx == l_idx)
                continue;

            int l = unassigned_loc[l_idx];
            int dist = D[k][l];

            if (dist < min1)
            {
                min2 = min1;
                min1 = dist;
                idx1 = l_idx;
            }
            else if (dist < min2)
            {
                min2 = dist;
            }
        }
        best[k_idx] = {min1, idx1, min2};
    }

    // Build reduced L-matrix
    for (int i_idx = 0; i_idx < u; ++i_idx)
    {
        int i = unassigned_fac[i_idx];

        for (int k_idx = 0; k_idx < r; ++k_idx)
        {
            int k = unassigned_loc[k_idx];
            int cost = 0;

            // Interaction with other unassigned facilities
            for (int j_idx = 0; j_idx < u; ++j_idx)
            {
                int j = unassigned_fac[j_idx];

                if (i == j)
                    continue;

                // Pick best or second-best distance if best is disallowed
                int d = (best[k_idx].idx1 == k_idx) ? best[k_idx].min2 : best[k_idx].min1;

                cost += F[i][j] * d;
            }

            // Interaction with assigned facilities
            for (int a_idx = 0; a_idx < (int) assigned_fac.size(); ++a_idx)
            {
                int j = assigned_fac[a_idx];
                int l = partial_mapping[j];

                cost += F[i][j] * D[k][l];
            }

            L[i_idx*r + k_idx] = cost;
        }
    }

    return L;
}


int Node::bound_GLB (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m)
{
    vector<int> partial_mapping = this->mapping;
    int dp = this->depth;

    vector<int> L = this->Assemble_LAP(D, F, n, m);

    int fixed_cost = Objective(partial_mapping, D, F, n);

    int remaining_lb = Hungarian(L, n - dp, m - dp);

    return fixed_cost + remaining_lb;
}


int Node::bound_GLB_shift (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m, float alpha)
{
    const vector<int>& partial_mapping = this->mapping;
    int dp = this->depth;
    int u = n - dp;
    int r = m - dp;

    // ----- 1. build the LAP matrix -----
    vector<int> L = this->Assemble_LAP(D, F, n, m);

    // ----- 2. compute baseline bound and get duals -----
    vector<int> rowDual, colDual;

    int fixed_cost = Objective(partial_mapping, D, F, n);
    int remaining_lb = HungarianWithDuals(L, u, r, rowDual, colDual);

    int baseBound = fixed_cost + remaining_lb;

    // ----- 3. build damped separable shift -----
    vector<int> L_shifted(L.size());
    for (int i = 0; i < u; ++i)
    {
        for (int k = 0; k < r; ++k)
        {
            int idx = i*r + k;
            double shift = - alpha * (static_cast<double>(rowDual[i]) + static_cast<double>(colDual[k]));

            L_shifted[idx] = static_cast<int>(std::round(L[idx] - shift));
        }
    }

    // ----- 4. recompute LAP lower bound on shifted matrix -----
    int remaining_lb_shifted = Hungarian(L_shifted, u, r);
    int improvedBound = fixed_cost + remaining_lb_shifted;

    // ----- 5. keep the better bound -----
    return std::max(baseBound, improvedBound);
}