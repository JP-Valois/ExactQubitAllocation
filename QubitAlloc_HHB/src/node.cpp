#include "../include/node.hpp"
#include "../include/hungarian.hpp"


Node Node::Root (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m)
{
    Solution sol;
    vector<int> map(n, -1);
    sol.mapping = map;
    sol.cost = -1;
    vector<bool> available(m, true);
    CostMatrix CM = CostMatrix::Assemble(D, F, m);

    return Node{sol, 0, 0, CM, available};
}


int localLogicalQubitIndex (const vector<int>& mapping, int i)
{
    int j{0}, k{0};

    while (j < i)
    {
        if (mapping[j++] == -1)
            ++k;  
    }

    return k;
}

int localPhysicalQubitIndex (const vector<bool>& av, int j)
{
    int l{0};

    for (int i = 0; i < j; ++i)
    {
        if (av[i])
        ++l;
    }

    return l;
}


int Node::bound (int it_max, int min_cost, int& it, bool early_stop, double& rt)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    int& lb = this->lower_bound;
    CostMatrix& CM = this->costMatrix;

    const int m = CM.get_size();

    assert(m > 0 && "Error: Cannot bound problem of size 0.");

    vector<int>& C = CM.get_costs();
    vector<int>& L = CM.get_leader();

    int i, j;
    int cost, incre, idx_submat;

    it = 0;

    while (it < it_max && lb <= min_cost)
    {
        ++ it;

        CM.distributeLeader();

        CM.halveComplementary();

        // apply Hungarian algorithm to each sub-matrix
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < m; ++j)
            {
                idx_submat = i*m + j;

                cost = Hungarian(C.data(), i, j, m);

                L[idx_submat] += cost;
            }
        }

        // apply Hungarian algorithm to the leader matrix
        incre = Hungarian(L.data(), 0, 0, m);

        if (early_stop && incre == 0)
            break;

        lb += incre;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delta = t1 - t0;
    rt += delta.count();

    return lb;
}


int Node::bound_OMP (int it_max, int min_cost, int& it, bool early_stop, double& rt)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    int& lb = this->lower_bound;
    CostMatrix& CM = this->costMatrix;

    const int m = CM.get_size();

    assert(m > 0 && "Error: Cannot bound problem of size 0.");

    vector<int>& C = CM.get_costs();
    vector<int>& L = CM.get_leader();

    int i, j;
    int cost, incre, idx_submat;
    int m2 = m*m;

    const int nb_threads = std::min(m2, omp_get_max_threads());

    it = 0;

    while (it < it_max && lb <= min_cost)
    {
        ++ it;

        CM.distributeLeader_OMP();

        CM.halveComplementary_OMP();

        #pragma omp parallel for collapse(2) default(none) private(i, j, cost, idx_submat) shared(C, L, m, m2) num_threads(nb_threads)
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < m; ++j)
            {
                idx_submat = i*m + j;

                cost = Hungarian(C.data(), i, j, m);

                L[idx_submat] += cost;
            }
        }

        incre = Hungarian(L.data(), 0, 0, m);

        if (early_stop && incre == 0)
            break;

        lb += incre;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delta = t1 - t0;
    rt += delta.count();

    return lb;
}


vector<Node> Node::decompose (const vector<int>& priority, int n, int m, int min_cost)
{
    vector<Node> children;

    Solution& sol = this->solution;
    const int sz = this->size;
    const int lb = this->lower_bound;
    const CostMatrix& CM = this->costMatrix;
    vector<bool>& av = this->available;

    // next logical qubit q_i to assign
    int i = priority[sz];

    // local index of q_i in the cost matrix
    int k = localLogicalQubitIndex(sol.mapping, i);

    // iterate over available physical qubits
    for (int j = m - 1; j >= 0; --j)
    {
        if (!av[j]) continue; // skip if not available

        // local index of P_j in the cost matrix
        int l = localPhysicalQubitIndex(av, j);
        
        // increment lower bound
        int incre = CM.get_leader()[k * (m - sz) + l];
        int lb_new = lb + incre;

        // prune
        if (lb_new > min_cost)
        {
            continue;
        }

        // assign q_i to P_j
        sol.mapping[i] = j;
        av[j] = false;

        // reduce cost matrix according to the new sub-problem
        CostMatrix CM_new = CM.reduce(k, l);

        // insert in children vector
        children.push_back(Node{sol, sz+1, lb_new, CM_new, av});

        // restore data
        sol.mapping[i] = -1;
        av[j] = true;
    }

    return children;
}
