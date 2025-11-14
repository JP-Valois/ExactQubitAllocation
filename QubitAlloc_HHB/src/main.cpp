#include "../include/node.hpp"
#include "../include/heuristics.hpp"
#include "../include/objective.hpp"


// ARGUMENTS:
// argv[1] = interaction matrix path
// argv[2] = distance matrix path
// argv[3] = it_max                           (optional, default: 5)
// argv[4] = sequential or parallel bounding  (optional, default: true)


int main (int argc, char** argv)
{
    
    //----------------------------------- Parameterization -----------------------------------//

    // assertion
    if (argc < 3)
    {
        std::cerr << "Error: two arguments are required\n";
        std::cerr << "(1): interaction matrix path\n";
        std::cerr << "(2): distance matrix path" << std::endl;
        exit(1);
    }

    // distance matrix & physical qubits
    string distance_matrix_path = argv[2];
    vector<vector<int>> D;
    int m = Load_distanceMatrix(D, distance_matrix_path);

    // interaction matrix & logical qubits
    string interaction_matrix_path = argv[1];
    vector<vector<int>> F;
    int n = Load_interactionMatrix(F, interaction_matrix_path, m);

    // maximal number of bounding iterations
    int it_max = 5;
    if (argc > 3)
    {
        it_max = std::stoi(argv[3]);
    }

    // sequential or parallel bounding
    bool bound_par = true;
    if (argc > 4)
    {
        bound_par = static_cast<bool>(std::stoi(argv[4]));
    }

    // displays
    std::cout << "\nnumber of physical qubits = " << m << "\n";
    std::cout << "number of logical qubits = " << n << std::endl;
    std::cout << "\nit_max = " << it_max << std::endl;
    if (bound_par)
        std::cout << "number of threads = " << omp_get_max_threads() << std::endl;

    // profiling variables
    double rt, rt_bound{0};


    //----------------------------------- Initialization -----------------------------------//

    // START TIMER //
    auto begin = std::chrono::high_resolution_clock::now();
    
    // prioritization queue of logical qubits
    vector<int> priority = Prioritization(F, n, m);

    // greedy heuristic allocation
    vector<int> map0;
    int min_cost = GreedyAllocation(D, F, priority, n, m, map0);

    // TEST //
    //std::cout << "\ngreedy = " << map0 << std::endl;
    //std::cout << "cost = " << min_cost << std::endl;
    //////////

    // variables
    vector<vector<int>> optimal_solutions;
    double nb_nodes_visited = 0;
    double bound_iter_sum = 0;
    double nb_nodes_bounded = 0;
    int it;


    //----------------------------------- B&B tree exploration -----------------------------------//

    std::cout << "\n...exploration\n" << std::endl;

    stack<Node> Pool;
    Pool.push(Node::Root(D, F, n, m));

    // TEST //
    /*
    std::cout << F << std::endl;
    std::cout << D << std::endl;
    Node root0 = Node::Root(D, F, n, m);
    std::cout << root0.get_costMatrix() << std::endl;
    */
    //////////

    while (!Pool.empty())
    {
        // extract next node
        Node node = Pool.top();
        Pool.pop();

        ++ nb_nodes_visited;

        // leaf node
        if (node.isLeaf(n))
        {
            const vector<int>& map = node.get_solution().mapping;
            int cost = ObjectiveFunction(map, D, F, n);

            if (cost < min_cost)
            {
                optimal_solutions = {map};
                min_cost = cost;
            }
            else if (cost == min_cost)
            {
                optimal_solutions.push_back(map);
            }

            continue;
        }
        
        // bounding
        int lb;

        if (bound_par)
        {
            lb = node.bound_OMP(it_max, min_cost, it, true, rt_bound);
        }
        else
        {
            lb = node.bound(it_max, min_cost, it, true, rt_bound);
        }

        ++ nb_nodes_bounded;
        bound_iter_sum += it;

        // pruning
        if (lb > min_cost)
        {
            continue;
        }

        // decomposition
        vector<Node> children = node.decompose(priority, n, m, min_cost);

        // append in Pool
        for (Node node : children)
            Pool.push(node);
    }


    //----------------------------------- Termination -----------------------------------//

    // STOP TIMER //
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - begin;
    rt = duration.count();

    // results display
    std::cout << "number of solutions = " << optimal_solutions.size() << "\n";

    assert(ValidSolutions(optimal_solutions, D, F, n, min_cost) && "Error: Invalid set of optimal solutions");
    std::cout << "cost = " << min_cost << std::endl;

    // visited nodes & meain bounding iterations display
    const double nb_nodes_total = ComputeTotalNumberOfNodes(n, m);
    std::cout << "\ntotal number of nodes = " << nb_nodes_total << "\n";
    std::cout << "number of visited nodes = " << nb_nodes_visited << "\n";
    std::cout << "bounding mean iter = " << bound_iter_sum / nb_nodes_bounded << std::endl;

    // runtime 
    std::cout << "\nruntime = " << rt << " s\n" << std::endl;
    std::cout << "bounding: " << (rt_bound / rt) * 100 << " %\n" << std::endl;


    return 0;
}
