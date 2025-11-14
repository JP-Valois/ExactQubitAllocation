#include "../include/utils.hpp"


double ComputeTotalNumberOfNodes (int n, int m)
{
    double nb_nodes{1}, nb_nodes_layer{1};
    
    for (int i = 0; i < n; ++i)
    {
        nb_nodes_layer *= (m - i);
        nb_nodes += nb_nodes_layer;
    }

    return nb_nodes;
}


vector<int> Arange (int n)
{
    vector<int> v(n, 0);

    for (int i = 1; i < n; ++i)
        v[i] = i;

    return v;
}


int idx4D (int i, int j, int k, int l, int n)
{
    return n*n*n*i + n*n*j + n*k + l;
}
