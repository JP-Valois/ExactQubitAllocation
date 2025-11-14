#include "../include/hungarian.hpp"


/**
 * Sets a = min(a, b)
 * @return true if b < a
 */
bool ckmin(int &a, const int &b) 
{ 
    return b < a ? (a = b, true) : false; 
}

/**
 * Hungarian algorithm implementation that modifies the input cost matrix and returns
 * the optimal assignment cost, while modifying the permutation vector.
 *
 * @param C A reference to a matrix of dimensions n×n such that C[w][j] is the cost
 *          of assigning the w-th worker to the j-th job (possibly negative).
 * @param i0 row index of the cost sub-matrix C
 * @param j0 column index of the cost sub-matrix C
 * @param n size of the input cost matrix C
 * @return The minimum assignment cost.
 */
int Hungarian (int* C, int i0, int j0, int n)
{
    int w, j, w_cur, j_cur, j_next;

    const int inf = INF / 2;

    // job[j] = worker assigned to job j, or -1 if unassigned
    vector<int> job(n + 1, -1);

    // yw[w] is the potential for worker w
    // yj[j] is the potential for job j
    vector<int> yw(n, 0), yj(n + 1, 0);

    // main Hungarian algorithm
    for (w_cur = 0; w_cur < n; ++w_cur)
    {
        j_cur = n;
        job[j_cur] = w_cur;

        vector<int> min_to(n + 1, inf);
        vector<int> prv(n + 1, -1);
        vector<bool> in_Z(n + 1, false);

        while (job[j_cur] != -1)
        {
            in_Z[j_cur] = true;
            w = job[j_cur];
            int delta = inf;
            j_next = 0;

            for (j = 0; j < n; ++j)
            {
                if (!in_Z[j])
                {
                    // reduced cost = C[w][j] - yw[w] - yj[j]
                    int cur_cost = C[idx4D(i0, j0, w, j, n)] - yw[w] - yj[j];

                    if (ckmin(min_to[j], cur_cost))
                        prv[j] = j_cur;
                    if (ckmin(delta, min_to[j])) 
                        j_next = j;
                }
            }

            // update potentials
            for (j = 0; j <= n; ++j)
            {
                if (in_Z[j])
                {
                    yw[job[j]] += delta;
                    yj[j] -= delta;
                }
                else
                {
                    min_to[j] -= delta;
                }
            }

            j_cur = j_next;
        }

        // update worker assignment along the found augmenting path
        for (; j_cur != n; j_cur = j)
        {
            j = prv[j_cur];
            job[j_cur] = job[j];
        }
    }

    // compute total cost
    int total_cost = 0;

    // for j in [0..n-1], job[j] is the worker assigned to job j
    for (j = 0; j < n; ++j)
    {
        if (job[j] != -1)
            total_cost += C[idx4D(i0, j0, job[j], j, n)];
    }

    // OPTIONAL: Reflecting the "reduced costs" after the Hungarian
    // algorithm by applying the final potentials:
    for (w = 0; w < n; ++w)
    {
        for (j = 0; j < n; ++j)
        {
            if (C[idx4D(i0, j0, w, j, n)] < inf)
            {
                // subtract the final potentials from the original cost
                C[idx4D(i0, j0, w, j, n)] = C[idx4D(i0, j0, w, j, n)] - yw[w] - yj[j];
            }
        }
    }

    return total_cost;
}
