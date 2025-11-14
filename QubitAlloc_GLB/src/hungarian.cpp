#include "../include/hungarian.hpp"


/**
 * Sets a = min(a, b)
 * @return true if b < a
 */
bool ckmin(int &a, const int &b) 
{ 
    return b < a ? (a = b, true) : false; 
}

int HungarianWithDuals(const std::vector<int>& L, int n, int m, std::vector<int>& rowDual, std::vector<int>& colDual)
{
    const int INF2 = INF / 2;

    std::vector<int> job(m + 1, -1);
    std::vector<int> yw(n, 0), yj(m + 1, 0);

    for (int w_cur = 0; w_cur < n; ++w_cur)
    {
        int j_cur = m;             // dummy job
        job[j_cur] = w_cur;

        std::vector<int> min_to(m + 1, INF2);
        std::vector<int> prv(m + 1, -1);
        std::vector<bool> in_Z(m + 1, false);

        while (job[j_cur] != -1)
        {
            in_Z[j_cur] = true;
            int w = job[j_cur];
            int delta = INF2, j_next = 0;

            for (int j = 0; j < m; ++j)
            {
                if (!in_Z[j])
                {
                    int cur_cost = L[w * m + j] - yw[w] - yj[j];

                    if (ckmin(min_to[j], cur_cost))
                        prv[j] = j_cur;
                    if (ckmin(delta, min_to[j]))
                        j_next = j;
                }
            }

            // update potentials
            for (int j = 0; j <= m; ++j)
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

        // augment along the found path
        for (int j; j_cur != m; j_cur = j)
        {
            j = prv[j_cur];
            job[j_cur] = job[j];
        }
    }

    // compute total cost
    int total_cost = 0;
    for (int j = 0; j < m; ++j)
        if (job[j] != -1)
            total_cost += L[job[j] * m + j];

    // ---- export duals (ignore dummy job m) ----
    rowDual = yw;                                // n entries
    colDual.assign(yj.begin(), yj.begin() + m);  // m entries (discard dummy)

    return total_cost;
}


int Hungarian (const std::vector<int>& L, int n, int m)
{
    std::vector<int> dummyRow, dummyCol;
    return HungarianWithDuals(L, n, m, dummyRow, dummyCol);
}