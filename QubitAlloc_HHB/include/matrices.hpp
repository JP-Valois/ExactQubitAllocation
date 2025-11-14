#ifndef __MATRICES__
    #define __MATRICES__ 

    #include "utils.hpp"

    
    class CostMatrix
    {
        private:

            vector<int> costs;
            vector<int> leader;
            int size;

        public:

            CostMatrix () = default;

            CostMatrix (const vector<int>& costs0, const vector<int>& leader0, int n) : costs(costs0), leader(leader0), size(n) {}

            ~ CostMatrix () = default;

            const vector<int>& get_costs () const { return costs; }
            vector<int>& get_costs () { return costs; }
            const vector<int>& get_leader () const { return leader; }
            vector<int>& get_leader () { return leader; }
            int get_size () const { return size; }

            static CostMatrix Assemble (const vector<vector<int>>& D, const vector<vector<int>>& F, int n);

            void distributeLeader ();
            void distributeLeader_OMP ();

            void halveComplementary ();
            void halveComplementary_OMP ();

            CostMatrix reduce (int i, int j) const;

            friend std::ostream& operator<<(std::ostream& os, const CostMatrix& CM)
            {
                const int n = CM.size;
                const auto& C = CM.costs;
                const auto& L = CM.leader;

                os << "Cost matrix (size = " << n << ")\n";

                for (int i = 0; i < n; ++i)
                {
                    for (int k = 0; k < n; ++k)
                    {
                        os << "[ ";
                        for (int j = 0; j < n; ++j)
                        {
                            for (int l = 0; l < n; ++l)
                            {
                                if ((k == i) ^ (l == j))
                                    os << "* ";
                                else
                                    os << C[idx4D(i, j, k, l, n)] << " ";
                            }
                            os << "| ";
                        }
                        os << "\b\b]\n";
                    }
                    if (i != n - 1)
                    {
                        os << "  ";
                        for (int k = 0; k < n*n + n - 1; ++k)
                            os << "- ";
                        os << "\n";
                    }
                }

                os << "Leader =\n";
                os << L;
                return os;
            }

            // Boost serialization function
            template <class Archive>
            void serialize (Archive &ar, const unsigned int version)
            {
                ar & costs;
                ar & leader;
                ar & size;
            }
    };

    int Load_distanceMatrix (vector<vector<int>>& D, const std::string& file_path);

    int Load_interactionMatrix (vector<vector<int>>& M, const std::string& file_path, int m);

    vector<vector<int>> LinearCouplingDistanceMatrix (int n);

#endif
