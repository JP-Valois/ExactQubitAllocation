#ifndef __NODE__
    #define __NODE__

    #include "matrices.hpp"

    struct Solution
    {
        vector<int> mapping;
        int cost;

        // Boost serialization support
        template<class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & mapping;
            ar & cost;
        }
    };

    class Node
    {
        private:

            Solution solution;
            //vector<int> solution;
            int size;
            int lower_bound;
            CostMatrix costMatrix;
            vector<bool> available;

        public:

            Node () = default;

            ~ Node () = default;

            Node (Solution solution0, int size0, int lower_bound0, CostMatrix& costMatrix0, vector<bool> available0) : solution(solution0),
                                                                                                                       size(size0),
                                                                                                                       lower_bound(lower_bound0),
                                                                                                                       costMatrix(costMatrix0),
                                                                                                                       available(available0) {}

            const Solution& get_solution () const { return solution; }
            Solution& get_solution () { return solution; }
            int get_size () const { return size; }
            int get_lowerBound () const { return lower_bound; }
            const CostMatrix& get_costMatrix () const { return costMatrix; }
            CostMatrix& get_costMatrix () { return costMatrix; }
            const vector<bool>& get_available () const { return available; }
            vector<bool>& get_available () { return available; }

            bool isLeaf (int n) const { return n == size; }

            int bound (int it_max, int min_cost, int& it, bool early_stop, double& rt);
            int bound_OMP (int it_max, int min_cost, int& it, bool early_stop, double& rt);

            vector<Node> decompose (const vector<int>& priority, int n, int m, int min_cost);

            // Boost serialization support
            template <class Archive>
            void serialize (Archive &ar, const unsigned int version)
            {
                ar & solution;
                ar & size;
                ar & lower_bound;
                ar & costMatrix;
                ar & available;
            }

        static Node Root (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m);
    };

    double ComputeTotalNumberOfNodes (int n, int m);

#endif
