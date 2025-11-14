#ifndef __NODE__
    #define __NODE__

    #include "preprocessing.hpp"
    #include "objective.hpp"

    class Node
    {
        private:

            vector<int> mapping;
            vector<bool> available;
            int depth;

        public:

            Node () = default;

            ~ Node () = default;

            Node (vector<int> mapping0, vector<bool> available0, int depth0) : mapping(mapping0),
                                                                               available(available0),
                                                                               depth(depth0) {}

            const vector<int>& get_mapping () const { return mapping; }
            vector<int>& get_mapping () { return mapping; }
            int get_depth () const { return depth; }
            int& get_depth () { return depth; }
            const vector<bool>& get_available () const { return available; }
            vector<bool>& get_available () { return available; }

            bool isLeaf (int n) const { return n == depth; }

            vector<int> Assemble_LAP (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m);

            int bound_GLB (const vector<vector<int>>& D,  const vector<vector<int>>& F, int n, int m);
            int bound_GLB_shift (const vector<vector<int>>& D, const vector<vector<int>>& F, int n, int m, float alpha);

            vector<Node> decompose (const vector<int>& priority, int n, int m);

        static Node Root (int n, int m);
    };

    double ComputeTotalNumberOfNodes (int n, int m);

#endif
