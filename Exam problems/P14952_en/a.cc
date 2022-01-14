
#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <queue>

using namespace std;



class Graph
{
private:
    struct pqcomp
    {
        bool operator()(size_t a, size_t b)
        {
            return a > b;
        }
    };

    vector<vector<size_t>>          m_adjacency_list;
    vector<size_t>                  m_topologicalSort;
    vector<char>                    m_mark;
    priority_queue<size_t
                 , vector<size_t>
                 , pqcomp>          m_stack;

    void m_topSort()
    {
        priority_queue<size_t, vector<size_t>, pqcomp> m_stack = priority_queue<size_t
                                                                              , vector<size_t>
                                                                              , pqcomp>();

        // We first need the set of all vertices without edges
        vector<size_t> node_inc_edges = vector<size_t>(m_adjacency_list.size(), 0);
        for (size_t node = 0; node < m_adjacency_list.size(); ++node)
            for (auto& node2 : m_adjacency_list[node])
                ++node_inc_edges[node2];
        // Now we search for those nodes without incoming edges
        for (size_t i = 0; i < node_inc_edges.size(); ++i)
            if (node_inc_edges[i] == 0) m_stack.push(i);

        // Now we run the algorithm.
        while(!m_stack.empty())
        {
            // The node is guaranteed to have no incoming edges, hence we insert it into the list.
            size_t node = m_stack.top(); m_stack.pop();
            m_topologicalSort.push_back(node);

            // Now we remove all the outgoing edges from our node.
            for (auto& conn_node : m_adjacency_list[node])
            {
                --node_inc_edges[conn_node];
                if (node_inc_edges[conn_node] == 0)
                    m_stack.push(conn_node);
            }
        }

        
    }

public:
    Graph(size_t MAX_NODES)
    {
        m_adjacency_list = vector<vector<size_t>>(MAX_NODES, vector<size_t>());
    }

    void link(size_t x, size_t y)
    {
        m_adjacency_list[x].push_back(y);
    }

    void topological_sort()
    {
        m_mark            = vector<char>(m_adjacency_list.size(), 0);
        m_topologicalSort = vector<size_t>();
        m_topSort();

    }

    void print_topSort()
    {
        int i;
        for (i = 0; i < m_topologicalSort.size() - 1; ++i)
            cout << m_topologicalSort[i] << ' ';
        cout << m_topologicalSort[i];
    }

    void print_links()
    {
        size_t i = 0;
        for (auto& node : m_adjacency_list)
        {
            for (auto& child : node)
                cout << i << ' ' << child << "   ";
            cout << endl;
            ++i;
        }
    }

};

int main()
{
    size_t n, m;
    while (cin >> n)
    {
        Graph myGraph(n);
        cin >> m;
        
        size_t x, y;
        for (size_t i = 0; i < m; ++i)
        {
            cin >> x >> y;
            myGraph.link(x, y);
        }
        myGraph.topological_sort();
        myGraph.print_topSort();
        cout << endl;
    }
}