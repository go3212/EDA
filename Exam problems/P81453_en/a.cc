#include <iostream>
#include <vector>
#include <queue>
#include <set>

using namespace std;

class Graph
{
private:
    struct Node 
    {
        size_t index;
        size_t weight;
        Node*  predecessor;

        Node(size_t x)
        : index(x)
        , predecessor(nullptr)
        , weight(0)
        {}
    };

    vector<set<size_t>> adj_list;
    vector<size_t>      backtrack;
    
public:
    Graph(size_t MAX_NODES)
    : adj_list(vector<set<size_t>>(MAX_NODES, set<size_t>()))
    {};

    void add_directed_edge(size_t x, size_t y)
    {
        adj_list[x].insert(y);
    }

    void pathfinder(size_t x, size_t y)
    {
        queue<Node*> myQueue;
        vector<bool> visited(adj_list.size(), 0);

        // We start at the x node.
        Node* current;
        myQueue.push(new Node(x));
        while(!myQueue.empty())
        {
            current = myQueue.front(); myQueue.pop();
            visited[current->index] = true;
            if (current->index == y) break;

            for (auto& e_node : adj_list[current->index])
            {
                if (!visited[e_node])
                {
                    Node* n_node = new Node(e_node);
                    n_node->predecessor = current;
                    myQueue.push(n_node);
                }
            }
        }
        backtrack = vector<size_t>();
        while(current != nullptr)
        {
            backtrack.push_back(current->index);
            current = current->predecessor;
        }
    }

    void print()
    {
        for (size_t i = backtrack.size() - 1; i > 0; --i)
        {
            cout << backtrack[i] << ' ';
        }
        cout << backtrack[0] << endl;
    }

};

int main()
{
    size_t n, m;
    while(cin >> n)
    {
        cin >> m;
        Graph graph(n);
        size_t x, y;
        for (size_t i = 0; i < m; ++i)
        {
            cin >> x >> y;
            graph.add_directed_edge(x, y);
        }
        graph.pathfinder(0, n - 1);
        graph.print();

    }
}