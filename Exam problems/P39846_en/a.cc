#include <iostream>
#include <vector>
#include <queue>

using namespace std;

class TreasureMap
{
private:
    struct Node
    {
        Node*  m_parent;
        size_t m_distance;
        size_t m_i, m_j;

        Node() 
        : m_parent(nullptr)
        {};

        Node(size_t i, size_t j)
        : m_i(i)
        , m_j(j)
        , m_parent(nullptr)
        {};
    };

    vector<vector<char>>  m_map;
    vector<vector<Node*>> m_nodes;
    vector<Node*>         m_treasures;
    queue<Node*>          m_queue;

    void visit(Node* parent, size_t i, size_t j)
    {
        if (i < 0 || i >= m_map.size())     return;
        if (j < 0 || j >= m_map[0].size())  return;

        if (m_map[i][j] == 'X')       return;
        if (m_nodes[i][j] != nullptr) return;
        Node* n_node = new Node(i, j);
        n_node->m_parent = parent;
        n_node->m_distance = parent->m_distance + 1;
        m_nodes[i][j] = n_node;
        m_queue.push(n_node);
        if (m_map[i][j] == 't')       
            m_treasures.push_back(n_node);
    }

    void shortest_path_treasures(size_t i, size_t j)
    {
        m_queue     = queue<Node*>();
        m_nodes     = vector<vector<Node*>>(m_map.size()
                    , vector<Node*>(m_map[0].size(), nullptr)); 
        m_treasures = vector<Node*>();

        Node* current = new Node(i, j);
        current->m_distance = 0;
        m_queue.push(current);
        m_nodes[i][j] = current;
        if (m_map[i][j] == 't')       
            m_treasures.push_back(current);

        while(!m_queue.empty())
        {
            current = m_queue.front();
            m_queue.pop();

            // Visitamos los nodos adyacentes.
            visit(current, current->m_i + 1, current->m_j    );
            visit(current, current->m_i - 1, current->m_j    );
            visit(current, current->m_i    , current->m_j + 1);
            visit(current, current->m_i    , current->m_j - 1);
        } 
    }

public:
    TreasureMap() {}
    void read(size_t i, size_t j)
    {
        m_map       = vector<vector<char>>(i, vector<char>(j)); 

        for (size_t _i = 0; _i < i; ++_i)
            for (size_t _j = 0; _j < j; ++_j)
                cin >> m_map[_i][_j];
    }

    int furthest_treasure(size_t i, size_t j)
    {
        shortest_path_treasures(i, j);
        if (m_treasures.size() == 0) return -1;
        return m_treasures[m_treasures.size() - 1]->m_distance;
    } 
};

int main()
{
    size_t i, j; cin >> i >> j;
    TreasureMap treasureMap;
    treasureMap.read(i, j);
    cin >> i >> j;
    int x = treasureMap.furthest_treasure(i - 1, j - 1);
    if (x == -1) cout << "no treasure can be reached" << endl;
    else cout << "maximum distance: " << x << endl; 
    
}