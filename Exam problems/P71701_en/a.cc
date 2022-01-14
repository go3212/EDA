#include <iostream>
#include <vector>
#include <queue>
using namespace std;

class NKings
{
private:
    struct Node
    {
        size_t n_kings;
        vector<vector<char>> m_map;

        Node() 
        : n_kings(0)
        {};
    };

    vector<vector<char>> m_solutions;
    queue<Node*> m_queue;

    void place_king(int i, int j, vector<vector<char>>& map)
    {
        map[i][j] = 'K';
        if (i - 1 >= 0)            map[i - 1][j] = 'X';
        if (i + 1 < map.size())    map[i + 1][j] = 'X';
        if (j - 1 >= 0)            map[i][j - 1] = 'X';
        if (j + 1 < map[0].size()) map[i - 1][j] = 'X';
    }

    void run_analysis()
    {

    }

public:
    NKings()
    {}

    void nkings(size_t n, size_t k)
    {
        // We'll only initialize kings for first row
        for (size_t i = 0; i < k; ++i)
        {
            Node* node   = new Node();
            node->m_map  = vector<vector<char>>(k, vector<char>(k, '.'));
            place_king(i, 0, node->m_map);
            m_queue.push(node);
        }
    }
};

int main()
{


}