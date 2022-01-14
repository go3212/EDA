#include <iostream>
#include <vector>
#include <string>
#include <stack>

using namespace std;

class Graph
{
private:
    struct Node
    {
        size_t id;
        string m_word;
        vector<Node*> m_childs;

        Node(const string& word)
        : m_word(word)
        , m_childs(vector<Node*>())
        {};
    };

    vector<Node*>  m_nodes;

    const void printIdWord(const vector<size_t>& word) const
    {
        for (size_t i = 0; i < word.size(); ++i)
        {
            cout << m_nodes[word[i]]->m_word;
        }
        cout << endl;
    }

    void printAllConcat(Node* node, vector<size_t>& word, vector<bool> visited)
    {
        word.push_back(node->id);
        visited[node->id] = true;
        if (word.size() == m_nodes.size()) printIdWord(word);
        for (auto& c_node : node->m_childs)
        {
            if (!visited[c_node->id])
                printAllConcat(c_node, word, visited);
        }
        word.pop_back();
    }
    
public:
    Graph()
    : m_nodes(vector<Node*>())
    {};

    void addWord(const string& str)
    {
        Node* node = new Node(str);
        node->id = m_nodes.size();
        m_nodes.push_back(node);
    }

    void getCombinations()
    {
        // We first link all the words with each other.
        for (size_t i = 0; i < m_nodes.size(); ++i)
        {
            for (size_t j = 1 + i; j < m_nodes.size(); ++j)
            {
                auto& wordi = m_nodes[i]->m_word;
                auto& wordj = m_nodes[j]->m_word;
                if (wordi.back() != wordj[0]) 
                    m_nodes[i]->m_childs.push_back(m_nodes[j]);
                if (wordj.back() != wordi[0]) 
                    m_nodes[j]->m_childs.push_back(m_nodes[i]);
            }
        }
        // Now we run DFS and we print when we go through all nodes.
        for (auto& start_node : m_nodes)
        {
            vector<bool> visited(m_nodes.size(), false);
            vector<size_t> word;
            printAllConcat(start_node, word, visited); 
        }

    }

};

int main()
{
    size_t n; cin >> n;
    Graph graph;
    for (size_t i = 0; i < n; ++i)
    {
        string word; cin >> word;
        graph.addWord(word);
    }
    graph.getCombinations();

}