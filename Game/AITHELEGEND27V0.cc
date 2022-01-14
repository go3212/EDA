#include "Player.hh"

#include <iostream>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <memory>
#include <stack>
#include <map>
#include <cstdlib>
#include <queue>

using namespace std;

enum COORD { x, y, z };

int sign(int num)
{
    if (num == 0) return 0;
    if (num > 0) return 1;
    if (num < 0) return -1;
}
size_t abs_m(int x)
{
    if (x < 0) return -1*x;
    return x;
}

// Normalized vector in squares
struct Vec2d
{
    int i, j;
    Vec2d()
    : i()
    , j()
    {}

    Vec2d(const int& i, const int& j)
    : i(i)
    , j(j)
    {}

    Vec2d(Dir dir) 
    : i(0)
    , j(0)
    {
            switch(dir)
            {
                case LEFT:  j = -1; break;
                case RIGHT: j = 1;  break;
                case TOP:   i = -1; break;
                case BOTTOM: i = 1; break;
            }
    }

    // Builds cell vector from origin (p) and destination (q)
    Vec2d(const Pos& p, const Pos& q)
    : i(sign(q.i - p.i))
    , j(sign(q.j - p.j))
    {}

    Dir toDir()
    {
        // We prioritize directions in i axis.
        if (i > 0) return BOTTOM;
        if (i < 0) return TOP;
        if (j < 0) return LEFT;
        if (j > 0) return RIGHT;
        return NONE;
    }
};

Vec2d operator+(const Vec2d& p, const Vec2d& q) {return Vec2d(p.i + q.i, q.j + p.j); };

Pos operator+(const Pos& p, const Vec2d& vec) { return Pos(p.i + vec.i, p.j + vec.j); }
Pos operator+(const Vec2d& vec, const Pos& p) { return Pos(p.i + vec.i, p.j + vec.j); }

size_t distance(const Pos& p, const Pos& q)  { return (q.i - p.i)*(q.i - p.i) + (q.j - p.j)*(q.j - p.j); };
struct distance_comparator
{
    Pos m_point;
    distance_comparator(const Pos& point)
    : m_point(point)
    {}
    bool operator() (const Pos& p, const Pos& q) { return distance(m_point, p) < distance(m_point, q); }
};

double operator^(const Vec2d& a, const Vec2d& b) { return a.j*b.j - a.i*b.i; }


vector<Pos>::iterator find_min_y (vector<Pos>& points)
{
    if (points.size() == 0) throw "Points set cannot be empty";
    vector<Pos>::iterator pmin  = points.begin();
    vector<Pos>::iterator begin;
    for (begin = ++points.begin(); begin != points.end(); ++begin)
    {
        if ((*begin).j < (*pmin).j)  
            pmin = begin;
        else if ((*begin).j == (*pmin).j && (*begin).i < (*pmin).i)  
            pmin = begin;
    }
    return pmin;
};

vector<Pos>::iterator find_max_y (vector<Pos>& points)
{
    if (points.size() == 0) throw "Points set cannot be empty";
    vector<Pos>::iterator pmax  = points.begin();
    vector<Pos>::iterator begin;
    for (begin = ++points.begin(); begin != points.end(); ++begin)
    {
        if ((*begin).j > (*pmax).j)  
            pmax = begin;
        else if ((*begin).j == (*pmax).j && (*begin).i > (*pmax).i)  
            pmax = begin;
    }
    return pmax;
}

/* To find orientation of ordered triplet (p, q, r).
 * The function returns following values
 * 0 --> p, q and r are collinear
 * 1 --> Clockwise
 * -1 --> Counterclockwise 
 */ 
int orientation (const Pos& p, const Pos& q, const Pos& r)
{
    int cross_product = (q.j - p.j)*(r.i - q.i) - (q.i - p.i)*(r.j - q.j);
    if (cross_product == 0) return 0;    // Colinear.
    return (cross_product > 0) ? 1 : -1; // CW / CCW

};

Pos pcomp;
static bool polar_comparator(const Pos& p, const Pos& q)
{
    int o = orientation(pcomp, p, q);
    if (o == 0)
        return (distance(pcomp, p) < distance(pcomp, q)) ? true : false;
 
   return (o == 1) ? false : true;
    
};

vector<Pos> convexHull (vector<Pos>& points)
{
    if (points.size() < 3) return vector<Pos>();

    // First we search for the bottom left point.
    vector<Pos>::iterator pmin = find_min_y(points);
    // We swap the points.
    Pos p_sw = points[0];
    points[0] = *pmin; *pmin = p_sw;
    pmin = points.begin();

    pcomp = *pmin;

    sort(points.begin() + 1, points.end(), polar_comparator);

    // We assume that a convex hull is possible
    vector<Pos> hull;
    hull.push_back(*pmin);
    hull.push_back(points[1]);
    hull.push_back(points[2]);
    for (size_t i = 3; i < points.size(); ++i)
    {
        while (hull.size() > 1 && orientation(*(hull.end() - 2), *(hull.end() - 1), points[i]) != -1)
            hull.pop_back();
        hull.push_back(points[i]);
    }

    return hull;
};


class kdtree
{
private:
    struct Node
    {
        Pos* m_point;
        Node* m_left;
        Node* m_right;
        //Node* m_parent;

        Node() 
            : m_point(nullptr)
            , m_left(nullptr)
            , m_right(nullptr)
            //, m_parent(nullptr)
        {}
        Node(Pos* point, Node* left, Node* right) 
            : m_point(point)
            , m_left(left)
            , m_right(right)
            //, m_parent(nullptr)
        {}
    };

    Node* m_root;
    vector<Pos> m_points;
    Node* m_best;
    double m_best_dist;
    size_t m_visited;

    struct point_cmp 
    {
        point_cmp(size_t index) : m_index(index) {}
        bool operator() (const Pos& n1, const Pos& n2) const 
        {
			switch (m_index)
			{
				case 0: return n1.i < n2.i;
				case 1: return n1.j < n2.j;
			}
            return false;
        }
        size_t m_index;
    };



    Node* make_tree(size_t begin, size_t end, size_t index)
    {
        if (begin >= end) return nullptr;
        size_t n = (begin + end)/2;
        nth_element(m_points.begin() + begin, m_points.begin() + n, m_points.begin() + end, point_cmp(index));
        index = (index + 1)%2;
        Node* m_left =  make_tree(begin, n, index);
        Node* m_right =  make_tree(n + 1, end, index);
        Node* _root = new Node(&m_points[n], m_left, m_right);
        //if(m_left != nullptr)  m_left->m_parent  = _root;
        //if(m_right != nullptr) m_right->m_parent = _root;
        return _root;
    }

	int indexpoint (const Pos& p, int index) const
	{
		switch (index)
		{
			case 0: return p.i;
			case 1: return p.j;
		}
	}

    void nearest(Node* root, const Pos& point, size_t index)
    {
        if (root == nullptr) return;
        ++m_visited;
        double dist = distance(*((*root).m_point), point);
        double dx = indexpoint(*(*root).m_point, index) - indexpoint(point, index);
		
        if (m_best == nullptr || dist < m_best_dist)
        {
            m_best = root;
            m_best_dist = dist;
        }
        if (m_best_dist == 0) return;

        index = (index + 1)%2;

        nearest((dx > 0) ? root->m_left : root->m_right, point, index);
        if (dx*dx > m_best_dist) return;
        nearest((dx > 0) ? root->m_right : root->m_left, point, index);
    }

    void nearest( Node* root
                , const Pos& point
                , size_t index
                , priority_queue<Pos, vector<Pos>, distance_comparator>& closest
                , int k)
    {
        if (root == nullptr) return;
        ++m_visited;
        double dist = distance(*((*root).m_point), point);
		
        double dx = indexpoint(*(*root).m_point, index) - indexpoint(point, index);
        if (closest.size() >= k && distance(point, *root->m_point) < distance(point, closest.top()))
        {
            closest.pop();
            closest.push(*root->m_point);
        }
        if (closest.size() < k) closest.push(*root->m_point);

        if (m_best == nullptr || dist < m_best_dist)
        {
            m_best = root;
            m_best_dist = dist;
        }
        //if (m_best_dist == 0) return;
        
        double dx2 = indexpoint(*(*root).m_point, index) - indexpoint(closest.top(), index);
        index = (index + 1)%2;

        nearest((dx2 > 0) ? root->m_left : root->m_right, point, index, closest, k);
        if (dx2*dx2 >= distance(point, closest.top())) return;
        nearest((dx2 > 0) ? root->m_right : root->m_left, point, index, closest, k);
    }

public:
	kdtree() = default;

    template<typename iterator>
    kdtree(const iterator begin, const iterator end) 
        : m_points(begin, end)
    {
        m_root = make_tree(0, m_points.size(), 0);
    }

    const Pos nearest(const Pos& point)
    {
        if (m_root == nullptr) throw std::logic_error("tree is empty");
        m_best = nullptr;
        m_visited = 0;
        m_best_dist = 0;
        nearest(m_root, point, 0);
        return *(*m_best).m_point;
    }

    vector<Pos> nearest(const Pos& point, int k)
    {
        if (m_root == nullptr) throw std::logic_error("tree is empty");
        m_best = nullptr;
        m_visited = 0;
        m_best_dist = 0;
        priority_queue<Pos, vector<Pos>, distance_comparator> closest(point);
        nearest(m_root, point, 0, closest, k);
        // cerr << m_visited << 'x' << endl;
        vector<Pos> r_vec(closest.size());
        int i = closest.size() - 1;
        while(!closest.empty()) r_vec[i] = closest.top(), closest.pop(), --i;
        return r_vec;
    }

    vector<Pos> getPoints()
    {
        return m_points;
    }

};

string pp (const Pos& p) {return "(" + to_string(p.i) + "," + to_string(p.j) + ")"; };

void print_points(const vector<Pos>& points, int i, int j)
{
    vector<vector<char>> map(i, vector<char>(j, '.'));
    for (auto& p : points)
    {
        map[p.i][p.j] = 'X'; 
    }
    for (auto& row : map)
    {
        for (int k = 0; k < row.size(); ++k)
        {
            cerr << row[k];
        }
        cerr << endl;
    }
}


class AStar
{
private:


    struct Node
    {
        Node* m_parent;

        Pos m_point;

        size_t m_g;
        size_t m_h;

        bool checked;

        Node()
        : m_parent(nullptr)
        , m_point()
        , m_g(0)
        , m_h(0)
        , checked(false)
        {};

        Node(const Pos& point)
        : m_parent(nullptr)
        , m_point(point)
        , m_g(0)
        , m_h(0)
        , checked(false)
        {};
        
        Node(const Pos& pos, const Pos& dest)
        : m_parent(nullptr)
        , m_point(pos)
        , m_g(0)
        , m_h(manhattanDistance(pos, dest))
        , checked(false)
        {};

        Node(Node* parent, const Pos& pos, const Pos& dest)
        : m_parent(parent)
        , m_point(pos)
        , m_g(parent->m_g + 1)
        , m_h(euclideanDistance(pos, dest))
        , checked(false)
        {};

        size_t getScore() const {return m_g + m_h; };
        const Pos pointDir(const Dir& dir) const 
        {
            switch(dir)
            {
                case LEFT:   return m_point + Vec2d(0 , -1);
                case RIGHT:  return m_point + Vec2d(0 , 1);
                case TOP:    return m_point + Vec2d(-1 , 0);
                case BOTTOM: return m_point + Vec2d(1 ,0);
            }
            return m_point;
        }
    };


    struct cmp_cost
    {
        bool operator()(const Node* n1, const Node* n2) const { return n1->getScore() > n2->getScore(); }
    };

    vector<vector<bool>> m_map;
    vector<vector<Node*>> m_nodes;
    priority_queue<Node*, vector<Node*>, cmp_cost> openSet;
    vector<Node*> closedSet;

    Node* visit(Node* root, const Pos p, const Pos dest)
    {
        if (p.i >= m_map.size() || p.i < 0) return nullptr;
        if (p.j >= m_map[0].size() || p.j < 0) return nullptr;
        if (m_nodes[p.i][p.j] != nullptr && m_nodes[p.i][p.j]->checked) return nullptr;
        // To vist a node, we first have to check if it has been visited. And whether it is a wall or not.
        if (m_map[p.i][p.j] == 1) return nullptr;
        if (m_nodes[p.i][p.j] == nullptr)
        {
            m_nodes[p.i][p.j] = new Node(root, p, dest);
            openSet.push(m_nodes[p.i][p.j]);
            return m_nodes[p.i][p.j];
        }
        // Now we check if the new path is better than the parent's, and we update accordingly.
        Node* d_node = m_nodes[p.i][p.j];
        if (root->m_g + 1 > d_node->m_g) return d_node;
        d_node->m_parent = root;
        d_node->m_g = root->m_g + 1;
        return d_node;
    }

    Node* visit(Node* root, const Pos p, const Pos dest, const vector<vector<bool>>& mask)
    {
        if (p.i >= m_map.size() || p.i < 0) return nullptr;
        if (p.j >= m_map[0].size() || p.j < 0) return nullptr;
        if (m_nodes[p.i][p.j] != nullptr && m_nodes[p.i][p.j]->checked) return nullptr;
        // If the node we're accessing is masked, and it is the target, we still return it.
        if (mask[p.i][p.j] && p == dest) 
        {
            m_nodes[p.i][p.j] = new Node(root, p, dest);
            openSet.push(m_nodes[p.i][p.j]);
            return m_nodes[p.i][p.j];
        }

        // To vist a node, we first have to check if it has been visited. And whether it is a wall or not.
        if (m_map[p.i][p.j] == 1 || mask[p.i][p.j] == 1) return nullptr;
        if (m_nodes[p.i][p.j] == nullptr)
        {
            m_nodes[p.i][p.j] = new Node(root, p, dest);
            openSet.push(m_nodes[p.i][p.j]);
            return m_nodes[p.i][p.j];
        }
        // Now we check if the new path is better than the parent's, and we update accordingly.
        Node* d_node = m_nodes[p.i][p.j];
        if (root->m_g + 1 > d_node->m_g) return d_node;
        d_node->m_parent = root;
        d_node->m_g = root->m_g + 1;
        return d_node;
    }

    vector<Pos> m_findPath(const Pos& source, const Pos& target)
    {
        m_nodes   = vector<vector<Node*>>(m_map.size(), vector<Node*>(m_map[0].size(), nullptr));
        closedSet = vector<Node*>();
        openSet   = priority_queue<Node*, vector<Node*>, cmp_cost>();
        Node* current = nullptr;
        vector<Pos> ret;

        openSet.push(new Node(source, target));
        while(!openSet.empty())
        {
            current = openSet.top();
            openSet.pop();
            closedSet.push_back(current);
            current->checked = true;
            if (current->m_point == target) break;

            visit(current, current->pointDir(LEFT), target);
            visit(current, current->pointDir(RIGHT), target);
            visit(current, current->pointDir(BOTTOM), target);
            visit(current, current->pointDir(TOP), target);
        }
        // Current must contain the node path. We backtrack and then we deallocate memory.
        while (current->m_point != source)
        {
            ret.push_back(current->m_point);
            current = current->m_parent;

        }
        for (auto& node : closedSet) delete node;
        while(!openSet.empty()) delete openSet.top(), openSet.pop();        
        return ret;
    }

    vector<Pos> m_findPath(const Pos& source, const Pos& target, const vector<vector<bool>>& mask)
    {
        m_nodes   = vector<vector<Node*>>(m_map.size(), vector<Node*>(m_map[0].size(), nullptr));
        closedSet = vector<Node*>();
        openSet   = priority_queue<Node*, vector<Node*>, cmp_cost>();
        Node* current = nullptr;
        vector<Pos> ret;

        openSet.push(new Node(source, target));
        while(!openSet.empty())
        {
            current = openSet.top();
            openSet.pop();
            closedSet.push_back(current);
            current->checked = true;
            if (current->m_point == target) break;

            visit(current, current->pointDir(LEFT), target, mask);
            visit(current, current->pointDir(RIGHT), target, mask);
            visit(current, current->pointDir(BOTTOM), target, mask);
            visit(current, current->pointDir(TOP), target, mask);
        }
        // Current must contain the node path. We backtrack and then we deallocate memory.
        while (current->m_point != source)
        {
            ret.push_back(current->m_point);
            current = current->m_parent;
        }
        for (auto& node : closedSet) delete node;
        while(!openSet.empty()) delete openSet.top(), openSet.pop();     
        // for (int i = 0; i < m_map.size(); ++i) 
        // {
        //     for (int j = 0; j < m_map[0].size(); ++j) cerr << m_map[i][j] ? 'x' : ' ';  
        //     cerr << endl;
        // }

        // print_points(ret, m_map.size(), m_map[0].size()); 
        return ret;
    }

public:
    static size_t manhattanDistance(const Pos p, const Pos q) { return abs(q.i - p.i) + abs(q.j - q.j); }
    static size_t euclideanDistance(const Pos p, const Pos q) { return distance(p, q); }

    AStar() {};

    AStar(vector<vector<bool>> m_map)
    : m_map(m_map)
    {};

    vector<Pos> findPath(const Pos& source, const Pos& target, const vector<vector<bool>>& mask)
    {
        return m_findPath(source, target, mask);
    }

    vector<Pos> findPath(const Pos& source, const Pos& target)
    {
        return m_findPath(source, target);
    }

    Vec2d findBestDir(const Pos& source, const Pos& target)
    {
        vector<Pos>m_path = m_findPath(source, target);
        return Vec2d(source, m_path.size() > 0 ? m_path[m_path.size() - 1] : source);
    }

    Vec2d findBestDir(const Pos& source, const Pos& target, const vector<vector<bool>>& mask)
    {
        vector<Pos>m_path = m_findPath(source, target, mask);
        return Vec2d(source, m_path.size() > 0 ? m_path[m_path.size() - 1] : source);
    }

};

/**
 * Write the name of your player and save this file
 * with the same name and .cc extension.
 */
#define PLAYER_NAME THELEGEND272

struct PLAYER_NAME : public Player {
	kdtree kdcity;
	kdtree kdpaths;
    kdtree kdalies;
    kdtree kdpath_entrance;
    AStar  astar;
    vector<vector<bool>> mask;
    

  	/**
  	 * Do not modify this function.
  	 */
  	static Player* factory () { return new PLAYER_NAME; };

    struct StructureInfo
    {
        size_t alies, enemies;
        StructureInfo() : alies(0), enemies(0) {}
    };

    struct CityInfo
    {
        size_t protecting;
        size_t cells;  // The number of cells it has.
        vector<size_t> paths_id;

        CityInfo() : protecting(0), cells(0)
        {}
    };

    vector<StructureInfo> cities, paths;
    vector<CityInfo> citiesinfo;

    enum Quest { GOTO_ALIE, GOTO_CLOSEST, CONQUER_CITY, GOTO_UNCONQUERED_CITY, CONQUER_PATH, GOTO_MASK, KILL_ENEMY, INFECT_ENEMIES, GOTO_REQUEST, CONQUER_CLOSEST, PROTECT_CITY, PROTECT_PATH};
    enum Params { AVOID_VIRUS, AVOID_ENEMIES };
    struct Memory
    {
        Pos              pos;
        int              health;
        Quest            quest;
        Pos              target;
        Pos              protecting;
        bool             has_target;

        Memory()
        : pos()
        , quest(GOTO_CLOSEST)
        , target()
        , has_target(false)
        {};

        void setTarget(const Pos& pos) { target = pos, has_target = true; }

        void clearTarget() { has_target = false; }
    };

    map<int, Memory> units;

    void scan_surrounding(const Pos& pos, vector<Pos>& alies, vector<Pos>& enemies, vector<Pos>& infected, vector<Pos>& cities, vector<Pos>& paths, const int RADIUS)
    {
        // We scan in a diamond-like way, because enemies cannot move diagonally, doesnt make sense to check further.
        for (int i = -1*RADIUS; i <= RADIUS; ++i)
        {
            for (int j = -1*RADIUS; j <=  RADIUS; ++j)
            {
                if (i == 0 && j == 0) continue;
                if (abs(i) + abs(j) > RADIUS) continue;
                Pos cell_ = Pos(pos.i+i, pos.j+j);
                if (cell_.i >= rows() || cell_.i < 0) break;
                if (cell_.j >= cols() || cell_.j < 0) break;

                Cell i_cell = cell(cell_);
                switch(i_cell.type)
                {
                    WALL: break;
                    CITY: if (i_cell.unit_id == -1) cities.push_back(cell_); break;
                    PATH: if (i_cell.unit_id == -1) cities.push_back(cell_); break;
                    GRASS: break;
                }
                if(i_cell.unit_id != -1) 
                {
                    if (unit(i_cell.unit_id).player != me()) enemies.push_back(cell_);
                    else alies.push_back(cell_);
                }
                if (i_cell.virus > 0) infected.push_back(cell_);
            }
        }     
    }   

    void scan_surrounding_fixed(const Pos& pos, vector<Pos>& alies, vector<Pos>& enemies, vector<Pos>& infected, vector<Pos>& cities, vector<Pos>& paths, const int RADIUS)
    {
        // We scan in a diamond-like way, because enemies cannot move diagonally, doesnt make sense to check further.
        for (int i = -1*RADIUS; i <= RADIUS; ++i)
        {
            for (int j = -1*RADIUS; j <=  RADIUS; ++j)
            {
                if (abs(i) + abs(j) != RADIUS) continue;
                Pos cell_ = Pos(pos.i+i, pos.j+j);
                if (cell_.i >= rows() || cell_.i < 0) break;
                if (cell_.j >= cols() || cell_.j < 0) break;

                Cell i_cell = cell(cell_);
                switch(i_cell.type)
                {
                    WALL: break;
                    CITY: if (i_cell.unit_id == -1) cities.push_back(cell_); break;
                    PATH: if (i_cell.unit_id == -1) paths.push_back(cell_); break;
                    GRASS: break;
                }
                if(i_cell.unit_id != -1) 
                {
                    if (unit(i_cell.unit_id).player != me()) enemies.push_back(cell_);
                    else alies.push_back(cell_);
                }
                if (i_cell.virus > 0) infected.push_back(cell_);
            }
        }
                
    }   

    // Computes the best move within a region to escape enemies' sight 
    Vec2d avoid_enemies(const Pos& pos, const vector<Pos>& enemies)
    {
        if (enemies.size() == 0) return Vec2d(NONE);
        vector<bool> dirs(4, 1);
        for (auto& enemy : enemies)
            dirs[Vec2d(pos, Pos(enemy.i, enemy.j)).toDir()] = 0;
        // We fetch all the available directions.
        vector<Vec2d> dir_;
        for (int i = 0; i < dirs.size(); ++i) if (dirs[i] && !mask[pos.i + Vec2d((Dir)i).i][pos.j + Vec2d((Dir)i).j]) dir_.push_back(Vec2d((Dir)i));
        if (dir_.size() == 0) return Vec2d(NONE);
        // Now we check if there is a direction in a path/city;
        for (auto& x : dir_) if (cell(pos + x).type == CITY || cell(pos + x).type == PATH) return x;
        return dir_[0];
    }

    // Returns the closest path that I dont own, it returns the closest path that I own if I own all the closest paths.
    Pos closest_path (const Pos& pos, const int k)
    {
        int id_path = cell(pos).path_id;
        vector<Pos> points = kdpaths.nearest(pos, k);

        if (points.size() == 0) return pos; 
        vector<Pos> extremes;
        vector<int> path_ids;
        for (int i = 0; i < points.size(); ++i) 
        {
            int p_id = cell(points[i]).path_id;
            if (p_id == -1) continue;
            if (find(path_ids.begin(), path_ids.end(), p_id) == path_ids.end()) 
            {
                extremes.push_back(points[i]);
                path_ids.push_back(p_id);
            }
        }
        if (extremes.size() == 0) return pos;
        // We have all the points that match our criteria. We sort them by distance.
        sort(extremes.begin(), extremes.end(), 
        [pos](const Pos& p, const Pos& q) -> bool 
        {
            return distance(pos, p) < distance(pos, q);
        });

        if (id_path != -1 && id_path == cell(extremes[0]).path_id) extremes.erase(extremes.begin());

        //print_points(extremes, rows(), cols());
        // We first check if either of the paths is not conquered
        for (auto& p : extremes) if (path_owner(cell(p).path_id) == -1) return p;
        // Now we check if either of the paths is not conquered by me.
        for (auto& p : extremes) if (path_owner(cell(p).path_id) != me()) return p;
        // If there is no path close, we return the closest that i own.
        return extremes[0];
        
    }
    
    // Returns the closest city that I dont own, it returns the closest city that I own if I own all the closest cities.
    // It doesnt return the city I'm in if i'm in a city. k is the number of cities to search, the one i'm in included if i'm.
    Pos closest_city (const Pos& pos, const int k)
    {
        int id_city = cell(pos).city_id;
        vector<Pos> nearest_city_vec = kdcity.nearest(pos, k);
        if (nearest_city_vec.size() == 0) return Pos(0, 0);
        
        if (id_city != -1 && id_city == cell(nearest_city_vec[0]).city_id) nearest_city_vec.erase(nearest_city_vec.begin());

        // We first check if either a city is not conquered.
        for (auto& city : nearest_city_vec) if (city_owner(cell(city).city_id) == -1) return city;
        // Now we check the city with the less enemies per alies.
        for (auto& city : nearest_city_vec) if (cities[cell(city).city_id].enemies - citiesinfo[cell(city).city_id].protecting == 0) return city;
        // Now I check if either one of the closer cities is conquered by an enemie.
        for (auto& city : nearest_city_vec) if (city_owner(cell(city).city_id) != me()) return city;
        return nearest_city_vec[0];
    }

    // We fetch the closest unconquered path of a city.
    // The cell must be close to the city in order to perform a good scan.
    Pos closest_city_path (const Pos& pos, const size_t& id_city)
    {   
        #define SCAN_STRENGHT 10
        // We first do a scan of the citie's paths.
        vector<int> u_paths;
        for (auto& p_id : citiesinfo[id_city].paths_id) 
            if (path_owner(p_id) != me()) u_paths.push_back(p_id);
        if (u_paths.size() == 0) return pos;
        // We now have all the unconquered paths. We scan for all the close paths.
        vector<Pos> available_paths = kdpath_entrance.nearest(pos, citiesinfo[id_city].paths_id.size());
        //print_points(kdpath_entrance.getPoints(), rows(), cols());
        //print_points(available_paths, rows(), cols());
        // We return the first that matches an unconquered one.
        for (auto& p : available_paths)
        {
            if (find(u_paths.begin(), u_paths.end(), cell(p).path_id) != u_paths.end())
                return p;
        }
        // if there's no match, we return the current pos.
        return pos;
    }

    StructureInfo players_in_city(size_t id_city)
    {
        vector<Pos> city_points = city(id_city);
        StructureInfo info;
        for (auto& point : city_points) 
        {
            if (cell(point).unit_id == -1) continue;
            if(unit(cell(point).unit_id).player != me()) ++info.enemies;
            else if (unit(cell(point).unit_id).player == me()) ++info.alies;
        }
        return info;
    }    

    StructureInfo players_in_path(size_t id_path)
    {
        vector<Pos> path_points = path(id_path).second;
        StructureInfo info;
        for (auto& point : path_points) 
        {
            if (cell(point).unit_id == -1) continue;
            if(unit(cell(point).unit_id).player != me()) ++info.enemies;
            else if (unit(cell(point).unit_id).player == me()) ++info.alies;
        }
        return info;
    }

	void round_special_tasks()
	{
		switch (round())
		{
			// En la ronda 1 definiremos el kdtree que contiene todos los puntos de interes.
			case 0:
            {
				// Construimos el kdtree de ciudades.
				//vector<Pos> city_points;
                vector<Pos> city_center_points;
				size_t n_cities = nb_cities();
                cities = vector<StructureInfo>(n_cities);
                citiesinfo = vector<CityInfo>(n_cities);
				for (size_t i = 0; i < n_cities; ++i)
				{
					City _city = city(i);
                    citiesinfo[i].cells = _city.size();
					vector<Pos> hull = convexHull(_city);
                    //for (auto& point : hull) city_center_points.push_back(point);
                    auto miny = find_min_y(hull);
                    auto maxy = find_max_y(hull);
                    city_center_points.push_back(Pos((maxy->i - miny->i) >> 1,(maxy->j - miny->j) >> 1) + *miny);
				}
				kdcity = kdtree(city_center_points.begin(), city_center_points.end());

				vector<Pos> path_points;
                vector<Pos> path_entrance_points;
				// Ahora inertamos los paths en nuestro set de puntos y su kd.
				size_t n_paths = nb_paths();
                paths = vector<StructureInfo>(n_paths);
				for (size_t i = 0; i < n_paths; ++i)
				{
					Path _path = path(i);
                    citiesinfo[_path.first.first].paths_id.push_back(i);
                    citiesinfo[_path.first.second].paths_id.push_back(i);
                    auto hull = convexHull(_path.second);
                    for (auto& p : hull) 
                    {
                        size_t count = 0;
                        // A point is a extrem of a line if there is only one other point surronding it.
                        if (cell(p + Vec2d(TOP)).path_id != -1) ++count;
                        if (cell(p + Vec2d(BOTTOM)).path_id != -1) ++count;
                        if (cell(p + Vec2d(LEFT)).path_id != -1) ++count;
                        if (cell(p + Vec2d(RIGHT)).path_id != -1) ++count;
                        if (count == 1) path_entrance_points.push_back(p);
                    }
					for (auto& pos : _path.second) path_points.push_back(pos);
				}
				kdpaths = kdtree(path_points.begin(), path_points.end());
                kdpath_entrance = kdtree(path_entrance_points.begin(), path_entrance_points.end());
                // Ahora barrimos el mapa para encontrar los muros.
                vector<vector<bool>> walls(rows(), vector<bool>(cols(), 0));
                for (int i = 0; i < walls.size(); ++i)
                    for (int j = 0; j < walls[0].size(); ++j)
                        if (cell(i, j).type == WALL) walls[i][j] = 1;

                astar = AStar(walls);
            }
                break;
            default:
                // Actualizamos los enemigos/aliados.
                for (int s = 0; s < nb_cities(); ++s) cities[s] = players_in_city(s);
                for (int s = 0; s < nb_paths(); ++s)  paths[s] = players_in_path(s);
                
                break;
		}
	}


    Dir best_move(const Unit& unit_)
    {
        #define ATTACK_HEALTH 5

        // We first see if the unit is in the system.
        auto unit__ = units.find(unit_.id);
        if (unit__ == units.end()) units[unit_.id] = Memory(), unit__ = units.find(unit_.id), unit__->second.pos = unit_.pos;

        auto& unit  = unit__->second;
        // The unit stores its previous position, we'd like to see whether it has died or not.
        // The new position is at most manhattan distance 1. If we're within a close range i.e 2, we dont have to reset the memory.
        if(AStar::manhattanDistance(unit_.pos, unit.pos) > 4) 
        {
            if (unit.quest == PROTECT_CITY)
            {
                // If the unit was protecting, we must remove the protecting state.
                citiesinfo[cell(unit.protecting).city_id].protecting -= 1;
            }

            units[unit_.id] = Memory();
        }
        unit.pos    = unit_.pos;
        unit.health = unit_.health;

        Cell cell_ = cell(unit_.pos);

        // HAY QUE AÑADIR LOGICA DE PROTECCIÓN CONTRA VIRUS

        // We first scan our surrounding at distance 1. This procedure must be done always, it is independent from other actions.
        vector<Pos> alies, enemies, infected, cities, paths;
        scan_surrounding(unit_.pos, alies, enemies, infected, cities, paths, 1);
        // We create the path mask, later we'll update the mechanism.
        for (int i = 0; i < alies.size()  ; ++i)  mask[alies[i].i][alies[i].j] = 1;
        for (int i = 0; i < enemies.size(); ++i)  mask[enemies[i].i][enemies[i].j] = 1;
        //for (int i = 0; i < infected.size(); ++i)  mask[infected[i].i][infected[i].j] = 1;
        //for (int i = 0; i < infected.size(); ++i) mask[infected[i].i][infected[i].j] = 1;
        // For now, we'll not avoid infected cells. If enemies are a distance 1, 
        // we'll attack if our health is over 40, when the enemy cell has lower health 
        // than ours or if we have an alie closer than 2 units manhattan distance.
        // Si tenemos enemigos a 1 unidad de distancia.
        if (enemies.size() > 0)
        {
            //if (unit.health > 45 && enemies.size() == 1)
            //    return Vec2d(unit.pos, enemies[0]).toDir();
            //if (unit.health < ATTACK_HEALTH && alies.size() == 0)
            //    return avoid_enemies(unit.pos, enemies).toDir();
            return Vec2d(unit.pos, enemies[0]).toDir();
            
        }
        

        // We check distance 2 just in case.
        scan_surrounding_fixed(unit.pos, alies, enemies, infected, cities, paths, 2);

        // Before we execute our task, we'll se if we can protect a city/path.
        int city_id = cell(unit.pos).city_id;
        int path_id = cell(unit.pos).path_id;

        if (city_id != -1) // we are in a city or a path, we consider doing something about it.
        {
            // If we are in a city, we check if it has the minimum cells to protect it.
            if (city_id != -1 && PLAYER_NAME::citiesinfo[city_id].protecting < PLAYER_NAME::citiesinfo[city_id].paths_id.size())
            {
                unit.clearTarget();
                unit.quest = PROTECT_CITY;
                unit.protecting = unit.pos;
                PLAYER_NAME::citiesinfo[city_id].protecting += 1;
            }
        }



        Pos nearest_city, nearest_path;
        vector<Pos> nearest_city_vec, nearest_path_vec;
        size_t change_quest = false;
        size_t i = 0;
        #define MAX_SEARCH_QUERY 2
        do 
        {
            change_quest = false;
            switch (unit.quest)
            {
                case GOTO_ALIE:     break;
                case GOTO_CLOSEST:  
                    // The closest structure is either a city or a path. 
                    // Once we get there, if conquered, we move on to the closest conquerable structure.
                    if (cell_.type == GRASS)
                    {
                        nearest_city = kdcity.nearest(unit.pos);
                        nearest_path = kdpaths.nearest(unit.pos);
                        unit.setTarget(distance(unit.pos, nearest_city) < distance(unit.pos, nearest_path) 
                                       ? cities.size() == 0 ? kdcity.nearest(unit.pos) : cities[0]
                                       : paths.size() == 0 ? kdpaths.nearest(unit.pos) : paths[0]);
                        break;
                    }
                    // We update the quest, if we've exhausted our quest targets.
                    if (cell_.type == CITY || cell_.type == PATH) 
                    {
                        unit.clearTarget();
                        // We search our closest objective and we update accordingly.
                        nearest_city = closest_city(unit.pos, 3 + i);
                        nearest_path = closest_path(unit.pos, 20);

                        unit.quest = distance(unit.pos, nearest_city) < distance(unit.pos, nearest_path) ? CONQUER_CITY : CONQUER_PATH;
                        change_quest = true;
                    }
                    break;
                case CONQUER_CITY:     
                    // When we are in this state, we'd like to go to an unconquered city first, if possible. Enemies are avoided if necessary.
                    if (!unit.has_target) 
                    {
                        unit.setTarget(closest_city(unit.pos, 3 + i));
                        //break;
                    }
                    // In this quest, we accomplish our target when we get to a city, doesnt matter if it is
                    // our target's one, since we are searching for the 2 closest, they could be collinear.
                    // If the city we're going to has been conquered by me.
                    if (city_owner(cell(unit.target).city_id) == me() && cell(unit.pos).type == GRASS) 
                    {
                        unit.clearTarget();
                        if (i <= MAX_SEARCH_QUERY) change_quest = true;
                        if (unit.quest == PROTECT_CITY) break;
                        unit.quest = CONQUER_CLOSEST;
                        break;
                    }
                    // When we get to our target, we clear it.
                    if (cell(unit.pos).city_id == cell(unit.target).city_id)
                    { 
                        unit.clearTarget();
                        // ** AÑADIR MAS LOGICA, ESTADO DE PROTECCION Y DEMÁS.
                        if (i <= MAX_SEARCH_QUERY) change_quest = true;
                        if (unit.quest == PROTECT_CITY) break;
                        unit.quest = CONQUER_CLOSEST;
                    }
                    break;
                case CONQUER_PATH:     
                    // In this state we do de same as with the latter one. 
                    // We search for the closest "two paths". We'll search for a total of 15 points, 
                    // then we search the ones that have different id's.
                    if (!unit.has_target) 
                    {
                        unit.setTarget(closest_path(unit.pos, 20 + 2*i));
                        //break;
                    }
                    // If the target has been conquered before our cell got there, we go to the closest structure.
                    if (path_owner(cell(unit.target).path_id) == me() && cell(unit.pos).type == GRASS) 
                    {
                        unit.clearTarget();
                        if (i <= MAX_SEARCH_QUERY) change_quest = true;
                        unit.quest = CONQUER_CLOSEST;
                        break;
                    }
                    // We complete our quest when we get to a path, it doesnt matter if its our target
                    if (cell(unit.pos).path_id == cell(unit.target).path_id) 
                    { 
                        unit.clearTarget();
                        // ** AÑADIR MAS LOGICA, ESTADO DE PROTECCION Y DEMÁS.
                        if (i <= MAX_SEARCH_QUERY) change_quest = true;
                        unit.quest = CONQUER_CLOSEST;
                        break;
                    }
                    break;
                case CONQUER_CLOSEST:
                    // We try to conquer the closest unconquered structure.
                    // We update the quest, if we've exhausted our quest targets.
                    unit.clearTarget();
                    // We search our closest objective and we update accordingly.
                    nearest_city = closest_city(unit.pos, 6 + i);
                    nearest_path = closest_path(unit.pos, 15 + 2*i);
                    if (city_owner(cell(nearest_city).city_id) == -1 && path_owner(cell(nearest_path).path_id)  == -1)
                        unit.quest = distance(unit.pos, nearest_city) < distance(unit.pos, nearest_path) ? CONQUER_CITY : CONQUER_PATH;
                    else if (city_owner(cell(nearest_city).city_id) != me() && path_owner(cell(nearest_path).path_id) != me())
                        unit.quest = distance(unit.pos, nearest_city) < distance(unit.pos, nearest_path) ? CONQUER_CITY : CONQUER_PATH;
                    else if (city_owner(cell(nearest_city).city_id) != me())
                        unit.quest =  CONQUER_CITY;
                    else unit.quest = CONQUER_PATH;
                        
                    // Now we must check if either the path or the city are conquered.
                    // For comunication, if everything's conquered we stay in a city.
                    change_quest = true;
                    break;
                case PROTECT_CITY:
                {
                    #define SEARCH_PATH_MULTIPLIER 10
                    // We must protect our city against intruders, whilst conquering a close path. If our city is endangered, we protect it at all cost.
                    // A safe city is efficient at conquering paths.
                    // A cell always goes to the closest unconquered path whithin its range. We search 5 cells*path.
                    // If we are in a path, or in grass, we must check wheter our city has been conquered by an enemy.
                    if (city_owner(cell(unit.protecting).city_id) != me()) 
                    {
                        unit.clearTarget();
                        unit.quest = CONQUER_CITY;
                        if(i <= MAX_SEARCH_QUERY) change_quest = true;
                        break;
                    }
                    // We always search for unconquered paths.
                    Pos p = closest_city_path(unit.pos, cell(unit.protecting).city_id);
                    if (p == unit.pos) break;

                    // If there is a path to conquer, we go to it.
                    unit.setTarget(p);

    
                    break;
                }
                case GOTO_MASK:     break;
                case KILL_ENEMY:    break;
                case GOTO_REQUEST:  break;
            }
            ++i;
        } while(change_quest);

        // We now check if we're stuck, being stuck means that if we move one position, chances are that we wont make it.
        // If we get stuck 

        if (unit.has_target) 
        {
            Vec2d move = astar.findBestDir(unit.pos, unit.target, mask);
            return move.toDir();
        }

        return NONE;

        
    }

  	/**
  	 * Play method, invoked once per each round.
  	 */
  	virtual void play () 
  	{
		round_special_tasks();

        // En cada ronda revisamos 

		vector<int> id_units = my_units(me());
		mask = vector<vector<bool>>(rows(), vector<bool>(rows(), 0));
        // We add our unit's positions to our mask.
        for (auto& id_unit : id_units)
        {
            mask[unit(id_unit).pos.i][unit(id_unit).pos.j] = 1;
        }

		for (auto& id_unit : id_units)
		{
			Unit my_unit = unit(id_unit);
			// La unidad mirará sus alrededores para celdas != GRASS, para así decidir que hacer.
		    //Pos nearest = kdmap.nearest(my_unit.pos);
            move(id_unit, best_move(my_unit));

            //move(id_unit, RIGHT);
			//// Ahora tenemos que hacer algo al respecto, sabemos si una celda tiene virus o si hay un enemigo.
			//Cell my_cell = cell(my_unit.pos.i, my_unit.pos.j);
			
		}
  	}

};



/**
 * Do not modify the following line.
 */
RegisterPlayer(PLAYER_NAME);


