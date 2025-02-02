#include "Player.hh"

#include <iostream>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <memory>
#include <stack>
#include <map>
#include <queue>

using namespace std;

enum COORD { x, y, z };

int sign(int num)
{
    if (num == 0) return 0;
    if (num > 0) return 1;
    if (num < 0) return -1;
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
#define PLAYER_NAME THELEGEND270

struct PLAYER_NAME : public Player {
	kdtree kdcity;
	kdtree kdpaths;
    kdtree kdalies;
    AStar  astar;
    vector<vector<bool>> mask;

  	/**
  	 * Do not modify this function.
  	 */
  	static Player* factory () { return new PLAYER_NAME; };
	
	void round_special_tasks()
	{
		switch (round())
		{
			// En la ronda 1 definiremos el kdtree que contiene todos los puntos de interes.
			case 0:
				// Construimos el kdtree de ciudades.
				//vector<Pos> city_points;
                vector<Pos> city_center_points;
				size_t n_cities = nb_cities();
				for (size_t i = 0; i < n_cities; ++i)
				{
					City _city = city(i);
					vector<Pos> hull = convexHull(_city);
                    auto miny = find_min_y(hull);
                    auto maxy = find_max_y(hull);
                    city_center_points.push_back(Pos((maxy->i - miny->i) >> 1,(maxy->j - miny->j) >> 1) + *miny);
				}

				vector<Pos> path_points;
				// Ahora inertamos los paths en nuestro set de puntos y su kd.
				size_t n_paths = nb_paths();
				for (size_t i = 0; i < n_paths; ++i)
				{
					Path _path = path(i);
					for (auto& pos : _path.second) path_points.push_back(pos);
				}
				kdcity = kdtree(city_center_points.begin(), city_center_points.end());
				kdpaths = kdtree(path_points.begin(), path_points.end());

                // Ahora barrimos el mapa para encontrar los muros.
                vector<vector<bool>> walls(rows(), vector<bool>(cols(), 0));
                for (int i = 0; i < walls.size(); ++i)
                    for (int j = 0; j < walls[0].size(); ++j)
                        if (cell(i, j).type == WALL) walls[i][j] = 1;

                astar = AStar(walls);

                break;

		}
	}

    enum Quest { GOTO_ALIE, GOTO_CLOSEST, GOTO_CITY, GOTO_UNCONQUERED_CITY, GOTO_PATH, GOTO_MASK, KILL_ENEMY, INFECT_ENEMIES, GOTO_REQUEST };
    enum Params { AVOID_VIRUS, AVOID_ENEMIES, PROTECT };
    struct Request {};
    struct Response {};
    struct Memory
    {
        Pos              last_pos;
        Quest            quest;
        vector<Pos>      next_steps;
        vector<Params>   args;
        Pos              target;
        int              path;
        int              path_count;
        bool             has_target;
        vector<Request>  requests;
        vector<Response> responses;

        Memory()
        : last_pos()
        , quest(GOTO_CLOSEST)
        , next_steps()
        , args()
        , target()
        , has_target(false)
        , requests()
        , responses()
        {};

        void setTarget(const Pos& pos) { target = pos, has_target = true; }

        void clearTarget() {has_target = false; }
    };

    map<int, Memory> units;

    void scan_surrounding(const Pos& pos, vector<Pos>& alies, vector<Pos>& enemies, vector<Pos>& infected, vector<Pos>& cities, vector<Pos>& paths, const int RADIUS)
    {
        // We scan in a diamond-like way, because enemies cannot move diagonally, doesnt make sense to check further.
        for (int i = -1*RADIUS; i <= RADIUS; ++i)
        {
            for (int j = -1*RADIUS; j <=  RADIUS; ++j)
            {
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


    Dir best_move(const Unit& unit_)
    {
        #define ATTACK_HEALTH 10

        // We first see if the unit is in the system.
        auto unit = units.find(unit_.id);
        if (unit == units.end()) units[unit_.id] = Memory(), unit = units.find(unit_.id);

        Cell cell_ = cell(unit_.pos);

        // HAY QUE AÑADIR LOGICA DE PROTECCIÓN CONTRA VIRUS

        // We first scan our surrounding at distance 1. This procedure must be done always, it is independent from other actions.
        vector<Pos> alies, enemies, infected, cities, paths;
        scan_surrounding(unit_.pos, alies, enemies, infected, cities, paths, 1);
        // We create the path mask, later we'll update the mechanism.
        for (int i = 0; i < alies.size()  ; ++i)  mask[alies[i].i][alies[i].j] = 1;
        for (int i = 0; i < enemies.size(); ++i)  mask[enemies[i].i][enemies[i].j] = 1;
        //for (int i = 0; i < infected.size(); ++i) mask[infected[i].i][infected[i].j] = 1;
        // For now, we'll not avoid infected cells. If enemies are a distance 1, 
        // we'll attack if our health is over 40, when the enemy cell has lower health 
        // than ours or if we have an alie closer than 2 units manhattan distance.
        // Si tenemos enemigos a 1 unidad de distancia.
        if (enemies.size() > 0)
        {
            if (unit_.health < ATTACK_HEALTH && alies.size() == 0)
                return avoid_enemies(unit_.pos, enemies).toDir();
            if (unit_.health > ATTACK_HEALTH && alies.size() > 0)
                return Vec2d(unit_.pos, enemies[0]).toDir();
            return avoid_enemies(unit_.pos, enemies).toDir();
        }
        // We check distance 2 just in case.
        scan_surrounding_fixed(unit_.pos, alies, enemies, infected, cities, paths, 2);


        // If no requests, we follow normal processing.
        //if (GOTO_CLOSEST &&  cell_.type == GRASS) unit->second.quest = GOTO_CLOSEST;
        Pos nearest_city, nearest_path;
        vector<Pos> nearest_city_vec, nearest_path_vec;

        switch (unit->second.quest)
        {
            case GOTO_ALIE:     break;
            case GOTO_CLOSEST:  
                // The closest structure is either a city or a path. 
                // Once we get there, if conquered, we move on to the closest conquerable structure.
                if (cell_.type == GRASS)
                {
                    nearest_city = kdcity.nearest(unit_.pos);
                    nearest_path = kdpaths.nearest(unit_.pos);
                    unit->second.setTarget(distance(unit_.pos, nearest_city) < distance(unit_.pos, nearest_path) 
                                          ? cities.size() == 0 ? kdcity.nearest(unit_.pos) : cities[0]
                                          : paths.size() == 0 ? kdpaths.nearest(unit_.pos) : paths[0]);
                    break;
                }
                // We update the quest, if we've exhausted our quest targets.
                if (cell_.type == CITY || cell_.type == PATH) unit->second.clearTarget();
                // We search our closest objective and we update accordingly.
                if (cell_.type == CITY) unit->second.quest = GOTO_PATH;
                if (cell_.type == PATH) unit->second.quest = GOTO_CITY;
                break;
            case GOTO_CITY:     
                // When we are in this state, we'd like to go to an unconquered city first, if possible. Enemies are avoided if necessary.
                if (!unit->second.has_target)
                {
                    nearest_city_vec = kdcity.nearest(unit_.pos, 2);
                    if (nearest_city_vec.size() == 1)
                        unit->second.setTarget(cities.size() == 0 ? nearest_city_vec[1] : cities[0]);
                    if (nearest_city_vec.size() >= 2)
                    {
                        // If we have two cities, we must search for the one that is closest and not conquered.
                        int closest_owner  = city_owner(cell(nearest_city_vec[0]).city_id);
                        int furthest_owner = city_owner(cell(nearest_city_vec[1]).city_id);
                        if (closest_owner == -1)    { unit->second.setTarget(cities.size() == 0 ? nearest_city_vec[0] : cities[0]); break; }
                        if (furthest_owner == -1)   { unit->second.setTarget(nearest_city_vec[1]); break; }
                        if (closest_owner != me())  { unit->second.setTarget(cities.size() == 0 ? nearest_city_vec[0] : cities[0]); break; }
                        if (furthest_owner != me()) { unit->second.setTarget(nearest_city_vec[0]); break; }
                        // If both cities are conquered by me, we go to the closest one.
                        unit->second.setTarget(nearest_city_vec[0]);
                        break;
                    }
                }
                // In this quest, we accomplish our target when we get to a city, doesnt matter if it is
                // our target's one, since we are searching for the 2 closest, they could be collinear.
                // When we get to our target, we clear it.
                if (cell(unit_.pos).type == CITY)
                { 
                    unit->second.clearTarget();
                    // ** AÑADIR MAS LOGICA, ESTADO DE PROTECCION Y DEMÁS.
                    unit->second.quest = GOTO_PATH;
                }

                // If we've gotten to a city, we update accordingly.

                break;
            case GOTO_PATH:     
                // In this state we do de same as with the latter one. 
                // We search for the closest "two paths". We'll search for a total of 15 points, 
                // then we search the ones that have different id's.
                if (!unit->second.has_target)
                {
                    vector<Pos> points = kdpaths.nearest(unit_.pos, 4);

                    if (points.size() == 0) break; 
                    vector<Pos> extremes;
                    int first_id = cell(points[0]).path_id;
                    extremes.push_back(points[0]);
                    for (int i = 1; i < points.size(); ++i) 
                        if (cell(points[i]).path_id != first_id && cell(points[i]).path_id != -1) 
                        {
                            extremes.push_back(points[i]);
                            break;
                        }

                    // Now we have two points that meet our requirements. We search for the closest one.
                    Pos closest = distance(unit_.pos, extremes[0]) < distance(unit_.pos, extremes[1])
                                  ? extremes[0] : extremes[1];
                    Pos furthest = distance(unit_.pos, extremes[0]) > distance(unit_.pos, extremes[1])
                                  ? extremes[0] : extremes[1];
                    
                    int closest_owner  = me();
                    int furthest_owner = me();
                    if (cell(closest).path_id != -1) closest_owner = path_owner(cell(closest).path_id);
                    if (cell(furthest).path_id != -1) furthest_owner = path_owner(cell(furthest).path_id);


                    if (closest_owner == -1)   { unit->second.setTarget(closest); break; }
                    if (furthest_owner == -1)  { unit->second.setTarget(furthest); break;}
                    if (closest_owner != me()) { unit->second.setTarget(closest); break; }
                    if (furthest_owner != me()){ unit->second.setTarget(furthest); break; }
                    break;
                }
                // We complete our quest when we get to a path, it doesnt matter if its our target
                if (cell(unit_.pos).type == PATH) 
                { 
                    unit->second.clearTarget();
                    // ** AÑADIR MAS LOGICA, ESTADO DE PROTECCION Y DEMÁS.
                    unit->second.quest = GOTO_CITY;
                }
                break;
            case GOTO_MASK:     break;
            case KILL_ENEMY:    break;
            case GOTO_REQUEST:  break;
        }
        Vec2d move = astar.findBestDir(unit_.pos, unit->second.target, mask);


        // Now we scan the surroinding of the next move. We'd like to avoid enemies if possible.
        Pos next = unit_.pos + move;
        vector<Pos> alies_, enemies_, infected_;
        vector<Pos> fill;
        scan_surrounding(next, alies_, enemies_, infected_, fill, fill, 1);
        if (enemies_.size() > 0) return avoid_enemies(next, enemies_).toDir();

        return move.toDir();
    }

  	/**
  	 * Play method, invoked once per each round.
  	 */
  	virtual void play () 
  	{
		vector<int> id_units = my_units(me());
		mask = vector<vector<bool>>(rows(), vector<bool>(rows(), 0));
		round_special_tasks();

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


