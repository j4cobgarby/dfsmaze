#include <vector>
#include <map>

#include <raylib.h>

class Graph;

class Node {
protected:
    Graph *graph;
    
    int x;
    int y;
    bool discovered;

    int maze_prev;

    std::vector<int> neighbours;
public:
    Node() {}
    Node(Graph *graph, int x, int y, std::vector<int> neighbours);

    std::vector<int> &get_neighbours() {return neighbours;}
    void add_neighbour(int n) {neighbours.push_back(n);}
    int getx() {return x;}
    int gety() {return y;}

    bool is_discovered() {return discovered;}
    void set_discovered(bool b) {discovered = b;}

    void set_prev(int i) {maze_prev = i;}
    int get_prev() {return maze_prev;}
};

class Graph {
protected:
    std::map<int, Node> nodes;
    int next_id;
public:
    Graph();

    int add_node(Node);
    void del_node(int);

    Node &get_node(int);
    std::map<int, Node> &get_nodes() {return nodes;}
};

