#include "graph.hpp"

#include <algorithm>

Node::Node(Graph *graph, int x, int y, std::vector<int> neighbours) :
    x(x), y(y), graph(graph), neighbours(neighbours) {
    discovered = false;
    maze_prev = -1;
}

Graph::Graph() {
    next_id = 0;
}

int Graph::add_node(Node n) {
    nodes[next_id] = n;
    return next_id++;
}

void Graph::del_node(int to_delete) {
    for (auto & [node_i, node] : nodes) {
        std::vector<int> &neighbs = node.get_neighbours();
        auto f = std::find(neighbs.begin(), neighbs.end(), to_delete);
        if (f != neighbs.end()) {
            neighbs.erase(f);
        }
    }
}

Node &Graph::get_node(int node_i) {
    return nodes[node_i];
}
