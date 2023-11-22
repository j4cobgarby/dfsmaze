#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stack>
#include <string>

#include <raylib.h>

#include "graph.hpp"

#define UPDATE_PERIOD 0.1
#define TILE_SPACING 5
#define LINE_THICKNESS 2
#define HALF_SIZE ((float)TILE_SPACING / 2)

/* Returns:
    0 if the A B C are colinear
    1 if A B C are clockwise
    2 if A B C are counterclockwise
*/
int orientation(int ax, int ay, int bx, int by, int cx, int cy) {
    int o = (cy-ay) * (bx-ax) - (by-ay) * (cx-ax);
    if (o == 0) return 0;
    if (o > 0) return 1;
    else return 2;
}

bool point_on_segment(int ax, int ay, int bx, int by, int px, int py) {
    if ((px - ax) * (by - ay) == (bx - ax) * (py - ay)) {
        if (std::min(ax, bx) <= px && px <= std::max(ax, bx) &&
            std::min(ay, by) <= py && py <= std::max(ay, by)) {
            return true;
        }
    }
    return false;
}

bool lines_cross(int ax, int ay, int bx, int by, int cx, int cy, int dx, int dy) {
    // return ccw(ax, ay, cx, cy, dx, dy) != ccw(bx, by, cx, cy, dx, dy) && ccw(ax, ay, bx, by, cx, cy) != ccw(ax, ay, bx, by, dx, dy);
    int o_abc = orientation(ax, ay, bx, by, cx, cy);
    int o_abd = orientation(ax, ay, bx, by, dx, dy);

    int o_cda = orientation(cx, cy, dx, dy, ax, ay);
    int o_cdb = orientation(cx, cy, dx, dy, bx, by);
    
    if (o_abc != o_abd && o_cda != o_cdb) return true; // No colinearity, simple intersection

    if (o_abc == 0 && point_on_segment(ax, ay, bx, by, cx, cy)) return true;
    if (o_abd == 0 && point_on_segment(ax, ay, bx, by, dx, dy)) return true;
    if (o_cda == 0 && point_on_segment(cx, cy, dx, dy, ax, ay)) return true;
    if (o_cdb == 0 && point_on_segment(cx, cy, dx, dy, bx, by)) return true;

    return false;
}

bool node2node_intersects(int node_a_i, int node_b_i, Graph &G) {
    Node &node_a = G.get_node(node_a_i);
    Node &node_b = G.get_node(node_b_i);
    
    bool overlaps = false;
    for (auto & [test_i, test_node] : G.get_nodes()) {
        if (test_node.is_discovered() && test_node.get_prev() >= 0) {
            if (test_i == node_a_i || test_node.get_prev() == node_a_i) continue; // Ignore simply extension of path
            Node &test_prev = G.get_node(test_node.get_prev());
            if (lines_cross(node_a.getx(), node_a.gety(), node_b.getx(), node_b.gety(), 
                test_node.getx(), test_node.gety(), test_prev.getx(), test_prev.gety())) {
                overlaps = true;
                break;
            }
        }
    }
    return overlaps;
}

int dfs_step(Graph &G, std::stack<int> &stack, bool check_intersections=true, bool just_one_step=true) {
    while (!stack.empty()) {
        int next_i = stack.top();
        stack.pop(); 

        Node &next = G.get_node(next_i);
        if (!next.is_discovered()) {
            next.set_discovered(true);
            for (int edge : next.get_neighbours()) {
                if (!check_intersections || !node2node_intersects(next_i, edge, G)) {
                    if (!G.get_node(edge).is_discovered())
                        G.get_node(edge).set_prev(next_i);
                    stack.push(edge);
                }
            }
            if (just_one_step) return next_i;
        }
    }

    return -1;
}

Graph make_grid_graph(int argc, char **argv, int *winwidth, int *winheight) {
    int rows, cols;

    // argv = a.out grid [rows] [cols]
    if (argc == 4) {
        rows = atoi(argv[2]);
        cols = atoi(argv[3]);
    } else {
        std::cout << "Usage: " << argv[0] << " grid [n rows] [n cols]" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Rows: " << rows << ", Columns: " << cols << std::endl;

    Graph G;

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            int i = G.add_node(Node(&G, c * TILE_SPACING + HALF_SIZE, r * TILE_SPACING + HALF_SIZE, {}));

            Node &node = G.get_node(i);
            if (r > 0) node.add_neighbour(i - cols);
            if (r < rows - 1) node.add_neighbour(i + cols);
            if (c > 0) node.add_neighbour(i - 1);
            if (c < cols - 1) node.add_neighbour(i + 1);

            // if (c > 1) node.add_neighbour(i - 2);
            // if (c < cols - 2) node.add_neighbour(i + 2);

            // if (c > 19) node.add_neighbour(i - 20);
            // if (c < cols - 20) node.add_neighbour(i + 20);

            // if (r > 0 && c > 0) node.add_neighbour(i - cols - 1);
            // if (r > 0 && c < cols - 1) node.add_neighbour(i - cols + 1);
            // if (r < rows - 1 && c > 0) node.add_neighbour(i + cols - 1);
            // if (r < rows - 1 && c < cols - 1) node.add_neighbour(i + cols + 1);

            // if (r % 20 < 10) {
            //     if (r > 0 && c > 2) node.add_neighbour(i - cols - 3);
            //     if (r < rows - 1 && c < cols - 3) node.add_neighbour(i + cols + 3);
            // } else {
            //     if (r > 0 && c < cols - 3) node.add_neighbour(i - cols + 3);
            //     if (r < rows - 1 && c > 2) node.add_neighbour(i + cols - 3);
            // }

            // if (r > 1 && c > 1) node.add_neighbour(i - 2 * cols - 2);

            std::shuffle(node.get_neighbours().begin(), node.get_neighbours().end(), gen);
        }
    }

    *winwidth = cols * TILE_SPACING;
    *winheight = rows * TILE_SPACING;

    return G;
}

float dist(int x1, int y1, int x2, int y2) {
    return sqrtf((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

Graph make_loose_graph(int argc, char **argv, int *winwidth, int *winheight) {
    int n_nodes, width, height;
    float threshold;

    // argv = a.out loose [nodes] [width] [height] [threshold]
    if (argc == 6) {
        n_nodes = atoi(argv[2]);
        width = atoi(argv[3]);
        height = atoi(argv[4]);
        threshold = atof(argv[5]);
    } else {
        std::cout << "Usage: " << argv[0] << " loose [nodes] [width] [height] [threshold (float)]" << std::endl;
        exit(EXIT_FAILURE);
    }

    Graph G;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::mt19937::result_type> width_dist(0, width-1);
    std::uniform_int_distribution<std::mt19937::result_type> height_dist(0, height-1);

    for (int i = 0; i < n_nodes; i++) {
        G.add_node(Node(&G, width_dist(gen), height_dist(gen), {}));
    }

    for (auto & [node_i, node] : G.get_nodes()) {
        for (auto & [to_i, to] : G.get_nodes()) {
            float dist_to = dist(node.getx(), node.gety(), to.getx(), to.gety());
            if (dist_to < threshold) {
                node.add_neighbour(to_i);
                // to.add_neighbour(node_i);
            }
        }

        std::sort(node.get_neighbours().begin(), node.get_neighbours().end(), [&](int &a, int &b){
            Node &na = G.get_node(a);
            Node &nb = G.get_node(b);
            return dist(node.getx(), node.gety(), na.getx(), na.gety()) > dist(node.getx(), node.gety(), nb.getx(), nb.gety());
        });        
    }

    *winwidth = width;
    *winheight = height;

    return G;
}

Graph make_radial_graph(int argc, char **argv, int *winwidth, int *winheight) {
    // argv = a.out radial [layers] [layer_height] [theta step] [angle between connecting nodes]

    int layers, layer_height;
    float theta_step, angle_between_connections;

    if (argc == 6) {
        layers = atoi(argv[2]);
        layer_height = atoi(argv[3]);
        theta_step = atof(argv[4]);
        angle_between_connections = atof(argv[5]);
    } else {
        exit(EXIT_FAILURE);
    }

    Graph G;
    std::random_device rd;
    std::mt19937 gen(rd());

    *winwidth = *winheight = layers * layer_height * 2 + 20; // 10px padding each side

    int middle_node = G.add_node(Node(&G, *winwidth/2, *winheight/2, {})); // Center

    for (int l = 0; l < layers; l++) {
        float layer_radius = (float)layer_height * l;
        int last_ni = -1;
        int first_in_layer = -1;
        
        for (float theta = 0; theta < 2 * 3.14159; theta += theta_step) {
            float xoff = std::sin(theta) * layer_radius;
            float yoff = std::cos(theta) * layer_radius;
            int new_ni = G.add_node(Node(&G, *winwidth/2 + (int)xoff, *winheight/2 + (int)yoff, {}));

            if (first_in_layer == -1) first_in_layer = new_ni;

            if (l == 1) G.get_node(middle_node).add_neighbour(new_ni);
            if (last_ni >= 0) {
                G.get_node(new_ni).add_neighbour(last_ni);
                G.get_node(last_ni).add_neighbour(new_ni);
            }
            if (theta + theta_step >= 2 * 3.14159) {
                G.get_node(first_in_layer).add_neighbour(new_ni);
                G.get_node(new_ni).add_neighbour(first_in_layer);
            }
            if (l > 1) {
                G.get_node(new_ni).add_neighbour(new_ni - (int)((2*3.14159)/theta_step) - 1);
                G.get_node(new_ni - (int)((2 * 3.14159)/theta_step) - 1).add_neighbour(new_ni);
            }
            // if (l > 1 && l < layers-1) G.get_node(new_ni + (int)((2*3.14159)/theta_step)).add_neighbour(new_ni);

            std::shuffle(G.get_node(new_ni).get_neighbours().begin(), G.get_node(new_ni).get_neighbours().end(), gen);

            last_ni = new_ni;
        }
    }

    return G;
}

int main(int argc, char **argv) {
    SetTraceLogLevel(LOG_ERROR);

    int winwidth, winheight;

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <graph type> [type options]" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string gtype = std::string(argv[1]);
    Graph G;

    if (gtype == "grid") {
        G = make_grid_graph(argc, argv, &winwidth, &winheight);
    } else if (gtype == "loose") {
        G = make_loose_graph(argc, argv, &winwidth, &winheight);
    } else if (gtype == "radial") {
        G = make_radial_graph(argc, argv, &winwidth, &winheight);
    }


    InitWindow(winwidth, winheight, "Depth-First Search Maze");

    std::stack<int> dfs_stack;
    dfs_stack.push(0);
    double last_update_time = GetTime();
    int dfs_head = -1;
    int nodes_found = 0;
    bool done = false;

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        // for (auto & [node_i, node] : G.get_nodes()) {
            // DrawCircle(node.getx(), node.gety(), 2, GRAY);
            // DrawPixel(node.getx(), node.gety(), WHITE);
        // }

        if (IsKeyDown(KEY_N)) {
            for (auto & [node_i, node] : G.get_nodes()) {
                for (int to_i : node.get_neighbours()) {
                    Node &to = G.get_node(to_i);
                    DrawLine(node.getx(), node.gety(), to.getx(), to.gety(), ColorAlpha(WHITE, 0.15));
                }
                DrawCircle(node.getx(), node.gety(), 2, RED);
            }
        }

        for (auto & [node_i, node] : G.get_nodes()) {
            if (node.is_discovered() && node.get_prev() >= 0) {
                Node &maze_prev = G.get_node(node.get_prev());
                if (maze_prev.is_discovered()) {
                    DrawLineEx(
                        {(float)maze_prev.getx(), (float)maze_prev.gety()},
                        {(float)node.getx(), (float)node.gety()},
                        (float)LINE_THICKNESS, ColorAlpha(WHITE, 1)
                    );
                    // DrawCircle(maze_prev.getx(), maze_prev.gety(), (float)LINE_THICKNESS/2, WHITE);
                    // DrawCircle(node.getx(), node.gety(), (float)LINE_THICKNESS/2, WHITE);
                }
            }
        }

        if (dfs_head >= 0) {
            Node &head = G.get_node(dfs_head);
            DrawCircle(head.getx(), head.gety(), 2.5, RED);

            bool first = true;
            for (auto it = head.get_neighbours().rbegin(); it != head.get_neighbours().rend(); it++) {
                int n = *it;
                Node &neighb = G.get_node(n);
                if (neighb.is_discovered()) continue;
                DrawLineEx(
                    {(float)head.getx(), (float)head.gety()},
                    {(float)neighb.getx(), (float)neighb.gety()},
                    first ? 2 : 1, ColorAlpha(first ? GREEN : WHITE, 0.4)
                );
                first = false;
            }
        }

        double now = GetTime();
        // if (!done) {
        if (!done && now - last_update_time > UPDATE_PERIOD) {
        // if (IsKeyPressed(KEY_SPACE)) {
            last_update_time = now;            
            dfs_head = dfs_step(G, dfs_stack, true, true);
            if (dfs_head != -1) nodes_found++;
            else {
                std::cout << "Finished calculating maze" << std::endl;
                done = true;
            }
        }
        
        EndDrawing();
    }

    CloseWindow();

    return EXIT_SUCCESS;
}
