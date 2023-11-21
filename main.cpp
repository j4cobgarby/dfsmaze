#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stack>
#include <string>

#include <raylib.h>

#include "graph.hpp"

#define TILE_SPACING 30
#define LINE_THICKNESS 7
#define TILE_SIZE 5
#define HALF_SIZE ((float)TILE_SPACING / 2)
#define TILE_OFFSET ((TILE_SPACING/2) - (TILE_SIZE/2))

bool ccw(int ax, int ay, int bx, int by, int cx, int cy) {
    return (cy-ay) * (bx-ax) > (by-ay) * (cx-ax);
}

bool lines_cross(int ax, int ay, int bx, int by, int cx, int cy, int dx, int dy) {
    return ccw(ax, ay, cx, cy, dx, dy) != ccw(bx, by, cx, cy, dx, dy) && ccw(ax, ay, bx, by, cx, cy) != ccw(ax, ay, bx, by, dx, dy);
}

int dfs_step(Graph &G, std::stack<int> &stack) {
    while (!stack.empty()) {
        int next_i = stack.top();
        stack.pop(); 

        Node &next = G.get_node(next_i);
        if (!next.is_discovered()) {
            next.set_discovered(true);
            for (int edge : next.get_neighbours()) {
                bool overlaps = false;

                
                for (auto & [test_i, test_node] : G.get_nodes()) {
                    if (test_node.is_discovered() && test_node.get_prev() >= 0) {
                        Node &test_prev = G.get_node(test_node.get_prev());
                        if (lines_cross(next.getx(), next.gety(), G.get_node(edge).getx(), G.get_node(edge).gety(), 
                            test_node.getx(), test_node.gety(), test_prev.getx(), test_prev.gety())) {
                            overlaps = true;
                            break;
                        }
                    }
                }

                if (!overlaps) {
                    if (!G.get_node(edge).is_discovered())
                        G.get_node(edge).set_prev(next_i);
                    stack.push(edge);
                }
            }
            return next_i;
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

            if (c > 1) node.add_neighbour(i - 2);
            if (c < cols - 2) node.add_neighbour(i + 2);

            if (c > 2) node.add_neighbour(i - 3);
            if (c < cols - 3) node.add_neighbour(i + 3);

            // if (r > 0 && c > 0) node.add_neighbour(i - cols - 1);
            // if (r > 0 && c < cols - 1) node.add_neighbour(i - cols + 1);
            // if (r < rows - 1 && c > 0) node.add_neighbour(i + cols - 1);
            // if (r < rows - 1 && c < cols - 1) node.add_neighbour(i + cols + 1);

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

int main(int argc, char **argv) {
    SetTraceLogLevel(LOG_ERROR);

    std::cout << "Perpendicular touching lines cross? " << lines_cross(0,0, 5,0, 0,0, 0,5) << std::endl;
    
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
    }


    InitWindow(winwidth, winheight, "Depth-First Search Maze");

    std::stack<int> dfs_stack;
    dfs_stack.push(0);
    double last_update_time = GetTime();
    int dfs_head = -1;
    int nodes_found = 0;

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        // for (auto & [node_i, node] : G.get_nodes()) {
        //     DrawCircle(node.getx(), node.gety(), 2, GRAY);
        // }

        for (auto & [node_i, node] : G.get_nodes()) {
            auto col = node.is_discovered() ? WHITE : GRAY;
            if (node_i == dfs_head) col = RED;

            if (node.is_discovered() && node.get_prev() >= 0) {
                Node &maze_prev = G.get_node(node.get_prev());
                if (maze_prev.is_discovered()) {
                    DrawLineEx(
                        {(float)maze_prev.getx(), (float)maze_prev.gety()},
                        {(float)node.getx(), (float)node.gety()},
                        (float)LINE_THICKNESS, ColorAlpha(WHITE, 1)
                    );
                    DrawCircle(maze_prev.getx(), maze_prev.gety(), (float)LINE_THICKNESS/2, WHITE);
                    DrawCircle(node.getx(), node.gety(), (float)LINE_THICKNESS/2, WHITE);
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
        // if (true) {
        if (now - last_update_time > 0.005) {
        // if (IsKeyPressed(KEY_SPACE)) {
            last_update_time = now;            
            dfs_head = dfs_step(G, dfs_stack);
            if (dfs_head != -1) nodes_found++;
        }
        
        EndDrawing();
    }

    CloseWindow();

    return EXIT_SUCCESS;
}
