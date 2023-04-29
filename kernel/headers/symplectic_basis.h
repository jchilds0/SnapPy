//
// Created by joshu on 20/03/2023.
//

#ifndef SNAPPEA_SYMPLECTIC_BASIS_H
#define SNAPPEA_SYMPLECTIC_BASIS_H

#endif //SNAPPEA_SYMPLECTIC_BASIS_H

#include "triangulation.h"
#include <stdbool.h>

/**
 * Graph
 */

struct PathEndPoint {
    int face;                       // face containg the short rectangle carrying the curve
    int vertex;                     // vertex we dive through the manifold along
    int graphVertex;                // vertex of the graph
    struct CuspRegion *region;
};

struct EdgeNode {
    int y;
    struct EdgeNode *next;
    struct EdgeNode *prev;
};

struct Graph {
    struct EdgeNode **edge_list_begin;
    struct EdgeNode **edge_list_end;
    int *pCuspRegionIndex;
    int *degree;
    int nvertices;
    int nedges;
    int directed;
};

/**
 * Cusp Triangulation
 */

struct CuspRegion {
    struct CuspTriangle *tri;
    int tetIndex;
    int tetVertex;
    int dist[3];        // Distance to cuspVertex
    int adjTri[3];      // Indicates the cusp triangle sides that can be reached
    int adjNodes[3];    // Index of adjacent triangles
};

struct CuspFace {
    int index;
};

struct CuspVertex {
    int edgeIndex;
    int vertexIndex;
    EdgeClass *edge;
    int v1;
    int v2;
};

struct CuspTriangle {
    Tetrahedron *tet;
    int tetVertex;
    struct CuspVertex vertices[3];
    struct CuspFace faces[3];
    struct CuspTriangle *neighbours[4];
    int orientVertices[4][4];
    bool oriented;
};

/**
 * Queue
 */

struct Queue {
    int front;      // First element of queue
    int rear;       // Last element of queue
    int len;        // num of elements
    int size;       // array size
    int *array;
};

/**
 * Stack
 */

struct Node {
    int item;
    struct Node *next;
};

extern void             find_standard_basepoint(Triangulation *manifold, Cusp *cusp);

// Graph
struct Graph            *init_graph(int maxVertices, bool directed);
void                    free_graph(struct Graph *);
int                     insert_edge(struct Graph *, int, int, bool);
void                    delete_edge(struct Graph *, int, int, bool);
int                     edge_exists(struct Graph *, int, int);

// Graph Splitting
struct Graph            *split_along_path(struct Graph *, struct CuspRegion **, int *, int);
void                    init_vertices(struct Graph *, struct Graph *, int *, int);
void                    add_non_path_edges(struct Graph *, struct Graph *, int *, int);
void                    add_path_edges(struct Graph *, struct Graph *, struct CuspRegion **, int *, int);
bool                    inclusion(int *, int, int);

// Dual Graph
struct CuspTriangle     **init_cusp_triangulation(Triangulation *);
void                    cusp_vertex_index(struct CuspTriangle **);
void                    vertex_orientation(struct CuspTriangle **);
void                    walk_around_vertex(struct CuspTriangle **, struct CuspTriangle *, int, int);
void                    free_cusp_triangulation(struct CuspTriangle **);
int                     **get_symplectic_equations(Triangulation *, int, int, int);
void                    find_intersection_triangle(struct CuspTriangle **pTriangle);
int                     num_cusp_regions(Triangulation *, struct CuspTriangle **);
void                    construct_dual_graph(Triangulation *, struct Graph *, struct CuspTriangle **, struct CuspRegion **);
void                    init_zero_vertex(struct CuspRegion *, struct CuspTriangle *);
int                     insert_triangle_edge(struct Graph *, int, int, struct CuspTriangle *, struct CuspTriangle *,
        struct CuspRegion **);
int                     is_equal(struct CuspRegion *, struct CuspRegion *, struct CuspTriangle *, struct CuspTriangle *,
                                 int, int, int, int, int, int);
void                    init_vertex(struct CuspRegion *, struct CuspRegion *, struct CuspTriangle *, struct CuspTriangle *,
                                    int, int, int, int, int, int);
int                     is_center_vertex(struct CuspRegion *);
int                     flow(struct CuspTriangle *, int);
struct CuspTriangle     *find_cusp_triangle(struct CuspTriangle **, int, int);
int                     find_cusp_triangle_index(struct CuspTriangle **, int, int);
void                    print_debug_info(struct CuspTriangle **, struct Graph *, struct CuspRegion **, int **, int *, int);
void                    label_triangulation_edges(Triangulation *);
void                    update_orientation(struct CuspTriangle *, struct CuspTriangle *, int);
struct Graph            *construct_dual_curves(struct Graph *, struct CuspTriangle **, struct CuspRegion **,
        int, struct PathEndPoint **, int **, int *);
void                    find_index(struct Graph *, struct CuspRegion **, int, struct PathEndPoint *, struct PathEndPoint *);
void                    find_holonomies(struct Graph *, struct CuspRegion **, struct PathEndPoint **, int **, int **, int *, int);
void                    find_path_holonomy(struct Graph *, struct CuspRegion **, struct PathEndPoint *, struct PathEndPoint *, int *, int *, int);
void                    inside_vertex(struct CuspRegion *, int, int, int, int *, int *);

// Queue
void                    initialise_queue(struct Queue *, int);
struct Queue            *enqueue(struct Queue *, int);
int                     dequeue(struct Queue *);
void                    resize_queue(struct Queue *);
int                     empty_queue(struct Queue *);
void                    free_queue(struct Queue *);

// Stack
void                    push(struct Node *, int);
int                     pop(struct Node *);
int                     is_empty(struct Node *);

/**
 * Graph for Breadth First Search
 */

void                    init_search(struct Graph *, bool *, bool *, int *);
void                    bfs(struct Graph *, int, bool *, bool *, int *);
void                    process_vertex_early(int);
void                    process_edge(int, int);
void                    process_vertex_late(int);
void                    find_path(int, int, int *, int *, int, int *);