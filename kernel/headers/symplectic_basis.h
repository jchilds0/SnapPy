//
// Created by joshu on 20/03/2023.
//

#ifndef SNAPPEA_SYMPLECTIC_BASIS_H
#define SNAPPEA_SYMPLECTIC_BASIS_H

#endif //SNAPPEA_SYMPLECTIC_BASIS_H

#include "triangulation.h"
#include <stdbool.h>

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

/**
* Cusp Triangulation
*/

struct CuspVertex {
    int edgeClass;
    int edgeIndex;
    EdgeClass *edge;
    int v1;
    int v2;
};

struct CuspTriangle {
    Tetrahedron *tet;
    int tetIndex;
    int tetVertex;
    int numCurves;
    struct CuspVertex vertices[4];
    struct CuspTriangle *neighbours[4];
    int curve[4];
};

struct CuspRegion {
    struct CuspTriangle *tri;
    int tetIndex;
    int tetVertex;
    int curve[4][4];            /** looking at face, number of curves between the region and vertex */
    int adjTri[4];              /** does the region meet this edge of the cusp triangle */
    int adjRegions[4];          /** index of the adjacent regions */
    int dive[4][4];             /** can we dive along the face into this vertex */
};

struct PathEndPoint {
    int face;                       // face containg the short rectangle carrying the curve
    int vertex;                     // vertex we dive through the manifold along
    int regionIndex;
    struct CuspRegion *region;
};

struct DualCurves {
    struct EdgeNode *curves[2][2];
    struct PathEndPoint *endpoints[2][2];
};

/**
 * Graph
 */

struct EdgeNode {
    int y;
    int face;
    int insideVertex;
    struct EdgeNode *next;
    struct EdgeNode *prev;
};

struct Graph {
    struct EdgeNode **edge_list_begin;
    struct EdgeNode **edge_list_end;
    struct CuspRegion **pRegion;
    int *degree;
    int nVertices;
    int directed;
};

// Graph
struct Graph            *init_graph(int maxVertices, bool directed);
void                    free_graph(struct Graph *);
void                    reduce_graph_size(struct Graph *);
int                     insert_edge(struct Graph *, int, int, bool);
void                    delete_edge(struct Graph *, int, int, bool);
int                     edge_exists(struct Graph *, int, int);

// Symplectic Basis
int                     **get_symplectic_equations(Triangulation *, int, int, int);
struct CuspTriangle     **init_cusp_triangulation(Triangulation *);
void                    cusp_vertex_index(struct CuspTriangle **);
void                    walk_around_vertex(struct CuspTriangle **, struct CuspTriangle *, int, int);
void                    free_cusp_triangulation(struct CuspTriangle **);
struct CuspRegion       **init_cusp_region(struct CuspTriangle **);
void                    free_cusp_region(struct CuspRegion **);
void                    init_region(struct CuspRegion *, struct CuspTriangle *, int [4], int [4]);
int                     find_adj_region_index(struct CuspRegion **pCuspRegion, struct CuspRegion *region, int face);
void                    update_adj_regions(struct CuspRegion **);
struct DualCurves       **init_oscillating_curves(int);
void                    free_oscillating_curves(struct DualCurves **);
void                    find_intersection_triangle(Triangulation *);
int                     num_cusp_regions(struct CuspTriangle **);
struct Graph *          construct_dual_graph(struct CuspRegion **);
int                     flow(struct CuspTriangle *, int);
int                     find_cusp_triangle_index(struct CuspTriangle **, int, int);
void                    print_debug_info(struct CuspTriangle **, struct Graph *, struct CuspRegion **, struct DualCurves **, int);
void                    label_triangulation_edges(Triangulation *);
struct CuspRegion       **construct_dual_curves(struct CuspTriangle **, struct CuspRegion **, struct DualCurves **, int);
void                    find_path_endpoints_e0(struct Graph *, struct CuspRegion **, struct DualCurves *, struct DualCurves *, int, int);
void                    find_path_endpoints(struct Graph *, struct CuspRegion **, struct DualCurves *, int, int, int);
void                    update_path_info(struct CuspRegion **, struct DualCurves *, int);
struct CuspRegion       **update_cusp_regions(struct CuspRegion **, struct DualCurves *, int);
struct CuspRegion       *update_cusp_region_node(struct CuspRegion *, struct EdgeNode *, struct PathEndPoint *, int);
void                    update_cusp_triangle(struct CuspRegion **, struct CuspRegion *, struct EdgeNode *);
void                    update_cusp_triangle_endpoints(struct CuspRegion **, struct CuspRegion *, struct PathEndPoint *, struct EdgeNode *);
void                    copy_region(struct CuspRegion *, struct CuspRegion *);
void                    find_holonomies(struct CuspRegion **, struct DualCurves **, int, int **);
void                    find_path_holonomy(struct CuspRegion **, struct DualCurves *, int, int *);
void                    inside_vertex(struct CuspRegion *, struct EdgeNode *, int *, int *);

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
void                    find_path(int, int, int *, struct EdgeNode *);