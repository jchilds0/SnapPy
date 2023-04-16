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

typedef struct edgenode {
    int y;
    struct edgenode *next;
} edgenode;

typedef struct cuspnode {
    struct CuspTriangle *tri;
    int tetIndex;
    int tetVertex;
    int dist[3];
    int adjTri[3];
} cuspnode;

typedef struct graph {
    struct edgenode **edges;
    struct cuspnode **vertexData;
    int *degree;
    int nvertices;
    int nedges;
    int directed;
} graph;

/**
 * Cusp Triangulation
 */

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
    struct CuspTriangle *neighbours[4];
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

// Graph
void init_graph(struct graph *g, int maxVertices, bool directed);
void free_graph(struct graph *);
int insert_edge(struct graph *, int, int, bool);
int edge_exists(struct graph *, int, int);

// Dual Graph
void init_cusp_triangulation(Triangulation *, struct CuspTriangle **);
void cusp_vertex_index(struct CuspTriangle **);
void walk_around_vertex(struct CuspTriangle **, struct CuspTriangle *, int, int);
void free_cusp_triangulation(struct CuspTriangle **);
int **get_symplectic_equations(Triangulation *, struct CuspTriangle **, int, int **);
void construct_dual_graph(struct graph *, struct CuspTriangle **);
int insert_triangle_edge(struct graph *, int, int, struct CuspTriangle *, struct CuspTriangle *, bool);
int is_equal(struct cuspnode *, struct cuspnode *, struct CuspTriangle *, struct CuspTriangle *, int, int, int, int, int, int);
void init_vertex(struct cuspnode *, struct cuspnode *, struct CuspTriangle *, struct CuspTriangle *, int, int, int, int, int, int);
int is_center_vertex(struct cuspnode *);
int flow(struct CuspTriangle *, int);
int visited(int **, int *, int, int);
struct CuspTriangle *find_cusp_triangle(struct CuspTriangle **, int, int);
void print_debug_info(struct CuspTriangle **, struct graph *, int);
void remove_extra_edges(struct graph *);
void add_misc_edges(struct graph *);
void label_triangulation_edges(Triangulation *);

// Queue
void initialise_queue(struct Queue *, int);
struct Queue *enqueue(struct Queue *, int);
int dequeue(struct Queue *);
void resize_queue(struct Queue *);
int empty_queue(struct Queue *);
void free_queue(struct Queue *);

// Stack
void push(struct Node *, int);
int pop(struct Node *);
int is_empty(struct Node *);

/**
 * Graph for Breadth First Search
 */

void init_search(struct graph *, bool *, bool *, int *);
void bfs(struct graph *, int, bool *, bool *, int *);
void process_vertex_early(int);
void process_edge(int, int);
void process_vertex_late(int);
void find_path(int, int, int *, int *, int);