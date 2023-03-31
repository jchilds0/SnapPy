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

typedef struct graph {
    struct edgenode **edges;
    int **vertexHomology;     // [tet index, tet vertex, dist. to v0, dist. to v1, dist. to v2]
    int *degree;
    int nvertices;
    int nedges;
    int directed;
} graph;

/**
 * Cusp Triangulation
 */

struct CuspVertex {
    int index;
    EdgeClass *edge;
    int v1;
    int v2;
};

struct CuspTriangle {
    Tetrahedron *tet;
    int tetVertex;
    struct CuspVertex vertices[3];
    int edgesThreeToFour[3];
    int edgesFourToThree[4];
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
void initialise_graph(struct graph *, int, bool);
void free_graph(struct graph *);
int insert_edge(struct graph *, int, int, struct CuspTriangle *, struct CuspTriangle *, bool);
void update_vertex_homology(struct graph *, int, int, struct CuspTriangle *, struct CuspTriangle *);

// Dual Graph
void init_cusp_triangulation(Triangulation *, struct CuspTriangle **);
void free_cusp_triangulation(Triangulation *, struct CuspTriangle **);
int **get_symplectic_equations(Triangulation *manifold, struct CuspTriangle **, int, int **);
void construct_dual_graph(struct graph *, Triangulation *, struct CuspTriangle **);
int minCuspDistance(struct graph *, int);
int flow(struct CuspTriangle *, int);
int visited(int **, int *, int, int);
struct CuspTriangle *findTriangle(Triangulation *, struct CuspTriangle **, int, int);
void printDebugInfo(Triangulation *, struct CuspTriangle **, struct graph *);
void remove_extra_edges(struct graph *);
void add_misc_edges(struct graph *);

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

void initialise_search(struct graph *, bool *, bool *, int *);
void bfs(struct graph *, int, bool *, bool *, int *);
void process_vertex_early(int);
void process_edge(int, int);
void process_vertex_late(int);
void find_path(int, int, int *, int *, int);