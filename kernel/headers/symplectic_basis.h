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
    int **vertexHomology;     // [tet vertex, cusp vertex, distance from cusp vertex]
    int *degree;
    int nvertices;
    int nedges;
    int directed;
} graph;

void initialise_graph(struct graph *, int, bool);
void free_graph(struct graph *);
void insert_edge(struct graph *, int, int, bool);

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
    int edges[3];
};

void init_cusp_triangulation(Triangulation *, struct CuspTriangle **);
int **get_symplectic_equations(Triangulation *manifold, struct CuspTriangle **, int, int **);
void construct_dual_graph(struct graph *, Triangulation *, struct CuspTriangle **);
int flow(struct CuspTriangle *, int, int);
void printTriangleInfo(Triangulation *, struct CuspTriangle **);
void remove_extra_edges(struct graph *);
void add_misc_edges(struct graph *);

/**
 * Queue
 */

struct queue {
    int front;      // First element of queue
    int rear;       // Last element of queue
    int len;        // num of elements
    int size;       // array size
    int *array;
};

void initialise_queue(struct queue *, int);
struct queue *enqueue(struct queue *, int);
int dequeue(struct queue *);
void resize_queue(struct queue *);
int empty_queue(struct queue *);
void free_queue(struct queue *);

/**
 * Graph for Breadth First Search
 */

void initialise_search(struct graph *, bool *, bool *, int *);
void bfs(struct graph *, int, bool *, bool *, int *);
void process_vertex_early(int);
void process_edge(int, int);
void process_vertex_late(int);
void find_path(int, int, int *, int *, int);