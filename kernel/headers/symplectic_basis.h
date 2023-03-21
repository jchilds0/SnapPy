//
// Created by joshu on 20/03/2023.
//

#ifndef SNAPPEA_SYMPLECTIC_BASIS_H
#define SNAPPEA_SYMPLECTIC_BASIS_H

#endif //SNAPPEA_SYMPLECTIC_BASIS_H

#include "triangulation.h"
#include <stdbool.h>

#define     MAXV    10

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

typedef struct edgenode {
    int y;
    int weight;
    struct edgenode *next;
} edgenode;

typedef struct {
    struct edgenode *edges[MAXV];
    int degree[MAXV];
    int nvertices;
    int nedges;
    int directed;
} graph;

void initialise_graph(graph *, int, int, bool);
void free_graph(graph *);
void initialise_search(graph *, bool *, bool *, int *);
void bfs(graph *, int, bool *, bool *, int *);
void process_vertex_early(int);
void process_edge(int, int);
void process_vertex_late(int);
void find_path(int, int, int *, int *, int);
void insert_edge(graph *, int, int, bool);