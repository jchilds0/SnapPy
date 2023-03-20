//
// Created by joshu on 20/03/2023.
//

#ifndef SNAPPEA_SYMPLECTIC_BASIS_H
#define SNAPPEA_SYMPLECTIC_BASIS_H

#endif //SNAPPEA_SYMPLECTIC_BASIS_H

#include "triangulation.h"
#include <stdbool.h>

struct queue {
    int head;   // First element of queue
    int tail;   // Last element of queue
    int len;
    struct Tetrahedron **array;
};

void initialise_queue(struct queue *, int);
struct queue *enqueue(struct queue *, Tetrahedron *);
Tetrahedron *dequeue(struct queue *);
void resize_queue(struct queue *);
int empty_queue(struct queue *);
void free_queue(struct queue *);

//#define MAXV   100
//
//bool processed[MAXV + 1];
//bool discovered[MAXV + 1];
//int parent[MAXV + 1];
//
//typedef struct edgenode {
//    int y;
//    int weight;
//    struct edgenode *next;
//} edgenode;
//
//typedef struct {
//    struct edgenode *edges[MAXV + 1];
//    int degree[MAXV + 1];
//    int nvertices;
//    int nedges;
//    int directed;
//} graph;
//
//void initialise_search(graph *);
//void process_vertex_early(int);
//void process_edge(int, int);
//void process_vertex_late(int);
//void find_path(int, int, int[]);
//void insert_edge(graph *, int, int, bool);