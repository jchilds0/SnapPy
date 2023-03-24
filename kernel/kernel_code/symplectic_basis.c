/*
 *  Symplectic Basis
 *
 *  Computes the symplectic basis of a cusped 3-manifold
*/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "kernel.h"
#include "kernel_namespace.h"
#include "symplectic_basis.h"
#include "../addl_code/addl_code.h"

int** get_symplectic_basis(Triangulation *manifold) {
    int i, **c_eqns, *num_rows, *num_cols, dual_rows;

    // Get Gluing Equations
    c_eqns = get_gluing_equations(manifold, num_rows, num_cols);
    dual_rows = *num_rows - 1;

    // Allocate Equations and Cusp Arrays

    int **eqns = NEW_ARRAY(dual_rows, int*);

    for (i = 0; i < dual_rows; i ++)
        eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

    struct Triangle **pTriangle = NEW_ARRAY(4 * manifold->num_tetrahedra, struct Triangle *);

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++)
        pTriangle[i] = malloc(sizeof( struct Triangle ));

    // Algorithm

    init_cusp_triangulation(manifold, pTriangle);

    eqns = construct_equations(manifold, dual_rows, eqns);

    // Free Cusp Triangulation Array

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++)
        free(pTriangle[i]);

    my_free(pTriangle);

    return eqns;
}

void init_cusp_triangulation(Triangulation *manifold, struct Triangle **pTriangle) {
    int i, j;
    Tetrahedron *tet= manifold->tet_list_begin.next;

    for (i = 0; i < manifold->num_tetrahedra; i++) {
        pTriangle[i]->tet = tet;
        pTriangle[i]->vertex = i % 4;

        switch (pTriangle[i]->vertex) {
            case 0:
                pTriangle[i]->edges[0] = 1;
                pTriangle[i]->edges[1] = 2;
                pTriangle[i]->edges[2] = 3;
                break;
            case 1:
                pTriangle[i]->edges[0] = 0;
                pTriangle[i]->edges[1] = 2;
                pTriangle[i]->edges[2] = 3;
                break;
            case 2:
                pTriangle[i]->edges[0] = 0;
                pTriangle[i]->edges[1] = 1;
                pTriangle[i]->edges[2] = 3;
                break;
            case 3:
                pTriangle[i]->edges[0] = 0;
                pTriangle[i]->edges[1] = 1;
                pTriangle[i]->edges[2] = 2;
                break;
            default:
                break;
        }

        for (j = 0; j < 3; j++) {
            // Edge between pTriangle[i]->vertex and pTriangle[i]->edges[j]
            pTriangle[i]->vertices[j].v1 = pTriangle[i]->vertex;
            pTriangle[i]->vertices[j].v2 = pTriangle[i]->edges[j];

            pTriangle[i]->vertices[j].edge = tet->edge_class[
                    edge_between_vertices[
                            pTriangle[i]->vertices[j].v1][
                            pTriangle[i]->vertices[j].v2]];

            pTriangle[i]->vertices[j].index = (int) pTriangle[i]->vertices[j].edge->incident_edge_index;
        }

        tet = tet->next;
    }
}

int **construct_equations(Triangulation *manifold, int dual_rows, int **eqns) {
    int i, j, T = manifold -> num_tetrahedra;

    // Dual Curve Equations
    for (i = 0; i < dual_rows; i ++) {
        for (j = 0; j < 3 * T; j ++)
            eqns[i][j] = i + j;
    }

    return eqns;
}

void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

// -----------------------------------------------------------

// BFS from Skiena Algorithm Design Manual

// Queue

void initialise_queue(struct queue *q, int size) {
    q->front = 0;
    q->rear = -1;
    q->len = 0;
    q->size = size;
    q->array = NEW_ARRAY(size, int);
}

struct queue *enqueue(struct queue *q, int i) {
    // Queue is full
    if ( q->size == q->len ) {
        resize_queue(q);
        q = enqueue(q, i);
    } else {
        q->rear = (q->rear + 1) % (q->size - 1);
        q->array[q->rear] = i;
        q->len++;
    }

    return q;
}

int dequeue(struct queue *q) {
    // User to verify queue is not empty
    int i = q->array[q->front];

    q->front = (q->front + 1) % (q->size - 1);
    q->len--;

    return i;
}

int empty_queue(struct queue *q) {
    return (!q->len);
}

void resize_queue(struct queue *q) {
    int i;
    struct queue p;

    initialise_queue(&p, 2 * q->size);

    // Copy elements to new array
    while (!empty_queue(q)) {
        i = dequeue(q);
        enqueue(&p, i);
    }

    free_queue(q);

    // Move p queue to q
    q->front = p.front;
    q->rear = p.rear;
    q->len = p.len;
    q->size = p.size;
    q->array = p.array;
}

void free_queue(struct queue *q) {
    free(q->array);
}

// Breadth First Search

void initialise_graph(graph *g, int maxVertices, int maxEdges, bool directed) {
    int i;

    g->nvertices = maxVertices;
    g->nedges = 0;
    g->directed = directed;

    for (i = 0; i < MAXV; i++) {
        g->degree[i] = 0;
    }
    for (i = 0; i < MAXV; i++) {
        g->edges[i] = NULL;
    }
//    g->edges = NEW_ARRAY(maxEdges, edgenode *);
//    g->degree = NEW_ARRAY(maxVertices, int);
}

void free_graph(graph *g) {
//    free(g->edges);
//    free(g->degree);
}

void insert_edge(graph *g, int x, int y, bool directed) {
    edgenode *p;

    p = malloc(sizeof( edgenode ));
    p->weight = 0;
    p->y = y;
    p->next = g->edges[x];

    g->edges[x] = p;

    g->degree[x]++;

    if (!directed) {
        insert_edge(g, y, x, true);
    }
}

//bool *processed = NEW_ARRAY(nvertices, bool);
//bool *discovered = NEW_ARRAY(nvertices, bool);
//int *parent = NEW_ARRAY(nvertices, int);

void initialise_search(graph *g, bool *processed, bool *discovered, int *parent) {
    int i;

    for (i = 0; i <= g->nvertices; i ++) {
        processed[i] = false;
        discovered[i] = false;
        parent[i] = -1;
    }
}

void bfs(graph *g, int start, bool *processed, bool *discovered, int *parent) {
    struct queue q;
    int v, y;
    edgenode *p;

    initialise_queue(&q, 10);
    enqueue(&q, start);
    discovered[start] = true;

    while (!empty_queue(&q)) {
        v = dequeue(&q);
        process_vertex_early(v);
        processed[v] = true;
        p = g->edges[v];
        while (p != NULL) {
            y = p->y;
            if ((!processed[y]) || g->directed) {
                process_edge(v, y);
            }
            if (!discovered[y]) {
                enqueue(&q, y);
                discovered[y] = true;
                parent[y] = v;
            }
            p = p->next;
        }

        process_vertex_late(v);
    }

    free_queue(&q);
}

void process_vertex_early(int v) {

}

void process_edge(int x, int y) {
    printf("    Processed edge (%d, %d)\n", x, y);
}

void process_vertex_late(int v) {

}

void find_path(int start, int end, int *parents, int *path, int index) {
    if ((start == end) || (end == -1)) {
        path[index] = start;
    } else {
        find_path(start, parents[end], parents, path, index + 1);
        path[index] = end;
    }
}