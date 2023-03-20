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

int** get_symplectic_basis(Triangulation *manifold, int* dual_rows)
{
    int i, j, T;
    int genus = 1;
    int** eqns;

    T = manifold -> num_tetrahedra;

    *dual_rows = 2;

    eqns = NEW_ARRAY(*dual_rows, int*);

    for (i = 0; i < *dual_rows; i ++) {
        eqns[i] = NEW_ARRAY(3*T, int);

        for (j = 0; j < 3 * T; j ++)
            eqns[i][j] = i + j;
    }

    return eqns;
}

// -----------------------------------------------------------

// BFS from Skiena Algorithm Design Manual

// Queue

void initialise_queue(struct queue *q, int size) {
    q->head = -1;
    q->tail = 0;
    q->len = size;

    q->array = malloc( sizeof(Tetrahedron *) * size);
}

struct queue *enqueue(struct queue *q, Tetrahedron *tet) {
    // Queue is full
    if ((q->head == 0 && q->tail == q->len - 1) || (q->head == q->tail - 1 && q->head != -1)) {
        resize_queue(q);
        enqueue(q, tet);
        return q;
    }

    if (q->head == q->len - 1) {
        q->head = 0;
    } else {
        q->head++;
    }

    q->array[q->head] = tet;

    return q;
}

Tetrahedron *dequeue(struct queue *q) {
    Tetrahedron *tet = q->array[q->tail];

    if (q->tail == q->len - 1) {
        q->tail = 0;
    } else {
        q->tail++;
    }

    return tet;
}

int empty_queue(struct queue *q) {
    return (q->head == q->len - 1 && q->tail == 0) || (q->head == q->tail - 1);
}

void resize_queue(struct queue *q) {
    Tetrahedron *tet;
    struct queue p;

    initialise_queue(&p, 2 * q->len);

    // Copy elements to new array
    while (!empty_queue(q)) {
        tet = dequeue(q);
        enqueue(&p, tet);
    }

    // Move p queue to q
    q->tail = p.tail;
    q->head = p.head;
    q->len = p.len;
    q->array = p.array;

    free_queue(&p);
}

void free_queue(struct queue *q) {
    free(q->array);
}

// Breadth First Search

//void insert_edge(graph *g, int x, int y, bool directed) {
//    edgenode *p;
//
//    p = malloc(sizeof( edgenode));
//    p->weight = 0;
//    p->y = y;
//    p->next = g->edges[x];
//
//    g->edges[x] = p;
//
//    g->degree[x]++;
//
//    if (!directed) {
//        insert_edge(g, y, x, true);
//    }
//}
//
//void initialise_search(graph *g) {
//    int i;
//
//    for (i = 0; i <= g->nvertices; i ++) {
//        processed[i] = false;
//        discovered[i] = false;
//        parent[i] = -1;
//    }
//}
//
//void bfs(graph *g, Triangulation *start) {
//    struct queue q;
//    int v, y;
//    edgenode *p;
//
//    initialise_queue(&q, 10);
//    enqueue(&q, (Tetrahedron *) start);
//    discovered[start] = true;
//
//    while (!empty_queue(&q)) {
//        v = dequeue(&q);
//        process_vertex_early(v);
//        processed[v] = true;
//        p = g->edges[v];
//        while (p != NULL) {
//            y = p->y;
//            if ((!processed[y]) || g->directed) {
//                process_edge(v, y);
//            }
//            if (!discovered[y]) {
//                enqueue(&q, y);
//                discovered[y] = true;
//                parent[y] = v;
//            }
//            p = p->next;
//        }
//
//        process_vertex_late(v);
//    }
//
//    free_queue(&q);
//}
//
//void process_vertex_early(int v) {
//
//}
//
//void process_edge(int x, int y) {
//    printf("processed edge (%d, %d)\n", x, y);
//}
//
//void process_vertex_late(int v) {
//
//}
//
//void find_path(int start, int end, int parents[]) {
//    if ((start == end) || (end == -1)) {
//        printf("\n%d", start);
//    } else {
//        find_path(start, parents[end], parents);
//        printf(" %d", end);
//    }
//}