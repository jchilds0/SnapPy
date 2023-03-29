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

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i;

    // Edge Curves C_i -> gluing equations
    int **edge_eqns, edge_num_rows;

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;

    // Get Gluing Equations
    edge_eqns = get_gluing_equations(manifold, &edge_num_rows, num_cols);
    symp_num_rows = edge_num_rows;

    // Allocate Symplectic Equations Array
    symp_eqns = NEW_ARRAY(symp_num_rows, int*);

    for (i = 0; i < symp_num_rows; i ++)
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

    struct CuspTriangle **pTriangle = NEW_ARRAY(4 * manifold->num_tetrahedra, struct Triangle *);

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++)
        pTriangle[i] = malloc(sizeof( struct CuspTriangle ));

    // Get Symplectic Equation
    symp_eqns = get_symplectic_equations(manifold, pTriangle, symp_num_rows, symp_eqns);

    // Free Cusp Triangulation Array
    for (i = 0; i < 4 * manifold->num_tetrahedra; i++)
        free(pTriangle[i]);

    my_free(pTriangle);

    // Construct return array
    int **eqns = NEW_ARRAY(edge_num_rows + symp_num_rows, int *);

    for (i = 0; i < edge_num_rows; i++) {
        eqns[2 * i] = edge_eqns[i];
        eqns[2 * i + 1] = symp_eqns[i];
    }

    *num_rows = 2 * edge_num_rows;
    return eqns;
}

void init_cusp_triangulation(Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int i, j, k;
    EdgeClass *edge;
    Tetrahedron *tet= manifold->tet_list_begin.next;

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++) {
        pTriangle[i]->tet = tet;
        pTriangle[i]->tetVertex = i % 4;

        switch (pTriangle[i]->tetVertex) {
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
            pTriangle[i]->vertices[j].v1 = pTriangle[i]->tetVertex;
            pTriangle[i]->vertices[j].v2 = pTriangle[i]->edges[j];

            pTriangle[i]->vertices[j].edge = tet->edge_class[
                    edge_between_vertices[
                            pTriangle[i]->vertices[j].v1][
                            pTriangle[i]->vertices[j].v2]];

            k = 0;
            for (edge = manifold->edge_list_begin.next; edge != pTriangle[i]->vertices[j].edge; edge = edge->next, k++);

            pTriangle[i]->vertices[j].index = k;
        }

        if (i % 4 == 3) {
            tet = tet->next;
        }
    }
}

int **get_symplectic_equations(Triangulation *manifold, struct CuspTriangle **pTriangle, int num_rows, int **eqns) {
    int i, j, T = manifold -> num_tetrahedra;
    int genus = 2 * manifold->num_tetrahedra - num_rows;
    struct graph graph1;

    init_cusp_triangulation(manifold, pTriangle);
    initialise_graph(&graph1, 4 * (genus + 2) * manifold->num_tetrahedra, FALSE);
    construct_dual_graph(&graph1, manifold, pTriangle);
    printTriangleInfo(manifold, pTriangle);

    // Dual Curve Equations
    for (i = 0; i < num_rows; i ++) {
        for (j = 0; j < 3 * T; j ++)
            eqns[i][j] = i + j;
    }

    free_graph(&graph1);

    return eqns;
}

void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

/* Construct Dual Graph
 *
 * Start in the corner of a triangle of the cusp triangulation
 * and walk around the boundary of the manifold, add edges to
 * the dual graph.
 */

void construct_dual_graph(struct graph *graph1, Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int index;
    struct CuspTriangle *tri;
    struct queue queue1;

    initialise_queue(&queue1, 3 * manifold->num_tetrahedra);

    // Start at the inside corner of triangle 1.
    enqueue(&queue1, 0);
    graph1->vertexHomology[0][0] = 0;          // tet vertex
    graph1->vertexHomology[0][1] = 0;          // cusp vertex
    graph1->vertexHomology[0][2] = 0;          // distance

    while (!empty_queue(&queue1)) {
        index = dequeue(&queue1);
        tri = pTriangle[index];


    }
}

void printTriangleInfo(Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int i, j;

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++) {
        for (j = 0; j < 3; j++) {
            printf("Cusp Edge %d of Tet Vertex %d of Tet %d glues to Cusp Edge %d of Tet Vertex %d on Tet %d\n",
                   pTriangle[i]->edges[j],      // Cusp Edge
                   pTriangle[i]->tetVertex,     // Tet Vertex
                   pTriangle[i]->tet->index,    // Tet Index
                   EVALUATE(pTriangle[i]->tet->gluing[pTriangle[i]->edges[j]], pTriangle[i]->edges[j]),    // Cusp Edge
                   EVALUATE(pTriangle[i]->tet->gluing[pTriangle[i]->edges[j]], pTriangle[i]->tetVertex),     // Tet Vertex
                   pTriangle[i]->tet->neighbor[pTriangle[i]->edges[j]]->index);                             // Tet Index
        }
    }
}


int flow(struct CuspTriangle *pTri, int i, int j) {
//    Tetrahedron *tet = pTri->tet;
//    int mflow = FLOW(tet->curve[0][right_handed][v][f], tet->curve[0][right_handed][v][ff]);
//    int lflow = FLOW(tet->curve[1][right_handed][v][f], tet->curve[1][right_handed][v][ff]);
//
//    return mflow + lflow;
}

/* Remove Extra Edges
 *
 * Remove edges associated to the curves that dive into the
 * manifold since we can only do this once.
 */
void remove_extra_edges(struct graph *graph1) {

}

/* Add Misc Edges
 *
 * Add edges associated to
 *  - Diving through the manifold
 *  -
 */
void add_misc_edges(struct graph *graph1) {

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

void initialise_graph(struct graph *g, int maxVertices, bool directed) {
    int i, j;

    g->nvertices = maxVertices;
    g->nedges = 0;
    g->directed = directed;

    g->edges = NEW_ARRAY(maxVertices, edgenode *);
    g->degree = NEW_ARRAY(maxVertices, int);
    g->vertexHomology = NEW_ARRAY(maxVertices, int *);

    for (i = 0; i < maxVertices; i++) {
        g->degree[i] = 0;
        g->edges[i] = NULL;

        g->vertexHomology[i] = NEW_ARRAY(7, int);
        for (j = 0; j < 6; j++)
            g->vertexHomology[i][j] = 0;
    }
}

void free_graph(struct graph *g) {
    int i;

    free(g->edges);
    free(g->degree);

    for (i = 0; i < g->nvertices; i++)
        free(g->vertexHomology[i]);
    free(g->vertexHomology);
}

void insert_edge(struct graph *g, int x, int y, bool directed) {
    edgenode *p;

    p = malloc(sizeof( edgenode ));
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

void initialise_search(struct graph *g, bool *processed, bool *discovered, int *parent) {
    int i;

    for (i = 0; i <= g->nvertices; i ++) {
        processed[i] = false;
        discovered[i] = false;
        parent[i] = -1;
    }
}

void bfs(struct graph *g, int start, bool *processed, bool *discovered, int *parent) {
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