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

/*
 * Get Symplectic Basis
 *
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the get_gluing_equations and get_symplectic_equations functions
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i, j;

    // Edge Curves C_i -> gluing equations
    int **edge_eqns, edge_num_rows;

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;

    // Get Gluing Equations
    edge_eqns = get_gluing_equations(manifold, &edge_num_rows, num_cols);
    symp_num_rows = edge_num_rows;

    // Allocate Symplectic Equations Array
    symp_eqns = NEW_ARRAY(symp_num_rows, int *);

    for (i = 0; i < symp_num_rows; i ++)
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

    // Allocate Cusp Triangulation Array
    struct CuspTriangle **pTriangle = NEW_ARRAY(4 * manifold->num_tetrahedra, struct CuspTriangle *);

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++)
        pTriangle[i] = malloc(sizeof( struct CuspTriangle ));

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, pTriangle, symp_num_rows, symp_eqns);

    // Construct return array
    int **eqns = NEW_ARRAY(edge_num_rows + symp_num_rows, int *);

    for (i = 0; i < edge_num_rows; i++) {
        eqns[2 * i] = edge_eqns[i];
        eqns[2 * i + 1] = symp_eqns[i];
    }

    free_cusp_triangulation(manifold, pTriangle);
    *num_rows = edge_num_rows + symp_num_rows;
    return eqns;
}

/*
 * Initialise Cusp Triangulation
 *
 * Construct the pTriangle array which consists of the triangles in the cusp triangulation
 */

void init_cusp_triangulation(Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int i, j, k;
    EdgeClass *edge;
    Tetrahedron *tet= manifold->tet_list_begin.next;

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++) {
        pTriangle[i]->tet = tet;
        pTriangle[i]->tetVertex = i % 4;

        switch (pTriangle[i]->tetVertex) {
            case 0:
                pTriangle[i]->edgesThreeToFour[0] = 1;
                pTriangle[i]->edgesThreeToFour[1] = 2;
                pTriangle[i]->edgesThreeToFour[2] = 3;
                break;
            case 1:
                pTriangle[i]->edgesThreeToFour[0] = 0;
                pTriangle[i]->edgesThreeToFour[1] = 2;
                pTriangle[i]->edgesThreeToFour[2] = 3;
                break;
            case 2:
                pTriangle[i]->edgesThreeToFour[0] = 0;
                pTriangle[i]->edgesThreeToFour[1] = 1;
                pTriangle[i]->edgesThreeToFour[2] = 3;
                break;
            case 3:
                pTriangle[i]->edgesThreeToFour[0] = 0;
                pTriangle[i]->edgesThreeToFour[1] = 1;
                pTriangle[i]->edgesThreeToFour[2] = 2;
                break;
            default:
                break;
        }

        switch (pTriangle[i]->tetVertex) {
            case 0:
                pTriangle[i]->edgesFourToThree[0] = -1;
                pTriangle[i]->edgesFourToThree[1] = 0;
                pTriangle[i]->edgesFourToThree[2] = 1;
                pTriangle[i]->edgesFourToThree[3] = 2;
                break;
            case 1:
                pTriangle[i]->edgesFourToThree[0] = 0;
                pTriangle[i]->edgesFourToThree[1] = -1;
                pTriangle[i]->edgesFourToThree[2] = 1;
                pTriangle[i]->edgesFourToThree[3] = 2;
                break;
            case 2:
                pTriangle[i]->edgesFourToThree[0] = 0;
                pTriangle[i]->edgesFourToThree[1] = 1;
                pTriangle[i]->edgesFourToThree[2] = -1;
                pTriangle[i]->edgesFourToThree[3] = 2;
                break;
            case 3:
                pTriangle[i]->edgesFourToThree[0] = 0;
                pTriangle[i]->edgesFourToThree[1] = 1;
                pTriangle[i]->edgesFourToThree[2] = 2;
                pTriangle[i]->edgesFourToThree[3] = -1;
                break;
            default:
                break;
        }

        for (j = 0; j < 3; j++) {
            // Edge between pTriangle[i]->vertex and pTriangle[i]->edges[j]
            pTriangle[i]->vertices[j].v1 = pTriangle[i]->tetVertex;
            pTriangle[i]->vertices[j].v2 = pTriangle[i]->edgesThreeToFour[j];

            pTriangle[i]->vertices[j].edge = tet->edge_class[
                    edge_between_vertices[
                            pTriangle[i]->vertices[j].v1][
                            pTriangle[i]->vertices[j].v2]];

            k = 0;
            for (edge = manifold->edge_list_begin.next; edge != pTriangle[i]->vertices[j].edge; edge = edge->next, k++);

            pTriangle[i]->vertices[j].index = k;
        }

        for (j = 0; j < 4; j++) {
            pTriangle[i]->neighbours[j] = pTriangle[4 * (pTriangle[i]->tet->neighbor[j]->index)
                                                    + EVALUATE(pTriangle[i]->tet->gluing[j], pTriangle[i]->tetVertex)];
        }


        if (i % 4 == 3) {
            tet = tet->next;
        }
    }
}

// Free memory for cusp triangulation data structure
void free_cusp_triangulation(Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int i;

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++)
        free(pTriangle[i]);

    my_free(pTriangle);
}

/*
 * Get Symplectic Equations
 *
 * Setup graph and cusp triangulation, and run construct dual graph.
 */

int **get_symplectic_equations(Triangulation *manifold, struct CuspTriangle **pTriangle, int num_rows, int **eqns) {
    int i, j, T = manifold -> num_tetrahedra;
    int genus = (2 * manifold->num_tetrahedra - num_rows) / 2;
    struct graph graph1;

    init_cusp_triangulation(manifold, pTriangle);
    initialise_graph(&graph1, 12 * (2 * genus + 1) * manifold->num_tetrahedra, FALSE);
    construct_dual_graph(&graph1, manifold, pTriangle);
    //printDebugInfo(manifold, pTriangle, NULL);

    // Dual Curve Equations
    for (i = 0; i < num_rows; i ++) {
        for (j = 0; j < 3 * T; j ++)
            eqns[i][j] = i + j;
    }

    free_graph(&graph1);

    return eqns;
}

// Free memory used to find the symplectic basis
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
    int i, index, index1, index2, cuspEdge, cuspVertex;
    struct CuspTriangle *tri, *adjTri;
//    struct Queue queue1;

//    initialise_queue(&queue1, 3 * manifold->num_tetrahedra);

    struct Node stack;
    stack.item = -1;

    int visited[graph1->nvertices];

    for (i = 0; i < graph1->nvertices; i++)
        visited[i] = 0;

    // Start at the inside corner of triangle 1.
    push(&stack, 0);
    graph1->vertexHomology[0][0] = 0;          // Tet Index
    graph1->vertexHomology[0][1] = 0;          // Tet Vertex
    graph1->vertexHomology[0][2] = 0;          // dist to v1
    graph1->vertexHomology[0][3] = flow(pTriangle[0], pTriangle[0]->edgesThreeToFour[0])
            + flow(pTriangle[0], pTriangle[0]->edgesThreeToFour[1]);                       // dist to v2
    graph1->vertexHomology[0][4] = flow(pTriangle[0], pTriangle[0]->edgesThreeToFour[0])
            + flow(pTriangle[0], pTriangle[0]->edgesThreeToFour[2]);                       // dist to v3

    while (!is_empty(&stack)) {
        index = pop(&stack);

        if (visited[index])
            continue;

        // Find the triangle that vertex (index) of the dual graph lies on.
        tri = findTriangle(manifold, pTriangle, graph1->vertexHomology[index][0], graph1->vertexHomology[index][1]);
        cuspVertex = tri->edgesThreeToFour[minCuspDistance(graph1, index)];

        // First two edges can always be reached
        cuspEdge = (int) remaining_face[tri->tetVertex][cuspVertex];
        adjTri = tri->neighbours[cuspEdge];
        index1 = insert_edge(graph1, index, cuspEdge, tri, adjTri, FALSE);

        if (index1 == -1) {
            printf("Inserting edge from triangle (%d, %d) across edge %d to triangle (%d, %d) failed\n",
                   tri->tet->index,
                   tri->tetVertex,
                   cuspVertex,
                   adjTri->tet->index,
                   adjTri->tetVertex
            );
            printDebugInfo(manifold, pTriangle, graph1);
        }

        cuspEdge = (int) remaining_face[cuspVertex][tri->tetVertex];
        adjTri = tri->neighbours[cuspEdge];
        index2 = insert_edge(graph1, index, cuspEdge, tri, adjTri, FALSE);

        if (index2 == -1) {
            printf("Inserting edge from triangle (%d, %d) across edge %d to triangle (%d, %d) failed\n",
                   tri->tet->index,
                   tri->tetVertex,
                   cuspVertex,
                   adjTri->tet->index,
                   adjTri->tetVertex
            );
            printDebugInfo(manifold, pTriangle, graph1);
        }

        // Add vertices to the queue
        if (!visited[index1])
            push(&stack, index1);
        if (!visited[index2])
            push(&stack, index2);

        // Add the third edge if possible
        if (flow(tri, cuspVertex) <= graph1->vertexHomology[index][tri->edgesFourToThree[cuspVertex] + 2]) {
            adjTri = tri->neighbours[cuspVertex];
            index1 = insert_edge(graph1, index, cuspVertex, tri, adjTri, FALSE);

            if (index1 == -1) {
                printf("Inserting edge from triangle (%d, %d) across edge %d to triangle (%d, %d) failed\n",
                       tri->tet->index,
                       tri->tetVertex,
                       cuspVertex,
                       adjTri->tet->index,
                       adjTri->tetVertex
                       );
                printDebugInfo(manifold, pTriangle, graph1);
            }

            push(&stack, index1);
        }

        visited[index] = 1;
    }
}

/*
 * Min Cusp Distance
 *
 * Calculate the cusp vertex which the graph vertex is closest to.
 */

int minCuspDistance(struct graph *g, int index) {
    int dist0, dist1, dist2, min_dist;

    dist0 = g->vertexHomology[index][2];
    dist1 = g->vertexHomology[index][3];
    dist2 = g->vertexHomology[index][4];

    min_dist = MIN(dist0, MIN(dist1, dist2));

    if (min_dist == dist0)
        return 0;
    else if (min_dist == dist1)
        return 1;
    else
        return 2;
}


/*
 * Print Triangle Information
 *
 * Debug function for printing the gluing information
 */

void printDebugInfo(Triangulation *manifold, struct CuspTriangle **pTriangle, struct graph *g) {
    int i, j;

    struct CuspTriangle *tri;

//    for (i = 0; i < 4 * manifold->num_tetrahedra; i++) {
//        for (j = 0; j < 3; j++) {
//            printf("Cusp Edge %d of Tet Vertex %d of Tet %d glues to Cusp Edge %d of Tet Vertex %d on Tet %d\n",
//                   pTriangle[i]->edges[j],      // Cusp Edge
//                   pTriangle[i]->tetVertex,     // Tet Vertex
//                   pTriangle[i]->tet->index,    // Tet Index
//                   EVALUATE(pTriangle[i]->tet->gluing[pTriangle[i]->edges[j]], pTriangle[i]->edges[j]),    // Cusp Edge
//                   EVALUATE(pTriangle[i]->tet->gluing[pTriangle[i]->edges[j]], pTriangle[i]->tetVertex),     // Tet Vertex
//                   pTriangle[i]->tet->neighbor[pTriangle[i]->edges[j]]->index);                             // Tet Index
//        }
//    }

    for (i = 0; i < g->nvertices; i++) {
        if (g->vertexHomology[i][0] == -1 || g->vertexHomology[i][1] == -1) {
            tri = pTriangle[0];
        } else {
            tri = findTriangle(manifold, pTriangle, g->vertexHomology[i][0], g->vertexHomology[i][1]);
        }

        printf("Tet: %d, Tet Vertex: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d.\n",
               g->vertexHomology[i][0],
               g->vertexHomology[i][1],
               tri->edgesThreeToFour[0],
               g->vertexHomology[i][2],
               tri->edgesThreeToFour[1],
               g->vertexHomology[i][3],
               tri->edgesThreeToFour[2],
               g->vertexHomology[i][4]
        );
    }

    printf("-------------------------------\n");
}

/*
 * Flow
 *
 * Calculate the number of curves passing around a vertex in the cusp triangulation.
 */

int flow(struct CuspTriangle *tri, int vertex) {
    int mflow, lflow, retval;

    // Contribution from meridian curves
    mflow = FLOW(tri->tet->curve[0][right_handed][tri->tetVertex][remaining_face[tri->tetVertex][vertex]],
                 tri->tet->curve[0][right_handed][tri->tetVertex][remaining_face[vertex][tri->tetVertex]]);

    // Contribution from longitudinal curves
    lflow = FLOW(tri->tet->curve[1][right_handed][tri->tetVertex][remaining_face[tri->tetVertex][vertex]],
                 tri->tet->curve[1][right_handed][tri->tetVertex][remaining_face[vertex][tri->tetVertex]]);

    retval = ABS(mflow) + ABS(lflow);
    return retval;
}

/*
 * Find Triangle
 *
 * Returns a pointer to the triangle with tetIndex and tetVertex.
 */
struct CuspTriangle *findTriangle(Triangulation *manifold, struct CuspTriangle **pTriangle, int tetIndex, int tetVertex) {
    int i;

    for (i = 0; i < 4 * manifold->num_tetrahedra; i++) {
        if (tetIndex == pTriangle[i]->tet->index && tetVertex == pTriangle[i]->tetVertex)
            return pTriangle[i];
    }

    return NULL;
}

int visited(int **array, int *point, int lenArray, int lenPoint) {
    int i, j, equal;

    for (i = 0; i < lenArray; i++) {
        equal = 1;

        for (j = 0; j < lenPoint; j++) {
            if (array[i][j] != point[j])
                equal = 0;
        }

        if (equal == 1)
            return 1;
    }

    return 0;
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

/*
 * Queue Data Structure
 */

void initialise_queue(struct Queue *q, int size) {
    q->front = 0;
    q->rear = -1;
    q->len = 0;
    q->size = size;
    q->array = NEW_ARRAY(size, int);
}

struct Queue *enqueue(struct Queue *q, int i) {
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

int dequeue(struct Queue *q) {
    // User to verify queue is not empty
    int i = q->array[q->front];

    q->front = (q->front + 1) % (q->size - 1);
    q->len--;

    return i;
}

int empty_queue(struct Queue *q) {
    return (!q->len);
}

void resize_queue(struct Queue *q) {
    int i;
    struct Queue p;

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

void free_queue(struct Queue *q) {
    free(q->array);
}

/*
 * Stack Data Structure
 */

void push(struct Node *stack, int item) {
    struct Node *node = malloc(sizeof( struct Node));
    node->item = stack->item;
    node->next = stack->next;
    stack->item = item;
    stack->next = node;
}

int pop(struct Node *stack) {
    int item = stack->item;
    struct Node *node = stack->next;
    stack->item = stack->next->item;
    stack->next = stack->next->next;
    free(node);

    return item;
}

int is_empty(struct Node *stack) {
    return stack->item == -1;
}

// Breadth First Search

/*
 * Initialise Graph
 *
 * Initialise the arrays of the graph 'g' to their default values
 */
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

        g->vertexHomology[i] = NEW_ARRAY(6, int);
        for (j = 0; j < 6; j++)
            g->vertexHomology[i][j] = -1;
    }
}

// Free the memory of the graph data structure
void free_graph(struct graph *g) {
    int i;

    free(g->edges);
    free(g->degree);

    for (i = 0; i < g->nvertices; i++)
        free(g->vertexHomology[i]);
    free(g->vertexHomology);
}

/*
 * Insert Edge
 *
 * Insert an edge into the graph 'g' from vertex 'index' on triangle 'x' to triangle 'y'
 */

int insert_edge(struct graph *g, int index, int face, struct CuspTriangle *x, struct CuspTriangle *y, bool directed) {
    edgenode *p;
    int i, x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, temp;

    // Vertices on triangle x which are glued to triangle y
    x_vertex1 = (int) remaining_face[x->tetVertex][face];
    x_vertex2 = (int) remaining_face[face][x->tetVertex];

    // Vertices on triangle y which correspond to the vertices v1 and v2
    y_vertex1 = EVALUATE(x->tet->gluing[face], x_vertex1);
    y_vertex2 = EVALUATE(x->tet->gluing[face], x_vertex2);
    y_face = EVALUATE(x->tet->gluing[face], face);

    // Find index of y
    for (i = 0; i < g->nvertices; i ++) {
        if (g->vertexHomology[i][0] == y->tet->index                                    // Tet Index
        && g->vertexHomology[i][1] == y->tetVertex                                      // Tet Vertex
        && g->vertexHomology[i][y->edgesFourToThree[y_vertex1] + 2] ==
            g->vertexHomology[index][x->edgesFourToThree[x_vertex1] + 2]                // Dist. to y_vertex1
        && g->vertexHomology[i][y->edgesFourToThree[y_vertex2] + 2] ==
            g->vertexHomology[index][x->edgesFourToThree[x_vertex2] + 2]                // Dist. to y_vertex2
        && g->vertexHomology[i][y->edgesFourToThree[y_face] + 2] == flow(y, y_face)            // Dist. to vertex 3
        ) {
            break;
        }
    }

    // Insert new vertex
    if (i == g->nvertices) {
        for (i = 0; i < g->nvertices && g->vertexHomology[i][0] != -1; i++);

        if (i == g->nvertices)
            return -1;

        g->vertexHomology[i][0] = y->tet->index;
        g->vertexHomology[i][1] = y->tetVertex;
        g->vertexHomology[i][y->edgesFourToThree[y_vertex1] + 2] = g->vertexHomology[index][x->edgesFourToThree[x_vertex1] + 2];
        g->vertexHomology[i][y->edgesFourToThree[y_vertex2] + 2] = g->vertexHomology[index][x->edgesFourToThree[x_vertex2] + 2];
        g->vertexHomology[i][y->edgesFourToThree[y_face] + 2] = flow(y, y_face);
    }

    // update_vertex_homology(g, index, i, x, y);

    p = malloc(sizeof( edgenode ));
    p->y = i;
    p->next = g->edges[index];

    g->edges[index] = p;
    g->degree[index]++;

    if (!directed) {
        insert_edge(g, i, EVALUATE(x->tet->gluing[face], face), y, x, true);
    }

    return i;
}

/*
 * Update Vertex Homology
 *
 * Updates the vertex homology array of the graph
 *
 * !! Deprecated !!
 */

void update_vertex_homology(struct graph *g, int index1, int index2, struct CuspTriangle *tri1, struct CuspTriangle *tri2) {
    int flow_current_vertex, flow_vertex1, flow_vertex2;

    flow_current_vertex = flow(tri2, g->vertexHomology[index2][2]);
    flow_vertex1 = flow(tri2, remaining_face[g->vertexHomology[index2][2]][tri2->tetVertex]);
    flow_vertex2 = flow(tri2, remaining_face[tri2->tetVertex][g->vertexHomology[index2][2]]);

    // Check if index2 lies in the "center" of tri2
    if (g->vertexHomology[index2][3] == flow_current_vertex) {
        g->vertexHomology[index2][3] = -1;
    } else if (flow_vertex2 < flow_current_vertex || flow_vertex1 < flow_current_vertex) {
        if (flow_vertex1 <= flow_vertex2) {
            g->vertexHomology[index2][2] = (int) remaining_face[g->vertexHomology[index2][2]][tri2->tetVertex];
            g->vertexHomology[index2][3] = flow_vertex1 + flow_current_vertex - g->vertexHomology[index2][3];
        } else {
            g->vertexHomology[index2][2] = (int) remaining_face[tri2->tetVertex][g->vertexHomology[index2][2]];
            g->vertexHomology[index2][3] = flow_vertex2 + flow_current_vertex - g->vertexHomology[index2][3];
        }
    }
}

/*
 * Initialise Search
 *
 * Initialise default values for bfs arrays
 */

void initialise_search(struct graph *g, bool *processed, bool *discovered, int *parent) {
    int i;

    for (i = 0; i <= g->nvertices; i ++) {
        processed[i] = false;
        discovered[i] = false;
        parent[i] = -1;
    }
}

/*
 * Breadth First Search
 *
 * Graph search algorithm starting at vertex 'start'.
 */

void bfs(struct graph *g, int start, bool *processed, bool *discovered, int *parent) {
    struct Queue q;
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

/*
 * Find Path
 *
 * Recover the path through the graph from the parents array
 */
void find_path(int start, int end, int *parents, int *path, int index) {
    if ((start == end) || (end == -1)) {
        path[index] = start;
    } else {
        find_path(start, parents[end], parents, path, index + 1);
        path[index] = end;
    }
}

