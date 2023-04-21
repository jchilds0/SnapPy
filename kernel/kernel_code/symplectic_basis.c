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

#define atleast_two(a, b, c)    (a && b) || (a && c) || (b && c)

int edgesThreeToFour[4][3] = {{1, 2, 3},
                              {0, 2, 3},
                              {0, 1, 3},
                              {0, 1, 2}};
int edgesFourToThree[4][4] = {{9, 0, 1, 2},
                              {0, 9, 1, 2},
                              {0, 1, 9, 2},
                              {0, 1, 2, 9}};

int numCuspTriangles;

/*
 * Get Symplectic Basis
 *
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the get_gluing_equations and get_symplectic_equations functions
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i, j, e0 = 0;
    numCuspTriangles = 4 * manifold->num_tetrahedra;

    // Edge Curves C_i -> gluing equations
    int **edge_eqns, edge_num_rows;

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;

    // Get Gluing Equations
    edge_eqns = get_gluing_equations(manifold, &edge_num_rows, num_cols);
    symp_num_rows = edge_num_rows;

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, symp_num_rows, *num_cols, e0);

    // Construct return array
    *num_rows = edge_num_rows + symp_num_rows - 2;
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < edge_num_rows; i++) {
        if (i == e0)
            continue;

        eqns[2 * j] = edge_eqns[i];
        eqns[2 * j + 1] = symp_eqns[i];
        j++;
    }

    my_free(symp_eqns);
    my_free(edge_eqns);
    return eqns;
}

/*
 * Initialise Cusp Triangulation
 *
 * Construct the pTriangle array which consists of the triangles in the cusp triangulation
 */

static int numEdgeClasses;

struct CuspTriangle **init_cusp_triangulation(Triangulation *manifold) {
    int i, j;
    Tetrahedron *tet= manifold->tet_list_begin.next;

    // Allocate Cusp Triangulation Array
    struct CuspTriangle **pTriangle = NEW_ARRAY(numCuspTriangles, struct CuspTriangle *);

    for (i = 0; i < numCuspTriangles; i++)
        pTriangle[i] = malloc(sizeof( struct CuspTriangle ));

    label_triangulation_edges(manifold);

    for (i = 0; i < numCuspTriangles; i++) {
        pTriangle[i]->tet = tet;
        pTriangle[i]->tetVertex = i % 4;

        for (j = 0; j < 4; j++) {
            if (pTriangle[i]->tetVertex == j) {
                pTriangle[i]->neighbours[j] = NULL;
            } else {
                pTriangle[i]->neighbours[j] = pTriangle[4 * (pTriangle[i]->tet->neighbor[j]->index)
                                                        + EVALUATE(pTriangle[i]->tet->gluing[j], pTriangle[i]->tetVertex)];
            }
        }

        for (j = 0; j < 3; j++) {
            pTriangle[i]->faces[j].index = -1;

            pTriangle[i]->vertices[j].v1 = pTriangle[i]->tetVertex;
            pTriangle[i]->vertices[j].v2 = edgesThreeToFour[pTriangle[i]->tetVertex][j];
            pTriangle[i]->vertices[j].edge = pTriangle[i]->tet->edge_class[
                    edge_between_vertices[pTriangle[i]->vertices[j].v1][pTriangle[i]->vertices[j].v2]];
            pTriangle[i]->vertices[j].edgeIndex = pTriangle[i]->vertices[j].edge->index;
            pTriangle[i]->vertices[j].vertexIndex = -1;
        }

        if (i % 4 == 3) {
            tet = tet->next;
        }
    }

    cusp_vertex_index(pTriangle);
    //label_cusp_faces(pTriangle);
    print_debug_info(pTriangle, NULL, NULL, NULL, NULL, 4);
    return pTriangle;
}

/*
 * Label edges
 *
 * Give each edge of the triangulation an index to identify the cusp vertices
 */

void label_triangulation_edges(Triangulation *manifold) {
    int i = 0;
    EdgeClass *edge = &manifold->edge_list_begin;

    while ((edge = edge->next)->next != NULL)
        edge->index = i++;
    
    numEdgeClasses = i;
}

/*
 * Find cusp vertex index
 *
 * Each edge class of the manifold appears as two vertices
 * in the cusp triangulation. We iterate over the cusp triangulation,
 * walking around each vertex to give it the same index.
 */

void cusp_vertex_index(struct CuspTriangle **pTriangle) {
    int i, j, currentIndex[numEdgeClasses];
    
    for (i = 0; i < numEdgeClasses; i++) 
        currentIndex[i] = 0;
    
    for (i = 0; i < numCuspTriangles; i++) {
        for (j = 0; j < 3; j++) {
            if (pTriangle[i]->vertices[j].vertexIndex != -1)
                continue;

            walk_around_vertex(pTriangle, pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][j],
                               currentIndex[pTriangle[i]->vertices[j].edgeIndex]);

            currentIndex[pTriangle[i]->vertices[j].edgeIndex]++;
        }
    }
}

/*
 * Walk around cusp vertex
 *
 * Walk around vertex cuspVertex of triangle pTriangle[cuspTriangle]
 * and set vertexIndex to index.
 */

void walk_around_vertex(struct CuspTriangle **pTriangle, struct CuspTriangle *tri, int cuspVertex, int index) {
    int gluing_vertex, outside_vertex, old_gluing_vertex, old_cusp_vertex, old_outside_vertex;
    gluing_vertex = (int) remaining_face[cuspVertex][tri->tetVertex];
    outside_vertex = (int) remaining_face[tri->tetVertex][cuspVertex];

    while (tri->vertices[edgesFourToThree[tri->tetVertex][cuspVertex]].vertexIndex == -1) {
        tri->vertices[edgesFourToThree[tri->tetVertex][cuspVertex]].vertexIndex = index;

        // Move to the next cusp triangle
        old_cusp_vertex = cuspVertex;
        old_gluing_vertex = gluing_vertex;
        old_outside_vertex = outside_vertex;

        cuspVertex = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
        gluing_vertex = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_outside_vertex);
        outside_vertex = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_gluing_vertex);
        tri = find_cusp_triangle(pTriangle, tri->neighbours[old_gluing_vertex]->tet->index,
                                          EVALUATE(tri->tet->gluing[old_gluing_vertex], tri->tetVertex));
    }
}

// Free memory for cusp triangulation data structure
void free_cusp_triangulation(struct CuspTriangle **pTriangle) {
    int i;

    for (i = 0; i < numCuspTriangles; i++)
        free(pTriangle[i]);

    my_free(pTriangle);
}

/*
 * Get Symplectic Equations
 *
 * Setup graph and cusp triangulation, and run construct dual graph.
 */
static int genus;
static int numDualCurves;

int **get_symplectic_equations(Triangulation *manifold, int num_rows, int numCols, int e0) {
    int i, j, numVertices, T = manifold -> num_tetrahedra;
    genus = (2 * manifold->num_tetrahedra - num_rows) / 2;
    numDualCurves = num_rows;
    struct Graph *graph1;
    struct CuspTriangle **pTriangle;

    int **dualCurve = NEW_ARRAY(2 * numDualCurves, int *);
    int *dualCurveLen = NEW_ARRAY(2 * numDualCurves, int);

    for (i = 0; i < 2 * numDualCurves; i++) {
        dualCurveLen[i] = 0;
        dualCurve[i] = NEW_ARRAY(numCols, int);
    }

    // Allocate Symplectic Equations Array
    int **symp_eqns = NEW_ARRAY(num_rows, int *);

    for (i = 0; i < num_rows; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    pTriangle = init_cusp_triangulation(manifold);
    numVertices = (2 * genus + 1) * numCuspTriangles + 4 * genus;
    graph1 = init_graph(numVertices, FALSE);
    struct CuspNode **vertexData = NEW_ARRAY(numVertices, struct CuspNode *);

    for (i = 0; i < numVertices; i++) {
        vertexData[i] = NEW_STRUCT(struct CuspNode);

        vertexData[i]->adjTri[0] = 0;
        vertexData[i]->adjTri[1] = 0;
        vertexData[i]->adjTri[2] = 0;
        vertexData[i]->dist[0] = -1;
        vertexData[i]->dist[1] = -1;
        vertexData[i]->dist[2] = -1;
        vertexData[i]->tetVertex = -1;
        vertexData[i]->tetIndex = -1;
    }

    construct_dual_graph(manifold, graph1, pTriangle, vertexData);
    print_debug_info(pTriangle, graph1, vertexData, NULL, NULL, 0);
    print_debug_info(pTriangle, graph1, vertexData, NULL, NULL, 2);
    print_debug_info(pTriangle, graph1, vertexData, NULL, NULL, 3);

    graph1 = construct_dual_curves(graph1, pTriangle, vertexData, e0, dualCurve, dualCurveLen);
    find_holonomies(graph1, symp_eqns, dualCurve, dualCurveLen);

    print_debug_info(pTriangle, graph1, vertexData, dualCurve, dualCurveLen, 5);

    // Dual Curve Equations
    for (i = 0; i < num_rows; i ++) {
        for (j = 0; j < 3 * T; j ++)
            symp_eqns[i][j] = i + j;
    }

    free_graph(graph1);
    free_cusp_triangulation(pTriangle);

    for (i = 0; i < numVertices; i++)
        my_free(vertexData[i]);
    my_free(vertexData);

    for (i = 0; i < 2 * numDualCurves; i++)
        my_free(dualCurve[i]);
    my_free(dualCurve);
    my_free(dualCurveLen);

    return symp_eqns;
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
static int intersectTetIndex, intersectTetVertex;

void construct_dual_graph(Triangulation *manifold, struct Graph *graph1, struct CuspTriangle **pTriangle, struct CuspNode **vertexData) {
    int i, index, cuspEdge;
    int indices[3] = { 0, 0, 0 };
    struct CuspTriangle *adjTri[3];
    struct CuspNode *node;

    // Call fundamental group to find standard basepoint of the cusp triangulation
    fundamental_group(manifold, FALSE, FALSE, FALSE, FALSE);
    intersectTetIndex = manifold->cusp_list_begin.next->basepoint_tet->index;
    intersectTetVertex = (int) manifold->cusp_list_begin.next->basepoint_vertex;

    struct Node stack;
    stack.item = -1;

    int visited[graph1->nvertices];

    for (i = 0; i < graph1->nvertices; i++)
        visited[i] = 0;

    // Start at the inside corner of triangle 1.
    push(&stack, 0);
    vertexData[0]->tri = pTriangle[0];
    vertexData[0]->tetVertex = 0;          // Tet Index
    vertexData[0]->tetIndex = 0;             // Tet Vertex
    vertexData[0]->dist[0] = 0;             // dist to v1
    vertexData[0]->dist[1] = flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][0])
            + flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][1]);                       // dist to v2
    vertexData[0]->dist[2] = flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][0])
            + flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][2]);                       // dist to v3
    vertexData[0]->adjTri[0] = (flow(pTriangle[0], 1) == 0);
    vertexData[0]->adjTri[1] = 1;
    vertexData[0]->adjTri[2] = 1;
    graph1->vertexData[0] = 0;

    // Walk around the cusp triangulation inserting edges
    while (!is_empty(&stack)) {
        index = pop(&stack);
        node = vertexData[index];

        if (visited[index])
            continue;

        for (i = 0; i < 3; i++) {
            if (vertexData[index]->adjTri[i]) {
                cuspEdge = edgesThreeToFour[node->tetVertex][i];
                adjTri[i] = node->tri->neighbours[cuspEdge];
                indices[i] = insert_triangle_edge(graph1, index, cuspEdge, node->tri, adjTri[i], vertexData);
            } else {
                indices[i] = -2;
            }
        }


        for (i = 0; i < 3; i++) {
            if (indices[i] == -1) {
                printf("Inserting edge from triangle (%d, %d) to triangle (%d, %d) failed\n",
                       node->tri->tet->index,
                       node->tetVertex,
                       adjTri[i]->tet->index,
                       adjTri[i]->tetVertex
                );
                print_debug_info(pTriangle, graph1, NULL, NULL, NULL, 2);
            }

            if (!visited[indices[i]] && indices[i] != -2)
                push(&stack, indices[i]);
        }

        visited[index] = 1;
    }
    
    //remove_extra_edges(graph1);
    //add_misc_edges(graph1, edges);
}

int insert_triangle_edge(struct Graph *g, int index, int face, struct CuspTriangle *x, struct CuspTriangle *y, struct CuspNode **vertexData) {
    int i, x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist;

    // Vertices on triangle x which are glued to triangle y
    x_vertex1 = (int) remaining_face[x->tetVertex][face];
    x_vertex2 = (int) remaining_face[face][x->tetVertex];

    // Vertices on triangle y which correspond to the vertices v1 and v2
    y_vertex1 = EVALUATE(x->tet->gluing[face], x_vertex1);
    y_vertex2 = EVALUATE(x->tet->gluing[face], x_vertex2);
    y_face = EVALUATE(x->tet->gluing[face], face);

    // Calculate distance to vertex y_face
    if (y->tetVertex == intersectTetVertex && y->tet->index == intersectTetIndex) {
        // Intersect vertex: cross to the closest cusp vertex and then to y_face vertex
        dist = MIN(flow(y, y_vertex1) + vertexData[index]->dist[edgesFourToThree[x->tetVertex][x_vertex1]],
                   flow(y, y_vertex2) + vertexData[index]->dist[edgesFourToThree[x->tetVertex][x_vertex2]])
                           + flow(y, y_face);
    } else {
        // Normal cusp triangle
        dist = flow(y, y_face);

        if (vertexData[index]->dist[edgesFourToThree[x->tetVertex][x_vertex1]] < flow(y, y_vertex1)) {
            // Inside the flows around y_vertex1
            dist += (flow(y, y_vertex1) - vertexData[index]->dist[edgesFourToThree[x->tetVertex][x_vertex1]]);
        } else {
            // Inside the flows around y_vertex2 or center
            dist += (flow(y, y_vertex2) - vertexData[index]->dist[edgesFourToThree[x->tetVertex][x_vertex2]]);
        }
    }

    // Find index of y
    for (i = 0; i < g->nvertices; i++) {
        if (is_equal(vertexData[index], vertexData[i], x, y,
                     x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist)) {
            break;
        }
    }

    // Insert new vertex
    if (i == g->nvertices) {
        for (i = 0; i < g->nvertices && vertexData[i]->tetIndex != -1; i++);

        if (i == g->nvertices)
            return -1;

        init_vertex(vertexData[index], vertexData[i], x, y,
                    x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist);
    }

    // Ignore edge if it already exists
    if (!edge_exists(g, index, i))
        insert_edge(g, index, i, g->directed);

    g->vertexData[i] = i;
    return i;
}

/*
 * Is Equal
 *
 * Check if two holonomies correspond to the same vertex
 */

int is_equal(struct CuspNode *xNode, struct CuspNode *yNode, struct CuspTriangle *x, struct CuspTriangle *y,
             int x_vertex1, int x_vertex2, int y_vertex1, int y_vertex2, int y_face, int dist) {
    int tetIndex, tetVertex, distTriVertex1, distTriVertex2, distTriVertex3, intersectFace;

    tetIndex = y->tet->index == yNode->tetIndex;
    tetVertex = y->tetVertex == yNode->tetVertex;

    distTriVertex1 = (yNode->dist[edgesFourToThree[yNode->tetVertex][y_vertex1]] ==
                      xNode->dist[edgesFourToThree[xNode->tetVertex][x_vertex1]]);
    distTriVertex2 = (yNode->dist[edgesFourToThree[yNode->tetVertex][y_vertex2]] ==
                      xNode->dist[edgesFourToThree[xNode->tetVertex][x_vertex2]]);
    distTriVertex3 = yNode->dist[edgesFourToThree[yNode->tetVertex][y_face]] == dist;

    if (yNode->tetIndex == intersectTetIndex && yNode->tetVertex == intersectTetVertex) {
        intersectFace = yNode->adjTri[edgesFourToThree[yNode->tetVertex][y_face]];
    } else {
        intersectFace = TRUE;
    }

    return tetVertex && tetIndex && distTriVertex1 && distTriVertex2 && distTriVertex3 && intersectFace;
}

/*
 * Initialise Vertex
 *
 * Initialise vertex on the dual graph
 */

void init_vertex(struct CuspNode *xNode, struct CuspNode *yNode, struct CuspTriangle *x, struct CuspTriangle *y,
                 int x_vertex1, int x_vertex2, int y_vertex1, int y_vertex2, int y_face, int dist) {
    int i, cuspVertex;

    yNode->tri = y;
    yNode->tetIndex = y->tet->index;
    yNode->tetVertex = y->tetVertex;
    yNode->dist[edgesFourToThree[y->tetVertex][y_vertex1]] = xNode->dist[edgesFourToThree[x->tetVertex][x_vertex1]];
    yNode->dist[edgesFourToThree[y->tetVertex][y_vertex2]] = xNode->dist[edgesFourToThree[x->tetVertex][x_vertex2]];
    yNode->dist[edgesFourToThree[y->tetVertex][y_face]] = dist;

    // add faces which can be reached
    if (y->tetVertex == intersectTetVertex && yNode->tetIndex == intersectTetIndex) {
        // Vertex lies on triangle with intersection
        if (atleast_two(!yNode->dist[0], !yNode->dist[1], !yNode->dist[2])) {
            // Center vertex, add to all adj tri's
            yNode->adjTri[0] = 1;
            yNode->adjTri[1] = 1;
            yNode->adjTri[2] = 1;
        } else if (!yNode->dist[0] || !yNode->dist[1] || !yNode->dist[2]){
            // Can reach two edges
            for (i = 0; i < 3; i++) {
                cuspVertex = edgesThreeToFour[yNode->tetVertex][i];
                if (!yNode->dist[i]) {
                    yNode->adjTri[edgesFourToThree[yNode->tetVertex][remaining_face[yNode->tetVertex][cuspVertex]]] = 1;
                    yNode->adjTri[edgesFourToThree[yNode->tetVertex][remaining_face[cuspVertex][yNode->tetVertex]]] = 1;
                }
            }
        } else {
            // Can reach one edge
            yNode->adjTri[edgesFourToThree[yNode->tetVertex][y_face]] = 1;
        }
    } else {
        // Normal cusp triangle
        if (is_center_vertex(yNode)) {
            // Center vertex, add to all adj tri's
            yNode->adjTri[0] = 1;
            yNode->adjTri[1] = 1;
            yNode->adjTri[2] = 1;
        } else {
            // Can reach two edges
            for (i = 0; i < 3; i++) {
                cuspVertex = edgesThreeToFour[yNode->tetVertex][i];
                if (yNode->dist[i] <= flow(yNode->tri, cuspVertex)) {
                    yNode->adjTri[edgesFourToThree[yNode->tetVertex][remaining_face[yNode->tetVertex][cuspVertex]]] = 1;
                    yNode->adjTri[edgesFourToThree[yNode->tetVertex][remaining_face[cuspVertex][yNode->tetVertex]]] = 1;
                }
            }
        }
    }

}

/*
 * Is Center Vertex
 *
 * Check if a vertex can reach all three edges of the cusp triangle
 */

int is_center_vertex(struct CuspNode *y) {
    int dist1, dist2, dist3;

    dist1 = (flow(y->tri, edgesThreeToFour[y->tetVertex][0]) <= y->dist[0]);
    dist2 = (flow(y->tri, edgesThreeToFour[y->tetVertex][1]) <= y->dist[1]);
    dist3 = (flow(y->tri, edgesThreeToFour[y->tetVertex][2]) <= y->dist[2]);

    return dist1 && dist2 && dist3;
}

/*
 * Print Triangle Information
 *
 * Debug function for printing the gluing information
 */

void print_debug_info(struct CuspTriangle **pTriangle, struct Graph *g, struct CuspNode **vertexData,
        int **dualCurve, int *dualCurveLen, int flag) {
    int i, j, x_vertex1, x_vertex2, y_vertex1, y_vertex2;

    struct CuspTriangle *tri;
    struct EdgeNode *edge ;

    if (!flag) {
        // Gluing Info
        for (i = 0; i < numCuspTriangles; i++) {
            for (j = 0; j < 3; j++) {
                tri = pTriangle[i];
                x_vertex1 = (int) remaining_face[tri->tetVertex][edgesThreeToFour[tri->tetVertex][j]];
                x_vertex2 = (int) remaining_face[edgesThreeToFour[tri->tetVertex][j]][tri->tetVertex];
                y_vertex1 = EVALUATE(tri->tet->gluing[edgesThreeToFour[tri->tetVertex][j]], x_vertex1);
                y_vertex2 = EVALUATE(tri->tet->gluing[edgesThreeToFour[tri->tetVertex][j]], x_vertex2);

                printf("(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                       "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                       tri->tet->index,               // Tet Index
                       tri->tetVertex,                // Tet Vertex
                       edgesThreeToFour[tri->tetVertex][j],      // Cusp Edge
                       tri->tet->neighbor[edgesThreeToFour[tri->tetVertex][j]]->index,                              // Tet Index
                       EVALUATE(tri->tet->gluing[edgesThreeToFour[tri->tetVertex][j]], tri->tetVertex),             // Tet Vertex
                       EVALUATE(tri->tet->gluing[edgesThreeToFour[tri->tetVertex][j]], edgesThreeToFour[tri->tetVertex][j]),   // Cusp Edge
                       x_vertex1, y_vertex1,
                       x_vertex2, y_vertex2
                       );
            }
        }
    } else if (flag == 1) {
        // Vertex Distance Info
        for (i = 0; i < g->nvertices; i++) {
            if (vertexData[i]->tetIndex == -1 || vertexData[g->vertexData[i]]->tetVertex == -1) {
                tri = pTriangle[0];
            } else {
                tri = find_cusp_triangle(pTriangle, vertexData[i]->tetIndex, vertexData[g->vertexData[i]]->tetVertex);
            }

            printf("Tet Index: %d, Tet Vertex: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d.\n",
                   vertexData[g->vertexData[i]]->tetIndex,
                   vertexData[g->vertexData[i]]->tetVertex,
                   edgesThreeToFour[tri->tetVertex][0],
                   vertexData[g->vertexData[i]]->dist[0],
                   edgesThreeToFour[tri->tetVertex][1],
                   vertexData[g->vertexData[i]]->dist[1],
                   edgesThreeToFour[tri->tetVertex][2],
                   vertexData[g->vertexData[i]]->dist[2]
            );
        }
    } else if (flag == 2) {
        // Graph Info
        for (i = 0; i < g->nvertices; i++) {
            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Dist: %d, %d, %d): ",
                   i,
                   vertexData[g->vertexData[i]]->tetIndex,
                   vertexData[g->vertexData[i]]->tetVertex,
                   vertexData[g->vertexData[i]]->adjTri[0],
                   vertexData[g->vertexData[i]]->adjTri[1],
                   vertexData[g->vertexData[i]]->adjTri[2],
                   vertexData[g->vertexData[i]]->dist[0],
                   vertexData[g->vertexData[i]]->dist[1],
                   vertexData[g->vertexData[i]]->dist[2]
                   );
            edge = g->edges[i];

            if (edge == NULL) {
                printf("no edges\n");
            } else {
                printf("%d", edge->y);
                while ((edge = edge->next) != NULL) {
                    printf(", %d", edge->y);
                }
                printf("\n");
            }
        }
    } else if (flag == 3) {
        // Homology Info
        printf("Longitudinal\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d: %d, Cusp Edge %d: %d, Cusp Edge %d: %d\n",
                   tri->tet->index,
                   tri->tetVertex,
                   edgesThreeToFour[tri->tetVertex][0],
                   tri->tet->curve[0][right_handed][tri->tetVertex][edgesThreeToFour[tri->tetVertex][0]],
                   edgesThreeToFour[tri->tetVertex][1],
                   tri->tet->curve[0][right_handed][tri->tetVertex][edgesThreeToFour[tri->tetVertex][1]],
                   edgesThreeToFour[tri->tetVertex][2],
                   tri->tet->curve[0][right_handed][tri->tetVertex][edgesThreeToFour[tri->tetVertex][2]]
                );
        }
        printf("Meridian\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d: %d, Cusp Edge %d: %d, Cusp Edge %d: %d\n",
                   tri->tet->index,
                   tri->tetVertex,
                   edgesThreeToFour[tri->tetVertex][0],
                   tri->tet->curve[1][right_handed][tri->tetVertex][edgesThreeToFour[tri->tetVertex][0]],
                   edgesThreeToFour[tri->tetVertex][1],
                   tri->tet->curve[1][right_handed][tri->tetVertex][edgesThreeToFour[tri->tetVertex][1]],
                   edgesThreeToFour[tri->tetVertex][2],
                   tri->tet->curve[1][right_handed][tri->tetVertex][edgesThreeToFour[tri->tetVertex][2]]
                );
        }
    } else if (flag == 4) {
        // Edge indices
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) Vertex %d: (Edge Class: %d, Edge Index: %d, Face: %d), "
                   "Vertex %d: (Edge Class: %d, Edge Index: %d, Face: %d), Vertex %d: (Edge Class: %d, Edge Index: %d, Face: %d)\n",
                   tri->tet->index,
                   tri->tetVertex,
                   edgesThreeToFour[tri->tetVertex][0],
                   tri->vertices[0].edgeIndex,
                   tri->vertices[0].vertexIndex,
                   tri->faces[0].index,
                   edgesThreeToFour[tri->tetVertex][1],
                   tri->vertices[1].edgeIndex,
                   tri->vertices[1].vertexIndex,
                   tri->faces[1].index,
                   edgesThreeToFour[tri->tetVertex][2],
                   tri->vertices[2].edgeIndex,
                   tri->vertices[2].vertexIndex,
                   tri->faces[2].index
            );
        }
    } else if (flag == 5) {
        // Dual Curve Paths
        for (i = 0; i < numDualCurves; i++) {
            printf("Vertex %d: ", i);

            for (j = 0; j < dualCurveLen[i]; j++) {
                printf("%d ", dualCurve[i][j]);
            }
            printf("\n");
        }
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
 * Returns a pointer to the triangle with a given tetIndex and tetVertex.
 */

struct CuspTriangle *find_cusp_triangle(struct CuspTriangle **pTriangle, int tetIndex, int tetVertex) {
    int i;

    for (i = 0; i < numCuspTriangles; i++) {
        if (tetIndex == pTriangle[i]->tet->index && tetVertex == pTriangle[i]->tetVertex)
            return pTriangle[i];
    }

    return NULL;
}

/*
 * Construct Dual Curves
 */

struct Graph *construct_dual_curves(struct Graph *g, struct CuspTriangle **pTriangle, struct CuspNode **vertexData,
        int e0, int **dualCurves, int *dualCurveLen) {
    int i, startIndex, endIndex, startDualIndex, endDualIndex, pathLen;
    bool *processed = NEW_ARRAY(g->nvertices, bool);
    bool *discovered = NEW_ARRAY(g->nvertices, bool);
    int *parent = NEW_ARRAY(g->nvertices, int);
    int *path = NEW_ARRAY(g->nvertices, int);

    for (i = 0; i < numDualCurves; i++) {
        /*
         * Remove edges associated to diving through the manifold,
         * store the start and end edges for each path in the endpoints
         * array, and store the edges removed in the edges array.
         */
        if (i == e0)
            continue;

        //remove_extra_edges(g, i, e0);

        find_index(g, vertexData, i, &startIndex, &startDualIndex);
        find_index(g, vertexData, e0, &endIndex, &endDualIndex);

        // Find path using bfs
        bfs(g, startIndex, processed, discovered, parent);
        find_path(startIndex, endIndex, parent, path, 0, &pathLen);
        dualCurves[2 * i] = path;
        dualCurveLen[2 * i] = pathLen;

        // Split graph
        g = split_along_path(g, path, pathLen);

        // Reallocate memory
        my_free(processed);
        my_free(discovered);
        my_free(parent);
        return g;
        processed = NEW_ARRAY(g->nvertices, bool);
        discovered = NEW_ARRAY(g->nvertices, bool);
        parent = NEW_ARRAY(g->nvertices, int);
        path = NEW_ARRAY(g->nvertices, int);

        // Repeat for the other half of the curve
        bfs(g, startDualIndex, processed, discovered, parent);
        find_path(startDualIndex, endDualIndex, parent, path, 0, &pathLen);
        dualCurves[2 * i + 1] = path;
        dualCurveLen[2 * i + 1] = pathLen;

        // Split graph
        g = split_along_path(g, path, pathLen);

        //add_edges_from_array(g, edges);
        // Re allocate memory
        my_free(processed);
        my_free(discovered);
        my_free(parent);
        processed = NEW_ARRAY(g->nvertices, bool);
        discovered = NEW_ARRAY(g->nvertices, bool);
        parent = NEW_ARRAY(g->nvertices, int);
        path = NEW_ARRAY(g->nvertices, int);
    }

    return g;
}

/*
 * Find index
 *
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass.
 */

void find_index(struct Graph *g, struct CuspNode **vertexData, int edgeClass, int *index, int *dualIndex) {
    int i, j, cuspVertex;
    bool found = FALSE;

    for (i = 0; i < g->nvertices; i++) {
        for (j = 0; j < 3; j++) {
            if (vertexData[g->vertexData[i]]->tri->vertices[j].edgeIndex == edgeClass &&
                vertexData[g->vertexData[i]]->tri->vertices[j].vertexIndex == 0 &&
                vertexData[g->vertexData[i]]->dist[j] == 0) {
                *index = i;
                cuspVertex = edgesThreeToFour[vertexData[g->vertexData[i]]->tetVertex][j];
                found = TRUE;
                break;
            }
        }

        if (found)
            break;
    }

    for (i = 0; i < g->nvertices; i++) {
        for (j = 0; j < 3; j++) {
            if (vertexData[g->vertexData[i]]->tri->vertices[j].edgeIndex == edgeClass &&
                vertexData[g->vertexData[i]]->tri->vertices[j].vertexIndex == 1 &&
                vertexData[g->vertexData[i]]->tetVertex == cuspVertex &&
                vertexData[g->vertexData[i]]->tetIndex == vertexData[g->vertexData[*index]]->tetIndex &&
                edgesThreeToFour[vertexData[g->vertexData[i]]->tetVertex][j] == vertexData[g->vertexData[*index]]->tetVertex) {
                *dualIndex = i;
                return;
            }
        }
    }
}

/* Remove Extra Edges
 *
 * Remove edges associated to the curves that dive into the
 * manifold since we can only do this once.
 */

void remove_extra_edges(struct Graph *g, struct CuspNode **vertexData, int **retval, int startEdgeClass, int endEdgeClass) {
    int i, j, startCuspVertex, endCuspVertex, startIndex, endIndex, startDualIndex, endDualIndex;

    // Find the indices of the start and end cusp vertices
    for (i = 0; i < g->nvertices; i++) {
        for (j = 0; j < 3; j++) {
            if (vertexData[g->vertexData[i]]->tri->vertices[j].edgeIndex == startEdgeClass &&
            vertexData[g->vertexData[i]]->tri->vertices[j].vertexIndex == 0) {
                startIndex = i;
                startCuspVertex = j;
            }

            if (vertexData[g->vertexData[i]]->tri->vertices[j].edgeIndex == endEdgeClass &&
            vertexData[g->vertexData[i]]->tri->vertices[j].vertexIndex == 0) {
                endIndex = i;
                endCuspVertex = j;
            }
        }
    }

    // Find the corresponding start and end vertices for the second half of the dual curve
    for (i = 0; i < g->nvertices; i++) {
        for (j = 0; j < 3; j++) {
            if (vertexData[g->vertexData[i]]->tri->vertices[j].edgeIndex == startEdgeClass &&
            vertexData[g->vertexData[i]]->tri->vertices[j].vertexIndex == 1 &&
            vertexData[g->vertexData[i]]->tetVertex == startCuspVertex &&
            j == vertexData[g->vertexData[startIndex]]->tetVertex) {
                startDualIndex = i;
            }

            if (vertexData[g->vertexData[i]]->tri->vertices[j].edgeIndex == endEdgeClass &&
                vertexData[g->vertexData[i]]->tri->vertices[j].vertexIndex == 1 &&
                vertexData[g->vertexData[i]]->tetVertex == endCuspVertex &&
                j == vertexData[g->vertexData[endIndex]]->tetVertex) {
                endDualIndex = j;
            }
        }
    }

    /*
     * The dual cusp triangle is linked to the start one as follows:
     * tri->tetVertex = dualTriCuspVertex
     * triCuspVertex = dualTri->tetVertex
     */

    int startVertex1, startVertex2, endVertex1, endVertex2;

    startVertex1 = (int) remaining_face[vertexData[g->vertexData[startIndex]]->tetVertex]
            [edgesThreeToFour[vertexData[g->vertexData[startIndex]]->tetVertex][startCuspVertex]];
    startVertex1 = (int) remaining_face[edgesThreeToFour[vertexData[g->vertexData[startIndex]]->tetVertex][startCuspVertex]]
            [vertexData[g->vertexData[startIndex]]->tetVertex];
    endVertex1 = (int) remaining_face[vertexData[g->vertexData[endIndex]]->tetVertex]
            [edgesThreeToFour[vertexData[g->vertexData[endIndex]]->tetVertex][endCuspVertex]];
    endVertex1 = (int) remaining_face[edgesThreeToFour[vertexData[g->vertexData[endIndex]]->tetVertex][endCuspVertex]]
            [vertexData[g->vertexData[endIndex]]->tetVertex];
}

/* Add Edges
 *
 * Add the edges from an array back to the graph
 * Used to undo remove extra edges.
 */

void add_edges_from_array(struct Graph *graph1, int **edges, int numEdges) {
    int i;

    for (i = 0; i < numEdges; i++) {
        insert_edge(graph1, edges[i][0], edges[i][1], graph1->directed);
    }
}

/*
 * Find Holonomies
 *
 * Construct the symplectic equations from the dual curves
 */

void find_holonomies(struct Graph *g, int **symp_eqns, int **dualCurves, int *dualCurveLen) {
    int i;

    for (i = 0; i < numDualCurves; i++) {
        find_path_holonomy(g, symp_eqns[i], dualCurves[2 * i], dualCurveLen[2 * i]);
        find_path_holonomy(g, symp_eqns[i], dualCurves[2 * i + 1], dualCurveLen[2 * i + 1]);
    }
}

/*
 * Find Path Holonomies
 *
 * Calculate holonomies along a path and update the row.
 */

void find_path_holonomy(struct Graph *g, int *row, int *path, int pathLen) {

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

// ---------------------------------------------------------

// Graph

/*
 * Initialise Graph
 *
 * Initialise the arrays of the graph 'g' to their default values
 */

struct Graph *init_graph(int maxVertices, bool directed) {
    int i;
    struct Graph *g = malloc(sizeof( struct Graph ));

    g->nvertices = maxVertices;
    g->nedges = 0;
    g->directed = directed;

    g->edges = NEW_ARRAY(maxVertices, struct EdgeNode *);
    g->degree = NEW_ARRAY(maxVertices, int);
    g->vertexData = NEW_ARRAY(maxVertices, int);

    for (i = 0; i < maxVertices; i++) {
        g->degree[i] = 0;
        g->edges[i] = NULL;
        g->vertexData[i] = i;
    }

    return g;
}

/*
 * Free Graph
 *
 * Free the memory of the graph data structure
 */

void free_graph(struct Graph *g) {
    int i;
    struct EdgeNode *edge;

    for (i = 0; i < g->nvertices; i++) {
        while (g->edges[i] != NULL) {
            if (g->edges[i]->next == NULL) {
                free(g->edges[i]);
                g->edges[i] = NULL;
            } else {
                edge = g->edges[i]->next;
                g->edges[i]->next = edge->next;
                free(edge);
            }
        }
    }

    my_free(g->edges);
    my_free(g->degree);
    my_free(g->vertexData);
    free(g);
}

/*
 * Insert Edge
 *
 * Insert an edge into the graph 'g' from vertex x to y.
 */

int insert_edge(struct Graph *g, int x, int y, bool directed) {
    struct EdgeNode *p;

    p = malloc(sizeof( struct EdgeNode ));
    p->y = y;
    p->next = g->edges[x];

    g->edges[x] = p;
    g->degree[x]++;

    if (!directed) {
        insert_edge(g, y, x, TRUE);
    }

    return x;
}

/*
 * Delete Edge
 *
 * Remove the edge from vertex x to vertex y
 */

void delete_edge(struct Graph *g, int vertex_x, int vertex_y, bool directed) {
    struct EdgeNode *node, *deleted_node;

    node = g->edges[vertex_x];

    if (node == NULL)
        return;

    if (node->y == vertex_y) {
        g->edges[vertex_x] = node->next;
    }

    while (node->next != NULL && node->next->y != vertex_y)
        node = node->next;

    if (node == NULL)
        return;

    deleted_node = node->next;
    node->next = node->next->next;
    free(deleted_node);

    if (!directed) {
        delete_edge(g, vertex_y, vertex_x, directed);
    }
}

/*
 * Edge Exists
 *
 * Check if an edge already exists in the graph
 */

int edge_exists(struct Graph *g, int v1, int v2) {
    struct EdgeNode *node = g->edges[v1];

    if (node == NULL) {
        return FALSE;
    }

    while (node != NULL) {
        if (node->y == v2) {
            return TRUE;
        }

        node = node->next;
    }

    return FALSE;
}

// Graph Splitting

/*
 * Split Along Path
 *
 * Given a graph graph1, add to each vertex 3 valence vertices and then join the
 * valence vertices such that there does not exist an edge crossing the path.
 * We identify valence vertices by the vertex of the cusp triangle they are
 * adjacent to.
 */

struct Graph *split_along_path(struct Graph *graph1, int *path, int pathLen) {
    struct Graph *graph2 = init_graph(graph1->nvertices + pathLen, graph1->directed);

    init_vertices(graph1, graph2, path, pathLen);
    add_non_path_edges(graph1, graph2, path, pathLen);
    add_path_edges(graph1, graph2, path, pathLen);

    free_graph(graph1);
    return graph2;
}

/*
 * Initialise Vertices
 *
 * Set the vertex data for the valence vertices.
 */

void init_vertices(struct Graph *graph1, struct Graph *graph2, int *array, int arrayLen) {
    int i, j;

    for (i = 0; i < graph1->nvertices; i++) {
        if (inclusion(array, arrayLen, i)) {
            graph2->vertexData[graph1->nvertices + i] = graph1->vertexData[i];
        }

        graph2->vertexData[i] = graph1->vertexData[i];
    }
}

/*
 * Add Non Path Edges
 *
 * For vertices not in the path, each of the valence vertices connects
 * to the central vertex and to the valence vertices of adjacent vertices
 */

void add_non_path_edges(struct Graph *graph1, struct Graph *graph2, int *path, int pathLen) {
    int i, j;
    struct EdgeNode *edge;

    for (i = 0; i < graph1->nvertices; i++) {
        edge = graph1->edges[i];

        if (inclusion(path, pathLen, i))
            continue;

        while (edge != NULL) {
            if (!inclusion(path, pathLen, edge->y))
                insert_edge(graph2, i, edge->y, graph2->directed);

            edge = edge->next;
        }
    }
}

/*
 * Add Path Edges
 *
 * For vertices on the path, we consider the two cases of the vertex connecting to
 * two vertices on the path (in the middle) or connecting to one vertex (end points).
 * In the latter case we treat them the same as vertices not on the path.
 */

void add_path_edges(struct Graph *graph1, struct Graph *graph2, int *path, int pathLen) {
    int i;

    for (i = 0; i < pathLen; i++) {
        if (i == 0 || i == pathLen - 1) {
            // Case 1: deg 1 in the path

        } else {
            // Case 2: deg 2 in the path

        }
    }
}

/*
 * Inclusion
 *
 * Check if the target is in the array
 */

bool inclusion(int *array, int arrayLen, int target) {
    int i;

    for (i = 0; i < arrayLen; i++) {
        if (array[i] == target) {
            return TRUE;
        }
    }

    return FALSE;
}

struct CuspNode *duplicate(struct CuspNode *node) {
    int i;
    struct CuspNode *newNode = NEW_STRUCT(struct CuspNode);

    newNode->tri = node->tri;
    newNode->tetIndex = node->tetIndex;
    newNode->tetVertex = node->tetVertex;

    for (i = 0; i < 3; i++) {
        newNode->dist[i] = node->dist[i];
        newNode->adjTri[i] = node->dist[i];
    }

    return newNode;
}


// ---------------------------------------------------------------

// Breadth First Search

/*
 * Initialise Search
 *
 * Initialise default values for bfs arrays
 */

void init_search(struct Graph *g, bool *processed, bool *discovered, int *parent) {
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

void bfs(struct Graph *g, int start, bool *processed, bool *discovered, int *parent) {
    struct Queue q;
    int v, y;
    struct EdgeNode *p;

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
void find_path(int start, int end, int *parents, int *path, int index, int *pathLen) {
    if ((start == end) || (end == -1)) {
        path[index] = start;
        *pathLen = index + 1;
    } else {
        find_path(start, parents[end], parents, path, index + 1, pathLen);
        path[index] = end;
    }
}
