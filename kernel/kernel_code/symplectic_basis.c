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
    int i;
    numCuspTriangles = 4 * manifold->num_tetrahedra;

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
    struct CuspTriangle **pTriangle = NEW_ARRAY(numCuspTriangles, struct CuspTriangle *);

    for (i = 0; i < numCuspTriangles; i++)
        pTriangle[i] = malloc(sizeof( struct CuspTriangle ));

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, pTriangle, symp_num_rows, symp_eqns);

    // Construct return array
    int **eqns = NEW_ARRAY(edge_num_rows + symp_num_rows, int *);

    for (i = 0; i < edge_num_rows; i++) {
        eqns[2 * i] = edge_eqns[i];
        eqns[2 * i + 1] = symp_eqns[i];
    }

    free_cusp_triangulation(pTriangle);
    my_free(symp_eqns);
    my_free(edge_eqns);
    *num_rows = edge_num_rows + symp_num_rows;
    return eqns;
}

/*
 * Initialise Cusp Triangulation
 *
 * Construct the pTriangle array which consists of the triangles in the cusp triangulation
 */

static int numEdgeClasses;

void init_cusp_triangulation(Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int i, j;
    Tetrahedron *tet= manifold->tet_list_begin.next;

    label_triangulation_edges(manifold);

    for (i = 0; i < numCuspTriangles; i++) {
        pTriangle[i]->tet = tet;
        pTriangle[i]->tetVertex = i % 4;

        for (j = 0; j < 4; j++) {
            pTriangle[i]->neighbours[j] = pTriangle[4 * (pTriangle[i]->tet->neighbor[j]->index)
                                                    + EVALUATE(pTriangle[i]->tet->gluing[j], pTriangle[i]->tetVertex)];
        }

        for (j = 0; j < 3; j++) {
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
    print_debug_info(pTriangle, NULL, 4);
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

int **get_symplectic_equations(Triangulation *manifold, struct CuspTriangle **pTriangle, int num_rows, int **eqns) {
    int i, j, numVertices, T = manifold -> num_tetrahedra;
    genus = (2 * manifold->num_tetrahedra - num_rows) / 2;
    struct graph graph1;

    init_cusp_triangulation(manifold, pTriangle);
    numVertices = (2 * genus + 1) * numCuspTriangles + 4 * genus;
    init_graph(&graph1, numVertices, FALSE);
    construct_dual_graph(&graph1, pTriangle);
    print_debug_info(pTriangle, &graph1, 0);
    print_debug_info(pTriangle, &graph1, 2);
    print_debug_info(pTriangle, &graph1, 3);

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
static int intersectTetIndex, intersectTetVertex;

void construct_dual_graph(struct graph *graph1, struct CuspTriangle **pTriangle) {
    int i, index, cuspEdge, directed = FALSE;
    int indices[3] = { 0, 0, 0 };
    struct CuspTriangle *tri, *adjTri[3];

    intersectTetIndex = 0;
    intersectTetVertex = 2;

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
    graph1->vertexHomology[0][3] = flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][0])
            + flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][1]);                       // dist to v2
    graph1->vertexHomology[0][4] = flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][0])
            + flow(pTriangle[0], edgesThreeToFour[pTriangle[0]->tetVertex][2]);                       // dist to v3

    // Walk around the cusp triangulation inserting edges
    while (!is_empty(&stack)) {
        index = pop(&stack);

        if (visited[index])
            continue;

        // Find the triangle that vertex (index) of the dual graph lies on.
        tri = find_cusp_triangle(pTriangle, graph1->vertexHomology[index][0],
                                 graph1->vertexHomology[index][1]);

        if (tri->tet->index == intersectTetIndex && tri->tetVertex == intersectTetVertex) {
            // 2 vertices of dist 0
            if ((!graph1->vertexHomology[index][2] && !graph1->vertexHomology[index][3]) ||
                    (!graph1->vertexHomology[index][3] && !graph1->vertexHomology[index][4]) ||
                    (!graph1->vertexHomology[index][2] && !graph1->vertexHomology[index][4])) {
                // Insert to all 3 edges
                for (i = 0; i < 3; i++) {
                    cuspEdge = edgesThreeToFour[tri->tetVertex][i];
                    adjTri[i] = tri->neighbours[cuspEdge];
                    indices[i] = insert_triangle_edge(graph1, index, cuspEdge, tri, adjTri[i], directed);
                }
            } else {
                // 1 vertex of dist 0
                for (i = 0; i < 3; i++) {
                    if (!graph1->vertexHomology[index][i + 2]) {
                        cuspEdge = (int) remaining_face[edgesThreeToFour[tri->tetVertex][i]][tri->tetVertex];
                        adjTri[0] = tri->neighbours[cuspEdge];
                        indices[0] = insert_triangle_edge(graph1, index, cuspEdge, tri, adjTri[0], directed);

                        cuspEdge = (int) remaining_face[tri->tetVertex][edgesThreeToFour[tri->tetVertex][i]];
                        adjTri[1] = tri->neighbours[cuspEdge];
                        indices[1] = insert_triangle_edge(graph1, index, cuspEdge, tri, adjTri[1], directed);
                    }
                }
            }
        } else {
            // Normal Triangle
            if (is_center_vertex(graph1, tri, index)) {
                // Vertex lies in the center, so we can add vertices across all edges
                for (i = 0; i < 3; i++) {
                    cuspEdge = edgesThreeToFour[tri->tetVertex][i];
                    adjTri[i] = tri->neighbours[cuspEdge];
                    indices[i] = insert_triangle_edge(graph1, index, cuspEdge, tri, adjTri[i], directed);
                }
            } else {
                // Add two vertices
                for (i = 0; i < 3; i++) {
                    if (graph1->vertexHomology[index][i + 2] <= flow(tri, edgesThreeToFour[tri->tetVertex][i])) {
                        cuspEdge = (int) remaining_face[edgesThreeToFour[tri->tetVertex][i]][tri->tetVertex];
                        adjTri[0] = tri->neighbours[cuspEdge];
                        indices[0] = insert_triangle_edge(graph1, index, cuspEdge, tri, adjTri[0], directed);

                        cuspEdge = (int) remaining_face[tri->tetVertex][edgesThreeToFour[tri->tetVertex][i]];
                        adjTri[1] = tri->neighbours[cuspEdge];
                        indices[1] = insert_triangle_edge(graph1, index, cuspEdge, tri, adjTri[1], directed);
                    }
                }
            }
        }

        for (i = 0; i < 3; i++) {
            if (indices[i] == -1) {
                printf("Inserting edge from triangle (%d, %d) to triangle (%d, %d) failed\n",
                       tri->tet->index,
                       tri->tetVertex,
                       adjTri[i]->tet->index,
                       adjTri[i]->tetVertex
                );
                print_debug_info(pTriangle, graph1, 1);
            }

            if (!visited[indices[i]])
                push(&stack, indices[i]);
        }

        visited[index] = 1;
    }
    
    remove_extra_edges(graph1);
    add_misc_edges(graph1);
}

int insert_triangle_edge(struct graph *g, int index, int face, struct CuspTriangle *x, struct CuspTriangle *y, bool directed) {
    int i, x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist;

    // Vertices on triangle x which are glued to triangle y
    x_vertex1 = (int) remaining_face[x->tetVertex][face];
    x_vertex2 = (int) remaining_face[face][x->tetVertex];

    // Vertices on triangle y which correspond to the vertices v1 and v2
    y_vertex1 = EVALUATE(x->tet->gluing[face], x_vertex1);
    y_vertex2 = EVALUATE(x->tet->gluing[face], x_vertex2);
    y_face = EVALUATE(x->tet->gluing[face], face);

    dist = flow(y, y_face);
    if (g->vertexHomology[index][edgesFourToThree[x->tetVertex][x_vertex1] + 2] < flow(y, y_vertex1)) {
        dist += (flow(y, y_vertex1) - g->vertexHomology[index][edgesFourToThree[x->tetVertex][x_vertex1] + 2]);
    } else {
        dist += (flow(y, y_vertex2) - g->vertexHomology[index][edgesFourToThree[x->tetVertex][x_vertex2] + 2]);
    }


    // Find index of y
    for (i = 0; i < g->nvertices; i++) {
        if (is_equal(g->vertexHomology[index], g->vertexHomology[i], x->tetVertex,
                     y->tet->index, y->tetVertex, x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist)) {
            break;
        }
    }

    // Insert new vertex
    if (i == g->nvertices) {
        for (i = 0; i < g->nvertices && g->vertexHomology[i][0] != -1; i++);

        if (i == g->nvertices)
            return -1;

        init_vertex(g->vertexHomology[index], g->vertexHomology[i], x->tetVertex,
                    y->tet->index, y->tetVertex, x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist);
    }

    // Ignore edge if it already exists
    if (!edge_exists(g, index, i))
        insert_edge(g, index, i, directed);

    return i;
}

/*
 * Is Equal
 *
 * Check if two holonomies correspond to the same vertex
 */

int is_equal(int *holonomyX, int *holonomyY, int x_tet_vertex, int y_tet_index, int y_tet_vertex,
             int x_vertex1, int x_vertex2, int y_vertex1, int y_vertex2, int y_face, int dist) {
    int i, tetIndex, tetVertex, distTriVertex1, distTriVertex2, distTriVertex3, intersectFace, other_face;
    bool found = FALSE;

    tetIndex = holonomyY[0] == y_tet_index;
    tetVertex = holonomyY[1] == y_tet_vertex;

    distTriVertex1 = (holonomyY[edgesFourToThree[y_tet_vertex][y_vertex1] + 2] ==
                      holonomyX[edgesFourToThree[x_tet_vertex][x_vertex1] + 2]);
    distTriVertex2 = (holonomyY[edgesFourToThree[y_tet_vertex][y_vertex2] + 2] ==
                      holonomyX[edgesFourToThree[x_tet_vertex][x_vertex2] + 2]);
    distTriVertex3 = holonomyY[edgesFourToThree[y_tet_vertex][y_face] + 2] == dist;

    if (y_tet_index == intersectTetIndex && y_tet_vertex == intersectTetVertex
    && !((!holonomyY[2] && !holonomyY[3]) || (!holonomyY[2] && !holonomyY[4]) || (!holonomyY[3] && !holonomyY[4]))) {
        // Check if dual graph vertex lies dist. 0 from a cusp vertex.
        for (i = 0; i< 3; i++) {
            if (y_face == remaining_face[y_tet_vertex][edgesThreeToFour[y_tet_vertex][i]]) {
                other_face = (int) remaining_face[edgesThreeToFour[y_tet_vertex][i]][y_tet_vertex];
            } else {
                other_face = (int) remaining_face[y_tet_vertex][edgesThreeToFour[y_tet_vertex][i]];
            }

            if (holonomyY[i + 2] == 0) {
                intersectFace = (holonomyY[5] == y_face || holonomyY[5] == other_face);
                found = TRUE;
            }
        }

        if (!found) {
            intersectFace = (holonomyY[5] == y_face);
        }
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

void init_vertex(int *holonomyX, int *holonomyY, int x_tet_vertex, int y_tet_index, int y_tet_vertex,
                 int x_vertex1, int x_vertex2, int y_vertex1, int y_vertex2, int y_face, int dist) {
    holonomyY[0] = y_tet_index;
    holonomyY[1] = y_tet_vertex;

    holonomyY[edgesFourToThree[y_tet_vertex][y_vertex1] + 2] = holonomyX[edgesFourToThree[x_tet_vertex][x_vertex1] + 2];
    holonomyY[edgesFourToThree[y_tet_vertex][y_vertex2] + 2] = holonomyX[edgesFourToThree[x_tet_vertex][x_vertex2] + 2];
    holonomyY[edgesFourToThree[y_tet_vertex][y_face] + 2] = dist;

    // Handle intersection vertex
    if (y_tet_index == intersectTetIndex && y_tet_vertex == intersectTetVertex) {
        holonomyY[5] = y_face;
    }
}

/*
 * Is Center Vertex
 *
 * Check if a vertex can reach all three edges of the cusp triangle
 */

int is_center_vertex(struct graph *g, struct CuspTriangle *tri, int index) {
    int dist1, dist2, dist3;

    dist1 = (flow(tri, edgesThreeToFour[tri->tetVertex][0]) <= g->vertexHomology[index][2]);
    dist2 = (flow(tri, edgesThreeToFour[tri->tetVertex][1]) <= g->vertexHomology[index][3]);
    dist3 = (flow(tri, edgesThreeToFour[tri->tetVertex][2]) <= g->vertexHomology[index][4]);

    return dist1 && dist2 && dist3;
}

/*
 * Print Triangle Information
 *
 * Debug function for printing the gluing information
 */

void print_debug_info(struct CuspTriangle **pTriangle, struct graph *g, int flag) {
    int i, j, x_vertex1, x_vertex2, y_vertex1, y_vertex2;

    struct CuspTriangle *tri;
    struct edgenode *edge ;

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
            if (g->vertexHomology[i][0] == -1 || g->vertexHomology[i][1] == -1) {
                tri = pTriangle[0];
            } else {
                tri = find_cusp_triangle(pTriangle, g->vertexHomology[i][0], g->vertexHomology[i][1]);
            }

            printf("Tet Index: %d, Tet Vertex: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d.\n",
                   g->vertexHomology[i][0],
                   g->vertexHomology[i][1],
                   edgesThreeToFour[tri->tetVertex][0],
                   g->vertexHomology[i][2],
                   edgesThreeToFour[tri->tetVertex][1],
                   g->vertexHomology[i][3],
                   edgesThreeToFour[tri->tetVertex][2],
                   g->vertexHomology[i][4]
            );
        }
    } else if (flag == 2) {
        // Graph Info
        for (i = 0; i < g->nvertices; i++) {
            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d) (Dist: %d, %d, %d, %d): ",
                   i,
                   g->vertexHomology[i][0],
                   g->vertexHomology[i][1],
                   g->vertexHomology[i][2],
                   g->vertexHomology[i][3],
                   g->vertexHomology[i][4],
                   g->vertexHomology[i][5]
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
            printf("(Tet Index: %d, Tet Vertex: %d) Cusp Vertex %d: (%d, %d), Cusp Vertex %d: (%d, %d), Cusp Vertex %d: (%d, %d)\n",
                   tri->tet->index,
                   tri->tetVertex,
                   edgesThreeToFour[tri->tetVertex][0],
                   tri->vertices[0].edgeIndex,
                   tri->vertices[0].vertexIndex,
                   edgesThreeToFour[tri->tetVertex][1],
                   tri->vertices[1].edgeIndex,
                   tri->vertices[1].vertexIndex,
                   edgesThreeToFour[tri->tetVertex][2],
                   tri->vertices[2].edgeIndex,
                   tri->vertices[2].vertexIndex
            );
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
 * Returns a pointer to the triangle with tetIndex and tetVertex.
 */

struct CuspTriangle *find_cusp_triangle(struct CuspTriangle **pTriangle, int tetIndex, int tetVertex) {
    int i;

    for (i = 0; i < numCuspTriangles; i++) {
        if (tetIndex == pTriangle[i]->tet->index && tetVertex == pTriangle[i]->tetVertex)
            return pTriangle[i];
    }

    return NULL;
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

void init_graph(struct graph *g, int maxVertices, bool directed) {
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
/*
 * Free Graph
 *
 * Free the memory of the graph data structure
 */

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
 * Insert an edge into the graph 'g' from vertex x to y.
 */

int insert_edge(struct graph *g, int x, int y, bool directed) {
    edgenode *p;

    p = malloc(sizeof( edgenode ));
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

void delete_edge(struct graph *g, int vertex_x, int vertex_y, bool directed) {
    struct edgenode *node;

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

    node->next = node->next->next;

    if (!directed) {
        delete_edge(g, vertex_y, vertex_x, directed);
    }
}

/*
 * Edge Exists
 *
 * Check if an edge already exists in the graph
 */

int edge_exists(struct graph *g, int v1, int v2) {
    struct edgenode *node = g->edges[v1];

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

/*
 * Initialise Search
 *
 * Initialise default values for bfs arrays
 */

void init_search(struct graph *g, bool *processed, bool *discovered, int *parent) {
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

