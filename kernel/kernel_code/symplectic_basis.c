/*
 *  Symplectic Basis
 *
 *  Computes the symplectic basis of a cusped 3-manifold
*/

#include <stdbool.h>
#include <stdio.h>
#include "kernel.h"
#include "kernel_namespace.h"
#include "symplectic_basis.h"
//#include "addl_code.h"                  // compile snappy
#include "../addl_code/addl_code.h"     // compile kernel

#define atleast_two(a, b, c)    ((a) && (b)) || ((a) && (c)) || ((b) && (c))

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
    int i, j, k;
    Tetrahedron *tet= manifold->tet_list_begin.next;

    // Allocate Cusp Triangulation Array
    struct CuspTriangle **pTriangle = NEW_ARRAY(numCuspTriangles, struct CuspTriangle *);

    for (i = 0; i < numCuspTriangles; i++)
        pTriangle[i] = NEW_STRUCT(struct CuspTriangle);

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

        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                pTriangle[i]->orientVertices[j][k] = 0;
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

        pTriangle[i]->oriented = FALSE;
    }

    cusp_vertex_index(pTriangle);
    vertex_orientation(pTriangle);
    //label_cusp_faces(pTriangle);
    //print_debug_info(pTriangle, NULL, NULL, NULL, NULL, 4);
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
    int i, j;

    int *currentIndex = NEW_ARRAY(numEdgeClasses, int);
    
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

    my_free(currentIndex);
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

/*
 * Vertex Orientation
 *
 *
 */

void vertex_orientation(struct CuspTriangle **pTriangle) {
    int i, j, index, adjIndex, face, vertex;
    struct Node indexStack, adjIndexStack, faceStack;
    indexStack.item = -1;
    adjIndexStack.item = -1;
    faceStack.item = -1;

    int zeroTri[4][4] = {
            {0, 0, 0, 0},
            {0, 0, -1, 1},
            {0, 1, 0, -1},
            {0, -1, 1, 0}};

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            pTriangle[0]->orientVertices[i][j] = zeroTri[i][j];

    for (i = 1; i < 4; i++) {
        index = find_cusp_triangle_index(pTriangle, pTriangle[0]->neighbours[i]->tet->index, pTriangle[0]->neighbours[i]->tetVertex);
        push(&indexStack, 0);
        push(&adjIndexStack, index);
        push(&faceStack, i);
    }

    pTriangle[0]->oriented = TRUE;

    while (!is_empty(&indexStack)) {
        index = pop(&indexStack);
        adjIndex = pop(&adjIndexStack);
        face = pop(&faceStack);

        update_orientation(pTriangle[index], pTriangle[adjIndex], face);
        pTriangle[index]->oriented = TRUE;

        for (i = 0; i < 3; i++) {
            vertex = edgesThreeToFour[pTriangle[adjIndex]->tetVertex][i];

            if (!pTriangle[adjIndex]->neighbours[vertex]->oriented) {
                push(&indexStack, adjIndex);
                push(&adjIndexStack, find_cusp_triangle_index(pTriangle,
                                                              pTriangle[adjIndex]->neighbours[vertex]->tet->index,
                                                              pTriangle[adjIndex]->neighbours[vertex]->tetVertex));
                push(&faceStack, vertex);
            }
        }
//        printf("Tri: (%d, %d), adj tri: (%d, %d), face: %d\n", pTriangle[index]->tet->index, pTriangle[index]->tetVertex,
//               pTriangle[adjIndex]->tet->index, pTriangle[adjIndex]->tetVertex, face);
//        print_debug_info(pTriangle, NULL, NULL, NULL, NULL, 6);
    }
}

/*
 * Update Orientation
 *
 * Cusp Triangle x glues to cusp triangle y across edge 'face' of x. Assume
 * x has an orientation, update the orientation on y to be consistent with x.
 */

void update_orientation(struct CuspTriangle *x, struct CuspTriangle *y, int face) {
    int x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face;

    // Vertices on triangle x which are glued to triangle y
    x_vertex1 = (int) remaining_face[x->tetVertex][face];
    x_vertex2 = (int) remaining_face[face][x->tetVertex];

    // Vertices on triangle y which correspond to the vertices v1 and v2
    y_vertex1 = EVALUATE(x->tet->gluing[face], x_vertex1);
    y_vertex2 = EVALUATE(x->tet->gluing[face], x_vertex2);
    y_face = EVALUATE(x->tet->gluing[face], face);

    y->orientVertices[y_vertex1][y_vertex2] = x->orientVertices[x_vertex1][face];
    y->orientVertices[y_vertex1][y_face] = x->orientVertices[x_vertex1][x_vertex2];

    y->orientVertices[y_vertex2][y_vertex1] = x->orientVertices[x_vertex2][face];
    y->orientVertices[y_vertex2][y_face] = x->orientVertices[x_vertex2][x_vertex1];

    y->orientVertices[y_face][y_vertex1] = y->orientVertices[y_vertex1][y_vertex2];
    y->orientVertices[y_face][y_vertex2] = y->orientVertices[y_vertex2][y_vertex1];
}

// Free memory for cusp triangulation data structure
void free_cusp_triangulation(struct CuspTriangle **pTriangle) {
    int i;

    for (i = 0; i < numCuspTriangles; i++)
        my_free(pTriangle[i]);

    my_free(pTriangle);
}

/*
 * Get Symplectic Equations
 *
 * Setup graph and cusp triangulation, and run construct dual graph.
 */
static int numDualCurves, genus;
static int intersectTetIndex, intersectTetVertex;


int **get_symplectic_equations(Triangulation *manifold, int num_rows, int numCols, int e0) {
    int i, j, numCuspRegions, T = manifold -> num_tetrahedra, face0;
    numDualCurves = num_rows;
    genus = (2 * manifold->num_tetrahedra - num_rows) / 2;
    struct Graph *graph1;
    struct CuspTriangle **pTriangle;

    // Call fundamental group to find standard basepoint of the cusp triangulation
    fundamental_group(manifold, FALSE, FALSE, FALSE, FALSE);
    intersectTetIndex = manifold->cusp_list_begin.next->basepoint_tet->index;
    intersectTetVertex = (int) manifold->cusp_list_begin.next->basepoint_vertex;

    int **dualCurve = NEW_ARRAY(2 * numDualCurves, int *);
    int *dualCurveLen = NEW_ARRAY(2 * numDualCurves, int);
    struct PathEndPoint **pPathEndPoint = NEW_ARRAY(2 * numDualCurves, struct PathEndPoint *);

    for (i = 0; i < 2 * numDualCurves; i++) {
        dualCurveLen[i] = 0;
        dualCurve[i] = NEW_ARRAY(numCols, int);
        pPathEndPoint[i] = NEW_STRUCT(struct PathEndPoint);
    }

    // Allocate Symplectic Equations Array
    int **symp_eqns = NEW_ARRAY(num_rows, int *);

    for (i = 0; i < num_rows; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    pTriangle = init_cusp_triangulation(manifold);
    numCuspRegions = num_cusp_regions(manifold, pTriangle);
    graph1 = init_graph(numCuspRegions, FALSE);
    struct CuspRegion **pCuspRegion = NEW_ARRAY(numCuspRegions, struct CuspRegion *);

    for (i = 0; i < numCuspRegions; i++) {
        pCuspRegion[i] = NEW_STRUCT(struct CuspRegion);

        pCuspRegion[i]->adjTri[0] = 0;
        pCuspRegion[i]->adjTri[1] = 0;
        pCuspRegion[i]->adjTri[2] = 0;
        pCuspRegion[i]->dist[0] = -1;
        pCuspRegion[i]->dist[1] = -1;
        pCuspRegion[i]->dist[2] = -1;
        pCuspRegion[i]->tetVertex = -1;
        pCuspRegion[i]->tetIndex = -1;
        pCuspRegion[i]->adjNodes[0] = -1;
        pCuspRegion[i]->adjNodes[1] = -1;
        pCuspRegion[i]->adjNodes[2] = -1;
    }

    construct_dual_graph(manifold, graph1, pTriangle, pCuspRegion);
//    print_debug_info(pTriangle, graph1, pCuspRegionIndex, NULL, NULL, 0);
    print_debug_info(pTriangle, graph1, pCuspRegion, NULL, NULL, 2);
//    print_debug_info(pTriangle, graph1, pCuspRegionIndex, NULL, NULL, 3);
    print_debug_info(pTriangle, NULL, NULL, NULL, NULL, 6);


    graph1 = construct_dual_curves(graph1, pTriangle, pCuspRegion, e0, pPathEndPoint, dualCurve, dualCurveLen);
    find_holonomies(graph1, pCuspRegion, pPathEndPoint, symp_eqns, dualCurve, dualCurveLen, e0);

    print_debug_info(pTriangle, graph1, pCuspRegion, dualCurve, dualCurveLen, 5);

    free_graph(graph1);
    free_cusp_triangulation(pTriangle);

    for (i = 0; i < numCuspRegions; i++)
        my_free(pCuspRegion[i]);
    my_free(pCuspRegion);

    for (i = 0; i < 2 * numDualCurves; i++) {
        my_free(dualCurve[i]);
        my_free(pPathEndPoint[i]);
    }

    my_free(pPathEndPoint);
    my_free(dualCurve);
    my_free(dualCurveLen);

    return symp_eqns;
}

/*
 * Number of Cusp Regions
 *
 * Count the number of cusp regions in the manifold to
 * allocate enough space for pCuspRegions
 */

int num_cusp_regions(Triangulation *manifold, struct CuspTriangle **pTriangle) {
    int i, j, numRegions = 0, vertex, intersectIndex;

    intersectIndex = find_cusp_triangle_index(pTriangle, intersectTetIndex, intersectTetVertex);

    // number of regions for triangle without intersection point
    for (i = 0; i < numCuspTriangles; i++) {
        for (j = 0; j < 3; j++) {
            vertex = edgesThreeToFour[pTriangle[i]->tetVertex][j];

            if (i == intersectIndex) {
                /*
                 * If we on the intersection triangle, each closed curve
                 * increases the number of regions by 2. This is because
                 * it splits two regions (one above the intersection point,
                 * and one below) into 4 regions.
                 */
                numRegions = numRegions + 2 * flow(pTriangle[i], vertex) + 1;

                if (flow(pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][0])
                + flow(pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][1])
                + flow(pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][2]) > 2 * genus)
                    // Too many curves on intersection vertex
                    uFatalError("num_cusp_regions", "symplectic_basis.c");
            } else {
                /*
                 * If we are not on the intersection triangle, we have flow
                 * number of curves around the vertex which turn 1 regions in
                 * flow plus one region, each curve adds one region.
                 */
                numRegions = numRegions + flow(pTriangle[i], vertex) + 1;
            }
        }
    }

    return numRegions;
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

void construct_dual_graph(Triangulation *manifold, struct Graph *graph1, struct CuspTriangle **pTriangle, struct CuspRegion **pCuspRegion) {
    int i, index, cuspEdge;
    int indices[3] = { 0, 0, 0 };
    struct CuspTriangle *adjTri[3];
    struct CuspRegion *node;

    struct Node stack;
    stack.item = -1;

    int *visited = NEW_ARRAY(graph1->nvertices, int);

    for (i = 0; i < graph1->nvertices; i++)
        visited[i] = 0;

    // Start at the inside corner of triangle 1.
    push(&stack, 0);
    init_zero_vertex(pCuspRegion[0], pTriangle[0]);
    graph1->pCuspRegionIndex[0] = 0;

    // Walk around the cusp triangulation inserting edges
    while (!is_empty(&stack)) {
        index = pop(&stack);
        node = pCuspRegion[index];

        if (visited[index])
            continue;

        for (i = 0; i < 3; i++) {
            if (pCuspRegion[index]->adjTri[i]) {
                cuspEdge = edgesThreeToFour[node->tetVertex][i];
                adjTri[i] = node->tri->neighbours[cuspEdge];
                indices[i] = insert_triangle_edge(graph1, index, cuspEdge, node->tri, adjTri[i], pCuspRegion);
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
                print_debug_info(pTriangle, graph1, pCuspRegion, NULL, NULL, 2);
            }

            if (!visited[indices[i]] && indices[i] != -2)
                push(&stack, indices[i]);
        }

        visited[index] = 1;
    }
    
    //remove_extra_edges(graph1);
    //add_misc_edges(graph1, edges);

    my_free(visited);
}

void init_zero_vertex(struct CuspRegion *vertex, struct CuspTriangle *pTri) {
    vertex->tri = pTri;
    vertex->tetVertex = 0;          // Tet Index
    vertex->tetIndex = 0;             // Tet Vertex
    vertex->dist[0] = 0;             // dist to v1
    vertex->dist[1] = flow(pTri, edgesThreeToFour[pTri->tetVertex][0])
            + flow(pTri, edgesThreeToFour[pTri->tetVertex][1]);                       // dist to v2
    vertex->dist[2] = flow(pTri, edgesThreeToFour[pTri->tetVertex][0])
            + flow(pTri, edgesThreeToFour[pTri->tetVertex][2]);                       // dist to v3
    vertex->adjTri[0] = (flow(pTri, 1) == 0);
    vertex->adjTri[1] = 1;
    vertex->adjTri[2] = 1;
}

/*
 * Insert Triangle Edge
 */

int insert_triangle_edge(struct Graph *g, int index, int face, struct CuspTriangle *x, struct CuspTriangle *y, struct CuspRegion **pCuspRegion) {
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
        dist = MIN(flow(y, y_vertex1) + pCuspRegion[index]->dist[edgesFourToThree[x->tetVertex][x_vertex1]],
                   flow(y, y_vertex2) + pCuspRegion[index]->dist[edgesFourToThree[x->tetVertex][x_vertex2]])
                           + flow(y, y_face);
    } else {
        // Normal cusp triangle
        dist = flow(y, y_face);

        if (pCuspRegion[index]->dist[edgesFourToThree[x->tetVertex][x_vertex1]] < flow(y, y_vertex1)) {
            // Inside the flows around y_vertex1
            dist += (flow(y, y_vertex1) - pCuspRegion[index]->dist[edgesFourToThree[x->tetVertex][x_vertex1]]);
        } else {
            // Inside the flows around y_vertex2 or center
            dist += (flow(y, y_vertex2) - pCuspRegion[index]->dist[edgesFourToThree[x->tetVertex][x_vertex2]]);
        }
    }

    // Find index of y
    for (i = 0; i < g->nvertices; i++) {
        if (is_equal(pCuspRegion[index], pCuspRegion[i], x, y,
                     x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist)) {
            break;
        }
    }

    // Insert new vertex
    if (i == g->nvertices) {
        for (i = 0; i < g->nvertices && pCuspRegion[i]->tetIndex != -1; i++);

        if (i == g->nvertices)
            return -1;

        init_vertex(pCuspRegion[index], pCuspRegion[i], x, y,
                    x_vertex1, x_vertex2, y_vertex1, y_vertex2, y_face, dist);
    }

    // Update pCuspRegionIndex
    pCuspRegion[index]->adjNodes[edgesFourToThree[x->tetVertex][face]] = i;

    insert_edge(g, index, i, g->directed);

    g->pCuspRegionIndex[i] = i;
    return i;
}

/*
 * Is Equal
 *
 * Check if two holonomies correspond to the same vertex
 */

int is_equal(struct CuspRegion *xNode, struct CuspRegion *yNode, struct CuspTriangle *x, struct CuspTriangle *y,
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

void init_vertex(struct CuspRegion *xNode, struct CuspRegion *yNode, struct CuspTriangle *x, struct CuspTriangle *y,
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

int is_center_vertex(struct CuspRegion *y) {
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

void print_debug_info(struct CuspTriangle **pTriangle, struct Graph *g, struct CuspRegion **pCuspRegion,
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
            if (pCuspRegion[i]->tetIndex == -1 || pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex == -1) {
                tri = pTriangle[0];
            } else {
                tri = find_cusp_triangle(pTriangle, pCuspRegion[i]->tetIndex, pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex);
            }

            printf("Tet Index: %d, Tet Vertex: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d, Cusp Vertex %d dist: %d.\n",
                   pCuspRegion[g->pCuspRegionIndex[i]]->tetIndex,
                   pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex,
                   edgesThreeToFour[tri->tetVertex][0],
                   pCuspRegion[g->pCuspRegionIndex[i]]->dist[0],
                   edgesThreeToFour[tri->tetVertex][1],
                   pCuspRegion[g->pCuspRegionIndex[i]]->dist[1],
                   edgesThreeToFour[tri->tetVertex][2],
                   pCuspRegion[g->pCuspRegionIndex[i]]->dist[2]
            );
        }
    } else if (flag == 2) {
        // Graph Info
        for (i = 0; i < g->nvertices; i++) {
            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d) Region %d (Adj Tri: %d, %d, %d) (Dist: %d, %d, %d): ",
                   i,
                   pCuspRegion[g->pCuspRegionIndex[i]]->tetIndex,
                   pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex,
                   g->pCuspRegionIndex[i],
                   pCuspRegion[g->pCuspRegionIndex[i]]->adjTri[0],
                   pCuspRegion[g->pCuspRegionIndex[i]]->adjTri[1],
                   pCuspRegion[g->pCuspRegionIndex[i]]->adjTri[2],
                   pCuspRegion[g->pCuspRegionIndex[i]]->dist[0],
                   pCuspRegion[g->pCuspRegionIndex[i]]->dist[1],
                   pCuspRegion[g->pCuspRegionIndex[i]]->dist[2]
                   );
            edge = g->edge_list_begin[i];
            while ((edge = edge->next)->next != NULL) {
                printf("%d ", edge->y);
            }
            printf("\n");
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
        for (i = 0; i < 2 * numDualCurves; i++) {
            printf("Curve %d: ", i);

            for (j = 0; j < dualCurveLen[i]; j++) {
                printf("%d ", dualCurve[i][j]);
            }
            printf("\n");
        }
    } else if (flag == 6) {
        // Inside Edge Info
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
                   tri->tet->index,               // Tet Index
                   tri->tetVertex,                // Tet Vertex
                   edge3_between_faces[edgesThreeToFour[tri->tetVertex][1]][edgesThreeToFour[tri->tetVertex][2]],
                   edge3_between_faces[edgesThreeToFour[tri->tetVertex][0]][edgesThreeToFour[tri->tetVertex][2]],
                   edge3_between_faces[edgesThreeToFour[tri->tetVertex][0]][edgesThreeToFour[tri->tetVertex][1]]
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
 * Find Triangle Index
 *
 * Returns a pointer to the triangle with a given tetIndex and tetVertex.
 */

int find_cusp_triangle_index(struct CuspTriangle **pTriangle, int tetIndex, int tetVertex) {
    int i;

    for (i = 0; i < numCuspTriangles; i++) {
        if (tetIndex == pTriangle[i]->tet->index && tetVertex == pTriangle[i]->tetVertex)
            return i;
    }

    return -1;
}

/*
 * Construct Dual Curves
 */

struct Graph *construct_dual_curves(struct Graph *g, struct CuspTriangle **pTriangle, struct CuspRegion **pCuspRegion,
        int e0, struct PathEndPoint **pEndPoint, int **dualCurves, int *dualCurveLen) {
    int i, startIndex, endIndex, startDualIndex, endDualIndex, pathLen;
    bool *processed = NEW_ARRAY(g->nvertices, bool);
    bool *discovered = NEW_ARRAY(g->nvertices, bool);
    int *parent = NEW_ARRAY(g->nvertices, int);
    int *path = NEW_ARRAY(g->nvertices, int);

    find_index(g, pCuspRegion, e0, pEndPoint[2 * e0], pEndPoint[2 * e0 + 1]);

    for (i = 0; i < numDualCurves; i++) {
        if (i == e0)
            continue;

        find_index(g, pCuspRegion, i, pEndPoint[2 * i], pEndPoint[2 * i + 1]);

        // Find path using bfs
        init_search(g, processed, discovered, parent);
        bfs(g, pEndPoint[2 * e0]->graphVertex, processed, discovered, parent);
        find_path(pEndPoint[2 * e0]->graphVertex, pEndPoint[2 * i]->graphVertex, parent, path, 0, &pathLen);
        dualCurves[2 * i] = path;
        dualCurveLen[2 * i] = pathLen;

        // Split graph
        g = split_along_path(g, pCuspRegion, path, pathLen);
        print_debug_info(pTriangle, g, pCuspRegion, dualCurves, dualCurveLen, 2);

        // Reallocate memory
        my_free(processed);
        my_free(discovered);
        my_free(parent);
        processed = NEW_ARRAY(g->nvertices, bool);
        discovered = NEW_ARRAY(g->nvertices, bool);
        parent = NEW_ARRAY(g->nvertices, int);
        path = NEW_ARRAY(g->nvertices, int);

        // Repeat for the other half of the curve
        init_search(g, processed, discovered, parent);
        bfs(g, pEndPoint[2 * e0 + 1]->graphVertex, processed, discovered, parent);
        find_path(pEndPoint[2 * e0 + 1]->graphVertex, pEndPoint[2 * i + 1]->graphVertex, parent, path, 0, &pathLen);
        dualCurves[2 * i + 1] = path;
        dualCurveLen[2 * i + 1] = pathLen;

        // Split graph
        g = split_along_path(g, pCuspRegion, path, pathLen);
        print_debug_info(pTriangle, g, pCuspRegion, dualCurves, dualCurveLen, 2);

        // Re allocate memory
        my_free(processed);
        my_free(discovered);
        my_free(parent);
        processed = NEW_ARRAY(g->nvertices, bool);
        discovered = NEW_ARRAY(g->nvertices, bool);
        parent = NEW_ARRAY(g->nvertices, int);
        path = NEW_ARRAY(g->nvertices, int);
    }

    my_free(processed);
    my_free(discovered);
    my_free(parent);

    return g;
}

/*
 * Find index
 *
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass.
 */

void find_index(struct Graph *g, struct CuspRegion **pCuspRegion, int edgeClass, struct PathEndPoint *endPoint0, struct PathEndPoint *endPoint1) {
    int i, j;
    bool found = FALSE;

    for (i = 0; i < g->nvertices; i++) {
        for (j = 0; j < 3; j++) {
            if (pCuspRegion[g->pCuspRegionIndex[i]]->tri->vertices[j].edgeIndex == edgeClass &&
                pCuspRegion[g->pCuspRegionIndex[i]]->tri->vertices[j].vertexIndex == 0 &&
                pCuspRegion[g->pCuspRegionIndex[i]]->dist[j] == 0) {
                endPoint0->region = pCuspRegion[g->pCuspRegionIndex[i]];
                endPoint0->vertex = edgesThreeToFour[pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex][j];
                endPoint0->face = (int) remaining_face[endPoint0->region->tetVertex][endPoint0->vertex];
                endPoint0->graphVertex = i;
                found = TRUE;
                break;
            }
        }

        if (found)
            break;
    }

    for (i = 0; i < g->nvertices; i++) {
        for (j = 0; j < 3; j++) {
            if (pCuspRegion[g->pCuspRegionIndex[i]]->tri->vertices[j].edgeIndex == edgeClass &&
                pCuspRegion[g->pCuspRegionIndex[i]]->tri->vertices[j].vertexIndex == 1 &&
                pCuspRegion[g->pCuspRegionIndex[i]]->tetIndex == endPoint0->region->tetIndex &&
                pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex == endPoint0->vertex &&
                edgesThreeToFour[pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex][j] == endPoint0->region->tetVertex &&
                pCuspRegion[g->pCuspRegionIndex[i]]->dist[j] == 0) {
                endPoint1->region = pCuspRegion[g->pCuspRegionIndex[i]];
                endPoint1->vertex = edgesThreeToFour[pCuspRegion[g->pCuspRegionIndex[i]]->tetVertex][j];
                endPoint1->face = (int) remaining_face[endPoint0->vertex][endPoint0->region->tetVertex];
                endPoint1->graphVertex = i;
                return;
            }
        }
    }
}

/*
 * Find Holonomies
 *
 * Construct the symplectic equations from the dual curves
 */

void find_holonomies(struct Graph *g, struct CuspRegion **pCuspRegion, struct PathEndPoint **pEndPoint,
        int **symp_eqns, int **dualCurves, int *dualCurveLen, int e0) {
    int i;

    for (i = 0; i < numDualCurves; i++) {
        find_path_holonomy(g, pCuspRegion, pEndPoint[2 * e0], pEndPoint[2 * i], symp_eqns[i],
                           dualCurves[2 * i], dualCurveLen[2 * i]);
        find_path_holonomy(g, pCuspRegion, pEndPoint[2 * e0 + 1], pEndPoint[2 * i + 1],symp_eqns[i],
                           dualCurves[2 * i + 1], dualCurveLen[2 * i + 1]);
        printf("%d: %d, %d: %d\n",
               pEndPoint[2 * e0]->graphVertex,
               pEndPoint[2 * e0]->face,
               pEndPoint[2 * i]->graphVertex,
               pEndPoint[2 * i]->face);
        printf("%d: %d, %d: %d\n",
               pEndPoint[2 * e0 + 1]->graphVertex,
               pEndPoint[2 * e0 + 1]->face,
               pEndPoint[2 * i + 1]->graphVertex,
               pEndPoint[2 * i + 1]->face);

    }
}

/*
 * Find Path Holonomies
 *
 * Calculate holonomies along a path and update the row.
 */

void find_path_holonomy(struct Graph *g, struct CuspRegion **pCuspRegion, struct PathEndPoint *endPoint0,
        struct PathEndPoint *endPoint1, int *row, int *path, int pathLen) {
    int i, j, index, insideVertex, face, dirFace;
    struct CuspRegion *midNode;
    struct PathEndPoint *endPoint;

    for (i = 0; i < pathLen; i++) {
        if (i == 0 || i == pathLen - 1) {
            // Path end points
            if (pathLen == 1) {
                // One vertex
                for (j = 0; j < 4; j++)
                    if (j != endPoint0->region->tetVertex && j != endPoint0->vertex && j != endPoint1->vertex)
                        break;

                index = 3 * endPoint0->region->tetIndex + edge3_between_faces[endPoint0->vertex][endPoint1->vertex];
                row[index] = row[index] + pCuspRegion[g->pCuspRegionIndex[path[i]]]->tri->orientVertices[j][endPoint0->vertex];
                continue;
            }

            if (i == 0) {
                midNode = pCuspRegion[g->pCuspRegionIndex[1]];
                endPoint = endPoint0;
            } else {
                midNode = pCuspRegion[g->pCuspRegionIndex[pathLen - 2]];
                endPoint = endPoint1;
            }

            for (j = 0; j < 3; j++) {
                if (midNode->adjNodes[j] == g->pCuspRegionIndex[path[i]]) {
                    face = EVALUATE(midNode->tri->tet->gluing[midNode->tetVertex], edgesThreeToFour[midNode->tetVertex][j]);
                    break;
                }
            }

            if (endPoint->vertex == face) {
                /*
                 * Curve passes through the face opposite the vertex it
                 * dives through. Picks up holonomy for the vertex between
                 * the face that carries the curve and the face the curve
                 * crosses
                 */
                index = 3 * endPoint->region->tetIndex + edge3_between_faces[endPoint->face][face];
                if (remaining_face[endPoint->vertex][endPoint->region->tetVertex] == face)
                    insideVertex = edgesFourToThree[endPoint->region->tetVertex][remaining_face[endPoint->region->tetVertex][endPoint->vertex]];
                else
                    insideVertex = edgesFourToThree[endPoint->region->tetVertex][face];
                row[index] = row[index] + pCuspRegion[g->pCuspRegionIndex[path[i]]]->tri->orientVertices[insideVertex][face];
            } else if (face == endPoint->face){
                 /*
                  * Curve passes through the face that carries it and
                  * thus picks up no holonomy.
                  */
                 ;
            } else {
                /*
                 * Curve picks up the holonomy for the vertex it dives through
                 */
                index = 3 * endPoint0->region->tetIndex + edge3_between_faces
                        [remaining_face[endPoint->vertex][endPoint->region->tetVertex]]
                        [remaining_face[endPoint->region->tetVertex][endPoint->vertex]];
                row[index] = row[index] + pCuspRegion[g->pCuspRegionIndex[path[i]]]->tri->orientVertices[endPoint->vertex][face];
            }

            continue;
        }

        midNode = pCuspRegion[g->pCuspRegionIndex[path[i]]];
        inside_vertex(midNode, g->pCuspRegionIndex[path[i - 1]], g->pCuspRegionIndex[path[i]], g->pCuspRegionIndex[path[i + 1]],
                      &insideVertex, &face);

        if (insideVertex == -1) {
            printf("didn't find inside vertex\n");
            continue;
        }

        index = 3 * midNode->tetIndex + edge3_between_faces[remaining_face[midNode->tetVertex][insideVertex]][remaining_face[insideVertex][midNode->tetVertex]];
        row[index] = row[index] + pCuspRegion[g->pCuspRegionIndex[path[i]]]->tri->orientVertices[insideVertex][face];
    }
}

/*
 * Inside Vertex
 *
 * Given three cusp nodes, which lie on adjacent cusp triangle,
 * find the inside vertex of the path first -> mid -> last
 */

void inside_vertex(struct CuspRegion *midNode, int first, int mid, int last, int *insideVertex, int *face) {
    int i, vertex, vertex1, vertex2;

    for (i = 0; i < 3; i++) {
        vertex = edgesThreeToFour[midNode->tetVertex][i];

        vertex1 = edgesFourToThree[midNode->tetVertex][(int) remaining_face[midNode->tetVertex][vertex]];
        vertex2 = edgesFourToThree[midNode->tetVertex][(int) remaining_face[vertex][midNode->tetVertex]];

        if (midNode->adjNodes[vertex1] == first && midNode->adjNodes[vertex2] == last) {
            *insideVertex = vertex;
            *face = (int) edgesThreeToFour[midNode->tetVertex][vertex1];
            return;
        } else if (midNode->adjNodes[vertex2] == first && midNode->adjNodes[vertex1] == last) {
            *insideVertex = vertex;
            *face = (int) edgesThreeToFour[midNode->tetVertex][vertex2];
            return;
        }
    }

    *insideVertex = -1;
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

struct Graph *split_along_path(struct Graph *graph1, struct CuspRegion **pCuspRegion, int *path, int pathLen) {
    struct Graph *graph2 = init_graph(graph1->nvertices + pathLen, graph1->directed);

    init_vertices(graph1, graph2, path, pathLen);
    add_non_path_edges(graph1, graph2, path, pathLen);
    add_path_edges(graph1, graph2, pCuspRegion, path, pathLen);

    free_graph(graph1);
    return graph2;
}

/*
 * Initialise Vertices
 *
 * Set the vertex data for the valence vertices.
 */

void init_vertices(struct Graph *graph1, struct Graph *graph2, int *array, int arrayLen) {
    int i, j = 0;

    for (i = 0; i < graph1->nvertices; i++) {
        if (inclusion(array, arrayLen, i)) {
            graph2->pCuspRegionIndex[graph1->nvertices + j] = graph1->pCuspRegionIndex[i];
            j++;
        }

        graph2->pCuspRegionIndex[i] = graph1->pCuspRegionIndex[i];
    }
}

/*
 * Add Non Path Edges
 *
 * For edges not incident to the path, copy the edge to graph2
 */

void add_non_path_edges(struct Graph *graph1, struct Graph *graph2, int *path, int pathLen) {
    int i;
    struct EdgeNode *edge;

    for (i = 0; i < graph1->nvertices; i++) {
        edge = graph1->edge_list_begin[i];

        if (inclusion(path, pathLen, i))
            continue;

        while ((edge = edge->next)->next != NULL) {
            if (!inclusion(path, pathLen, edge->y))
                insert_edge(graph2, i, edge->y, graph2->directed);
        }
    }
}

/*
 * Add Path Edges
 *
 * For vertices on the path, we consider the two cases of the vertex connecting to
 * two vertices on the path (in the middle) or connecting to one vertex (end points).
 * In the first case we split the vertex into two, left and right of the path.
 * In the second case we treat them the same as vertices not on the path.
 */

void add_path_edges(struct Graph *graph1, struct Graph *graph2, struct CuspRegion **pCuspRegion, int *path, int pathLen) {
    int i, j, xLeftVertex, xRightVertex, yLeftVertex, yRightVertex, face;
    struct EdgeNode *edge;
    struct CuspRegion *xRegion, *yRegion;
    bool left = FALSE;

    if (pathLen == 0) {
        return;
    }

    for (i = 0; i < pathLen; i++) {
        edge = graph1->edge_list_begin[path[i]];

        if (i != pathLen - 1) {
            xRegion = pCuspRegion[graph1->pCuspRegionIndex[path[i]]];
            yRegion = pCuspRegion[graph1->pCuspRegionIndex[path[i + 1]]];

            for (j = 0; j < 3; j++) {
                if (xRegion->adjNodes[j] == graph1->pCuspRegionIndex[path[i + 1]]) {
                    face = edgesThreeToFour[xRegion->tetVertex][j];
                    break;
                }
            }

            xRightVertex = (int) remaining_face[xRegion->tetVertex][face];
            xLeftVertex = (int) remaining_face[face][xRegion->tetVertex];

            if (i == 0) {
                yRightVertex = xRightVertex;
                yLeftVertex = xLeftVertex;
            }

            if (xRightVertex == yRightVertex) {
                // The path goes to the right so branches go left
                left = TRUE;
            } else if (xLeftVertex == yLeftVertex) {
                // The path goes to the left so branches go right
                left = FALSE;
            } else {
                // Shouldn't reach this point
                uFatalError("add_path_edges", "symplectic_basis.c");
            }

            yRightVertex = EVALUATE(xRegion->tri->tet->gluing[face], xRightVertex);
            yLeftVertex = EVALUATE(xRegion->tri->tet->gluing[face], xLeftVertex);

        }

        // Case 1: interior points of the path
        if (i != 0 && i != pathLen - 1) {
            // Insert to prev and next path vertex on left and right
            insert_edge(graph2, path[i - 1], path[i], graph2->directed);
            insert_edge(graph2, path[i], path[i + 1], graph2->directed);
            insert_edge(graph2, graph1->nvertices + i - 1, graph1->nvertices + i, graph2->directed);
            insert_edge(graph2, graph1->nvertices + i, graph1->nvertices + i + 1, graph2->directed);

            // Add remaining edges
            while ((edge = edge->next)->next != NULL) {
                if (!inclusion(path, pathLen, edge->y)) {
                    continue;
                }

                if (left) {
                    // Left
                    insert_edge(graph2, path[i], edge->y, graph2->directed);
                } else {
                    // Right
                    insert_edge(graph2, graph1->nvertices + i, edge->y, graph2->directed);
                }
            }
            continue;
        }

        // Case 2: Path End Points
        while ((edge = edge->next)->next != NULL)
            insert_edge(graph2, i, edge->y, graph2->directed);

        if (i == 0) {
            insert_edge(graph2, path[0], graph1->nvertices, graph2->directed);
        } else {
            insert_edge(graph2, path[0], graph1->nvertices + pathLen - 1, graph2->directed);
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


// ---------------------------------------------------------------

// Breadth First Search from Skiena Algorithm Design Manual

/*
 * Initialise Search
 *
 * Initialise default values for bfs arrays
 */

void init_search(struct Graph *g, bool *processed, bool *discovered, int *parent) {
    int i;

    for (i = 0; i < g->nvertices; i ++) {
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
        p = g->edge_list_begin[v];
        while ((p = p->next)->next != NULL) {
            y = p->y;
            if ((!processed[y]) || g->directed) {
                process_edge(v, y);
            }
            if (!discovered[y]) {
                enqueue(&q, y);
                discovered[y] = true;
                parent[y] = v;
            }
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

// -------------------------------------------------

// Data Structures

// Queue Data Structure

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
    my_free(q->array);
}

// Stack Data Structure

void push(struct Node *stack, int item) {
    struct Node *node = NEW_STRUCT(struct Node);
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
    my_free(node);

    return item;
}

int is_empty(struct Node *stack) {
    return stack->item == -1;
}

// Graph Data Structure

/*
 * Initialise Graph
 *
 * Initialise the arrays of the graph 'g' to their default values
 */

struct Graph *init_graph(int maxVertices, bool directed) {
    int i;
    struct Graph *g = NEW_STRUCT(struct Graph);

    g->nvertices = maxVertices;
    g->nedges = 0;
    g->directed = directed;

    g->edge_list_begin = NEW_ARRAY(maxVertices, struct EdgeNode *);
    g->edge_list_end = NEW_ARRAY(maxVertices, struct EdgeNode *);
    g->degree = NEW_ARRAY(maxVertices, int);
    g->pCuspRegionIndex = NEW_ARRAY(maxVertices, int);

    for (i = 0; i < maxVertices; i++) {
        g->degree[i] = 0;
        g->pCuspRegionIndex[i] = 0;

        g->edge_list_begin[i] = NEW_STRUCT(struct EdgeNode);
        g->edge_list_end[i] = NEW_STRUCT(struct EdgeNode);
        g->edge_list_begin[i]->next = g->edge_list_end[i];
        g->edge_list_begin[i]->prev = NULL;
        g->edge_list_end[i]->next = NULL;
        g->edge_list_end[i]->prev = g->edge_list_begin[i];
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
        while (g->edge_list_begin[i]->next != g->edge_list_end[i]) {
            edge = g->edge_list_begin[i]->next;
            REMOVE_NODE(edge);
            my_free(edge);
        }
    }

    my_free(g->edge_list_begin);
    my_free(g->edge_list_end);
    my_free(g->degree);
    my_free(g->pCuspRegionIndex);
    my_free(g);
}

/*
 * Insert Edge
 *
 * Insert an edge into the graph 'g' from vertex x to y.
 */

int insert_edge(struct Graph *g, int x, int y, bool directed) {
    // Ignore edge if it already exists
    if (edge_exists(g, x, y))
        return x;

    struct EdgeNode *p = NEW_STRUCT( struct EdgeNode);
    p->y = y;
    p->next = g->edge_list_begin[x]->next;
    p->prev = g->edge_list_begin[x];
    g->edge_list_begin[x]->next->prev = p;

    g->edge_list_begin[x]->next = p;
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

    node = g->edge_list_begin[vertex_x];

    while (node->next->y != vertex_y  && (node = node->next) != NULL);

    if (node == NULL)
        return;

    deleted_node = node->next;
    node->next = node->next->next;
    node->next->next->prev = node;
    my_free(deleted_node);

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
    struct EdgeNode *node = g->edge_list_begin[v1];

    while ((node = node->next)->next != NULL) {
        if (node->y == v2) {
            return TRUE;
        }
    }

    return FALSE;
}
