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
#include "addl_code.h"

#define FIRST                   0
#define SECOND                  1
#define START                   0
#define FINISH                  1

int edgesThreeToFour[4][3] = {{1, 2, 3},
                              {0, 2, 3},
                              {0, 1, 3},
                              {0, 1, 2}};
int edgesFourToThree[4][4] = {{9, 0, 1, 2},
                              {0, 9, 1, 2},
                              {0, 1, 9, 2},
                              {0, 1, 2, 9}};

/*
 * Get Symplectic Basis
 *
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the get_gluing_equations and get_symplectic_equations functions
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i, j, e0 = 0;

//    peripheral_curves(manifold);

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
        if (i == e0) {
            continue;
        }

        eqns[2 * j] = edge_eqns[i];
        eqns[2 * j + 1] = symp_eqns[i];
        j++;
    }

    my_free(symp_eqns);
    my_free(edge_eqns);
    return eqns;
}

/*
 * Get Symplectic Equations
 *
 * Setup graph and cusp triangulation, and run construct dual graph.
 */
static int numCuspTriangles, numCuspRegions, numDualCurves, numEdgeClasses, genus;
static int intersectTetIndex, intersectTetVertex;

int **get_symplectic_equations(Triangulation *manifold, int num_rows, int numCols, int e0) {
    int i, j, T = manifold -> num_tetrahedra;

    numCuspTriangles    = 4 * manifold->num_tetrahedra;
    numDualCurves       = num_rows;
    genus               = (2 * manifold->num_tetrahedra - num_rows) / 2;

    find_intersection_triangle(manifold);

    struct CuspTriangle **pTriangle     = init_cusp_triangulation(manifold);
    numCuspRegions                      = num_cusp_regions(pTriangle);
    struct CuspRegion **pCuspRegion     = init_cusp_region(pTriangle);
    struct Graph *graph1                = init_graph(numCuspRegions, FALSE);
    struct DualCurves **pDualCurves     = init_oscillating_curves(e0);

    print_debug_info(pTriangle, graph1, pCuspRegion, NULL, 0);
    print_debug_info(pTriangle, graph1, pCuspRegion, NULL, 2);
    print_debug_info(pTriangle, graph1, pCuspRegion, NULL, 3);

    // Allocate Symplectic Equations Array
    int **symp_eqns = NEW_ARRAY(num_rows, int *);

    for (i = 0; i < num_rows; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    construct_dual_graph(graph1, pCuspRegion);
    print_debug_info(pTriangle, NULL, NULL, NULL, 6);

    reduce_graph_size(graph1);
    graph1 = construct_dual_curves(graph1, pTriangle, pCuspRegion, pDualCurves, e0);
    find_holonomies(graph1, pCuspRegion, pDualCurves, e0, symp_eqns);

    print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 5);

    free_graph(graph1);
    free_cusp_triangulation(pTriangle);
    free_cusp_region(pCuspRegion);
    free_oscillating_curves(pDualCurves);
    return symp_eqns;
}

// Free memory used to find the symplectic basis
void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

/*
 * Find Intersection Vertex
 *
 * There exists a cusp triangle containing M and L
 * such that M and L cross the same cusp edge twice
 * exactly once, and we set the intersection of M and
 * L to be the interior of this triangle.
 *
 * To find this cusp triangle, we start where
 * periphal_curves.c places the intersection point, the base
 * tet and vertex. If this triangle is not suitable, we
 * continue forward until we find a suitable triangle.
 */

void find_intersection_triangle(Triangulation *manifold) {
    fundamental_group(manifold, FALSE, FALSE, FALSE, FALSE);

    intersectTetIndex = (int) manifold->cusp_list_begin.next->basepoint_tet->index;
    intersectTetVertex = (int) manifold->cusp_list_begin.next->basepoint_vertex;
}

/*
 * Initialise Cusp Triangulation
 *
 * Construct the pTriangle array which consists of the triangles in the cusp triangulation
 */

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
        pTriangle[i]->tetIndex = pTriangle[i]->tet->index;
        pTriangle[i]->tetVertex = i % 4;

        pTriangle[i]->numCurves = flow(pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][0])
                                  + flow(pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][1])
                                  + flow(pTriangle[i], edgesThreeToFour[pTriangle[i]->tetVertex][2]);

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

        for (j = 0; j < 4; j++) {
            if (j == pTriangle[i]->tetVertex)
                continue;

            pTriangle[i]->faces[j].index            = -1;

            pTriangle[i]->vertices[j].v1            = pTriangle[i]->tetVertex;
            pTriangle[i]->vertices[j].v2            = j;
            pTriangle[i]->vertices[j].edge          = pTriangle[i]->tet->edge_class[
                    edge_between_vertices[pTriangle[i]->vertices[j].v1][pTriangle[i]->vertices[j].v2]];
            pTriangle[i]->vertices[j].edgeIndex     = pTriangle[i]->vertices[j].edge->index;
            pTriangle[i]->vertices[j].vertexIndex   = -1;
        }

        if (i % 4 == 3) {
            tet = tet->next;
        }

        pTriangle[i]->oriented = FALSE;
    }

    cusp_vertex_index(pTriangle);
    vertex_orientation(pTriangle);
    print_debug_info(pTriangle, NULL, NULL, NULL, 4);
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
        for (j = 0; j < 4; j++) {
            if (j == pTriangle[i]->tetVertex || pTriangle[i]->vertices[j].vertexIndex != -1)
                continue;

            walk_around_vertex(pTriangle, pTriangle[i], j,currentIndex[pTriangle[i]->vertices[j].edgeIndex]);

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

    while (tri->vertices[cuspVertex].vertexIndex == -1) {
        tri->vertices[cuspVertex].vertexIndex = index;

        // Move to the next cusp triangle
        old_cusp_vertex         = cuspVertex;
        old_gluing_vertex       = gluing_vertex;
        old_outside_vertex      = outside_vertex;

        cuspVertex          = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
        gluing_vertex       = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_outside_vertex);
        outside_vertex      = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_gluing_vertex);
        tri                 = find_cusp_triangle(pTriangle, tri->neighbours[old_gluing_vertex]->tetIndex,
                                                 EVALUATE(tri->tet->gluing[old_gluing_vertex], tri->tetVertex));
    }
}

/*
 * Vertex Orientation
 *
 * Chooses an orientation for the zero -th
 * pTriangle and then makes every other tri
 * consistent with this one
 */

void vertex_orientation(struct CuspTriangle **pTriangle) {
    int i, j, index, adjIndex, face;
    struct Node indexStack, adjIndexStack, faceStack;
    indexStack.item = -1;
    adjIndexStack.item = -1;
    faceStack.item = -1;

    int zeroTri[4][4] = {
            {0, 0, 0, 0},
            {0, 0, 1, -1},
            {0, -1, 0, 1},
            {0, 1, -1, 0}};

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            pTriangle[0]->orientVertices[i][j] = zeroTri[i][j];

    for (i = 1; i < 4; i++) {
        index = find_cusp_triangle_index(pTriangle, pTriangle[0]->neighbours[i]->tetIndex, pTriangle[0]->neighbours[i]->tetVertex);
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

        for (i = 0; i < 4; i++) {
            if (pTriangle[adjIndex]->tetVertex == i || pTriangle[adjIndex]->neighbours[i]->oriented)
                continue;

            push(&indexStack, adjIndex);
            push(&adjIndexStack, find_cusp_triangle_index(pTriangle,
                                                          pTriangle[adjIndex]->neighbours[i]->tetIndex,
                                                          pTriangle[adjIndex]->neighbours[i]->tetVertex));
            push(&faceStack, i);
        }
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
    x_vertex1   = (int) remaining_face[x->tetVertex][face];
    x_vertex2   = (int) remaining_face[face][x->tetVertex];

    // Vertices on triangle y which correspond to the vertices v1 and v2
    y_vertex1   = EVALUATE(x->tet->gluing[face], x_vertex1);
    y_vertex2   = EVALUATE(x->tet->gluing[face], x_vertex2);
    y_face      = EVALUATE(x->tet->gluing[face], face);

    y->orientVertices[y_vertex1][y_vertex2]     = x->orientVertices[x_vertex1][face];
    y->orientVertices[y_vertex1][y_face]        = x->orientVertices[x_vertex1][x_vertex2];

    y->orientVertices[y_vertex2][y_vertex1]     = x->orientVertices[x_vertex2][face];
    y->orientVertices[y_vertex2][y_face]        = x->orientVertices[x_vertex2][x_vertex1];

    y->orientVertices[y_face][y_vertex1]        = y->orientVertices[y_vertex1][y_vertex2];
    y->orientVertices[y_face][y_vertex2]        = y->orientVertices[y_vertex2][y_vertex1];
}

// Free memory for cusp triangulation data structure
void free_cusp_triangulation(struct CuspTriangle **pTriangle) {
    int i;

    for (i = 0; i < numCuspTriangles; i++)
        my_free(pTriangle[i]);

    my_free(pTriangle);
}

/*
 * Number of Cusp Regions
 *
 * Count the number of cusp regions in the manifold to
 * allocate enough space for pCuspRegions
 */

int num_cusp_regions(struct CuspTriangle **pTriangle) {
    int i, j, numRegions = 0, intersectIndex;

    intersectIndex = find_cusp_triangle_index(pTriangle, intersectTetIndex, intersectTetVertex);

    // number of regions for triangle without intersection point
    for (i = 0; i < numCuspTriangles; i++) {
        if (i == intersectIndex) {
            /*
             * If we on the intersection triangle, each closed curve
             * increases the number of regions by 2. This is because
             * it splits two regions (one above the intersection point,
             * and one below) into 4 regions.
             */
            numRegions = numRegions + 2 * pTriangle[i]->numCurves + 1;

            if (pTriangle[i]->numCurves > 2 * genus)
                // Too many curves on intersection vertex
                uFatalError("num_cusp_regions", "symplectic_basis.c");
        } else {
            /*
             * If we are not on the intersection triangle, we have flow
             * number of curves around the vertex which turn 1 regions in
             * flow plus one region, each curve adds one region.
             */
            numRegions = numRegions + pTriangle[i]->numCurves + 1;
        }
    }

    return numRegions;
}

struct CuspRegion **init_cusp_region(struct CuspTriangle **pTriangle) {
    int i, j, k, l, index, vertex1, vertex2;
    struct CuspTriangle *tri;
    struct CuspRegion *region, **pCuspRegion = NEW_ARRAY(numCuspRegions, struct CuspRegion *);

    for (i = 0; i < numCuspRegions; i++) {
        pCuspRegion[i] = NEW_STRUCT(struct CuspRegion);

        pCuspRegion[i]->tetVertex   = -1;
        pCuspRegion[i]->tetIndex    = -1;

        for (j = 0; j < 4; j++) {
            pCuspRegion[i]->adjTri[j]           = 0;
            pCuspRegion[i]->adjRegions[j]       = -1;
            pCuspRegion[i]->dist[j]             = -1;
        }
    }

    index = 0;
    for (i = 0; i < numCuspTriangles; i++) {
        tri = pTriangle[i];

        // Intersection vertex doesn't have a center
        if (tri->tetIndex == intersectTetIndex && tri->tetVertex == intersectTetVertex) {
            // which vertex are we inside the flow of
            for (j = 0; j < 4; j++) {
                vertex1 = (int) remaining_face[tri->tetVertex][j];
                vertex2 = (int) remaining_face[j][tri->tetVertex];

                if (j == tri->tetVertex) {
                    continue;
                }

                for (k = 1; k < flow(tri, j); k++) {
                    for (l = 0; l < 2; l++) {
                        region = pCuspRegion[index];

                        region->tri                 = tri;
                        region->tetIndex            = region->tri->tetIndex;
                        region->tetVertex           = region->tri->tetVertex;
                        region->dist[vertex1]       = k;
                        region->dist[vertex2]       = MIN(region->dist[vertex1], 2 * flow(tri, j) - region->dist[vertex1])
                                                        + flow(tri, vertex2) + flow(tri, vertex1);
                        region->dist[j]             = flow(tri, j) - region->dist[vertex1] + flow(tri, vertex1);
                        region->adjTri[vertex1]     = 1;
                        index++;

                        // Swap vertices
                        vertex1 = (int) remaining_face[j][tri->tetVertex];
                        vertex2 = (int) remaining_face[tri->tetVertex][j];
                    }
                }

                // Region of dist 0 to j
                region = pCuspRegion[index];

                region->tri                 = tri;
                region->tetIndex            = region->tri->tetIndex;
                region->tetVertex           = region->tri->tetVertex;
                region->dist[j]             = 0;
                region->dist[vertex1]       = flow(tri, j) + flow(tri, vertex1);
                region->dist[vertex2]       = flow(tri, j) + flow(tri, vertex2);
                region->adjTri[vertex1]     = 1;
                region->adjTri[vertex2]     = 1;
                index++;

                // Region in the middle of face j
                if (flow(tri, vertex1) > 0 && flow(tri, vertex2) > 0) {
                    region = pCuspRegion[index];

                    region->tri             = tri;
                    region->tetIndex        = region->tri->tetIndex;
                    region->tetVertex       = region->tri->tetVertex;
                    region->dist[vertex1]   = flow(tri, vertex2);
                    region->dist[vertex2]   = flow(tri, vertex1);
                    region->dist[j]         = MIN(flow(tri, vertex1) + region->dist[vertex1],
                                                  flow(tri, vertex2) + region->dist[vertex2]) + flow(tri, j);
                    region->adjTri[j]       = 1;
                    index++;
                }
            }
            continue;
        }

        // which vertex are we inside the flow of
        for (j = 0; j < 4; j++) {
            vertex1 = (int) remaining_face[tri->tetVertex][j];
            vertex2 = (int) remaining_face[j][tri->tetVertex];

            if (j == tri->tetVertex) {
                continue;
            }

            for (k = 0; k < flow(tri, j); k++) {
                region = pCuspRegion[index];

                region->tri                 = tri;
                region->tetIndex            = region->tri->tetIndex;
                region->tetVertex           = region->tri->tetVertex;

                region->dist[j] = k;
                region->dist[vertex1]       = flow(region->tri, vertex1) + (flow(tri, j) - region->dist[j]);
                region->dist[vertex2]       = flow(region->tri, vertex2) + (flow(tri, j) - region->dist[j]);

                region->adjTri[vertex1]     = 1;
                region->adjTri[vertex2]     = 1;
                index++;
            }

        }

        // Center vertex
        region = pCuspRegion[index];

        region->tri = tri;
        region->tetIndex = region->tri->tetIndex;
        region->tetVertex = region->tri->tetVertex;

        for (j = 0; j < 4; j++) {
            if (j == tri->tetVertex)
                continue;

            region->dist[j] = flow(region->tri, j);
            region->adjTri[j] = 1;
        }
        index++;
    }

    // Add adjacent region info
    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i]->tetIndex == -1)
            break;

        for (j = 0; j < 4; j++) {
            if (!pCuspRegion[i]->adjTri[j]) {
                continue;
            }

            pCuspRegion[i]->adjRegions[j] = find_cusp_region_index(pCuspRegion, pCuspRegion[i], j);
        }
    }

    return pCuspRegion;
}

void free_cusp_region(struct CuspRegion **pCuspRegion) {
    int i;

    for (i = 0; i < numCuspRegions; i++)
        my_free(pCuspRegion[i]);
    my_free(pCuspRegion);
}

int find_cusp_region_index(struct CuspRegion **pCuspRegion, struct CuspRegion *region, int face) {
    int i, vertex1, vertex2, yVertex1, yVertex2, yFace, distV1, distV2, distV3;
    struct CuspTriangle *tri = region->tri;

    vertex1 = (int) remaining_face[tri->tetVertex][face];
    vertex2 = (int) remaining_face[face][tri->tetVertex];

    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i]->tri == tri)
            continue;

        yVertex1 = EVALUATE(tri->tet->gluing[face], vertex1);
        yVertex2 = EVALUATE(tri->tet->gluing[face], vertex2);
        yFace = EVALUATE(tri->tet->gluing[face], face);

        distV1 = (region->dist[vertex1] == pCuspRegion[i]->dist[yVertex1]);
        distV2 = (region->dist[vertex2] == pCuspRegion[i]->dist[yVertex2]);
        distV3 = (region->dist[yFace] >= flow(pCuspRegion[i]->tri, yFace));

        if (distV1 && distV2 && distV3)
            return i;
    }

    // We didn't find a cusp region
    //uFatalError("find_cusp_region", "symplectic_basis");
    return -1;
}

struct DualCurves **init_oscillating_curves(int e0) {
    int i, j;
    struct DualCurves **pDualCurves = NEW_ARRAY(numDualCurves, struct DualCurves *);

    // which curve
    for (i = 0; i < numDualCurves; i++) {
        pDualCurves[i] = NEW_STRUCT(struct DualCurves);

        // which half of the curve
        for (j = 0; j < 2; j++) {
            pDualCurves[i]->curves[j][START]            = NEW_STRUCT(struct EdgeNode);
            pDualCurves[i]->curves[j][START]->next      = pDualCurves[i]->curves[j][FINISH];
            pDualCurves[i]->curves[j][START]->prev      = NULL;
            pDualCurves[i]->endpoints[j][START]         = NEW_STRUCT(struct PathEndPoint);

            pDualCurves[i]->curves[j][FINISH]           = NEW_STRUCT(struct EdgeNode);
            pDualCurves[i]->curves[j][FINISH]->next     = NULL;
            pDualCurves[i]->curves[j][FINISH]->prev     = pDualCurves[i]->curves[j][START];
            pDualCurves[i]->endpoints[j][FINISH]        = pDualCurves[e0]->endpoints[j][START];
        }
    }

    return pDualCurves;
}

void free_oscillating_curves(struct DualCurves **pDualCurves) {
    int i, j;
    struct EdgeNode *edge;

    for (i = 0; i < numDualCurves; i++) {
        for (j = 0; j < 2; j++) {
            while (pDualCurves[i]->curves[j][START]->next != pDualCurves[i]->curves[j][FINISH]) {
                edge = pDualCurves[i]->curves[j][START]->next;
                REMOVE_NODE(edge);
                my_free(edge);
            }
            my_free(pDualCurves[i]->curves[j][START]);
            my_free(pDualCurves[i]->curves[j][FINISH]);
            my_free(pDualCurves[i]->endpoints[j][START]);
        }
        my_free(pDualCurves[i]);
    }

    my_free(pDualCurves);
}

/* Construct Dual Graph
 *
 * Start in the corner of a triangle of the cusp triangulation
 * and walk around the boundary of the manifold, add edges to
 * the dual graph.
 */

void construct_dual_graph(struct Graph *graph1, struct CuspRegion **pCuspRegion) {
    int i, index, face;
    struct CuspRegion *region;

    struct Node stack;
    stack.item = -1;

    int *visited = NEW_ARRAY(numCuspRegions, int);

    for (i = 0; i < graph1->nVertices; i++)
        visited[i] = 0;

    // Start at the inside corner of triangle 1.
    push(&stack, 0);
    graph1->pRegionIndex[0] = 0;

    // Walk around the cusp triangulation inserting edges
    while (!is_empty(&stack)) {
        index = pop(&stack);
        region = pCuspRegion[index];

        if (visited[index])
            continue;

        for (face = 0; face < 4; face++) {
            if (!region->adjTri[face])
                continue;

            insert_edge(graph1, index, region->adjTri[face], graph1->directed);
            graph1->pRegion[index] = region;
            graph1->pRegionIndex[index] = index;

            if (!visited[index])
                push(&stack, region->adjTri[face]);
        }

        visited[index] = 1;
    }

    my_free(visited);
}

/*
 * Print Triangle Information
 *
 * Debug function for printing the gluing information
 */

void print_debug_info(struct CuspTriangle **pTriangle, struct Graph *g, struct CuspRegion **pCuspRegion,
                      struct DualCurves **pDualCurves, int flag) {
    int i, j, k, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

    struct CuspTriangle *tri;
    struct EdgeNode *edge ;
    struct CuspRegion *region;

    if (!flag) {
        // Gluing Info
        printf("Triangle gluing info\n");
        for (i = 0; i < numCuspTriangles; i++) {
            for (j = 0; j < 4; j++) {
                if (j == pTriangle[i]->tetVertex)
                    continue;

                tri = pTriangle[i];
                x_vertex1 = (int) remaining_face[tri->tetVertex][j];
                x_vertex2 = (int) remaining_face[j][tri->tetVertex];
                y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                printf("(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                       "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                       tri->tetIndex,               // Tet Index
                       tri->tetVertex,                // Tet Vertex
                       j,      // Cusp Edge
                       tri->tet->neighbor[j]->index,                              // Tet Index
                       EVALUATE(tri->tet->gluing[j], tri->tetVertex),             // Tet Vertex
                       EVALUATE(tri->tet->gluing[j], j),   // Cusp Edge
                       x_vertex1, y_vertex1,
                       x_vertex2, y_vertex2
                       );
            }
        }
    } else if (flag == 1) {
        // Vertex Distance Info
        printf("Empty\n");
    } else if (flag == 2) {
        // Region Info
        printf("Cusp Region info\n");
        for (i = 0; i < numCuspRegions; i++) {
            region = pCuspRegion[i];
            printf("Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Dist: %d, %d, %d)\n",
                   i,
                   region->tetIndex,
                   region->tetVertex,
                   region->adjTri[edgesThreeToFour[region->tetVertex][0]],
                   region->adjTri[edgesThreeToFour[region->tetVertex][1]],
                   region->adjTri[edgesThreeToFour[region->tetVertex][2]],
                   region->dist[edgesThreeToFour[region->tetVertex][0]],
                   region->dist[edgesThreeToFour[region->tetVertex][1]],
                   region->dist[edgesThreeToFour[region->tetVertex][2]]
                   );
        }
    } else if (flag == 3) {
        // Homology Info
        printf("Homology info\n");
        printf("Meridian\n");

        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];

            printf("(Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tetIndex,
                   tri->tetVertex,
                   tri->tet->curve[M][right_handed][tri->tetVertex][0],
                   tri->tet->curve[M][right_handed][tri->tetVertex][1],
                   tri->tet->curve[M][right_handed][tri->tetVertex][2],
                   tri->tet->curve[M][right_handed][tri->tetVertex][3]
                );
        }
        printf("Longitude\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tetIndex,
                   tri->tetVertex,
                   tri->tet->curve[L][right_handed][tri->tetVertex][0],
                   tri->tet->curve[L][right_handed][tri->tetVertex][1],
                   tri->tet->curve[L][right_handed][tri->tetVertex][2],
                   tri->tet->curve[L][right_handed][tri->tetVertex][3]
                );
        }
        printf("Left handed\n");
        printf("Meridian\n");

        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];

            printf("(Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tetIndex,
                   tri->tetVertex,
                   tri->tet->curve[M][left_handed][tri->tetVertex][0],
                   tri->tet->curve[M][left_handed][tri->tetVertex][1],
                   tri->tet->curve[M][left_handed][tri->tetVertex][2],
                   tri->tet->curve[M][left_handed][tri->tetVertex][3]
            );
        }
        printf("Longitude\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tetIndex,
                   tri->tetVertex,
                   tri->tet->curve[L][left_handed][tri->tetVertex][0],
                   tri->tet->curve[L][left_handed][tri->tetVertex][1],
                   tri->tet->curve[L][left_handed][tri->tetVertex][2],
                   tri->tet->curve[L][left_handed][tri->tetVertex][3]
            );
        }
    } else if (flag == 4) {
        // Edge indices
        printf("Edge classes\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            v1 = edgesThreeToFour[tri->tetVertex][0];
            v2 = edgesThreeToFour[tri->tetVertex][1];
            v3 = edgesThreeToFour[tri->tetVertex][2];

            printf("(Tet Index: %d, Tet Vertex: %d) Vertex %d: (Edge Class: %d, Edge Index: %d, Face: %d), "
                   "Vertex %d: (Edge Class: %d, Edge Index: %d, Face: %d), Vertex %d: (Edge Class: %d, Edge Index: %d, Face: %d)\n",
                   tri->tetIndex, tri->tetVertex,
                   v1, tri->vertices[v1].edgeIndex, tri->vertices[v1].vertexIndex, tri->faces[v1].index,
                   v2, tri->vertices[v2].edgeIndex, tri->vertices[v2].vertexIndex, tri->faces[v2].index,
                   v3, tri->vertices[v3].edgeIndex, tri->vertices[v3].vertexIndex, tri->faces[v3].index
            );
        }
    } else if (flag == 5) {
        // Dual Curve Paths
        printf("Oscillating curve paths\n");
        for (i = 0; i < numDualCurves; i++) {
            for (j = 0; j < 2; j++) {
                edge = pDualCurves[i]->curves[j][START];

                while ((edge = edge->next)->next != NULL) {
                    printf("%d ", edge->y);
                }

                printf("\n");
            }
        }
    } else if (flag == 6) {
        // Inside Edge Info
        printf("Inside edge info\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            printf("(Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
                   tri->tetIndex,               // Tet Index
                   tri->tetVertex,                // Tet Vertex
                   edge3_between_faces[edgesThreeToFour[tri->tetVertex][1]][edgesThreeToFour[tri->tetVertex][2]],
                   edge3_between_faces[edgesThreeToFour[tri->tetVertex][0]][edgesThreeToFour[tri->tetVertex][2]],
                   edge3_between_faces[edgesThreeToFour[tri->tetVertex][0]][edgesThreeToFour[tri->tetVertex][1]]
            );
        }
    } else if (flag == 7) {
        printf("Graph info\n");
        for (i = 0; i < g->nVertices; i++) {\
            region = pCuspRegion[g->pRegionIndex[i]];
            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d) (Dist: %d, %d, %d): ",
                   i,
                   region->tetIndex,
                   region->tetVertex,
                   region->dist[edgesThreeToFour[tri->tetVertex][0]],
                   region->dist[edgesThreeToFour[tri->tetVertex][1]],
                   region->dist[edgesThreeToFour[tri->tetVertex][2]]
            );
            edge = g->edge_list_begin[i];
            while ((edge = edge->next)->next != NULL) {
                printf("%d ", edge->y);
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
    mflow = FLOW(tri->tet->curve[M][right_handed][tri->tetVertex][remaining_face[tri->tetVertex][vertex]],
                 tri->tet->curve[M][right_handed][tri->tetVertex][remaining_face[vertex][tri->tetVertex]]);

    // Contribution from longitudinal curves
    lflow = FLOW(tri->tet->curve[L][right_handed][tri->tetVertex][remaining_face[tri->tetVertex][vertex]],
                 tri->tet->curve[L][right_handed][tri->tetVertex][remaining_face[vertex][tri->tetVertex]]);

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
        if (tetIndex == pTriangle[i]->tetIndex && tetVertex == pTriangle[i]->tetVertex)
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
        if (tetIndex == pTriangle[i]->tetIndex && tetVertex == pTriangle[i]->tetVertex)
            return i;
    }

    return -1;
}

/*
 * Construct Dual Curves
 */

struct Graph *construct_dual_curves(struct Graph *g, struct CuspTriangle **pTriangle, struct CuspRegion **pCuspRegion,
        struct DualCurves **pDualCurves, int e0) {
    int i, j, *parent;
    bool *processed, *discovered;

    find_index(g, pCuspRegion, e0, pDualCurves[e0]->endpoints[FIRST][START],pDualCurves[e0]->endpoints[SECOND][START]);

    for (i = 0; i < numDualCurves; i++) {
        if (i == 0)
            continue;

        find_index(g, pCuspRegion, i, pDualCurves[i]->endpoints[FIRST][START],pDualCurves[i]->endpoints[SECOND][START]);

        // which half of the curve
        for (j = 0; j < 2; j++) {
            processed = NEW_ARRAY(g->nVertices, bool);
            discovered = NEW_ARRAY(g->nVertices, bool);
            parent = NEW_ARRAY(g->nVertices, int);

            // Find path using bfs
            init_search(g, processed, discovered, parent);
            bfs(g, pDualCurves[i]->endpoints[j][START]->graphVertex, processed, discovered, parent);
            find_path(
                    pDualCurves[i]->endpoints[j][START]->graphVertex,
                    pDualCurves[i]->endpoints[j][FINISH]->graphVertex, parent, pDualCurves[i]->curves[j][START]);

            // Split graph
            g = split_along_path(g, pCuspRegion, pDualCurves[i], j);
            print_debug_info(pTriangle, g, pCuspRegion, pDualCurves, 7);

            // Reallocate memory
            my_free(processed);
            my_free(discovered);
            my_free(parent);
        }
    }

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
    struct CuspRegion *pRegion;
    bool found = FALSE;

    for (i = 0; i < g->nVertices; i++) {
        pRegion = g->pRegion[i];
        for (j = 0; j < 4; j++) {
            if (j == pRegion->tetVertex)
                continue;

            if (pRegion->tri->vertices[j].edgeIndex == edgeClass &&
                pRegion->tri->vertices[j].vertexIndex == 0 &&
                pRegion->dist[j] == 0) {
                endPoint0->region = g->pRegion[i];
                endPoint0->regionIndex = g->pRegionIndex[i];
                endPoint0->vertex = j;
                endPoint0->face = (int) remaining_face[endPoint0->region->tetVertex][endPoint0->vertex];
                endPoint0->graphVertex = i;
                found = TRUE;
                break;
            }
        }

        if (found)
            break;
    }

    for (i = 0; i < g->nVertices; i++) {
        pRegion = g->pRegion[i];

        for (j = 0; j < 4; j++) {
            if (j == pRegion->tetVertex)
                continue;

            if (pRegion->tri->vertices[j].edgeIndex == edgeClass &&
                pRegion->tri->vertices[j].vertexIndex == 1 &&
                pRegion->tetIndex == endPoint0->region->tetIndex &&
                pRegion->tetVertex == endPoint0->vertex &&
                j == endPoint0->region->tetVertex &&
                pRegion->dist[j] == 0) {
                endPoint1->region = pRegion;
                endPoint1->regionIndex = g->pRegionIndex[i];
                endPoint1->vertex = j;
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

void find_holonomies(struct Graph *g, struct CuspRegion **pCuspRegion, struct DualCurves **pDualCurves, int e0, int **symp_eqns) {
    int i;

    for (i = 0; i < numDualCurves; i++) {
        if (i == e0)
            continue;

        find_path_holonomy(g, pCuspRegion, pDualCurves[i], FIRST, symp_eqns[i]);
        find_path_holonomy(g, pCuspRegion, pDualCurves[i], SECOND, symp_eqns[i]);
    }
}

/*
 * Find Path Holonomies
 *
 * Calculate holonomies along a path and update the row.
 */

void find_path_holonomy(struct Graph *g, struct CuspRegion **pCuspRegion, struct DualCurves *curve, int curveNum, int *row) {
    int i, j, index, insideVertex, face, dirFace, midNodeIndex;
    struct CuspRegion *midNode;
    struct EdgeNode *insideNode, *node = curve->curves[curveNum][START];
    struct PathEndPoint *endPoint, *pathStartPoint = curve->endpoints[curveNum][START], *pathEndPoint = curve->endpoints[curveNum][FINISH];

    // Edge cases
    if (node->next->next == NULL)
        return;
    else if (node->next->next->next == NULL) {
        // One vertex
        for (j = 0; j < 4; j++)
            if (j != pathStartPoint->region->tetVertex && j != pathStartPoint->vertex && j != pathEndPoint->vertex)
                break;

        index = 3 * pathStartPoint->region->tetIndex + edge3_between_faces[pathStartPoint->vertex][pathEndPoint->vertex];
        row[index] = row[index] + g->pRegion[node->next->y]->tri->orientVertices[j][pathStartPoint->vertex];
        return;
    }

    while((node = node->next)->next != NULL) {
        // Path end points
        if (node->prev->prev == NULL || node->next->next == NULL) {
            // End point index
            if (node->prev->prev == NULL) {
                insideNode = node->next;
                dirFace = 1;
            } else {
                insideNode = node->prev;
                dirFace = -1;
            }

            if (g->pRegionIndex[node->y] == pathStartPoint->region)
                endPoint = pathStartPoint;
            else
                endPoint = pathEndPoint;

            // Find the face the next vertex lies across
            for (j = 0; j < 4; j++) {
                if (j == endPoint->region->tetVertex)
                    continue;

                if (pCuspRegion[endPoint->region->adjRegions[j]] == g->pRegionIndex[insideNode->y]) {
                    face = j;
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

                if (remaining_face[endPoint->vertex][endPoint->region->tetVertex] == endPoint->vertex)
                    insideVertex = (int) remaining_face[endPoint->region->tetVertex][endPoint->vertex];
                else
                    insideVertex = (int) remaining_face[endPoint->vertex][endPoint->region->tetVertex];
            } else if (face == endPoint->face){
                 /*
                  * Curve passes through the face that carries it and
                  * thus picks up no holonomy.
                  */
                 continue;
            } else {
                /*
                 * Curve picks up the holonomy for the vertex it dives through
                 */
                index = 3 * pathStartPoint->region->tetIndex + edge3_between_faces
                        [remaining_face[endPoint->vertex][endPoint->region->tetVertex]]
                        [remaining_face[endPoint->region->tetVertex][endPoint->vertex]];

                insideVertex = endPoint->vertex;
            }

            row[index] = row[index] + dirFace * g->pRegion[node->y]->tri->orientVertices[insideVertex][face];

            continue;
        }

        midNode = g->pRegion[node->y];
        inside_vertex(midNode, g->pRegionIndex[node->prev->y], g->pRegionIndex[node->y],
                      g->pRegionIndex[node->next->y], &insideVertex, &face);

        if (insideVertex == -1) {
            printf("didn't find inside vertex\n");
            continue;
        }

        index = 3 * midNode->tetIndex + edge3_between_faces[remaining_face[midNode->tetVertex][insideVertex]][remaining_face[insideVertex][midNode->tetVertex]];
        row[index] = row[index] + g->pRegion[node->y]->tri->orientVertices[insideVertex][face];
    }
}

/*
 * Inside Vertex
 *
 * Given three cusp nodes, which lie on adjacent cusp triangle,
 * find the inside vertex of the path first -> mid -> last
 */

void inside_vertex(struct CuspRegion *midNode, int last, int mid, int first, int *insideVertex, int *face) {
    int i, vertex1, vertex2;

    for (i = 0; i < 4; i++) {
        if (i == midNode->tetVertex)
            continue;

        vertex1 = (int) remaining_face[midNode->tetVertex][i];
        vertex2 = (int) remaining_face[i][midNode->tetVertex];

        if (midNode->adjRegions[vertex1] == first && midNode->adjRegions[vertex2] == last) {
            *insideVertex = i;
            *face = vertex1;
            return;
        } else if (midNode->adjRegions[vertex2] == first && midNode->adjRegions[vertex1] == last) {
            *insideVertex = i;
            *face = vertex2;
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

struct Graph *split_along_path(struct Graph *graph1, struct CuspRegion **pCuspRegion, struct DualCurves *path, int curveNum) {
    int i = 1;
    struct EdgeNode *node = path->curves[curveNum][START];
    while ((node = node->next)->next != NULL)
        i++;

    struct Graph *graph2 = init_graph(graph1->nVertices + i, graph1->directed);

    init_vertices(graph1, graph2, path->curves[curveNum][START]);
    add_non_path_edges(graph1, graph2, path->curves[curveNum][START]);
    add_path_edges(graph1, graph2, pCuspRegion, path->curves[curveNum][START]);

    free_graph(graph1);
    return graph2;
}

/*
 * Initialise Vertices
 *
 * Set the vertex data for the valence vertices.
 */

void init_vertices(struct Graph *graph1, struct Graph *graph2, struct EdgeNode *node) {
    int i, j = 0;

    for (i = 0; i < graph1->nVertices; i++) {
        if (inclusion(node->next, i)) {
            graph2->pRegionIndex[graph1->nVertices + j] = graph1->pRegionIndex[i];
            j++;
        }

        graph2->pRegionIndex[i] = graph1->pRegionIndex[i];
    }
}

/*
 * Add Non Path Edges
 *
 * For edges not incident to the path, copy the edge to graph2
 */

void add_non_path_edges(struct Graph *graph1, struct Graph *graph2, struct EdgeNode *node) {
    int i;
    struct EdgeNode *edge;

    for (i = 0; i < graph1->nVertices; i++) {
        edge = graph1->edge_list_begin[i];

        if (inclusion(node->next, i))
            continue;

        while ((edge = edge->next)->next != NULL) {
            if (!inclusion(node->next, edge->y))
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

void add_path_edges(struct Graph *graph1, struct Graph *graph2, struct CuspRegion **pCuspRegion, struct EdgeNode *node) {
    int i = 0, j, xLeftVertex, xRightVertex, yLeftVertex, yRightVertex, face;
    struct EdgeNode *edge;
    struct CuspRegion *xRegion;
    bool left = FALSE;

    while ((node = node->next)->next != NULL) {
        edge = graph1->edge_list_begin[node->y];

        if (node->next->next != NULL) {
            if (graph1->pRegionIndex[node->y] == graph1->pRegionIndex[node->next->y]) {
                i++;
                continue;
            }

            xRegion = graph1->pRegion[node->y];

            for (j = 0; j < 4; j++) {
                if (j == xRegion->tetVertex)
                    continue;

                if (xRegion->adjRegions[j] == graph1->pRegionIndex[node->next->y]) {
                    face = j;
                    break;
                }
            }

            xRightVertex = (int) remaining_face[xRegion->tetVertex][face];
            xLeftVertex = (int) remaining_face[face][xRegion->tetVertex];

            if (node->prev->prev == NULL) {
                yRightVertex = xRightVertex;
                yLeftVertex = xLeftVertex;
            }

            if (xRightVertex == yRightVertex || xRightVertex == yLeftVertex) {
                // The path goes to the right so branches go left
                left = TRUE;
            } else if (xLeftVertex == yLeftVertex || xLeftVertex == xRightVertex) {
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
        if (node->prev->prev != NULL && node->next->next != NULL) {
            // Insert to prev and next path vertex on left and right
            insert_edge(graph2, node->prev->y, node->y, graph2->directed);
            insert_edge(graph2, node->y, node->next->y, graph2->directed);
            insert_edge(graph2, graph1->nVertices + i - 1, graph1->nVertices + i, graph2->directed);
            insert_edge(graph2, graph1->nVertices + i, graph1->nVertices + i + 1, graph2->directed);

            // Add remaining edges
            while ((edge = edge->next)->next != NULL) {
                if (!inclusion(node, edge->y)) {
                    continue;
                }

                if (left) {
                    // Left
                    insert_edge(graph2, node->y, edge->y, graph2->directed);
                } else {
                    // Right
                    insert_edge(graph2, graph1->nVertices + i, edge->y, graph2->directed);
                }
            }
            i++;
            continue;
        }

        // Case 2: Path End Points (should be covered by other cases)
        while ((edge = edge->next)->next != NULL)
            insert_edge(graph2, i, edge->y, graph2->directed);
    }
}

/*
 * Inclusion
 *
 * Check if the target is in the array
 */

bool inclusion(struct EdgeNode *node, int target) {
    if (node->next == NULL)
        return FALSE;
    else if (node->y == target)
        return TRUE;
    else
        return inclusion(node->next, target);
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

    for (i = 0; i < g->nVertices; i ++) {
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
void find_path(int start, int end, int *parents, struct EdgeNode *node) {
    struct EdgeNode *new_node = NEW_STRUCT(struct EdgeNode);

    INSERT_AFTER(new_node, node);

    if ((start == end) || (end == -1)) {
        new_node->y = start;
    } else {
        find_path(start, parents[end], parents, node);
        new_node->y = end;
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

    g->nVertices = maxVertices;
    g->directed = directed;

    g->edge_list_begin = NEW_ARRAY(maxVertices, struct EdgeNode *);
    g->edge_list_end = NEW_ARRAY(maxVertices, struct EdgeNode *);
    g->degree = NEW_ARRAY(maxVertices, int);
    g->pRegion = NEW_ARRAY(maxVertices, struct CuspRegion *);
    g->pRegionIndex = NEW_ARRAY(maxVertices, int);

    for (i = 0; i < maxVertices; i++) {
        g->degree[i] = 0;
        g->pRegionIndex[i] = 0;

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

    for (i = 0; i < g->nVertices; i++) {
        while (g->edge_list_begin[i]->next != g->edge_list_end[i]) {
            edge = g->edge_list_begin[i]->next;
            REMOVE_NODE(edge);
            my_free(edge);
        }
    }

    my_free(g->edge_list_begin);
    my_free(g->edge_list_end);
    my_free(g->degree);
    my_free(g->pRegion);
    my_free(g->pRegionIndex);
    my_free(g);
}

void reduce_graph_size(struct Graph *g) {
    return;
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
    INSERT_AFTER(p, g->edge_list_begin[x]);
    p->y = y;
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
