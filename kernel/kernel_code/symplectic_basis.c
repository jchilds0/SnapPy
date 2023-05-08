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

#define atleast_two(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))
#define orientate_vertex(tri, vertex, face)       ((face) == remaining_face[vertex][(tri)->tetVertex] ? 1 : -1)

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
static int debug = FALSE;

int **get_symplectic_equations(Triangulation *manifold, int num_rows, int numCols, int e0) {
    int i, j, T = manifold -> num_tetrahedra;

    numCuspTriangles    = 4 * manifold->num_tetrahedra;
    numDualCurves       = num_rows;
    genus               = (2 * manifold->num_tetrahedra - num_rows) / 2;

    find_intersection_triangle(manifold);

    struct CuspTriangle **pTriangle     = init_cusp_triangulation(manifold);
    numCuspRegions                      = num_cusp_regions(pTriangle);
    struct CuspRegion **pCuspRegion     = init_cusp_region(pTriangle);
    struct DualCurves **pDualCurves     = init_oscillating_curves(e0);

    print_debug_info(pTriangle, NULL, pCuspRegion, pDualCurves, 0);
    print_debug_info(pTriangle, NULL, pCuspRegion, pDualCurves, 3);
    print_debug_info(pTriangle, NULL, pCuspRegion, pDualCurves, 4);
    print_debug_info(pTriangle, NULL, pCuspRegion, pDualCurves, 6);
    print_debug_info(pTriangle, NULL, pCuspRegion, pDualCurves, 2);

    // Allocate Symplectic Equations Array
    int **symp_eqns = NEW_ARRAY(num_rows, int *);

    for (i = 0; i < num_rows; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    pCuspRegion = construct_dual_curves(pTriangle, pCuspRegion, pDualCurves, e0);
    find_holonomies(pCuspRegion, pDualCurves, e0, symp_eqns);

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
    int i, j;
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
            if (j == pTriangle[i]->tetVertex)
                continue;

            pTriangle[i]->vertices[j].v1            = pTriangle[i]->tetVertex;
            pTriangle[i]->vertices[j].v2            = j;
            pTriangle[i]->vertices[j].edge          = pTriangle[i]->tet->edge_class[
                    edge_between_vertices[pTriangle[i]->vertices[j].v1][pTriangle[i]->vertices[j].v2]];
            pTriangle[i]->vertices[j].edgeClass     = pTriangle[i]->vertices[j].edge->index;
            pTriangle[i]->vertices[j].edgeIndex   = -1;
        }

        if (i % 4 == 3) {
            tet = tet->next;
        }
    }

    cusp_vertex_index(pTriangle);
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
            if (j == pTriangle[i]->tetVertex || pTriangle[i]->vertices[j].edgeIndex != -1)
                continue;

            walk_around_vertex(pTriangle, pTriangle[i], j,currentIndex[pTriangle[i]->vertices[j].edgeClass]);

            currentIndex[pTriangle[i]->vertices[j].edgeClass]++;
        }
    }

    my_free(currentIndex);
}

/*
 * Walk around cusp vertex
 *
 * Walk around vertex cuspVertex of triangle pTriangle[cuspTriangle]
 * and set edgeIndex to index.
 */

void walk_around_vertex(struct CuspTriangle **pTriangle, struct CuspTriangle *tri, int cuspVertex, int index) {
    int gluing_vertex, outside_vertex, old_gluing_vertex, old_cusp_vertex, old_outside_vertex;
    gluing_vertex = (int) remaining_face[cuspVertex][tri->tetVertex];
    outside_vertex = (int) remaining_face[tri->tetVertex][cuspVertex];

    while (tri->vertices[cuspVertex].edgeIndex == -1) {
        tri->vertices[cuspVertex].edgeIndex = index;

        // Move to the next cusp triangle
        old_cusp_vertex         = cuspVertex;
        old_gluing_vertex       = gluing_vertex;
        old_outside_vertex      = outside_vertex;

        cuspVertex          = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
        gluing_vertex       = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_outside_vertex);
        outside_vertex      = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_gluing_vertex);
        tri                 = pTriangle[find_cusp_triangle_index(pTriangle, tri->neighbours[old_gluing_vertex]->tetIndex,
                                                 EVALUATE(tri->tet->gluing[old_gluing_vertex], tri->tetVertex))];
    }
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
    int i, numRegions = 0, intersectIndex;

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
    int i, j, k, l, index, vertex1, vertex2, vertex3, dist[4], adjTri[4];
    struct CuspTriangle *tri;
    struct CuspRegion **pCuspRegion = NEW_ARRAY(numCuspRegions, struct CuspRegion *);

    for (i = 0; i < numCuspRegions; i++) {
        pCuspRegion[i] = NEW_STRUCT(struct CuspRegion);

        pCuspRegion[i]->tetVertex   = -1;
        pCuspRegion[i]->tetIndex    = -1;

        for (j = 0; j < 4; j++) {
            pCuspRegion[i]->adjTri[j]           = 0;
            pCuspRegion[i]->adjRegions[j]       = -1;

            for (k = 0; k < 4; k ++)
                pCuspRegion[i]->curve[j][k] = -1;
        }
    }

    index = 0;
    for (i = 0; i < numCuspTriangles; i++) {
        tri = pTriangle[i];

        // Intersection vertex doesn't have a center
        if (tri->tetIndex == intersectTetIndex && tri->tetVertex == intersectTetVertex) {
            // which vertex are we inside the flow of
            for (j = 0; j < 4; j++) {
                if (j == tri->tetVertex) {
                    continue;
                }

                vertex1 = (int) remaining_face[tri->tetVertex][j];
                vertex2 = (int) remaining_face[j][tri->tetVertex];

                for (k = 1; k < flow(tri, j); k++) {
                    for (l = 0; l < 2; l++) {
                        dist[vertex1] = k;
                        dist[vertex2] = MIN(dist[vertex1], 2 * flow(tri, j) - dist[vertex1])
                                + flow(tri, vertex2) + flow(tri, vertex1);
                        dist[j] = flow(tri, j) - dist[vertex1] + flow(tri, vertex1);
                        dist[tri->tetVertex] = -1;

                        adjTri[vertex1] = 1;
                        adjTri[vertex2] = 0;
                        adjTri[j] = 0;
                        adjTri[tri->tetVertex] = -1;

                        init_region(pCuspRegion[index], tri, dist, adjTri);
                        index++;

                        // Swap vertices
                        vertex1 = (int) remaining_face[j][tri->tetVertex];
                        vertex2 = (int) remaining_face[tri->tetVertex][j];
                    }
                }

                // Region in the middle of face j
                if (flow(tri, vertex1) && flow(tri, vertex2)) {
                    dist[vertex1] = flow(tri, vertex2);
                    dist[vertex2] = flow(tri, vertex1);
                    dist[j] = MIN(flow(tri, vertex1) + dist[vertex1],
                                                  flow(tri, vertex2) + dist[vertex2]) + flow(tri, j);
                    dist[tri->tetVertex] = -1;

                    adjTri[vertex1] = 0;
                    adjTri[vertex2] = 0;
                    adjTri[j] = 1;
                    adjTri[tri->tetVertex] = -1;

                    init_region(pCuspRegion[index], tri, dist, adjTri);
                    index++;
                }
            }

            // Region of dist 0 to j
            vertex1 = edgesThreeToFour[tri->tetVertex][0];
            vertex2 = edgesThreeToFour[tri->tetVertex][1];
            vertex3 = edgesThreeToFour[tri->tetVertex][2];

            // Case 1; two regions
            // !!Double check!!
            if (atleast_two(!flow(tri, vertex1), !flow(tri, vertex2), !flow(tri, vertex3))) {
                dist[vertex1] = flow(tri, vertex1);
                dist[vertex2] = flow(tri, vertex2);
                dist[vertex3] = flow(tri, vertex3);
                dist[tri->tetVertex] = -1;

                adjTri[vertex1]     = 1;
                adjTri[vertex2]     = 1;
                adjTri[vertex3]     = 1;
                adjTri[tri->tetVertex] = -1;

                init_region(pCuspRegion[index], tri, dist, adjTri);
                index++;

                // Find vertex with non zero flow
                for (j = 0; j < 4; j++) {
                    if (j == tri->tetVertex)
                        continue;

                    if (flow(tri, j)) {
                        vertex1 = j;
                        vertex2 = (int) remaining_face[tri->tetVertex][vertex1];
                        vertex3 = (int) remaining_face[vertex1][tri->tetVertex];
                        break;
                    }
                }
                dist[vertex1] = 0;
                dist[vertex2] = flow(tri, vertex1);
                dist[vertex3] = flow(tri, vertex1);
                dist[tri->tetVertex] = -1;

                adjTri[vertex1]         = 0;
                adjTri[vertex2]         = 1;
                adjTri[vertex3]         = 1;
                adjTri[tri->tetVertex]  = 0;

                init_region(pCuspRegion[index], tri, dist, adjTri);
                index++;

            } else {
                // Case 2: three regions
                for (j = 0; j < 4; j++) {
                    if (j == tri->tetVertex)
                        continue;

                    vertex1 = (int) remaining_face[tri->tetVertex][j];
                    vertex2 = (int) remaining_face[j][tri->tetVertex];

                    dist[j] = 0;
                    dist[vertex1] = flow(tri, j) + flow(tri, vertex1);
                    dist[vertex2] = flow(tri, j) + flow(tri, vertex2);
                    dist[tri->tetVertex] = -1;

                    adjTri[j] = 0;
                    adjTri[vertex1]     = 1;
                    adjTri[vertex2]     = 1;
                    adjTri[tri->tetVertex] = 0;

                    init_region(pCuspRegion[index], tri, dist, adjTri);
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
                dist[j] = k;
                dist[vertex1] = flow(tri, vertex1) + (flow(tri, j) - dist[j]);
                dist[vertex2] = flow(tri, vertex2) + (flow(tri, j) - dist[j]);
                dist[tri->tetVertex] = -1;

                adjTri[j] = 0;
                adjTri[vertex1] = 1;
                adjTri[vertex2] = 1;
                adjTri[tri->tetVertex] = 0;

                init_region(pCuspRegion[index], tri, dist, adjTri);
                index++;
            }

        }

        // Center vertex
        for (j = 0; j < 4; j++) {
            if (j == tri->tetVertex)
                continue;

            dist[j] = flow(tri, j);
            adjTri[j] = 1;
        }

        dist[tri->tetVertex] = -1;
        adjTri[tri->tetVertex] = 0;

        init_region(pCuspRegion[index], tri, dist, adjTri);
        index++;
    }

    update_adj_regions(pCuspRegion);
    return pCuspRegion;

}

void free_cusp_region(struct CuspRegion **pCuspRegion) {
    int i;

    for (i = 0; i < numCuspRegions; i++)
        my_free(pCuspRegion[i]);
    my_free(pCuspRegion);
}

void init_region(struct CuspRegion *region, struct CuspTriangle *tri, int dist[4], int adjTri[4]) {
    int i, v1, v2, v3;

    region->tri                 = tri;
    region->tetIndex            = region->tri->tetIndex;
    region->tetVertex           = region->tri->tetVertex;

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tetVertex][i];
        v2 = edgesThreeToFour[tri->tetVertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tetVertex][(i + 2) % 3];

        region->curve[v2][v1] = dist[v1];
        region->curve[v3][v1] = dist[v1];
        region->dive[v2][v1] = (dist[v1] ? 0 : 1);
        region->dive[v3][v1] = (dist[v1] ? 0 : 1);

        region->adjTri[v1] = adjTri[v1];
    }
}

void update_adj_regions(struct CuspRegion **pCuspRegion) {
    int i, j;

    // Add adjacent region info
    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i] == NULL || pCuspRegion[i]->tetIndex == -1)
            continue;

        for (j = 0; j < 4; j++) {
            if (!pCuspRegion[i]->adjTri[j] || pCuspRegion[i]->tetVertex == j) {
                pCuspRegion[i]->adjRegions[j] = -1;
                continue;
            }

            pCuspRegion[i]->adjRegions[j] = find_adj_region_index(pCuspRegion, pCuspRegion[i], j);
        }
    }
}

int find_adj_region_index(struct CuspRegion **pCuspRegion, struct CuspRegion *region, int face) {
    int i, vertex1, vertex2, yVertex1, yVertex2, yFace, distV1, distV2, distV3, tetIndex, tetVertex;
    struct CuspTriangle *tri = region->tri;

    vertex1 = (int) remaining_face[tri->tetVertex][face];
    vertex2 = (int) remaining_face[face][tri->tetVertex];

    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i] == NULL || pCuspRegion[i]->tri == tri)
            continue;

        yVertex1    = EVALUATE(tri->tet->gluing[face], vertex1);
        yVertex2    = EVALUATE(tri->tet->gluing[face], vertex2);
        yFace       = EVALUATE(tri->tet->gluing[face], face);

        tetIndex    = (tri->neighbours[face]->tetIndex == pCuspRegion[i]->tetIndex);
        tetVertex   = (tri->neighbours[face]->tetVertex == pCuspRegion[i]->tetVertex);
        distV1      = (region->curve[face][vertex1] == pCuspRegion[i]->curve[yFace][yVertex1]);
        distV2      = (region->curve[face][vertex2] == pCuspRegion[i]->curve[yFace][yVertex2]);
        distV3      = pCuspRegion[i]->adjTri[yFace];

        // missing distance
        if (region->curve[face][vertex1] == -1 ||
        region->curve[face][vertex2] == -1)
            uFatalError("find_adj_region_index", "symplectic_basis");

        if (tetIndex && tetVertex && distV1 && distV2 && distV3)
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
            pDualCurves[i]->curves[j][FINISH]           = NEW_STRUCT(struct EdgeNode);
            pDualCurves[i]->curves[j][START]->next      = pDualCurves[i]->curves[j][FINISH];
            pDualCurves[i]->curves[j][START]->prev      = NULL;
            pDualCurves[i]->curves[j][FINISH]->next     = NULL;
            pDualCurves[i]->curves[j][FINISH]->prev     = pDualCurves[i]->curves[j][START];

            pDualCurves[i]->endpoints[j][START]         = NEW_STRUCT(struct PathEndPoint);
            pDualCurves[i]->endpoints[j][FINISH]        = NEW_STRUCT(struct PathEndPoint);
            pDualCurves[i]->endpoints[j][START]->region = NULL;
            pDualCurves[i]->endpoints[j][FINISH]->region= NULL;
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
            my_free(pDualCurves[i]->endpoints[j][FINISH]);
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

struct Graph *construct_dual_graph(struct CuspRegion **pCuspRegion) {
    int i, index, face;
    struct CuspRegion *region;

    struct Node stack;
    stack.item = -1;

    struct Graph *graph1 = init_graph(numCuspRegions, FALSE);

    int *visited = NEW_ARRAY(graph1->nVertices, int);

    for (i = 0; i < graph1->nVertices; i++)
        visited[i] = 0;

    // Start at the inside corner of triangle 1.
    push(&stack, 0);

    // Walk around the cusp triangulation inserting edges
    while (!is_empty(&stack)) {
        index = pop(&stack);
        region = pCuspRegion[index];

        if (visited[index])
            continue;

        for (face = 0; face < 4; face++) {
            if (!region->adjTri[face])
                continue;

            // Missing adj region data
            if (region->adjRegions[face] == -1)
                uFatalError("construct_dual_graph", "symplectic_basis");

            insert_edge(graph1, index, region->adjRegions[face], graph1->directed);
            graph1->pRegion[index] = region;

            if (!visited[region->adjRegions[face]])
                push(&stack, region->adjRegions[face]);
        }

        visited[index] = 1;
    }

    my_free(visited);
    return graph1;
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
    struct DualCurves *path;

    if (!debug)
        return;

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
            if (pCuspRegion[i] == NULL)
                continue;

            region = pCuspRegion[i];

            v1 = edgesThreeToFour[region->tetVertex][0];
            v2 = edgesThreeToFour[region->tetVertex][1];
            v3 = edgesThreeToFour[region->tetVertex][2];

            printf("Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
                   "(Curves: [[%d %d] [%d %d] [%d %d]]) (Dive: [[%d %d] [%d %d] [%d %d]])\n",
                   i, region->tetIndex, region->tetVertex,
                   region->adjTri[v1], region->adjTri[v2], region->adjTri[v3],
                   region->adjRegions[v1],
                   region->adjRegions[v2],
                   region->adjRegions[v3],
                   region->curve[v2][v1], region->curve[v3][v1],
                   region->curve[v1][v2], region->curve[v3][v2],
                   region->curve[v1][v3], region->curve[v2][v3],
                   region->dive[v2][v1], region->dive[v3][v1],
                   region->dive[v1][v2], region->dive[v3][v2],
                   region->dive[v1][v3], region->dive[v2][v3]
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
    } else if (flag == 4) {
        // Edge indices
        printf("Edge classes\n");
        for (i = 0; i < numCuspTriangles; i++) {
            tri = pTriangle[i];
            v1 = edgesThreeToFour[tri->tetVertex][0];
            v2 = edgesThreeToFour[tri->tetVertex][1];
            v3 = edgesThreeToFour[tri->tetVertex][2];

            printf("(Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
                   "Vertex %d: (%d %d), Vertex %d: (%d %d)\n",
                   tri->tetIndex, tri->tetVertex,
                   v1, tri->vertices[v1].edgeClass, tri->vertices[v1].edgeIndex,
                   v2, tri->vertices[v2].edgeClass, tri->vertices[v2].edgeIndex,
                   v3, tri->vertices[v3].edgeClass, tri->vertices[v3].edgeIndex
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
            if (g->pRegion[i] == NULL)
                continue;

            region = pCuspRegion[i];
            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d): ",
                   i,
                   region->tetIndex,
                   region->tetVertex
            );
            edge = g->edge_list_begin[i];
            while ((edge = edge->next)->next != NULL) {
                printf("%d ", edge->y);
            }
            printf("\n");
        }
    } else if (flag == 8) {
        // end point info
        printf("EndPoint Info\n");
        for (i = 0; i < numDualCurves; i++) {
            path = pDualCurves[i];
            for (j = 0; j < 2; j++) {
                for (k = 0; k < 2; k++) {
                    if (path->endpoints[j][k]->region == NULL)
                        continue;

                    printf("Region: %d (Tet Index: %d, Tet Vertex: %d) Curve num: %d Pos: %d Face: %d Vertex: %d\n",
                    path->endpoints[j][k]->regionIndex, path->endpoints[j][k]->region->tetIndex,
                    path->endpoints[j][k]->region->tetVertex, j, k,
                    path->endpoints[j][k]->face, path->endpoints[j][k]->vertex);
                }
            }
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

struct CuspRegion **construct_dual_curves(struct CuspTriangle **pTriangle, struct CuspRegion **pCuspRegion, struct DualCurves **pDualCurves, int e0) {
    int i, j, *parent;
    bool *processed, *discovered;
    struct Graph *graph1 = construct_dual_graph(pCuspRegion);

    print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 7);

    for (i = 0; i < numDualCurves; i++) {
        if (i == e0)
            continue;


        // which half of the curve
        for (j = 0; j < 2; j++) {
            find_path_endpoints(graph1, pCuspRegion, pDualCurves[i], e0, j, START);
            find_path_endpoints(graph1, pCuspRegion, pDualCurves[i], i, j, FINISH);
            processed = NEW_ARRAY(graph1->nVertices, bool);
            discovered = NEW_ARRAY(graph1->nVertices, bool);
            parent = NEW_ARRAY(graph1->nVertices, int);

            // Find path using bfs
            init_search(graph1, processed, discovered, parent);
            bfs(graph1, pDualCurves[i]->endpoints[j][START]->regionIndex, processed, discovered, parent);
            find_path(
                    pDualCurves[i]->endpoints[j][START]->regionIndex,
                    pDualCurves[i]->endpoints[j][FINISH]->regionIndex, parent, pDualCurves[i]->curves[j][START]);
            update_path_info(pCuspRegion, pDualCurves[i], j);

            print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 8);
            print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 5);

            // Reallocate memory
            my_free(processed);
            my_free(discovered);
            my_free(parent);
            free_graph(graph1);

            // Split the regions along the path
            pCuspRegion = update_cusp_regions(pCuspRegion, pDualCurves[i], j);
            print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 2);

            graph1 = construct_dual_graph(pCuspRegion);
            print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 7);
        }
    }

    my_free(graph1);
    return pCuspRegion;
}

/*
 * Find Path Endpoints
 *
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass.
 */

void find_path_endpoints(struct Graph *g, struct CuspRegion **pCuspRegion, struct DualCurves *path, int edgeClass, int curveNum, int pos) {
    int i, vertex, face1, face2, face;
    struct CuspRegion *pRegion;

    // which cusp region
    for (i = 0; i < g->nVertices; i++) {
        if (g->pRegion[i] == NULL)
            continue;

        pRegion = g->pRegion[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == pRegion->tetVertex)
                continue;

            // does the vertex belong to the correct class
            if (pRegion->tri->vertices[vertex].edgeClass != edgeClass ||
            pRegion->tri->vertices[vertex].edgeIndex != curveNum) {
                continue;
            }

            if (curveNum == 0) {
                // find a valid region and face to dive through
                face1 = (int) remaining_face[pRegion->tetVertex][vertex];
                face2 = (int) remaining_face[vertex][pRegion->tetVertex];

                if (pRegion->dive[face1][vertex])
                    face = face1;
                else if (pRegion->dive[face2][vertex])
                    face = face2;
                else
                    continue;
            } else if (curveNum == 1) {
                // find the region corresponding to curve 0
                if (pRegion->tetIndex != path->endpoints[0][pos]->region->tetIndex ||
                pRegion->tetVertex != path->endpoints[0][pos]->vertex ||
                vertex != path->endpoints[0][pos]->region->tetVertex)
                    continue;

                face = path->endpoints[0][pos]->face;

                if (!pRegion->dive[face][vertex])
                    continue;
            } else {
                // invalid curveNum
                uFatalError("find_path_endpoints", "symplectic_basis");
            }

            path->endpoints[curveNum][pos]->region = pRegion;
            path->endpoints[curveNum][pos]->vertex = vertex;
            path->endpoints[curveNum][pos]->face = face;
            path->endpoints[curveNum][pos]->regionIndex = i;
            return;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_path_endpoints", "symplectic_basis");
}

void update_path_info(struct CuspRegion **pCuspRegion, struct DualCurves *path, int curveNum) {
    int i, insideVertex, face, vertex1, vertex2;
    struct EdgeNode *node = path->curves[curveNum][START];

    // path len 0
    if (node->next->next == NULL)
        return;

    node = node->next;
    // path len 1
    if (node->next->next == NULL) {
        for (face = 0; face < 4; face++)
            if (pCuspRegion[node->y]->tetVertex != face &&
                path->endpoints[curveNum][START]->vertex != face &&
                path->endpoints[curveNum][FINISH]->vertex != face)
                break;

        node->insideVertex = face;
        return;
    }

    // Set Header node
    for (i = 0; i < 4; i++) {
        if (i == pCuspRegion[node->y]->tetVertex)
            continue;

        if (pCuspRegion[node->y]->adjRegions[i] == node->next->y)
            node->face = i;
    }

    vertex1 = (int) remaining_face[pCuspRegion[node->y]->tetVertex][path->endpoints[curveNum][START]->vertex];
    vertex2 = (int) remaining_face[path->endpoints[curveNum][START]->vertex][pCuspRegion[node->y]->tetVertex];

    if (node->face == path->endpoints[curveNum][START]->vertex) {
        node->insideVertex = path->endpoints[curveNum][START]->face == vertex1 ? vertex2 : vertex1;
    } else if (node->face == path->endpoints[curveNum][START]->face) {
        node->insideVertex = -1;
    } else {
        node->insideVertex = path->endpoints[curveNum][START]->vertex;
    }

    while ((node = node->next)->next->next != NULL) {
        inside_vertex(pCuspRegion[node->y], node, &insideVertex, &face);
        node->insideVertex = insideVertex == -1 ? face : insideVertex;
        node->face = face;
    }

    // Set Tail node
    for (i = 0; i < 4; i++) {
        if (i == pCuspRegion[node->y]->tetVertex)
            continue;

        if (pCuspRegion[node->y]->adjRegions[i] == node->prev->y)
            node->face = i;
    }

    vertex1 = (int) remaining_face[pCuspRegion[node->y]->tetVertex][path->endpoints[curveNum][FINISH]->vertex];
    vertex2 = (int) remaining_face[path->endpoints[curveNum][FINISH]->vertex][pCuspRegion[node->y]->tetVertex];

    if (node->face == path->endpoints[curveNum][FINISH]->vertex) {
        node->insideVertex = path->endpoints[curveNum][FINISH]->face == vertex1 ? vertex2 : vertex1;
    } else if (node->face == path->endpoints[curveNum][FINISH]->face) {
        node->insideVertex = -1;
    } else {
        node->insideVertex = path->endpoints[curveNum][FINISH]->vertex;
    }
}

/*
 * Find Edge Holonomies
 *
 */

struct CuspRegion **update_cusp_regions(struct CuspRegion **pCuspRegion, struct DualCurves *path, int curveNum) {
    int i, j, face, pathLen = 0, index;
    struct EdgeNode *node = path->curves[curveNum][START];
    struct CuspRegion *region, *newRegion;

    while ((node = node->next)->next != NULL)
        pathLen++;

    for (index = 0; index < numCuspRegions && pCuspRegion[index]->tetIndex != -1; index++);

    numCuspRegions = index + pathLen;
    struct CuspRegion **newCuspRegion = NEW_ARRAY(numCuspRegions, struct CuspRegion *);

    // Copy pointers to new array
    for (i = 0; i < index; i++)
        newCuspRegion[i] = pCuspRegion[i];

    for (i = index; i < numCuspRegions; i++)
        newCuspRegion[i] = NULL;

    my_free(pCuspRegion);

    node = path->curves[curveNum][START]->next;

    // empty path
    if (pathLen == 0)
        return pCuspRegion;

    // path of len 1
    if (pathLen == 1) {
        newCuspRegion[index] = NEW_STRUCT(struct CuspRegion);
        region = newCuspRegion[node->y];
        newRegion = newCuspRegion[index];
        copy_region(region, newRegion);

        face = node->insideVertex;

        newRegion->adjTri[path->endpoints[curveNum][START]->vertex]  = 0;
        newRegion->adjTri[path->endpoints[curveNum][FINISH]->vertex] = 0;
        newRegion->dive[path->endpoints[curveNum][START]->vertex][path->endpoints[curveNum][FINISH]->vertex] =
                (face != path->endpoints[curveNum][FINISH]->face);
        newRegion->dive[path->endpoints[curveNum][FINISH]->vertex][path->endpoints[curveNum][START]->vertex] =
                (face != path->endpoints[curveNum][START]->face);

        region->adjTri[face] = 0;
        region->dive[face][path->endpoints[curveNum][START]->vertex]  = (face == path->endpoints[curveNum][START]->face);
        region->dive[face][path->endpoints[curveNum][FINISH]->vertex] = (face == path->endpoints[curveNum][FINISH]->face);

        update_adj_regions(newCuspRegion);
        return newCuspRegion;
    }

    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and newRegion to the left of the curve.
     */
    update_cusp_triangle_endpoints(newCuspRegion, newCuspRegion[node->y], path->endpoints[curveNum][START], node);
    newCuspRegion[index] = update_cusp_region_node(newCuspRegion[node->y], node,
                                                   path->endpoints[curveNum][START],-1);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        update_cusp_triangle(newCuspRegion, newCuspRegion[node->y], node);
        newCuspRegion[index] = update_cusp_region_node(newCuspRegion[node->y], node,
                                                       path->endpoints[curveNum][START], 0);
        index++;
    }

    // update last region
    update_cusp_triangle_endpoints(newCuspRegion, newCuspRegion[node->y], path->endpoints[curveNum][FINISH], node);
    newCuspRegion[index] = update_cusp_region_node(newCuspRegion[node->y], node,
                                                   path->endpoints[curveNum][FINISH],1);

    update_adj_regions(newCuspRegion);
    return newCuspRegion;
}

struct CuspRegion *update_cusp_region_node(struct CuspRegion *region, struct EdgeNode *node, struct PathEndPoint *endPoint, int flag) {
    int vertex1, vertex2;
    struct CuspRegion *newRegion = NEW_STRUCT(struct CuspRegion);

    vertex1 = (int) remaining_face[region->tetVertex][node->insideVertex];
    vertex2 = (int) remaining_face[node->insideVertex][region->tetVertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * newRegion becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, newRegion);

    if (flag == 0) {
        // Update new region
        newRegion->curve[vertex1][node->insideVertex]++;
        newRegion->curve[vertex2][node->insideVertex]++;
        newRegion->dive[vertex1][node->insideVertex] = 0;
        newRegion->dive[vertex2][node->insideVertex] = 0;

        // Update region
        region->curve[vertex2][vertex1]++;
        region->curve[vertex1][vertex2]++;
        region->dive[vertex2][vertex1] = 0;
        region->dive[vertex1][vertex2] = 0;
        region->adjTri[node->insideVertex]   = 0;

        return newRegion;
    }

    vertex1 = (int) remaining_face[region->tetVertex][endPoint->vertex];
    vertex2 = (int) remaining_face[endPoint->vertex][region->tetVertex];

    if (node->face == endPoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        newRegion->curve[endPoint->vertex][vertex2]++;
        newRegion->dive[endPoint->vertex][vertex2] = 0;
        newRegion->adjTri[vertex1] = 0;

        region->curve[endPoint->vertex][vertex1]++;
        region->dive[endPoint->vertex][vertex1] = 0;
        region->adjTri[vertex2] = 0;
    } else if (node->face == endPoint->face) {
        // curve passes through the face that carries it
        newRegion->curve[endPoint->face][endPoint->face == vertex1 ? vertex2 : vertex1]++;
        newRegion->dive[endPoint->face][endPoint->face == vertex1 ? vertex2 : vertex1] = 0;
        newRegion->adjTri[endPoint->vertex] = 0;
        newRegion->adjTri[endPoint->face == vertex1 ? vertex2 : vertex1] = 0;

        region->curve[endPoint->face][endPoint->vertex]++;
    } else {
        // Curve goes around the vertex
        newRegion->curve[node->face][endPoint->face]++;
        newRegion->dive[node->face][endPoint->face] = 0;
        newRegion->adjTri[endPoint->face]   = 0;
        newRegion->adjTri[endPoint->vertex] = 0;

        region->curve[node->face][endPoint->vertex]++;
        region->dive[node->face][endPoint->vertex] = 0;
    }

    return newRegion;
}

void copy_region(struct CuspRegion *region1, struct CuspRegion *region2) {
    int i, j;

    if (region1 == NULL || region2 == NULL || region1->tri == NULL)
        uFatalError("copy_region", "symplectic_basis");

    region2->tri            = region1->tri;
    region2->tetIndex       = region1->tetIndex;
    region2->tetVertex      = region1->tetVertex;

    for (i = 0; i < 4; i++) {
        region2->adjTri[i]          = region1->adjTri[i];
        region2->adjRegions[i]      = -1;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]    = region1->curve[i][j];
            region2->dive[i][j]     = region1->dive[i][j];
        }
    }
}

void update_cusp_triangle(struct CuspRegion **pCuspRegion, struct CuspRegion *region, struct EdgeNode *node) {
    int i, face1, face2;

    face1 = (int) remaining_face[region->tetVertex][node->insideVertex];
    face2 = (int) remaining_face[node->insideVertex][region->tetVertex];

    for (i = 0; i < numCuspRegions; i++) {
        // is the region initialised?
        if (pCuspRegion[i] == NULL || pCuspRegion[i]->tetIndex == -1)
            continue;

        // which triangle are we in?
        if (pCuspRegion[i]->tetIndex != region->tetIndex || pCuspRegion[i]->tetVertex != region->tetVertex)
            continue;

        if (pCuspRegion[i]->curve[face1][node->insideVertex] > region->curve[face1][node->insideVertex])
            pCuspRegion[i]->curve[face1][node->insideVertex]++;
        else if (pCuspRegion[i]->curve[face1][node->insideVertex] < region->curve[face1][node->insideVertex])
            pCuspRegion[i]->curve[face1][face2]++;

        if (pCuspRegion[i]->curve[face2][node->insideVertex] > region->curve[face2][node->insideVertex])
            pCuspRegion[i]->curve[face2][node->insideVertex]++;
        else if (pCuspRegion[i]->curve[face2][node->insideVertex] < region->curve[face2][node->insideVertex])
            pCuspRegion[i]->curve[face2][face1]++;
    }
}

void update_cusp_triangle_endpoints(struct CuspRegion **pCuspRegion, struct CuspRegion *region, struct PathEndPoint *endPoint, struct EdgeNode *node) {
    int i, face1, face2;

    face1 = (int) remaining_face[region->tetVertex][endPoint->vertex];
    face2 = (int) remaining_face[endPoint->vertex][region->tetVertex];

    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i] == NULL || pCuspRegion[i]->tetIndex == -1)
            continue;

        // which triangle are we in?
        if (pCuspRegion[i]->tetIndex != region->tetIndex || pCuspRegion[i]->tetVertex != region->tetVertex)
            continue;

        if (node->face == endPoint->vertex) {
            // curve passes through the face opposite the vertex it dives through
            if (pCuspRegion[i]->curve[endPoint->vertex][face1] > region->curve[endPoint->vertex][face1]) {
                pCuspRegion[i]->curve[node->face][face1]++;
                pCuspRegion[i]->dive[node->face][face1] = 0;
            } else if (pCuspRegion[i]->curve[endPoint->vertex][face1] < region->curve[endPoint->vertex][face1]) {
                pCuspRegion[i]->curve[node->face][face2]++;
                pCuspRegion[i]->dive[node->face][face2] = 0;
            }
            continue;
        }

        // Curve goes around the vertex or passes through the face that carries it
        if (pCuspRegion[i]->curve[node->face][endPoint->vertex] > region->curve[node->face][endPoint->vertex]) {
            pCuspRegion[i]->curve[node->face][endPoint->vertex]++;
            pCuspRegion[i]->dive[node->face][endPoint->vertex] = 0;
        } else if (pCuspRegion[i]->curve[node->face][endPoint->vertex] < region->curve[node->face][endPoint->vertex]) {
            pCuspRegion[i]->curve[node->face][node->face == face1 ? face2 : face1]++;
        }
    }

}

int dist(struct CuspRegion *region, int vertex) {
    int v1, v2;

    v1 = (int) remaining_face[region->tetVertex][vertex];
    v2 = (int) remaining_face[vertex][region->tetVertex];

    return MAX(region->curve[v1][vertex], region->curve[v2][vertex]);
}

// ----------------------------------

/*
 * Find Holonomies
 *
 * Construct the symplectic equations from the dual curves
 */

void find_holonomies(struct CuspRegion **pCuspRegion, struct DualCurves **pDualCurves, int e0, int **symp_eqns) {
    int i;

    for (i = 0; i < numDualCurves; i++) {
        if (i == e0)
            continue;

        find_path_holonomy(pCuspRegion, pDualCurves[i], FIRST, symp_eqns[i]);
        find_path_holonomy(pCuspRegion, pDualCurves[i], SECOND, symp_eqns[i]);
    }
}

/*
 * Find Path Holonomies
 *
 * Calculate holonomies along a path and update the row.
 */

void find_path_holonomy(struct CuspRegion **pCuspRegion, struct DualCurves *curve, int curveNum, int *row) {
    int i, index, dirFace;
    struct CuspRegion *midNode;
    struct EdgeNode *node = curve->curves[curveNum][START];
    struct PathEndPoint *endPoint, *pathStartPoint = curve->endpoints[curveNum][START], *pathEndPoint = curve->endpoints[curveNum][FINISH];

    // Edge cases
    if (node->next->next == NULL)
        return;
    else if (node->next->next->next == NULL) {
        // One vertex
        node = node->next;

        for (i = 0; i < 4; i++)
            if (i != pathStartPoint->region->tetVertex && i != pathStartPoint->vertex && i != pathEndPoint->vertex)
                break;

        index = 3 * pathStartPoint->region->tetIndex + edge3_between_faces[pathStartPoint->vertex][pathEndPoint->vertex];
        row[index] = row[index] + orientate_vertex(pCuspRegion[node->y]->tri, i, pathStartPoint->vertex);
        return;
    }

    while((node = node->next)->next != NULL) {
        // Path end points
        if (node->prev->prev == NULL || node->next->next == NULL) {
            // End point index
            if (node->prev->prev == NULL) {
                dirFace = 1;
            } else {
                dirFace = -1;
            }

            if (node->y == pathStartPoint->regionIndex)
                endPoint = pathStartPoint;
            else
                endPoint = pathEndPoint;

            if (endPoint->vertex == node->face) {
                /*
                 * Curve passes through the face opposite the vertex it
                 * dives through. Picks up holonomy for the vertex between
                 * the face that carries the curve and the face the curve
                 * crosses
                 */

                index = 3 * pathStartPoint->region->tetIndex + edge3_between_faces
                [remaining_face[endPoint->region->tetVertex][node->insideVertex]]
                [remaining_face[node->insideVertex][endPoint->region->tetVertex]];

            } else if (endPoint->face == node->face){
                /*
                 * Curve passes through the face that carries it and
                 * thus picks up no holonomy.
                 */
                continue;
            } else {
                /*
                 * Curve picks up the holonomy for the vertex it dives through
                 */
                index = 3 * endPoint->region->tetIndex + edge3_between_faces[endPoint->face][node->face];
            }

            row[index] = row[index] + dirFace * orientate_vertex(pCuspRegion[node->y]->tri, node->insideVertex, node->face);

            continue;
        }

        midNode = pCuspRegion[node->y];

        if (node->insideVertex == node->face) {
            printf("didn't find inside vertex\n");
            continue;
        }

        index = 3 * midNode->tetIndex + edge3_between_faces[remaining_face[midNode->tetVertex][node->insideVertex]]
                [remaining_face[node->insideVertex][midNode->tetVertex]];
        row[index] = row[index] + orientate_vertex(pCuspRegion[node->y]->tri, node->insideVertex, node->face);
    }
}

/*
 * Inside Vertex
 *
 * Given three cusp nodes, which lie on adjacent cusp triangle,
 * find the inside vertex of the path first -> mid -> last
 */

void inside_vertex(struct CuspRegion *midNode, struct EdgeNode *node, int *insideVertex, int *face) {
    int i, vertex1, vertex2;

    for (i = 0; i < 4; i++) {
        if (i == midNode->tetVertex)
            continue;

        vertex1 = (int) remaining_face[midNode->tetVertex][i];
        vertex2 = (int) remaining_face[i][midNode->tetVertex];

        if (midNode->adjRegions[vertex1] == node->next->y && midNode->adjRegions[vertex2] == node->prev->y) {
            *insideVertex = i;
            *face = vertex1;
            return;
        } else if (midNode->adjRegions[vertex2] == node->next->y && midNode->adjRegions[vertex1] == node->prev->y) {
            *insideVertex = i;
            *face = vertex2;
            return;
        }
    }

    for (i = 0; i < 4; i++) {
        if (i == midNode->tetVertex)
            continue;

        if (midNode->adjRegions[i] == node->next->y || midNode->adjRegions[i] == node->prev->y) {
            *insideVertex = -1;
            *face = i;
            return;
        }
    }

    // where does the next node go?
    uFatalError("inside_vertex", "symplectic_basis");
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
//    printf("    Processed edge (%d, %d)\n", x, y);
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

    for (i = 0; i < maxVertices; i++) {
        g->degree[i] = 0;

        g->edge_list_begin[i] = NEW_STRUCT(struct EdgeNode);
        g->edge_list_end[i] = NEW_STRUCT(struct EdgeNode);
        g->edge_list_begin[i]->next = g->edge_list_end[i];
        g->edge_list_begin[i]->prev = NULL;
        g->edge_list_end[i]->next = NULL;
        g->edge_list_end[i]->prev = g->edge_list_begin[i];
        g->pRegion[i] = NULL;
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
