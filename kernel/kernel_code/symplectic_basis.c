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

    construct_dual_curves(pTriangle, pCuspRegion, pDualCurves, e0);
    find_holonomies(pCuspRegion, pDualCurves, e0, symp_eqns);

    print_debug_info(pTriangle, NULL, pCuspRegion, pDualCurves, 5);

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
            pTriangle[i]->vertices[j].edgeIndex     = pTriangle[i]->vertices[j].edge->index;
            pTriangle[i]->vertices[j].vertexIndex   = -1;
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
    int i, j, k, l, index, vertex1, vertex2, vertex3;
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

            for (k = 0; k < 4; k ++)
                pCuspRegion[i]->dive[j][k] = 0;
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

                // Region in the middle of face j
                if (flow(tri, vertex1) && flow(tri, vertex2)) {
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

            // Region of dist 0 to j
            vertex1 = edgesThreeToFour[tri->tetVertex][0];
            vertex2 = edgesThreeToFour[tri->tetVertex][1];
            vertex3 = edgesThreeToFour[tri->tetVertex][2];

            // Case 1; two regions
            if (atleast_two(!flow(tri, vertex1), !flow(tri, vertex2), !flow(tri, vertex3))) {
                region                      = pCuspRegion[index];
                region->tri                 = tri;
                region->tetIndex            = region->tri->tetIndex;
                region->tetVertex           = region->tri->tetVertex;
                region->dist[vertex1]       = flow(tri, vertex1);
                region->dist[vertex2]       = flow(tri, vertex2);
                region->dist[vertex3]       = flow(tri, vertex3);
                region->adjTri[vertex1]     = 1;
                region->adjTri[vertex2]     = 1;
                region->adjTri[vertex3]     = 1;
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

                region                      = pCuspRegion[index];
                region->tri                 = tri;
                region->tetIndex            = region->tri->tetIndex;
                region->tetVertex           = region->tri->tetVertex;
                region->dist[vertex1]       = 0;
                region->dist[vertex2]       = flow(tri, vertex1);
                region->dist[vertex3]       = flow(tri, vertex1);
                region->adjTri[vertex2]     = 1;
                region->adjTri[vertex3]     = 1;
                index++;
            } else {
                // Case 2: three regions
                for (j = 0; j < 4; j++) {
                    if (j == tri->tetVertex)
                        continue;

                    vertex1 = (int) remaining_face[tri->tetVertex][j];
                    vertex2 = (int) remaining_face[j][tri->tetVertex];

                    region                      = pCuspRegion[index];
                    region->tri                 = tri;
                    region->tetIndex            = region->tri->tetIndex;
                    region->tetVertex           = region->tri->tetVertex;
                    region->dist[j]             = 0;
                    region->dist[vertex1]       = flow(tri, j) + flow(tri, vertex1);
                    region->dist[vertex2]       = flow(tri, j) + flow(tri, vertex2);
                    region->adjTri[vertex1]     = 1;
                    region->adjTri[vertex2]     = 1;
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

    // update dive info
    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i]->tetIndex == -1)
            continue;

        for (j = 0; j < 4; j++) {
            if (pCuspRegion[i]->tetVertex == j) {
                continue;
            }

            vertex1 = (int) remaining_face[pCuspRegion[i]->tetVertex][j];
            vertex2 = (int) remaining_face[j][pCuspRegion[i]->tetVertex];

            if (!pCuspRegion[i]->dist[j])
                pCuspRegion[i]->dive[vertex1][j] = 1;
        }
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

void update_adj_regions(struct CuspRegion **pCuspRegion) {
    int i, j;

    // Add adjacent region info
    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i] == NULL || pCuspRegion[i]->tetIndex == -1)
            continue;

        for (j = 0; j < 4; j++) {
            if (!pCuspRegion[i]->adjTri[j] || pCuspRegion[i]->tetVertex == j) {
                continue;
            }

            pCuspRegion[i]->adjRegions[j] = find_adj_region_index(pCuspRegion, pCuspRegion[i], j);
        }
    }
}

int find_adj_region_index(struct CuspRegion **pCuspRegion, struct CuspRegion *region, int face) {
    int i, vertex1, vertex2, yVertex1, yVertex2, yFace, distV1, distV2, distV3, tetIndex, tetVertex, adj;
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
        distV1      = (region->dist[vertex1] * !region->dive[face][vertex1] == pCuspRegion[i]->dist[yVertex1] * !pCuspRegion[i]->dive[yFace][yVertex1]);
        distV2      = (region->dist[vertex2] * !region->dive[face][vertex2] == pCuspRegion[i]->dist[yVertex2] * !pCuspRegion[i]->dive[yFace][yVertex2]);
        distV3      = pCuspRegion[i]->adjTri[yFace];

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

    // init e0 first
    pDualCurves[e0] = NEW_STRUCT(struct DualCurves);

    // which half of the curve
    for (j = 0; j < 2; j++) {
        pDualCurves[e0]->curves[j][START]           = NEW_STRUCT(struct EdgeNode);
        pDualCurves[e0]->curves[j][FINISH]          = NEW_STRUCT(struct EdgeNode);
        pDualCurves[e0]->curves[j][START]->next     = pDualCurves[e0]->curves[j][FINISH];
        pDualCurves[e0]->curves[j][START]->prev     = NULL;
        pDualCurves[e0]->curves[j][FINISH]->next    = NULL;
        pDualCurves[e0]->curves[j][FINISH]->prev    = pDualCurves[e0]->curves[j][START];

        pDualCurves[e0]->endpoints[j][START]        = NEW_STRUCT(struct PathEndPoint);
        pDualCurves[e0]->endpoints[j][FINISH]       = pDualCurves[e0]->endpoints[j][START];
    }

    // which curve
    for (i = 0; i < numDualCurves; i++) {
        if (i == e0)
            continue;

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
    int i, j, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

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
            if (pCuspRegion[i] == NULL)
                continue;

            region = pCuspRegion[i];

            v1 = edgesThreeToFour[region->tetVertex][0];
            v2 = edgesThreeToFour[region->tetVertex][1];
            v3 = edgesThreeToFour[region->tetVertex][2];

            printf("Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) (Dist: %d, %d, %d)\n",
                   i, region->tetIndex, region->tetVertex,
                   region->adjTri[v1], region->adjTri[v2], region->adjTri[v3],
                   region->adjRegions[v1], region->adjRegions[v2], region->adjRegions[v3],
                   region->dist[v1], region->dist[v2], region->dist[v3]
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
                   v1, tri->vertices[v1].edgeIndex, tri->vertices[v1].vertexIndex,
                   v2, tri->vertices[v2].edgeIndex, tri->vertices[v2].vertexIndex,
                   v3, tri->vertices[v3].edgeIndex, tri->vertices[v3].vertexIndex
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
            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d) (Dist: %d, %d, %d): ",
                   i,
                   region->tetIndex,
                   region->tetVertex,
                   region->dist[edgesThreeToFour[region->tetVertex][0]],
                   region->dist[edgesThreeToFour[region->tetVertex][1]],
                   region->dist[edgesThreeToFour[region->tetVertex][2]]
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

void construct_dual_curves(struct CuspTriangle **pTriangle, struct CuspRegion **pCuspRegion, struct DualCurves **pDualCurves, int e0) {
    int i, j, *parent;
    bool *processed, *discovered;
    struct Graph *graph1 = construct_dual_graph(pCuspRegion);

    find_path_endpoints(graph1, pCuspRegion, e0, pDualCurves[e0]->endpoints[FIRST][START],
                        pDualCurves[e0]->endpoints[SECOND][START]);

    print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 7);

    for (i = 0; i < numDualCurves; i++) {
        if (i == 0)
            continue;

        find_path_endpoints(graph1, pCuspRegion, i, pDualCurves[i]->endpoints[FIRST][START],
                            pDualCurves[i]->endpoints[SECOND][START]);

        // which half of the curve
        for (j = 0; j < 2; j++) {
            processed = NEW_ARRAY(graph1->nVertices, bool);
            discovered = NEW_ARRAY(graph1->nVertices, bool);
            parent = NEW_ARRAY(graph1->nVertices, int);

            // Find path using bfs
            init_search(graph1, processed, discovered, parent);
            bfs(graph1, pDualCurves[i]->endpoints[j][START]->graphVertex, processed, discovered, parent);
            find_path(
                    pDualCurves[i]->endpoints[j][START]->graphVertex,
                    pDualCurves[i]->endpoints[j][FINISH]->graphVertex, parent, pDualCurves[i]->curves[j][START]);

            // Reallocate memory
            my_free(processed);
            my_free(discovered);
            my_free(parent);
            my_free(graph1);

            // Split the regions along the path
            pCuspRegion = update_cusp_regions(pCuspRegion, pDualCurves[i], j);
            print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 2);

            graph1 = construct_dual_graph(pCuspRegion);
            print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 7);
        }
    }

    my_free(graph1);
}

/*
 * Find Path Endpoints
 *
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass.
 */

void find_path_endpoints(struct Graph *g, struct CuspRegion **pCuspRegion, int edgeClass, struct PathEndPoint *endPoint0, struct PathEndPoint *endPoint1) {
    int i, vertex, face1, face2, face;
    struct CuspRegion *pRegion;
    bool found = FALSE;

    for (i = 0; i < g->nVertices; i++) {
        if (g->pRegion[i] == NULL)
            continue;

        pRegion = g->pRegion[i];
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == pRegion->tetVertex)
                continue;

            face1 = (int) remaining_face[pRegion->tetVertex][vertex];
            face2 = (int) remaining_face[vertex][pRegion->tetVertex];

            if (pRegion->dive[face1][vertex])
                face = face1;
            else if (pRegion->dive[face2][vertex])
                face = face2;
            else
                continue;

            if (pRegion->tri->vertices[vertex].edgeIndex == edgeClass &&
                pRegion->tri->vertices[vertex].vertexIndex == 0) {
                endPoint0->region = g->pRegion[i];
                endPoint0->vertex = vertex;
                endPoint0->face = face;
                endPoint0->graphVertex = i;
                found = TRUE;
                break;
            }
        }

        if (found)
            break;
    }

    for (i = 0; i < g->nVertices; i++) {
        if (g->pRegion[i] == NULL)
            continue;

        pRegion = g->pRegion[i];
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == pRegion->tetVertex)
                continue;

            if (pRegion->tri->vertices[vertex].edgeIndex == edgeClass &&
                pRegion->tri->vertices[vertex].vertexIndex == 1 &&
                pRegion->tetIndex == endPoint0->region->tetIndex &&
                pRegion->tetVertex == endPoint0->vertex &&
                vertex == endPoint0->region->tetVertex &&
                pRegion->dist[vertex] == 0) {
                endPoint1->region = pRegion;
                endPoint1->vertex = vertex;
                endPoint1->face = endPoint0->face;
                endPoint1->graphVertex = i;
                return;
            }
        }
    }
}

/*
 * Find Edge Holonomies
 *
 */

struct CuspRegion **update_cusp_regions(struct CuspRegion **pCuspRegion, struct DualCurves *path, int curveNum) {
    int i, j, face, pathLen = 0, index, insideVertex, vertex1, vertex2;
    struct EdgeNode *node = path->curves[curveNum][START];
    struct PathEndPoint *endPoint;
    struct CuspRegion *region, *newRegion;

    while ((node = node->next)->next != NULL)
        pathLen++;

    struct CuspRegion **newCuspRegion = NEW_ARRAY(numCuspRegions + pathLen, struct CuspRegion *);

    // Copy pointers to new array
    for (i = 0; i < numCuspRegions; i++)
        newCuspRegion[i] = pCuspRegion[i];

    for (i = numCuspRegions; i < numCuspRegions + pathLen; i++)
        newCuspRegion[i] = NULL;

    my_free(pCuspRegion);

    node = path->curves[curveNum][START]->next;

    // empty path
    if (pathLen == 0)
        return pCuspRegion;

    // path of len 1
    if (pathLen == 1) {
        newCuspRegion[numCuspRegions] = NEW_STRUCT(struct CuspRegion);
        region = newCuspRegion[node->y];
        newRegion = newCuspRegion[numCuspRegions];
        copy_region(region, newRegion);

        for (face = 0; face < 4; face++)
            if (region->tetVertex != face &&
            path->endpoints[curveNum][START]->vertex != face &&
            path->endpoints[curveNum][FINISH]->vertex != face)
                break;

        update_cusp_triangle(newCuspRegion, region, face);

        newRegion->dist[path->endpoints[curveNum][START]->vertex] += (face == path->endpoints[curveNum][START]->face);
        newRegion->dist[path->endpoints[curveNum][FINISH]->vertex] += (face == path->endpoints[curveNum][FINISH]->face);
        newRegion->adjTri[path->endpoints[curveNum][START]->vertex]  = 0;
        newRegion->adjTri[path->endpoints[curveNum][FINISH]->vertex] = 0;

        region->dist[path->endpoints[curveNum][START]->vertex] += (face != path->endpoints[curveNum][START]->face);
        region->dist[path->endpoints[curveNum][FINISH]->vertex] += (face != path->endpoints[curveNum][START]->face);
        region->dive[face][path->endpoints[curveNum][START]->vertex] = (face == path->endpoints[curveNum][START]->face);
        region->dive[face][path->endpoints[curveNum][FINISH]->vertex] = (face == path->endpoints[curveNum][FINISH]->face);

        update_adj_regions(newCuspRegion);
        numCuspRegions++;
        return newCuspRegion;
    }

    i = numCuspRegions;
    numCuspRegions += pathLen;
    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and newRegion to the left of the curve.
     */
    newCuspRegion[i] = update_cusp_region_node(newCuspRegion[node->y], node,
                                               path->endpoints[curveNum][START],
                                               path->endpoints[curveNum][START]->vertex,-1);
    i++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        inside_vertex(newCuspRegion[node->y], node, &insideVertex, &face);
        update_cusp_triangle(newCuspRegion, newCuspRegion[node->y], insideVertex);
        newCuspRegion[i] = update_cusp_region_node(newCuspRegion[node->y], node,
                                                   path->endpoints[curveNum][START],
                                                   insideVertex, 0);
    }

    // update last region
    newCuspRegion[i] = update_cusp_region_node(newCuspRegion[node->y], node,
                                               path->endpoints[curveNum][FINISH],
                                               path->endpoints[curveNum][FINISH]->vertex,1);

    update_adj_regions(newCuspRegion);
    return newCuspRegion;
}

struct CuspRegion *update_cusp_region_node(struct CuspRegion *region, struct EdgeNode *node, struct PathEndPoint *endPoint, int insideVertex, int flag) {
    int face, vertex1, vertex2;
    struct CuspRegion *newRegion = NEW_STRUCT(struct CuspRegion);
    struct EdgeNode *next;

    vertex1 = (int) remaining_face[region->tetVertex][insideVertex];
    vertex2 = (int) remaining_face[insideVertex][region->tetVertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * newRegion becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, newRegion);

    if (flag == 0) {
        // Update new region
        newRegion->dist[insideVertex]++;
        newRegion->dive[vertex1][insideVertex]  = 0;
        newRegion->dive[vertex2][insideVertex]  = 0;

        // Update region
        region->dist[vertex1]++;
        region->dist[vertex2]++;
        region->adjTri[insideVertex]       = 0;
        region->dive[vertex1][vertex2]     = 0;
        region->dive[vertex2][vertex1]     = 0;

        return newRegion;
    }

    if (flag == -1)
        next = node->next;
    else
        next = node->prev;

    for (face = 0; face < 4; face ++) {
        if (face == region->tetVertex)
            continue;

        if (region->adjRegions[face] == next->y)
            break;
    }

    if (face == endPoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        newRegion->dist[vertex2]++;
        newRegion->adjTri[vertex1]                  = 0;
        newRegion->dive[endPoint->vertex][vertex2]  = 0;

        region->dist[vertex1]++;
        region->adjTri[vertex2]                     = 0;
        region->dive[endPoint->vertex][vertex1]     = 0;

    } else if (face == endPoint->face) {
        // curve passes through the face that carries it
        // dont allow this to happen
        uFatalError("update_cusp_regions", "symplectic_basis");

    } else {
        // Curve goes around the vertex
        newRegion->dist[vertex1]++;
        newRegion->dist[vertex2]++;
        newRegion->adjTri[endPoint->face]   = 0;
        newRegion->adjTri[endPoint->vertex] = 0;

        region->dive[endPoint->face == vertex1 ? vertex2 : vertex1][endPoint->vertex] = 0;
        region->dist[endPoint->vertex]++;
    }

    return newRegion;
}

void copy_region(struct CuspRegion *region1, struct CuspRegion *region2) {
    int i, j;

    region2->tri            = region1->tri;
    region2->tetIndex       = region1->tetIndex;
    region2->tetVertex      = region1->tetVertex;

    for (i = 0; i < 4; i++) {
        region2->dist[i]            = region1->dist[i];
        region2->adjTri[i]          = region1->dist[i];
        region2->adjRegions[i]      = region1->dist[i];

        for (j = 0; j < 4; j++)
            region2->dive[i][j]     = region1->dive[i][j];
    }
}

void update_cusp_triangle(struct CuspRegion **pCuspRegion, struct CuspRegion *region, int insideVertex) {
    int i, vertex1, vertex2;

    vertex1 = (int) remaining_face[region->tetVertex][insideVertex];
    vertex2 = (int) remaining_face[insideVertex][region->tetVertex];

    for (i = 0; i < numCuspRegions; i++) {
        if (pCuspRegion[i] == NULL || pCuspRegion[i]->tetIndex == -1)
            continue;

        if (pCuspRegion[i]->tetIndex != region->tetIndex || pCuspRegion[i]->tetVertex != region->tetVertex)
            continue;

        if (pCuspRegion[i]->dist[insideVertex] > region->dist[insideVertex]) {
            pCuspRegion[i]->dist[insideVertex]++;
        } else if (pCuspRegion[i]->dist[insideVertex] == region->dist[insideVertex])
            continue;
        else {
            pCuspRegion[i]->dist[vertex1]++;
            pCuspRegion[i]->dist[vertex2]++;
        }
    }
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
    int i, index, insideVertex, face, dirFace;
    struct CuspRegion *midNode;
    struct EdgeNode *insideNode, *node = curve->curves[curveNum][START];
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
                insideNode = node->next;
                dirFace = 1;
            } else {
                insideNode = node->prev;
                dirFace = -1;
            }

            if (node->y == pathStartPoint->regionIndex)
                endPoint = pathStartPoint;
            else
                endPoint = pathEndPoint;

            // Find the face the next vertex lies across
            for (i = 0; i < 4; i++) {
                if (i == endPoint->region->tetVertex)
                    continue;

                if (endPoint->region->adjRegions[i] == insideNode->y) {
                    face = i;
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

                if (remaining_face[endPoint->vertex][endPoint->region->tetVertex] == endPoint->vertex) {
                    insideVertex = (int) remaining_face[endPoint->region->tetVertex][endPoint->vertex];
                }
                else {
                    insideVertex = (int) remaining_face[endPoint->vertex][endPoint->region->tetVertex];
                }

                index = 3 * pathStartPoint->region->tetIndex + edge3_between_faces
                [remaining_face[endPoint->region->tetVertex][insideVertex]]
                [remaining_face[insideVertex][endPoint->region->tetVertex]];

            } else if (endPoint->face == face){
                /*
                 * Curve passes through the face that carries it and
                 * thus picks up no holonomy.
                 */
                continue;
            } else {
                /*
                 * Curve picks up the holonomy for the vertex it dives through
                 */
                insideVertex = endPoint->vertex;
                index = 3 * endPoint->region->tetIndex + edge3_between_faces[endPoint->face][face];
            }

            row[index] = row[index] + dirFace * orientate_vertex(pCuspRegion[node->y]->tri, insideVertex, face);

            continue;
        }

        midNode = pCuspRegion[node->y];
        inside_vertex(midNode, node, &insideVertex, &face);

        if (insideVertex == -1) {
            printf("didn't find inside vertex\n");
            continue;
        }

        index = 3 * midNode->tetIndex + edge3_between_faces[remaining_face[midNode->tetVertex][insideVertex]][remaining_face[insideVertex][midNode->tetVertex]];
        row[index] = row[index] + orientate_vertex(pCuspRegion[node->y]->tri, insideVertex, face);
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

    *insideVertex = -1;
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
