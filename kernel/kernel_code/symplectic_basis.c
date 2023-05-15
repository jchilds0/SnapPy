/*
 *  Symplectic Basis
 *
 *  Computes the symplectic basis of a cusped 3-manifold
*/

#include <stdbool.h>
#include <stdio.h>
#include "kernel.h"
#include "kernel_namespace.h"
#include "addl_code.h"

#define atleast_two(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))
#define orientate_vertex(tri, vertex, face)       ((face) == remaining_face[vertex][(tri)->tetVertex] ? 1 : -1)

#define FIRST                   0
#define SECOND                  1
#define START                   0
#define FINISH                  1


/*
 * Queue
 */

struct Queue {
    int     front;      // First element of queue
    int     rear;       // Last element of queue
    int     len;        // num of elements
    int     size;       // array size
    int     *array;
};


/**
 * Cusp Triangulation
 *
 * CuspTriangle stores information about a triangle in the
 * cusp triangulation. The homology curves bound a fundamental
 * domain, and cusp regions store the information for
 * intersection of this domain with each cusp triangle. When
 * we add oscillating curves, these regions are divided further.
*/

struct CuspVertex {
    int                     edgeClass;
    int                     edgeIndex;
    EdgeClass               *edge;
    int                     v1;
    int                     v2;
};

struct CuspTriangle {
    Tetrahedron             *tet;                   /** tetrahedron the triangle comes from */
    int                     tetIndex;               /** tet->index */
    int                     tetVertex;              /** vertex the triangle comes from */
    int                     numCurves;              /** number of homology curves on the triangle */
    struct CuspVertex       vertices[4];            /** information about each vertex */
    struct CuspTriangle     *neighbours[4];         /** triangle neighbouring a face */
    struct CuspTriangle     *next;                  /** next cusp triangle on doubly linked list */
    struct CuspTriangle     *prev;                  /** prev cusp triangle on doubly linkled list */
};

struct CuspRegion {
    struct CuspTriangle     *tri;                   /** cusp triangle the region lies on */
    int                     tetIndex;               /** tri->tetIndex */
    int                     tetVertex;              /** tri->tetVertex */
    int                     index;                  /** index of the cusp region */
    int                     curve[4][4];            /** looking at face, number of curves between the region and vertex */
    int                     adjTri[4];              /** does the region meet this edge of the cusp triangle */
    struct CuspRegion       *adjRegions[4];         /** index of the adjacent regions */
    int                     dive[4][4];             /** can we dive along the face into this vertex */
    struct CuspRegion       *next;                  /** next cusp region in doubly linked list */
    struct CuspRegion       *prev;                  /** prev cusp region in doubly linked list */
};

struct ManifoldBoundary {
    int                     intersectTetIndex;      /** index of the intersection triangle */
    int                     intersectTetVertex;     /** vertex of the intersection triangle */
    int                     numEdgeClasses;         /** number of edge classes in the boundary */
    int                     numCuspTriangles;       /** number of cusp triangle in the boundary */
    int                     numCuspRegions;         /** number of cusp regions in the boundary */
    int                     numDualCurves;          /** number of dual curves in the boundary */
    Cusp                    *cusp;                  /** which cusp is the boundary in */
    struct Graph            *dual_graph;            /** dual graph of the cusp region */
    struct CuspTriangle     *cusp_triangle_begin;   /** header node of doubly linked list of cusp triangles */
    struct CuspTriangle     *cusp_triangle_end;     /** tail node of doubly linked list of cusp triangles */
    struct CuspRegion       *cusp_region_begin;     /** header node of doubly linked list of cusp regions */
    struct CuspRegion       *cusp_region_end;       /** tail node of doubly linked list of cusp regions */
    struct DualCurves       *dual_curve_begin;      /** header node of doubly linked list of dual curves */
    struct DualCurves       *dual_curve_end;        /** tail node of doubly linked list of dual curves */
};

/**
 * Dual Curves
 *
 * Each oscillating curve is made up of two components, accessed using macros
 * FIRST and SECOND. For each component we have the starting endpoint,
 * *endpoints[curveNum][START] and finish endpoint *endpoint[curveNum][FINISH].
 * The path of the curve is stored as a double linked list with header and tail nodes,
 * the header is curves[curveNum][START] and the tail is curve[curveNum][FINISH].
 */

struct PathEndPoint {
    int                     face;                   /** face containg the short rectangle carrying the curve */
    int                     vertex;                 /** vertex we dive through the manifold along */
    int                     regionIndex;            /** index of the region the endpoint lies in */
    struct CuspRegion       *region;                /** pointer to the region the endpoint lies in */
};

struct DualCurves {
    int                     edgeClass;              /** which edge class does the curve dive through */
    struct EdgeNode         *curves[2][2];          /** matrix of curves */
    struct PathEndPoint     *endpoints[2][2];       /** matrix of endpoints */
    struct DualCurves       *next;                  /** next dual curve in doubly linked list */
    struct DualCurves       *prev;                  /** prev dual curve in doubly linked list */
};

/**
 * Graph
 */

struct EdgeNode {
    int                     y;                      /** cusp region index */
    int                     nextFace;               /** face the path crosses to the next node */
    int                     prevFace;               /** face the path crosses to the prev node */
    int                     insideVertex;           /** inside vertex of the path */
    int                     intermediate;           /** (end multi graph) which vertex lies between two vertices */
    struct EdgeNode         *next;                  /** next node in doubly linked list */
    struct EdgeNode         *prev;                  /** prev node in doubly linked list */
};

struct Graph {
    struct EdgeNode         **edge_list_begin;      /** header node of doubly linked list */
    struct EdgeNode         **edge_list_end;        /** tail node of doubly linked list */
    struct CuspRegion       **pRegion;              /** list of regions in the graph */
    int                     *degree;                /** degree of each vertex */
    int                     nVertices;              /** number of vertices in the graph */
    int                     directed;               /** is the graph directed */
};

struct EndMultiGraph {
    int                     e0;                     /** base edge class */
    struct Graph            *multi_graph;           /** tree with extra edge of cusps */
    struct Graph            *double_graph;          /** double of the multi graph for finding paths of even length */
};

// Graph
struct Graph            *init_graph(int, bool);
void                    free_graph(struct Graph *);
int                     insert_edge(struct Graph *, int, int, bool);
void                    delete_edge(struct Graph *, int, int, bool);
int                     edge_exists(struct Graph *, int, int);

// Symplectic Basis
int                     **get_symplectic_equations(Triangulation *, int, int);

struct ManifoldBoundary *init_boundary(Triangulation *, Cusp *, int);
void                    free_boundary(struct ManifoldBoundary **, int);
void                    init_cusp_triangulation(Triangulation *, struct ManifoldBoundary *);
int                     flow(struct CuspTriangle *, int);
void                    label_triangulation_edges(Triangulation *);
struct CuspTriangle     *find_cusp_triangle(struct CuspTriangle *, struct CuspTriangle *, struct CuspTriangle *, int);
void                    label_cusp_vertex_indices(struct CuspTriangle *, struct CuspTriangle *, int);
void                    walk_around_cusp_vertex(struct CuspTriangle *, int, int);
void                    init_cusp_region(struct ManifoldBoundary *);
int                     init_intersect_cusp_region(struct ManifoldBoundary *, struct CuspTriangle *, int);
int                     init_normal_cusp_region(struct ManifoldBoundary *, struct CuspTriangle *, int);
void                    set_cusp_region_data(struct CuspRegion *, struct CuspTriangle *, int [4], int [4], int);
void                    update_adj_region_data(struct CuspRegion *, struct CuspRegion *);
struct CuspRegion       *find_adj_region(struct CuspRegion *, struct CuspRegion *, struct CuspRegion *, int);
void                    init_oscillating_curves(struct ManifoldBoundary *);
void                    free_oscillating_curves(struct DualCurves *);
void                    find_intersection_triangle(Triangulation *, struct ManifoldBoundary *);

/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    construct_oscillating_curves(Triangulation *, struct ManifoldBoundary **, int);
struct Graph *          construct_cusp_region_dual_graph(struct ManifoldBoundary *);
void                    print_debug_info(struct ManifoldBoundary **, int, int);
void                    find_path_endpoints(struct Graph *, struct PathEndPoint *, struct PathEndPoint *, int, int, bool);
void                    update_path_info(struct Graph *g, struct DualCurves *, int);
void                    split_cusp_regions_along_path(struct ManifoldBoundary *, struct DualCurves *, int);
struct CuspRegion       *update_cusp_region(struct CuspRegion *region, struct EdgeNode *, struct PathEndPoint *, int, int);
void                    update_cusp_triangle(struct CuspRegion *, struct CuspRegion *, struct CuspRegion *, struct EdgeNode *);
void                    update_cusp_triangle_endpoints(struct CuspRegion *, struct CuspRegion *, struct CuspRegion *, struct PathEndPoint *, struct EdgeNode *);
void                    copy_region(struct CuspRegion *, struct CuspRegion *);
void                    calculate_holonomy(Triangulation *, struct ManifoldBoundary **, int, int **);
void                    find_path_holonomy(struct Graph *, struct DualCurves *, int, int *);
void                    inside_vertex(struct CuspRegion *, struct EdgeNode *);

/**
 * Queue Data Structure
 */

void                    initialise_queue(struct Queue *, int);
struct Queue            *enqueue(struct Queue *, int);
int                     dequeue(struct Queue *);
void                    resize_queue(struct Queue *);
int                     empty_queue(struct Queue *);
void                    free_queue(struct Queue *);

/**
 * Graph for Breadth First Search
 */

void                    init_search(struct Graph *, bool *, bool *, int *);
void                    bfs(struct Graph *, int, bool *, bool *, int *);
void                    find_path(int, int, int *, struct EdgeNode *);

/**
 * Spanning Tree for End Multi Graph
 */

struct EndMultiGraph    *end_multi_graph(Triangulation *);
void                    spanning_tree(struct Graph *, struct Graph *, int, int *);
void                    cusp_graph(Triangulation *, struct Graph *);
void                    add_odd_cycle_edge(struct Graph *, int *);
int                     find_path_len(int, int, int *, int);
void                    construct_double_graph(struct Graph *, struct Graph *);

int edgesThreeToFour[4][3] = {{1, 2, 3},
                              {0, 2, 3},
                              {0, 1, 3},
                              {0, 1, 2}};

/*
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the get_gluing_equations and get_symplectic_equations functions
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i, j;

    peripheral_curves(manifold);

    // Edge Curves C_i -> gluing equations
    int **edge_eqns, edge_num_rows;

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;

    // Get Gluing Equations
    edge_eqns = get_gluing_equations(manifold, &edge_num_rows, num_cols);
    symp_num_rows = edge_num_rows;

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, symp_num_rows, *num_cols);

    // Construct return array
    *num_rows = edge_num_rows + symp_num_rows - 2;
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < edge_num_rows; i++) {
        if (i == 0) {
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
 * Setup graph and cusp triangulation, and run construct dual curves.
 */
static int debug = FALSE;

int **get_symplectic_equations(Triangulation *manifold, int num_rows, int numCols) {
    int i, j, e0;
    struct ManifoldBoundary **cusps     = NEW_ARRAY(manifold->num_cusps, struct ManifoldBoundary *);
    struct Cusp *cusp                   = manifold->cusp_list_begin.next;
    struct Graph *end_graph             = init_graph(manifold->num_cusps, FALSE);
    struct EndMultiGraph *multiGraph    = end_multi_graph(manifold);

    label_triangulation_edges(manifold);
    for (i = 0; i < manifold->num_cusps; i++) {
        cusps[i] = init_boundary(manifold, cusp, e0);
        cusp = cusp->next;
    }

    print_debug_info(cusps, e0, 0);       // Gluing
    print_debug_info(cusps, e0, 3);       // Homology
    print_debug_info(cusps, e0, 4);       // Edge classes
    print_debug_info(cusps, e0, 6);       // Inside Edge
    print_debug_info(cusps, e0, 2);       // Regions

    // Allocate Symplectic Equations Array
    int **symp_eqns = NEW_ARRAY(num_rows, int *);

    for (i = 0; i < num_rows; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    construct_oscillating_curves(manifold, cusps, e0);
    calculate_holonomy(manifold, cusps, e0, symp_eqns);

    print_debug_info(cusps, e0, 1);
    print_debug_info(cusps, e0, 5);

    free_boundary(cusps, manifold->num_cusps);
    my_free(cusps);
    free_graph(end_graph);
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
 * For each cusp, initalise the data structures
 * to describe the cusp
 */

struct ManifoldBoundary *init_boundary(Triangulation *manifold, Cusp *cusp, int e0) {
    struct ManifoldBoundary *boundary = NEW_STRUCT(struct ManifoldBoundary);

    // Invalid cusp topology
    if (cusp->topology == Klein_cusp)
        uFatalError("init_boundary", "symplectic_basis");

    boundary->cusp = cusp;
    boundary->numEdgeClasses = manifold->num_tetrahedra;
    boundary->numCuspTriangles = 0;
    boundary->numCuspRegions = 0;
    boundary->numDualCurves = manifold->num_tetrahedra;

    find_intersection_triangle(manifold, boundary);
    init_cusp_triangulation(manifold, boundary);
    init_cusp_region(boundary);
    init_oscillating_curves(boundary);

    return boundary;
}

void free_boundary(struct ManifoldBoundary **cusps, int numCusps) {
    int i;
    struct CuspTriangle *tri;
    struct CuspRegion *region;
    struct DualCurves *path;
    struct ManifoldBoundary *boundary;

    for (i = 0; i < numCusps; i++) {
        boundary = cusps[i];
        // free graph
        free_graph(boundary->dual_graph);

        // free dual curves
        while (boundary->dual_curve_begin->next != boundary->dual_curve_end) {
            path = boundary->dual_curve_begin->next;
            REMOVE_NODE(path);
            free_oscillating_curves(path);
        }

        // free cusp regions
        while (boundary->cusp_region_begin->next != boundary->cusp_region_end) {
            region = boundary->cusp_region_begin->next;
            REMOVE_NODE(region);
            my_free(region);
        }

        // free cusp triangle
        while (boundary->cusp_triangle_begin->next != boundary->cusp_triangle_end) {
            tri = boundary->cusp_triangle_begin->next;
            REMOVE_NODE(tri);
            my_free(tri);
        }

        my_free(boundary->cusp_triangle_begin);
        my_free(boundary->cusp_triangle_end);
        my_free(boundary->cusp_region_begin);
        my_free(boundary->cusp_region_end);
        my_free(boundary->dual_curve_begin);
        my_free(boundary->dual_curve_end);
        my_free(boundary);
    }
}

/*
 * peripheral_curves.c places a meridian and
 * longitude curve on each cusp. It starts at
 * a base triangle, the intersection point, and
 * searches outwards. We need information about
 * the intersection triangle for creating the
 * cusp regions, so we assume the intersection
 * triangle only contains eactly the curves which
 * intersect
 */

void find_intersection_triangle(Triangulation *manifold, struct ManifoldBoundary *boundary) {
    FaceIndex   face;
    Cusp *cusp = boundary->cusp;
    int n;

    for (cusp->basepoint_tet = manifold->tet_list_begin.next;
         cusp->basepoint_tet != &manifold->tet_list_end;
         cusp->basepoint_tet = cusp->basepoint_tet->next)

        for (cusp->basepoint_vertex = 0;
             cusp->basepoint_vertex < 4;
             cusp->basepoint_vertex++)
        {
            if (cusp->basepoint_tet->cusp[cusp->basepoint_vertex] != cusp)
                continue;

            for (face = 0; face < 4; face++)
            {
                if (face == cusp->basepoint_vertex)
                    continue;

                for (n = 0; n < 2; n++) {
                    cusp->basepoint_orientation = ORIENTATION(n);

                    if (cusp->basepoint_tet->curve
                        [M]
                        [cusp->basepoint_orientation]
                        [cusp->basepoint_vertex]
                        [face] != 0
                        && cusp->basepoint_tet->curve
                           [L]
                           [cusp->basepoint_orientation]
                           [cusp->basepoint_vertex]
                           [face] != 0) {
                        /*
                         *  We found the basepoint!
                         */

                        boundary->intersectTetIndex = cusp->basepoint_tet->index;
                        boundary->intersectTetVertex = (int) cusp->basepoint_vertex;
                        return;
                    }


                }
            }
        }
}

/*
 * Construct the cusp triangle doubly
 * linked list which consists of the
 * triangles in the cusp triangulation
 */

void init_cusp_triangulation(Triangulation *manifold, struct ManifoldBoundary *boundary) {
    int vertex, face, index = 0;
    Tetrahedron *tet;
    struct CuspTriangle *tri;

    // Allocate Cusp Triangulation Header and Tail Null nodes
    boundary->cusp_triangle_begin = NEW_STRUCT( struct CuspTriangle );
    boundary->cusp_triangle_end = NEW_STRUCT( struct CuspTriangle );

    boundary->cusp_triangle_begin->next = boundary->cusp_triangle_end;
    boundary->cusp_triangle_begin->prev = NULL;
    boundary->cusp_triangle_end->next = NULL;
    boundary->cusp_triangle_end->prev = boundary->cusp_triangle_begin;

    // which tetrahedron are we on
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // while vertex are we on
        for (vertex = 0; vertex < 4; vertex++) {
            // is this vertex on the right cusp
            if (tet->cusp[vertex] != boundary->cusp) {
                continue;
            }

            tri = NEW_STRUCT( struct CuspTriangle );
            INSERT_BEFORE(tri, boundary->cusp_triangle_end);
            index++;

            tri->tet = tet;
            tri->tetIndex = tri->tet->index;
            tri->tetVertex = vertex;

            tri->numCurves = flow(tri, edgesThreeToFour[tri->tetVertex][0])
                             + flow(tri, edgesThreeToFour[tri->tetVertex][1])
                             + flow(tri, edgesThreeToFour[tri->tetVertex][2]);

            for (face = 0; face < 4; face ++) {
                if (tri->tetVertex == face)
                    continue;

                tri->vertices[face].v1 = tri->tetVertex;
                tri->vertices[face].v2 = face;
                tri->vertices[face].edge = tri->tet->edge_class[edge_between_vertices[tri->vertices[face].v1][tri->vertices[face].v2]];
                tri->vertices[face].edgeClass = tri->vertices[face].edge->index;
                tri->vertices[face].edgeIndex = -1;
            }
        }
    }

    // which cusp triangle
    for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
        // which vertex
        for (face = 0; face < 4; face++) {
            if (face == tri->tetVertex)
                continue;

            tri->neighbours[face] = find_cusp_triangle(boundary->cusp_triangle_begin, boundary->cusp_triangle_end, tri, face);
        }
    }

    label_cusp_vertex_indices(boundary->cusp_triangle_begin, boundary->cusp_triangle_end, boundary->numEdgeClasses);
    boundary->numCuspTriangles = index;
}

/*
 * Give each edge of the triangulation an index to
 * identify the cusp vertices
 */

void label_triangulation_edges(Triangulation *manifold) {
    int i = 0;
    EdgeClass *edge = &manifold->edge_list_begin;

    while ((edge = edge->next)->next != NULL)
        edge->index = i++;

    // incorrect number of edge classes
    if (i != manifold->num_tetrahedra)
        uFatalError("label_triangulation_edges", "symplectic_basis");
}

/*
 * Each edge class of the manifold appears as two vertices
 * in the cusp triangulation. We iterate over the cusp triangulation,
 * walking around each vertex to give it the same index.
 */

void label_cusp_vertex_indices(struct CuspTriangle *cusp_triangle_begin, struct CuspTriangle *cusp_triangle_end, int numEdgeClasses) {
    int i, vertex;
    struct CuspTriangle *tri;

    int *currentIndex = NEW_ARRAY(numEdgeClasses, int);

    for (i = 0; i < numEdgeClasses; i++)
        currentIndex[i] = 0;

    for (tri = cusp_triangle_begin->next; tri != cusp_triangle_end; tri = tri->next) {
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tetVertex || tri->vertices[vertex].edgeIndex != -1)
                continue;

            walk_around_cusp_vertex(tri, vertex, currentIndex[tri->vertices[vertex].edgeClass]);
            currentIndex[tri->vertices[vertex].edgeClass]++;
        }
    }

    my_free(currentIndex);
}

/*
 * Walk around vertex cuspVertex of triangle *tri
 * and set edgeIndex to index.
 */

void walk_around_cusp_vertex(struct CuspTriangle *tri, int cuspVertex, int index) {
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
        tri                 = tri->neighbours[old_gluing_vertex];
    }
}

/*
 * Initialise the cusp region doubly linked list
 * to cotain the regions bounded by the meridian
 * and longitude curves.
 */

void init_cusp_region(struct ManifoldBoundary *boundary) {
    int index = 0;
    struct CuspTriangle *tri;

    // Header and tailer nodes.
    boundary->cusp_region_begin         = NEW_STRUCT( struct CuspRegion );
    boundary->cusp_region_end           = NEW_STRUCT( struct CuspRegion );
    boundary->cusp_region_begin->next   = boundary->cusp_region_end;
    boundary->cusp_region_begin->prev   = NULL;
    boundary->cusp_region_end->next     = NULL;
    boundary->cusp_region_end->prev     = boundary->cusp_region_begin;

    for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
        // Intersection vertex doesn't have a center
        if (tri->tetIndex == boundary->intersectTetIndex && tri->tetVertex == boundary->intersectTetVertex) {
            index = init_intersect_cusp_region(boundary, tri, index);
            continue;
        }

        index = init_normal_cusp_region(boundary, tri, index);
    }

    update_adj_region_data(boundary->cusp_region_begin, boundary->cusp_region_end);
    boundary->numCuspRegions = index;
}

int init_intersect_cusp_region(struct ManifoldBoundary *boundary, struct CuspTriangle *tri, int index) {
    int i, curveNum, vertex, vertex1, vertex2, vertex3;
    int adjTri[4], dist[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tetVertex) {
            continue;
        }

        vertex1 = (int) remaining_face[tri->tetVertex][vertex];
        vertex2 = (int) remaining_face[vertex][tri->tetVertex];

        for (i = 1; i < flow(tri, vertex); i++) {
            for (curveNum = 0; curveNum < 2; curveNum++) {
                dist[vertex1]           = i;
                dist[vertex2]           = MIN(dist[vertex1], 2 * flow(tri, vertex) - dist[vertex1])
                                          + flow(tri, vertex2) + flow(tri, vertex1);
                dist[vertex]            = flow(tri, vertex) - dist[vertex1] + flow(tri, vertex1);
                dist[tri->tetVertex]    = -1;

                adjTri[vertex1]         = 1;
                adjTri[vertex2]         = 0;
                adjTri[vertex]          = 0;
                adjTri[tri->tetVertex]  = -1;

                set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
                index++;

                // Swap vertices
                vertex1 = (int) remaining_face[vertex][tri->tetVertex];
                vertex2 = (int) remaining_face[tri->tetVertex][vertex];
            }
        }

        // Region in the middle of face vertex
        if (flow(tri, vertex1) && flow(tri, vertex2)) {
            dist[vertex1]           = flow(tri, vertex2);
            dist[vertex2]           = flow(tri, vertex1);
            dist[vertex]            = MIN(flow(tri, vertex1) + dist[vertex1],
                                          flow(tri, vertex2) + dist[vertex2]) + flow(tri, vertex);
            dist[tri->tetVertex]    = -1;

            adjTri[vertex1]         = 0;
            adjTri[vertex2]         = 0;
            adjTri[vertex]          = 1;
            adjTri[tri->tetVertex]  = -1;

            set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }
    }

    // Region of dist 0 to vertex
    vertex1 = edgesThreeToFour[tri->tetVertex][0];
    vertex2 = edgesThreeToFour[tri->tetVertex][1];
    vertex3 = edgesThreeToFour[tri->tetVertex][2];

    // Case 1; two regions
    if (atleast_two(!flow(tri, vertex1), !flow(tri, vertex2), !flow(tri, vertex3))) {
        dist[vertex1]           = flow(tri, vertex1);
        dist[vertex2]           = flow(tri, vertex2);
        dist[vertex3]           = flow(tri, vertex3);
        dist[tri->tetVertex]    = -1;

        adjTri[vertex1]         = 1;
        adjTri[vertex2]         = 1;
        adjTri[vertex3]         = 1;
        adjTri[tri->tetVertex]  = -1;

        set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
        index++;

        // Find vertex with non zero flow
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tetVertex)
                continue;

            if (flow(tri, vertex)) {
                vertex1 = vertex;
                vertex2 = (int) remaining_face[tri->tetVertex][vertex1];
                vertex3 = (int) remaining_face[vertex1][tri->tetVertex];
                break;
            }
        }
        dist[vertex1]           = 0;
        dist[vertex2]           = flow(tri, vertex1);
        dist[vertex3]           = flow(tri, vertex1);
        dist[tri->tetVertex]    = -1;

        adjTri[vertex1]         = 0;
        adjTri[vertex2]         = 1;
        adjTri[vertex3]         = 1;
        adjTri[tri->tetVertex]  = 0;

        set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
        index++;
    } else {
        // Case 2: three regions
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tetVertex)
                continue;

            vertex1 = (int) remaining_face[tri->tetVertex][vertex];
            vertex2 = (int) remaining_face[vertex][tri->tetVertex];

            dist[vertex]            = 0;
            dist[vertex1]           = flow(tri, vertex) + flow(tri, vertex1);
            dist[vertex2]           = flow(tri, vertex) + flow(tri, vertex2);
            dist[tri->tetVertex]    = -1;

            adjTri[vertex]          = 0;
            adjTri[vertex1]         = 1;
            adjTri[vertex2]         = 1;
            adjTri[tri->tetVertex]  = 0;

            set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }
    }

    return index;
}

int init_normal_cusp_region(struct ManifoldBoundary *boundary, struct CuspTriangle *tri, int index) {
    int i, vertex, vertex1, vertex2;
    int dist[4], adjTri[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tetVertex) {
            continue;
        }

        vertex1 = (int) remaining_face[tri->tetVertex][vertex];
        vertex2 = (int) remaining_face[vertex][tri->tetVertex];

        for (i = 0; i < flow(tri, vertex); i++) {
            dist[vertex] = i;
            dist[vertex1] = flow(tri, vertex1) + (flow(tri, vertex) - dist[vertex]);
            dist[vertex2] = flow(tri, vertex2) + (flow(tri, vertex) - dist[vertex]);
            dist[tri->tetVertex] = -1;

            adjTri[vertex] = 0;
            adjTri[vertex1] = 1;
            adjTri[vertex2] = 1;
            adjTri[tri->tetVertex] = 0;

            set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }

    }

    // Center vertex
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tetVertex)
            continue;

        dist[vertex] = flow(tri, vertex);
        adjTri[vertex] = 1;
    }

    dist[tri->tetVertex] = -1;
    adjTri[tri->tetVertex] = 0;

    set_cusp_region_data(boundary->cusp_region_end, tri, dist, adjTri, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which
 * allocates the attributes of the cusp region
 */

void set_cusp_region_data(struct CuspRegion *cusp_region_end, struct CuspTriangle *tri, int dist[4], int adjTri[4], int index) {
    int i, j, v1, v2, v3;
    struct CuspRegion *pRegion = NEW_STRUCT( struct CuspRegion );
    INSERT_BEFORE(pRegion, cusp_region_end);

    pRegion->tri                 = tri;
    pRegion->tetIndex            = pRegion->tri->tetIndex;
    pRegion->tetVertex           = pRegion->tri->tetVertex;
    pRegion->index               = index;

    // default values
    for (i = 0; i < 4; i++) {
        pRegion->adjTri[i] = 0;

        for (j = 0; j < 4; j++) {
            pRegion->curve[i][j] = -1;
            pRegion->dive[i][j] = 0;
        }
    }

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tetVertex][i];
        v2 = edgesThreeToFour[tri->tetVertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tetVertex][(i + 2) % 3];

        pRegion->curve[v2][v1] = dist[v1];
        pRegion->curve[v3][v1] = dist[v1];
        pRegion->dive[v2][v1] = (dist[v1] ? 0 : 1);
        pRegion->dive[v3][v1] = (dist[v1] ? 0 : 1);

        pRegion->adjTri[v1] = adjTri[v1];
    }
}

/*
 * Calculate which regions are located across cusp edges
 * and store the result in the adjRegions attribute
 */

void update_adj_region_data(struct CuspRegion *cusp_region_begin, struct CuspRegion *cusp_region_end) {
    int j;
    struct CuspRegion *pRegion;

    // Add adjacent region info
    for (pRegion = cusp_region_begin->next; pRegion != cusp_region_end; pRegion = pRegion->next) {
        for (j = 0; j < 4; j++) {
            if (!pRegion->adjTri[j] || pRegion->tetVertex == j) {
                pRegion->adjRegions[j] = NULL;
                continue;
            }

            pRegion->adjRegions[j] = find_adj_region(cusp_region_begin, cusp_region_end, pRegion, j);
        }
    }
}

/*
 * Find the cusp region which is adjacent
 * to xRegion across face.
 */

struct CuspRegion *find_adj_region(struct CuspRegion *cusp_region_begin, struct CuspRegion *cusp_region_end,
        struct CuspRegion *xRegion, int face) {
    int vertex1, vertex2, yVertex1, yVertex2, yFace, distV1, distV2, distV3, tetIndex, tetVertex;
    struct CuspTriangle *tri = xRegion->tri;
    struct CuspRegion *yRegion;

    vertex1 = (int) remaining_face[tri->tetVertex][face];
    vertex2 = (int) remaining_face[face][tri->tetVertex];

    for (yRegion = cusp_region_begin->next; yRegion != cusp_region_end; yRegion = yRegion->next) {
        yVertex1    = EVALUATE(tri->tet->gluing[face], vertex1);
        yVertex2    = EVALUATE(tri->tet->gluing[face], vertex2);
        yFace       = EVALUATE(tri->tet->gluing[face], face);

        tetIndex    = (tri->neighbours[face]->tetIndex == yRegion->tetIndex);
        tetVertex   = (tri->neighbours[face]->tetVertex == yRegion->tetVertex);

        if (!tetIndex || !tetVertex)
            continue;

        distV1      = (xRegion->curve[face][vertex1] == yRegion->curve[yFace][yVertex1]);
        distV2      = (xRegion->curve[face][vertex2] == yRegion->curve[yFace][yVertex2]);
        distV3      = yRegion->adjTri[yFace];

        // missing distance
        if (yRegion->curve[yFace][yVertex1] == -1 || yRegion->curve[yFace][yVertex2] == -1)
            uFatalError("find_adj_region", "symplectic_basis");

        if (distV1 && distV2 && distV3)
            return yRegion;
    }

    // We didn't find a cusp yRegion
    //uFatalError("find_cusp_region", "symplectic_basis");
    return NULL;
}

/*
 * Initialise dual curve doubly linked list which
 * stores the oscillating curves on the cusp
 */

void init_oscillating_curves(struct ManifoldBoundary *boundary) {
    int i, j;
    struct DualCurves *path;

    boundary->dual_curve_begin          = NEW_STRUCT( struct DualCurves );
    boundary->dual_curve_end            = NEW_STRUCT( struct DualCurves );
    boundary->dual_curve_begin->next    = boundary->dual_curve_end;
    boundary->dual_curve_begin->prev    = NULL;
    boundary->dual_curve_end->next      = NULL;
    boundary->dual_curve_end->prev      = boundary->dual_curve_begin;

    // which curve
    for (i = 0; i < boundary->numDualCurves; i++) {
        path = NEW_STRUCT(struct DualCurves);
        INSERT_AFTER(path, boundary->dual_curve_begin);

        path->edgeClass = i;

        // which half of the curve
        for (j = 0; j < 2; j++) {
            path->curves[j][START]              = NEW_STRUCT(struct EdgeNode);
            path->curves[j][FINISH]             = NEW_STRUCT(struct EdgeNode);
            path->curves[j][START]->next        = path->curves[j][FINISH];
            path->curves[j][START]->prev        = NULL;
            path->curves[j][FINISH]->next       = NULL;
            path->curves[j][FINISH]->prev       = path->curves[j][START];

            path->endpoints[j][START]           = NEW_STRUCT(struct PathEndPoint);
            path->endpoints[j][FINISH]          = NEW_STRUCT(struct PathEndPoint);
            path->endpoints[j][START]->region   = NULL;
            path->endpoints[j][FINISH]->region  = NULL;
        }
    }
}

void free_oscillating_curves(struct DualCurves *path) {
    int i;
    struct EdgeNode *edge;

    for (i = 0; i < 2; i++) {
        while (path->curves[i][START]->next != path->curves[i][FINISH]) {
            edge = path->curves[i][START]->next;
            REMOVE_NODE(edge);
            my_free(edge);
        }

        my_free(path->curves[i][START]);
        my_free(path->curves[i][FINISH]);
        my_free(path->endpoints[i][START]);
        my_free(path->endpoints[i][FINISH]);
    }

    my_free(path);
}

/*
 * Construct the graph dual to the cusp regions,
 * using region->index to label each vertex, and
 * adding edges using region->adjRegions[].
 */

struct Graph *construct_cusp_region_dual_graph(struct ManifoldBoundary *boundary) {
    int i, face;
    struct CuspRegion *region;

    struct Graph *graph1 = init_graph(boundary->numCuspRegions, FALSE);

    int *visited = NEW_ARRAY(graph1->nVertices, int);

    for (i = 0; i < graph1->nVertices; i++)
        visited[i] = FALSE;

    // Walk around the cusp triangulation inserting edges
    for (region = boundary->cusp_region_begin->next; region != boundary->cusp_region_end; region = region->next) {
        if (visited[region->index])
            continue;

        for (face = 0; face < 4; face++) {
            if (!region->adjTri[face])
                continue;

            // Missing adj region data
            if (region->adjRegions[face] == NULL)
                uFatalError("construct_cusp_region_dual_graph", "symplectic_basis");

            insert_edge(graph1, region->index, region->adjRegions[face]->index, graph1->directed);
            graph1->pRegion[region->index] = region;
        }

        visited[region->index] = 1;
    }

    my_free(visited);
    return graph1;
}

void print_debug_info(struct ManifoldBoundary **cusps, int e0, int flag) {
    int i, j, k, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

    struct CuspTriangle *tri;
    struct EdgeNode *edge ;
    struct CuspRegion *region;
    struct DualCurves *path;
    struct Graph *g;
    struct ManifoldBoundary *boundary = cusps[0];

    if (!debug)
        return;

    if (!flag) {
        // Gluing Info
        printf("Triangle gluing info\n");
        for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
            for (j = 0; j < 4; j++) {
                if (j == tri->tetVertex)
                    continue;

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
        // Dual Curve endpoints
        printf("Dual Curve Endpoint info\n");
        for (path = boundary->dual_curve_begin->next; path != boundary->dual_curve_end; path = path->next) {
            if (path->edgeClass == e0)
                continue;

            for (j = 0; j < 2; j++) {
                edge = path->curves[j][START];
                printf("Tet Vertex: %d (%d %d %d) Tet Vertex: %d (%d %d %d)\n",
                       path->endpoints[j][START]->region->tetVertex,
                       path->endpoints[j][START]->vertex,
                       path->endpoints[j][START]->face,
                       path->curves[j][START]->next->nextFace,
                       path->endpoints[j][FINISH]->region->tetVertex,
                       path->endpoints[j][FINISH]->vertex,
                       path->endpoints[j][FINISH]->face,
                       path->curves[j][FINISH]->prev->prevFace);
            }
        }
    } else if (flag == 2) {
        // Region Info
        printf("Cusp Region info\n");
        for (region = boundary->cusp_region_begin->next; region != boundary->cusp_region_end; region = region->next) {
            v1 = edgesThreeToFour[region->tetVertex][0];
            v2 = edgesThreeToFour[region->tetVertex][1];
            v3 = edgesThreeToFour[region->tetVertex][2];

            printf("Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
                   "(Curves: [[%d %d] [%d %d] [%d %d]]) (Dive: [[%d %d] [%d %d] [%d %d]])\n",
                   region->index, region->tetIndex, region->tetVertex,
                   region->adjTri[v1], region->adjTri[v2], region->adjTri[v3],
                   region->adjRegions[v1] == NULL ? -1 : region->adjRegions[v1]->index,
                   region->adjRegions[v2] == NULL ? -1 : region->adjRegions[v2]->index,
                   region->adjRegions[v3] == NULL ? -1 : region->adjRegions[v3]->index,
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

        for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
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
        for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
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
        for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
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
        for (path = boundary->dual_curve_begin->next; path != boundary->dual_curve_end; path = path->next) {
            if (path->edgeClass == e0)
                continue;

            for (j = 0; j < 2; j++) {
                edge = path->curves[j][START];

                while ((edge = edge->next)->next != NULL) {
                    printf("%d ", edge->y);
                }

                printf("\n");
            }
        }
    } else if (flag == 6) {
        // Inside Edge Info
        printf("Inside edge info\n");
        for (tri = boundary->cusp_triangle_begin->next; tri != boundary->cusp_triangle_end; tri = tri->next) {
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
        g = boundary->dual_graph;

        for (i = 0; i < g->nVertices; i++) {\
            if (g->pRegion[i] == NULL)
                continue;

            printf("Vertex %d (Tet Index: %d, Tet Vertex: %d): ",
                   i,
                   g->pRegion[i]->tetIndex,
                   g->pRegion[i]->tetVertex
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
        for (path = boundary->dual_curve_begin; path != boundary->dual_curve_end; path = path->next) {
            if (e0 == 0)
                continue;

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
 * Returns a pointer to the triangle with a given tetIndex and tetVertex.
 */

struct CuspTriangle *find_cusp_triangle(struct CuspTriangle *cusp_triangle_begin, struct CuspTriangle *cusp_triangle_end,
        struct CuspTriangle *tri, int face) {
    int tetIndex, tetVertex;
    struct CuspTriangle *pTri;

    tetIndex = tri->tet->neighbor[face]->index;
    tetVertex = EVALUATE(tri->tet->gluing[face], tri->tetVertex);

    for (pTri = cusp_triangle_begin->next; pTri != cusp_triangle_end; pTri = pTri->next) {
        if (pTri->tetIndex == tetIndex && pTri->tetVertex == tetVertex)
            return pTri;
    }

    // Didn't find a neighbour
    return NULL;
}

/*
 * Construct oscillating curves on the boundary components
 */

void construct_oscillating_curves(Triangulation *manifold, struct ManifoldBoundary **cusps, int e0) {
    int i, j, *parent;
    bool *processed, *discovered;
    struct ManifoldBoundary *boundary = cusps[0];
    struct DualCurves *path, *pathE0;

    // find e0 curve
    for (path = boundary->dual_curve_begin->next; path != boundary->dual_curve_end; path = path->next) {
        if (path->edgeClass == e0) {
            pathE0 = path;
            break;
        }
    }

    for (i = 0; i < manifold->num_cusps; i++) {
        cusps[i]->dual_graph = construct_cusp_region_dual_graph(cusps[i]);
    }

//    print_debug_info(pTriangle, graph1, pCuspRegion, pDualCurves, 7);
    find_path_endpoints(boundary->dual_graph, NULL,
                        pathE0->endpoints[FIRST][START], e0, FIRST, FALSE);
    find_path_endpoints(boundary->dual_graph, pathE0->endpoints[FIRST][START],
                        pathE0->endpoints[SECOND][START], e0, SECOND, TRUE);

    for (path = boundary->dual_curve_begin->next; path != boundary->dual_curve_end; path = path->next) {
        if (path->edgeClass == e0)
            continue;

        // which half of the curve
        for (j = 0; j < 2; j++) {
            // Find current endpoints
            find_path_endpoints(boundary->dual_graph, pathE0->endpoints[(j + 1) % 2][START],
                                path->endpoints[j][START], e0, j, TRUE);
            find_path_endpoints(boundary->dual_graph, path->endpoints[FIRST][FINISH],
                                path->endpoints[j][FINISH], path->edgeClass, j, j);

            processed = NEW_ARRAY(boundary->dual_graph->nVertices, bool);
            discovered = NEW_ARRAY(boundary->dual_graph->nVertices, bool);
            parent = NEW_ARRAY(boundary->dual_graph->nVertices, int);

            // Find path using bfs
            init_search(boundary->dual_graph, processed, discovered, parent);
            bfs(boundary->dual_graph, path->endpoints[j][START]->regionIndex, processed, discovered, parent);
            find_path(
                    path->endpoints[j][START]->regionIndex,
                    path->endpoints[j][FINISH]->regionIndex, parent, path->curves[j][START]);
            update_path_info(boundary->dual_graph, path, j);

            print_debug_info(cusps, e0, 8);
            print_debug_info(cusps, e0, 5);

            // Reallocate memory
            my_free(processed);
            my_free(discovered);
            my_free(parent);

            // Split the regions along the path
            split_cusp_regions_along_path(boundary, path, j);
            print_debug_info(cusps, e0, 2);

            free_graph(boundary->dual_graph);
            boundary->dual_graph = construct_cusp_region_dual_graph(boundary);
            print_debug_info(cusps, e0, 7);
        }
    }
}

/*
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass. If copy is true, find the path end point
 * corresponding to path1 and store it in path2. Else find a valid path end
 * point and store in path2
 */

void find_path_endpoints(struct Graph *g, struct PathEndPoint *path1, struct PathEndPoint *path2, int edgeClass, int edgeIndex, bool copy) {
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
            pRegion->tri->vertices[vertex].edgeIndex != edgeIndex) {
                continue;
            }

            face1 = (int) remaining_face[pRegion->tetVertex][vertex];
            face2 = (int) remaining_face[vertex][pRegion->tetVertex];

            // are we in the correct region for copy
            if (copy && (pRegion->tetIndex != path1->region->tetIndex ||
                pRegion->tetVertex != path1->vertex ||
                vertex != path1->region->tetVertex ||
                !pRegion->dive[path1->face][vertex]))
                continue;

            if (copy) {
                face = path1->face;
            } else if (pRegion->dive[face1][vertex]) {
                face = face1;
            } else if (pRegion->dive[face2][vertex]) {
                face = face2;
            } else
                continue;

            path2->region = pRegion;
            path2->vertex = vertex;
            path2->face = face;
            path2->regionIndex = i;
            return;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_path_endpoints", "symplectic_basis");
}

/*
 * After finding a path, each node contains the index of the
 * region it lies in. Update path info calculates the face
 * the path crosses to get to the next node and the vertex
 * it cuts off to simplify combinatorial holonomy calculation.
 */

void update_path_info(struct Graph *g, struct DualCurves *path, int curveNum) {
    int i, face, vertex1, vertex2;
    struct EdgeNode *node = path->curves[curveNum][START];

    // path len 0
    if (node->next->next == NULL)
        return;

    node = node->next;
    // path len 1
    if (node->next->next == NULL) {
        for (face = 0; face < 4; face++)
            if (g->pRegion[node->y]->tetVertex != face &&
                path->endpoints[curveNum][START]->vertex != face &&
                path->endpoints[curveNum][FINISH]->vertex != face)
                break;

        node->insideVertex = face;
        return;
    }

    // Set Header node
    for (i = 0; i < 4; i++) {
        if (i == g->pRegion[node->y]->tetVertex || !g->pRegion[node->y]->adjTri[i])
            continue;

        if (g->pRegion[node->y]->adjRegions[i]->index == node->next->y) {
            node->nextFace = i;
            node->prevFace = -1;
        }
    }

    vertex1 = (int) remaining_face[g->pRegion[node->y]->tetVertex][path->endpoints[curveNum][START]->vertex];
    vertex2 = (int) remaining_face[path->endpoints[curveNum][START]->vertex][g->pRegion[node->y]->tetVertex];

    if (node->nextFace == path->endpoints[curveNum][START]->vertex) {
        node->insideVertex = path->endpoints[curveNum][START]->face == vertex1 ? vertex2 : vertex1;
    } else if (node->nextFace == path->endpoints[curveNum][START]->face) {
        node->insideVertex = -1;
    } else {
        node->insideVertex = path->endpoints[curveNum][START]->vertex;
    }

    for (node = node->next; node->next->next != NULL; node = node->next)
        inside_vertex(g->pRegion[node->y], node);

    // Set Tail node
    for (i = 0; i < 4; i++) {
        if (i == g->pRegion[node->y]->tetVertex || !g->pRegion[node->y]->adjTri[i])
            continue;

        if (g->pRegion[node->y]->adjRegions[i]->index == node->prev->y) {
            node->prevFace = i;
            node->nextFace = -1;
        }
    }

    vertex1 = (int) remaining_face[g->pRegion[node->y]->tetVertex][path->endpoints[curveNum][FINISH]->vertex];
    vertex2 = (int) remaining_face[path->endpoints[curveNum][FINISH]->vertex][g->pRegion[node->y]->tetVertex];

    if (node->prevFace == path->endpoints[curveNum][FINISH]->vertex) {
        node->insideVertex = path->endpoints[curveNum][FINISH]->face == vertex1 ? vertex2 : vertex1;
    } else if (node->prevFace == path->endpoints[curveNum][FINISH]->face) {
        node->insideVertex = -1;
    } else {
        node->insideVertex = path->endpoints[curveNum][FINISH]->vertex;
    }
}

/*
 * The oscillating curve splits the region it passes through
 * into two regions. Split each region in two and update
 * attributes
 */

void split_cusp_regions_along_path(struct ManifoldBoundary *boundary, struct DualCurves *path, int curveNum) {
    int face, index = boundary->numCuspRegions;
    struct EdgeNode *node = path->curves[curveNum][START];
    struct CuspRegion *newRegion;
    struct Graph *g = boundary->dual_graph;

    node = path->curves[curveNum][START]->next;

    // empty path
    if (node->next == NULL)
        return ;

    // path of len 1
    if (node->next->next == NULL) {
        newRegion = NEW_STRUCT(struct CuspRegion);
        INSERT_BEFORE(newRegion, boundary->cusp_region_end)
        copy_region(g->pRegion[node->y], newRegion);

        face = node->insideVertex;

        newRegion->index = index;
        newRegion->adjTri[path->endpoints[curveNum][START]->vertex]  = 0;
        newRegion->adjTri[path->endpoints[curveNum][FINISH]->vertex] = 0;
        newRegion->dive[path->endpoints[curveNum][START]->vertex][path->endpoints[curveNum][FINISH]->vertex] =
                (face != path->endpoints[curveNum][FINISH]->face);
        newRegion->dive[path->endpoints[curveNum][FINISH]->vertex][path->endpoints[curveNum][START]->vertex] =
                (face != path->endpoints[curveNum][START]->face);

        g->pRegion[node->y]->adjTri[face] = 0;
        g->pRegion[node->y]->dive[face][path->endpoints[curveNum][START]->vertex]  = (face == path->endpoints[curveNum][START]->face);
        g->pRegion[node->y]->dive[face][path->endpoints[curveNum][FINISH]->vertex] = (face == path->endpoints[curveNum][FINISH]->face);

        update_adj_region_data(boundary->cusp_region_begin, boundary->cusp_region_end);
        boundary->numCuspRegions++;
        return;
    }

    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and newRegion to the left of the curve.
     */
    update_cusp_triangle_endpoints(boundary->cusp_region_begin, boundary->cusp_region_end,
                                   g->pRegion[node->y], path->endpoints[curveNum][START], node);
    newRegion = update_cusp_region(g->pRegion[node->y], node,
                                   path->endpoints[curveNum][START], index, -1);
    INSERT_BEFORE(newRegion, boundary->cusp_region_end);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        update_cusp_triangle(boundary->cusp_region_begin, boundary->cusp_region_end,
                             g->pRegion[node->y], node);
        newRegion = update_cusp_region(g->pRegion[node->y], node,
                                       path->endpoints[curveNum][START], index, 0);
        INSERT_BEFORE(newRegion, boundary->cusp_region_end);
        index++;
    }

    // update last region
    update_cusp_triangle_endpoints(boundary->cusp_region_begin, boundary->cusp_region_end,
                                   g->pRegion[node->y], path->endpoints[curveNum][FINISH], node);
    newRegion = update_cusp_region(g->pRegion[node->y], node,
                                   path->endpoints[curveNum][FINISH], index, 1);
    INSERT_BEFORE(newRegion, boundary->cusp_region_end);
    index++;

    update_adj_region_data(boundary->cusp_region_begin, boundary->cusp_region_end);
    boundary->numCuspRegions = index;
}

/*
 * Set the new and old region data. Draw a picture to see how
 * the attributes change in each case
 *   - flag = -1 : starting node
 *   - flag = 0  : interior node
 *   - flag = 1  : end node
 */

struct CuspRegion *update_cusp_region(struct CuspRegion *region, struct EdgeNode *node,
                                      struct PathEndPoint *endPoint, int index, int flag) {
    int face, vertex1, vertex2;
    struct CuspRegion *newRegion = NEW_STRUCT(struct CuspRegion);

    vertex1 = (int) remaining_face[region->tetVertex][node->insideVertex];
    vertex2 = (int) remaining_face[node->insideVertex][region->tetVertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * newRegion becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, newRegion);
    newRegion->index = index;

    if (flag == 0) {
        // Update new region
        newRegion->curve[vertex1][node->insideVertex]++;
        newRegion->curve[vertex2][node->insideVertex]++;
        newRegion->dive[node->insideVertex][vertex1]    = region->dive[node->insideVertex][vertex1];
        newRegion->dive[vertex2][vertex1]               = region->dive[vertex2][vertex1];
        newRegion->dive[node->insideVertex][vertex2]    = region->dive[node->insideVertex][vertex2];
        newRegion->dive[vertex1][vertex2]               = region->dive[vertex1][vertex2];

        // Update region
        region->curve[vertex2][vertex1]++;
        region->curve[vertex1][vertex2]++;
        region->dive[vertex2][vertex1]                  = 0;
        region->dive[node->insideVertex][vertex1]       = 0;
        region->dive[vertex1][vertex2]                  = 0;
        region->dive[node->insideVertex][vertex2]       = 0;
        region->adjTri[node->insideVertex]              = 0;

        return newRegion;
    }

    vertex1 = (int) remaining_face[region->tetVertex][endPoint->vertex];
    vertex2 = (int) remaining_face[endPoint->vertex][region->tetVertex];

    if (flag == -1) {
        face = node->nextFace;
    } else {
        face = node->prevFace;
    }

    if (face == endPoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        newRegion->curve[endPoint->vertex][vertex2]++;

        newRegion->dive[endPoint->vertex][vertex1]      = region->dive[endPoint->vertex][vertex1];
        newRegion->dive[vertex2][vertex1]               = region->dive[vertex2][vertex1];
        newRegion->dive[vertex1][endPoint->vertex]      = (endPoint->face == vertex1);
        newRegion->adjTri[vertex1]                      = 0;

        region->curve[endPoint->vertex][vertex1]++;
        region->dive[endPoint->vertex][vertex1]         = 0;
        region->dive[vertex2][endPoint->vertex]         = (endPoint->face == vertex2);
        region->dive[vertex2][vertex1]                  = 0;
        region->adjTri[vertex2]                         = 0;
    } else if (face == endPoint->face) {
        // curve passes through the face that carries it
        newRegion->curve[endPoint->face][endPoint->face == vertex1 ? vertex2 : vertex1]++;

        newRegion->dive[endPoint->face][endPoint->vertex]                   = 1;
        newRegion->adjTri[endPoint->vertex]                                 = 0;
        newRegion->adjTri[endPoint->face == vertex1 ? vertex2 : vertex1]    = 0;

        region->curve[endPoint->face][endPoint->vertex]++;
    } else {
        // Curve goes around the vertex
        newRegion->curve[face][endPoint->face]++;

        newRegion->dive[face][endPoint->vertex]                 = 1;
        newRegion->dive[endPoint->face][endPoint->vertex]       = 1;
        newRegion->adjTri[endPoint->face]                       = 0;
        newRegion->adjTri[endPoint->vertex]                     = 0;

        region->curve[face][endPoint->vertex]++;
        region->dive[face][endPoint->vertex]                    = 0;
    }

    return newRegion;
}

/*
 * region1 splits into region1 and region2, set them up
 * to be split
 */

void copy_region(struct CuspRegion *region1, struct CuspRegion *region2) {
    int i, j;

    if (region1 == NULL || region2 == NULL || region1->tri == NULL)
        uFatalError("copy_region", "symplectic_basis");

    region2->tri            = region1->tri;
    region2->tetIndex       = region1->tetIndex;
    region2->tetVertex      = region1->tetVertex;

    for (i = 0; i < 4; i++) {
        region2->adjTri[i]          = region1->adjTri[i];
        region2->adjRegions[i]      = NULL;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]    = region1->curve[i][j];
            region2->dive[i][j]     = 0;
        }
    }
}

/*
 * After splitting each region the path travels through,
 * the attributes for other regions in the same cusp
 * triangle is now out of date. Update cusp triangles
 * for nodes in the interior of the path.
 */

void update_cusp_triangle(struct CuspRegion *cusp_region_start, struct CuspRegion *cusp_region_end,
        struct CuspRegion *region, struct EdgeNode *node) {
    int face1, face2;
    struct CuspRegion *pRegion;

    face1 = (int) remaining_face[region->tetVertex][node->insideVertex];
    face2 = (int) remaining_face[node->insideVertex][region->tetVertex];

    for (pRegion = cusp_region_start->next; pRegion != cusp_region_end; pRegion = pRegion->next) {
        // is the region initialised?
        if (pRegion == NULL || pRegion->tetIndex == -1)
            continue;

        // which triangle are we in?
        if (pRegion->tetIndex != region->tetIndex || pRegion->tetVertex != region->tetVertex)
            continue;

        if (pRegion->curve[face1][node->insideVertex] > region->curve[face1][node->insideVertex]) {
            pRegion->curve[face1][node->insideVertex]++;
            pRegion->dive[face1][node->insideVertex] = 0;
        }
        else if (pRegion->curve[face1][node->insideVertex] < region->curve[face1][node->insideVertex]) {
            pRegion->curve[face1][face2]++;
        }

        if (pRegion->curve[face2][node->insideVertex] > region->curve[face2][node->insideVertex]) {
            pRegion->curve[face2][node->insideVertex]++;
            pRegion->dive[face2][node->insideVertex] = 0;
        }
        else if (pRegion->curve[face2][node->insideVertex] < region->curve[face2][node->insideVertex]) {
            pRegion->curve[face2][face1]++;
        }
    }
}

/*
 * After splitting each region the path travels through,
 * the attributes for other regions in the same cusp
 * triangle is now out of date. Update cusp triangles
 * for nodes at the end of the path.
 */

void update_cusp_triangle_endpoints(struct CuspRegion *cusp_region_start, struct CuspRegion *cusp_region_end, struct CuspRegion *region,
        struct PathEndPoint *endPoint, struct EdgeNode *node) {
    int face, face1, face2;
    struct CuspRegion *pRegion;

    face1 = (int) remaining_face[region->tetVertex][endPoint->vertex];
    face2 = (int) remaining_face[endPoint->vertex][region->tetVertex];

    if (node->prevFace == -1) {
        face = node->nextFace;
    } else {
        face = node->prevFace;
    }

    for (pRegion = cusp_region_start->next; pRegion != cusp_region_end; pRegion = pRegion->next) {
        if (pRegion == NULL || pRegion->tetIndex == -1)
            continue;

        // which triangle are we in?
        if (pRegion->tetIndex != region->tetIndex || pRegion->tetVertex != region->tetVertex)
            continue;

        if (face == endPoint->vertex) {
            // curve passes through the face opposite the vertex it dives through
            if (pRegion->curve[endPoint->vertex][face1] > region->curve[endPoint->vertex][face1]) {
                pRegion->curve[face][face1]++;
                pRegion->dive[face][face1] = 0;
            } else if (pRegion->curve[endPoint->vertex][face1] < region->curve[endPoint->vertex][face1]) {
                pRegion->curve[face][face2]++;
                pRegion->dive[face][face2] = 0;
            }
            continue;
        }

        // Curve goes around the vertex or passes through the face that carries it
        if (pRegion->curve[face][endPoint->vertex] > region->curve[face][endPoint->vertex]) {
            pRegion->curve[face][endPoint->vertex]++;
            pRegion->dive[face][endPoint->vertex] = 0;
        } else if (pRegion->curve[face][endPoint->vertex] < region->curve[face][endPoint->vertex]) {
            pRegion->curve[face][face == face1 ? face2 : face1]++;
            pRegion->dive[face][face == face1 ? face2 : face1] = 0;
        }
    }

}

// ----------------------------------

/*
 * Construct the symplectic equations from the dual curves
 */

void calculate_holonomy(Triangulation *manifold, struct ManifoldBoundary **cusps, int e0, int **symp_eqns) {
    int i;
    struct DualCurves *path;

    // which cusp are we looking at
    for (i = 0; i < manifold->num_cusps; i++) {

        // which curve on this cusp
        for (path = cusps[i]->dual_curve_begin->next; path != cusps[i]->dual_curve_end; path = path->next) {
            if (path->edgeClass == e0)
                continue;

            find_path_holonomy(cusps[i]->dual_graph, path, FIRST, symp_eqns[path->edgeClass]);
            find_path_holonomy(cusps[i]->dual_graph, path, SECOND, symp_eqns[path->edgeClass]);
        }
    }
}

/*
 * Calculate holonomies along a path and update the row.
 */

void find_path_holonomy(struct Graph *g, struct DualCurves *curve, int curveNum, int *row) {
    int i, index, dirFace, face;
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
        row[index] = row[index] + orientate_vertex(g->pRegion[node->y]->tri, i, pathStartPoint->vertex);
        return;
    }

    while((node = node->next)->next != NULL) {
        // Path end points
        if (node->prev->prev == NULL || node->next->next == NULL) {
            // End point index
            if (node->prev->prev == NULL) {
                face = node->nextFace;
                dirFace = 1;
            } else {
                face = node->prevFace;
                dirFace = -1;
            }

            if (node->y == pathStartPoint->regionIndex)
                endPoint = pathStartPoint;
            else
                endPoint = pathEndPoint;

            if (endPoint->vertex == face) {
                /*
                 * Curve passes through the face opposite the vertex it
                 * dives through. Picks up holonomy for the vertex between
                 * the face that carries the curve and the face the curve
                 * crosses
                 */

                index = 3 * pathStartPoint->region->tetIndex + edge3_between_faces
                [remaining_face[endPoint->region->tetVertex][node->insideVertex]]
                [remaining_face[node->insideVertex][endPoint->region->tetVertex]];

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
                index = 3 * endPoint->region->tetIndex + edge3_between_faces[endPoint->face][face];
            }

            row[index] = row[index] + dirFace * orientate_vertex(g->pRegion[node->y]->tri, node->insideVertex, face);

            continue;
        }

        index = 3 * g->pRegion[node->y]->tetIndex + edge3_between_faces[remaining_face[g->pRegion[node->y]->tetVertex][node->insideVertex]]
                [remaining_face[node->insideVertex][g->pRegion[node->y]->tetVertex]];
        row[index] = row[index] + orientate_vertex(g->pRegion[node->y]->tri, node->insideVertex, node->nextFace);
    }
}

/*
 * node lies in region midNode, find the vertex which the subpath
 * node->prev->y --> node->y --> node->next->y cuts off of the
 * cusp triangle midNode->tri.
 */

void inside_vertex(struct CuspRegion *midNode, struct EdgeNode *node) {
    int i, vertex1, vertex2;

    for (i = 0; i < 4; i++) {
        if (i == midNode->tetVertex)
            continue;

        vertex1 = (int) remaining_face[midNode->tetVertex][i];
        vertex2 = (int) remaining_face[i][midNode->tetVertex];

        if (!midNode->adjTri[vertex1] || !midNode->adjTri[vertex2])
            continue;

        if (midNode->adjRegions[vertex1]->index == node->next->y
        && midNode->adjRegions[vertex2]->index == node->prev->y) {
            node->insideVertex = i;
            node->nextFace = vertex1;
            node->prevFace = vertex2;
            return;
        } else if (midNode->adjRegions[vertex2]->index == node->next->y && midNode->adjRegions[vertex1]->index == node->prev->y) {
            node->insideVertex = i;
            node->nextFace = vertex2;
            node->prevFace = vertex1;
            return;
        }
    }

    // where does the next node go?
    uFatalError("inside_vertex", "symplectic_basis");
}

// ---------------------------------------------------------------

// Breadth First Search from Skiena Algorithm Design Manual

/*
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
        processed[v] = true;
        p = g->edge_list_begin[v];

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (!discovered[y]) {
                enqueue(&q, y);
                discovered[y] = true;
                parent[y] = v;
            }
        }
    }

    free_queue(&q);
}

/*
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
    if ( q->size - 1 == q->len ) {
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

// -------------------------------------------------------

// End Multi Graph

struct EndMultiGraph *end_multi_graph(Triangulation *manifold) {
    int start = 0, i;

    struct EndMultiGraph *multiGraph = NEW_STRUCT( struct EndMultiGraph );
    struct Graph *g = init_graph(manifold->num_cusps, FALSE);
    cusp_graph(manifold, g);


//    printf("--------------------------\n");
//    printf("Cusp Graph Info\n");
//
//    for (i = 0; i < graph1->nVertices; i++) {
//        printf("(Cusp Index: %d) ", i);
//        edge = graph1->edge_list_begin[i];
//        while ((edge = edge->next)->next != NULL) {
//            printf("%d ", edge->y);
//        }
//        printf("\n");
//    }
//
//    printf("--------------------------\n");

    multiGraph->multi_graph = init_graph(g->nVertices, g->directed);
    int *parent = NEW_ARRAY(multiGraph->multi_graph->nVertices, int);
    spanning_tree(g, multiGraph->multi_graph, start, parent);

//    printf("Cusp graph spanning tree\n");
//
//    for (i = 0; i < graph1->nVertices; i++) {
//        printf("(Cusp Index: %d) ", i);
//        edge = graph1->edge_list_begin[i];
//        while ((edge = edge->next)->next != NULL) {
//            printf("%d ", edge->y);
//        }
//        printf("\n");
//    }
//
//    printf("--------------------------\n");

    add_odd_cycle_edge(multiGraph->multi_graph, parent);

//    printf("Cusp graph with cycle\n");
//
//    for (i = 0; i < graph1->nVertices; i++) {
//        printf("(Cusp Index: %d) ", i);
//        edge = graph1->edge_list_begin[i];
//        while ((edge = edge->next)->next != NULL) {
//            printf("%d ", edge->y);
//        }
//        printf("\n");
//    }
//
//    printf("--------------------------\n");

    multiGraph->double_graph = init_graph(multiGraph->multi_graph->nVertices, multiGraph->multi_graph->directed);
    construct_double_graph(multiGraph->multi_graph, multiGraph->double_graph);
}

void cusp_graph(Triangulation *manifold, struct Graph *g) {
    int vertex1, vertex2;
    Tetrahedron *tet;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (vertex1 = 0; vertex1 < 4; vertex1++) {
            // which vertex of the cusp triangle at vertex1
            for (vertex2 = 0; vertex2 < 4; vertex2++) {
                if (vertex1 == vertex2)
                    continue;

                insert_edge(g, tet->cusp[vertex1]->index, tet->cusp[vertex2]->index, g->directed);
            }
        }
    }
}

/*
 * Find a spanning tree of graph1
 */

void spanning_tree(struct Graph *graph1, struct Graph *graph2, int start, int *parent) {
    int i;
    bool *processed = NEW_ARRAY(graph1->nVertices, bool);
    bool *discovered = NEW_ARRAY(graph1->nVertices, bool);

    // Find path using bfs
    init_search(graph1, processed, discovered, parent);
    bfs(graph1, start, processed, discovered, parent);

    for (i = 0; i < graph1->nVertices; i++) {
        if (parent[i] == -1)
            continue;

        insert_edge(graph2, i, parent[i], graph2->directed);
    }
}

/*
 * Add an edge to the tree *g which creates a cycle of odd length
 */

void add_odd_cycle_edge(struct Graph *g, int *parent) {

}

/*
 * Find the length of a path between start and end
 */

int find_path_len(int start, int end, int *parents, int pathLen) {
    if ((start == end) || (end == -1)) {
        return pathLen++;
    } else {
        return find_path_len(start, parents[end], parents, pathLen++);
    }
}

/*
 * Create a new graph g2 which consists of the
 * vertices of g1, and an edge uv from u to v iff
 * u and v share a neighbour
 */

void construct_double_graph(struct Graph *g1, struct Graph *g2) {
    int v1, v2;
    struct EdgeNode *node;

    for (v1 = 0; v1 < g1->nVertices; v1++) {
        for (v2 = 0; v2 < g1->nVertices; v2++) {
            if (v1 == v2)
                continue;

            for (node = g1->edge_list_begin[v1]->next; node != g1->edge_list_end[v1]; node = node->next) {
                if (edge_exists(g1, node->y, v2)) {
                    insert_edge(g2, v1, v2, g2->directed);
                    g2->edge_list_begin[v1]->next->intermediate = node->y;
                    g2->edge_list_begin[v2]->next->intermediate = node->y;
                }
            }
        }
    }
}

/*
 *
 */

void find_even_len_path(struct EndMultiGraph *multiGraph, int start, int end, struct EdgeNode *node) {
    struct EdgeNode *new_node = NEW_STRUCT(struct EdgeNode);

    INSERT_AFTER(new_node, node);

    if ((start == end) || (end == -1)) {
        new_node->y = start;
    } else {
        find_path(start, parents[end], parents, node);
        new_node->y = end;
    }
}
