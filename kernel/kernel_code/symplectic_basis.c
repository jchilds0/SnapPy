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

#define FIRST               0
#define MIDDLE              1
#define LAST                2
#define START               0
#define FINISH              1

/*
 * Queue
 */

typedef struct Queue {
    int                     front;      // First element of queue
    int                     rear;       // Last element of queue
    int                     len;        // num of elements
    int                     size;       // array size
    int                     *array;
} Queue ;

/**
 * Graph
 */

typedef struct EdgeNode {
    int                     y;                      /** cusp region index */
    int                     edgeClass;              /** edgeClass of the edge for edges in the end multi graph */
    int                     nextFace;               /** face the path crosses to the next node */
    int                     prevFace;               /** face the path crosses to the prev node */
    int                     insideVertex;           /** inside vertex of the path */
    struct EdgeNode         *next;                  /** next node in doubly linked list */
    struct EdgeNode         *prev;                  /** prev node in doubly linked list */
} EdgeNode;

typedef struct Graph {
    EdgeNode                *edge_list_begin;       /** header node of doubly linked list */
    EdgeNode                *edge_list_end;         /** tail node of doubly linked list */
    struct CuspRegion       **pRegion;              /** list of regions in the graph */
    int                     *degree;                /** degree of each vertex */
    int                     *color;                 /** color a tree bipartite */
    int                     nVertices;              /** number of vertices in the graph */
    int                     directed;               /** is the graph directed */
} Graph;

typedef struct CuspEndPoint {
    int                     cuspIndex;
    int                     edgeClass1;
    int                     edgeClass2;
    int                     pos;
} CuspEndPoint;

typedef struct EndMultiGraph {
    int                     e0;                     /** edge connecting vertices of the same color */
    Graph                   *multi_graph;           /** tree with extra edge of cusps */
} EndMultiGraph;

/**
 * Dual Curves
 *
 * Each oscillating curve is made up of two components, accessed using macros
 * FIRST and SECOND. For each component we have the starting endpoint,
 * *endpoints[curveNum][START] and finish endpoint *endpoint[curveNum][FINISH].
 * The path of the curve is stored as a double linked list with header and tail nodes,
 * the header is curves[curveNum][START] and the tail is curve[curveNum][FINISH].
 */

struct extra {
    int                     curve[4][4];            /** oscillating curve holonomy for a cusp triangle */
};

typedef struct PathEndPoint {
    int                     face;                   /** face containg the short rectangle carrying the curve */
    int                     vertex;                 /** vertex we dive through the manifold along */
    int                     regionIndex;            /** index of the region the endpoint lies in */
    struct CuspRegion       *region;                /** pointer to the region the endpoint lies in */
} PathEndPoint;

typedef struct DualCurves {
    int                     edgeClass[2];
    EdgeNode                curves_begin;           /** header node of doubbly linked list */
    EdgeNode                curves_end;             /** tailer node of doubbly linked list */
    PathEndPoint            endpoints[2];           /** path end points */
    struct DualCurves       *next;                  /** next dual curve in doubly linked list */
    struct DualCurves       *prev;                  /** prev dual curve in doubly linked list */
} DualCurves;

typedef struct OscillatingCurves {
    int                     numCurves;
    int                     *edgeClass;
    DualCurves              *dual_curve_begin;      /** array of doubly linked lists of dual curves */
    DualCurves              *dual_curve_end;        /** array of doubly linkek lists of dual curves */
} OscillatingCurves;

/**
 * Cusp Triangulation
 *
 * CuspTriangle stores information about a triangle in the
 * cusp triangulation. The homology curves bound a fundamental
 * domain, and cusp regions store the information for
 * intersection of this domain with each cusp triangle. When
 * we add oscillating curves, these regions are divided further.
*/

typedef struct CuspVertex {
    int                     edgeClass;
    int                     edgeIndex;
    EdgeClass               *edge;
    int                     v1;
    int                     v2;
} CuspVertex;

typedef struct CuspTriangle {
    Tetrahedron             *tet;                   /** tetrahedron the triangle comes from */
    int                     tetIndex;               /** tet->index */
    int                     tetVertex;              /** vertex the triangle comes from */
    int                     numCurves;              /** number of homology curves on the triangle */
    CuspVertex              vertices[4];            /** information about each vertex */
    struct CuspTriangle     *neighbours[4];         /** triangle neighbouring a face */
    struct CuspTriangle     *next;                  /** next cusp triangle on doubly linked list */
    struct CuspTriangle     *prev;                  /** prev cusp triangle on doubly linkled list */
} CuspTriangle;

typedef struct CuspRegion {
    CuspTriangle            *tri;                   /** cusp triangle the region lies on */
    int                     tetIndex;               /** tri->tetIndex */
    int                     tetVertex;              /** tri->tetVertex */
    int                     index;                  /** index of the cusp region */
    int                     curve[4][4];            /** looking at face, number of curves between the region and vertex */
    int                     adjTri[4];              /** does the region meet this edge of the cusp triangle */
    struct CuspRegion       *adjRegions[4];         /** index of the adjacent regions */
    int                     dive[4][4];             /** can we dive along the face into this vertex */
    struct CuspRegion       *next;                  /** next cusp region in doubly linked list */
    struct CuspRegion       *prev;                  /** prev cusp region in doubly linked list */
} CuspRegion;

typedef struct ManifoldBoundary {
    int                     intersectTetIndex;      /** index of the intersection triangle */
    int                     intersectTetVertex;     /** vertex of the intersection triangle */
    int                     numEdgeClasses;         /** number of edge classes in the boundary */
    int                     numCuspTriangles;       /** number of cusp triangle in the boundary */
    int                     numCuspRegions;         /** number of cusp regions in the boundary */
    int                     numDualCurves;          /** number of dual curves in the boundary */
    Triangulation           *manifold;              /** manifold */
    Cusp                    *cusp;                  /** which cusp is the boundary in */
    Graph                   *dual_graph;            /** dual graph of the cusp region */
    CuspTriangle            cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    CuspTriangle            cusp_triangle_end;      /** tail node of doubly linked list of cusp triangles */
    CuspRegion              cusp_region_begin;      /** header node of doubly linked list of cusp regions */
    CuspRegion              cusp_region_end;        /** tail node of doubly linked list of cusp regions */
} ManifoldBoundary;

// Graph
Graph                   *init_graph(int, bool);
void                    free_graph(Graph *);
int                     insert_edge(Graph *, int, int, bool);
void                    delete_edge(Graph *, int, int, bool);
int                     edge_exists(Graph *, int, int);

// Symplectic Basis
int                     *gluing_equations_for_edge_class(Triangulation *, int);
int                     **get_symplectic_equations(Triangulation *, int *, int *, int);

ManifoldBoundary        *init_boundary(Triangulation *, Cusp *);
void                    free_boundary(ManifoldBoundary **, int);
void                    init_cusp_triangulation(Triangulation *, ManifoldBoundary *);
int                     flow(CuspTriangle *, int);
void                    label_triangulation_edges(Triangulation *);
CuspTriangle            *find_cusp_triangle(CuspTriangle *, CuspTriangle *, CuspTriangle *, int);
void                    label_cusp_vertex_indices(CuspTriangle *, CuspTriangle *, int);
void                    walk_around_cusp_vertex(CuspTriangle *, int, int);
void                    init_cusp_region(ManifoldBoundary *);
int                     init_intersect_cusp_region(ManifoldBoundary *, CuspTriangle *, int);
int                     init_normal_cusp_region(ManifoldBoundary *, CuspTriangle *, int);
void                    set_cusp_region_data(CuspRegion *, CuspTriangle *, int [4], int [4], int);
void                    update_adj_region_data(CuspRegion *, CuspRegion *);
CuspRegion              *find_adj_region(CuspRegion *, CuspRegion *, CuspRegion *, int);
OscillatingCurves       *init_oscillating_curves(Triangulation *, int *);
void                    free_oscillating_curves(OscillatingCurves *);
void                    find_intersection_triangle(Triangulation *, ManifoldBoundary *);

/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    do_oscillating_curves(ManifoldBoundary **, OscillatingCurves *, EndMultiGraph *);
void                    do_one_dual_curve(ManifoldBoundary **, DualCurves *, DualCurves *, EndMultiGraph *, int);
void                    do_one_cusp(ManifoldBoundary *, DualCurves *, int);
Graph *                 construct_cusp_region_dual_graph(ManifoldBoundary *);
void                    print_debug_info(Triangulation *, ManifoldBoundary **, OscillatingCurves *, int);
DualCurves              *init_dual_curve(DualCurves *, int, int);
void                    find_path_endpoints(Graph *, DualCurves *, DualCurves *, int, int);
void                    find_single_endpoints(Graph *, PathEndPoint *, PathEndPoint *, int, int, bool);
void                    update_path_info(Graph *g, DualCurves *, int);
void                    update_path_endpoint_info(CuspRegion *, EdgeNode *, PathEndPoint *, int, int);
void                    split_cusp_regions_along_path(ManifoldBoundary *, DualCurves *);
CuspRegion              *update_cusp_region(CuspRegion *region, EdgeNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle(CuspRegion *, CuspRegion *, CuspRegion *, EdgeNode *);
void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, EdgeNode *, int);
void                    copy_region(CuspRegion *, CuspRegion *);
void                    calculate_holonomy(Triangulation *, int **, int);
void                    inside_vertex(CuspRegion *, EdgeNode *);

/**
 * Queue Data Structure
 */

Queue                   *init_queue(int);
Queue                   *enqueue(Queue *, int);
int                     dequeue(Queue *);
Queue                   *resize_queue(Queue *);
int                     empty_queue(Queue *);
void                    free_queue(Queue *);

/**
 * Graph for Breadth First Search
 */

void                    init_search(Graph *, bool *, bool *, int *);
void                    bfs(Graph *, int, bool *, bool *, int *);
void                    find_path(int, int, int *, EdgeNode *);


/**
 * Spanning Tree for End Multi Graph
 */

EndMultiGraph           *init_end_multi_graph(Triangulation *, int *);
void                    free_end_multi_graph(EndMultiGraph *);
int                     insert_edge_end_multi_graph(Graph *, int, int, int, bool);
void                    spanning_tree(Graph *, Graph *, int, int *);
void                    cusp_graph(Triangulation *, Graph *);
void                    color_graph(Graph *);
int                     *find_tree_edges(Graph *, int);
int                     find_same_color_edge(Graph *, Graph *, int *);
int                     find_path_len(int, int, int *, int);
CuspEndPoint            *find_multi_graph_path(Graph *, Triangulation *, int, int, int *);
CuspEndPoint            *find_cusp_endpoint_edge_classes(Graph *, EdgeNode *, EdgeNode *, int, int *);
void                    find_edge_ends(Graph *, Triangulation *, int, int *, int *);
void                    print_graph(Graph *, int);

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
    int *edgeClasses = NEW_ARRAY(manifold->num_tetrahedra, int);
    peripheral_curves(manifold);

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, edgeClasses, &symp_num_rows, *num_cols);

    // Construct return array
    *num_rows = 2 * (manifold->num_tetrahedra - manifold->num_cusps);
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edgeClasses[i]) {
            my_free(symp_eqns[i]);
            continue;
        }

        eqns[2 * j] = gluing_equations_for_edge_class(manifold, i);
        eqns[2 * j + 1] = symp_eqns[i];
        j++;
    }

    my_free(edgeClasses);
    my_free(symp_eqns);

    *num_cols = 3 * manifold->num_tetrahedra;
    return eqns;
}

/*
 * Copy of get_gluings_equations.c get_gluing_equations which finds 
 * the edge gluings equations for a given edge index. Used instead 
 * of get_gluing_equations to ensure we have the correct edge index 
 * and reduce the chance of memory leak since we dont need all of 
 * rows from gluing equations 
 */

int *gluing_equations_for_edge_class(Triangulation *manifold, int edgeClass) {
    int *eqns, i, T;
    EdgeClass *edge;
    PositionedTet ptet0, ptet;

    T = manifold->num_tetrahedra;
    eqns = NEW_ARRAY(3 * T, int);

    for (i = 0; i < 3 * T; i++)
        eqns[i] = 0;

    /*
     *  Build edge equations.
     */

    for (edge = manifold->edge_list_begin.next; edge != &manifold->edge_list_end; edge = edge->next) {
        if (edge->index == edgeClass)
            break;
    }

    set_left_edge(edge, &ptet0);
    ptet = ptet0;
    do {
        eqns[3 * ptet.tet->index + edge3_between_faces[ptet.near_face][ptet.left_face]]++;
        veer_left(&ptet);
    } while (same_positioned_tet(&ptet, &ptet0) == FALSE);

    return eqns;
}

/*
 * Setup graph and cusp triangulation, and run construct dual curves.
 */

static int debug = FALSE;

int **get_symplectic_equations(Triangulation *manifold, int *edgeClasses, int *num_rows, int numCols) {
    int i, j, k;
    label_triangulation_edges(manifold);

    ManifoldBoundary **cusps     = NEW_ARRAY(manifold->num_cusps, ManifoldBoundary *);
    Cusp *cusp                   = manifold->cusp_list_begin.next;
    EndMultiGraph *multiGraph    = init_end_multi_graph(manifold, edgeClasses);
    OscillatingCurves *oscCurv   = init_oscillating_curves(manifold, edgeClasses);

    for (i = 0; i < manifold->num_cusps; i++) {
        cusps[i] = init_boundary(manifold, cusp);
        cusp = cusp->next;
    }

    Tetrahedron *tet;
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        if (tet->extra != NULL)
            uFatalError("get_symplectic_equations", "symplectic_basis");

        tet->extra = NEW_ARRAY(manifold->num_tetrahedra, Extra);

        for (i = 0; i < manifold->num_tetrahedra; i++)
            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    tet->extra[i].curve[j][k] = 0;
    }

    print_debug_info(manifold, cusps, oscCurv, 0);       // Gluing
    print_debug_info(manifold, cusps, oscCurv, 3);       // Homology
    print_debug_info(manifold, cusps, oscCurv, 4);       // Edge classes
    print_debug_info(manifold, cusps, oscCurv, 6);       // Inside Edge
    print_debug_info(manifold, cusps, oscCurv, 2);       // Regions

    // Allocate Symplectic Equations Array
    *num_rows = oscCurv->numCurves;
    int **symp_eqns = NEW_ARRAY(manifold->num_tetrahedra, int *);

    for (i = 0; i < manifold->num_tetrahedra; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    do_oscillating_curves(cusps, oscCurv, multiGraph);
    calculate_holonomy(manifold, symp_eqns, manifold->num_tetrahedra);

    free_end_multi_graph(multiGraph);
    free_oscillating_curves(oscCurv);
    free_boundary(cusps, manifold->num_cusps);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        my_free(tet->extra);
        tet->extra = NULL;
    }

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

ManifoldBoundary *init_boundary(Triangulation *manifold, Cusp *cusp) {
    ManifoldBoundary *boundary = NEW_STRUCT(ManifoldBoundary);

    // Invalid cusp topology
    if (cusp->topology == Klein_cusp)
        uFatalError("init_boundary", "symplectic_basis");

    boundary->manifold = manifold;
    boundary->cusp = cusp;
    boundary->numEdgeClasses = manifold->num_tetrahedra;
    boundary->numCuspTriangles = 0;
    boundary->numCuspRegions = 0;
    boundary->numDualCurves = manifold->num_tetrahedra;

    find_intersection_triangle(manifold, boundary);
    init_cusp_triangulation(manifold, boundary);
    init_cusp_region(boundary);

    boundary->dual_graph = construct_cusp_region_dual_graph(boundary);

    return boundary;
}

void free_boundary(ManifoldBoundary **cusps, int numCusps) {
    int i;
    CuspTriangle *tri;
    CuspRegion *region;
    ManifoldBoundary *boundary;

    for (i = 0; i < numCusps; i++) {
        boundary = cusps[i];
        // free graph
        free_graph(boundary->dual_graph);

        // free cusp regions
        while (boundary->cusp_region_begin.next != &boundary->cusp_region_end) {
            region = boundary->cusp_region_begin.next;
            REMOVE_NODE(region);
            my_free(region);
        }

        // free cusp triangle
        while (boundary->cusp_triangle_begin.next != &boundary->cusp_triangle_end) {
            tri = boundary->cusp_triangle_begin.next;
            REMOVE_NODE(tri);
            my_free(tri);
        }

        my_free(boundary);
    }

    my_free(cusps);
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

void find_intersection_triangle(Triangulation *manifold, ManifoldBoundary *boundary) {
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

void init_cusp_triangulation(Triangulation *manifold, ManifoldBoundary *boundary) {
    int vertex, face, index = 0;
    Tetrahedron *tet;
    CuspTriangle *tri;

    // Allocate Cusp Triangulation Header and Tail Null nodes
    boundary->cusp_triangle_begin.next = &boundary->cusp_triangle_end;
    boundary->cusp_triangle_begin.prev = NULL;
    boundary->cusp_triangle_end.next = NULL;
    boundary->cusp_triangle_end.prev = &boundary->cusp_triangle_begin;

    // which tetrahedron are we on
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // while vertex are we on
        for (vertex = 0; vertex < 4; vertex++) {
            // is this vertex on the right cusp
            if (tet->cusp[vertex] != boundary->cusp) {
                continue;
            }

            tri = NEW_STRUCT( CuspTriangle );
            INSERT_BEFORE(tri, &boundary->cusp_triangle_end);
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
    for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
        // which vertex
        for (face = 0; face < 4; face++) {
            if (face == tri->tetVertex)
                continue;

            tri->neighbours[face] = find_cusp_triangle(&boundary->cusp_triangle_begin, &boundary->cusp_triangle_end, tri, face);
        }
    }

    label_cusp_vertex_indices(&boundary->cusp_triangle_begin, &boundary->cusp_triangle_end, boundary->numEdgeClasses);
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

void label_cusp_vertex_indices(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end, int numEdgeClasses) {
    int i, vertex;
    CuspTriangle *tri;

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

void walk_around_cusp_vertex(CuspTriangle *tri, int cuspVertex, int index) {
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

void init_cusp_region(ManifoldBoundary *boundary) {
    int index = 0;
    CuspTriangle *tri;

    // Header and tailer nodes.
    boundary->cusp_region_begin.next    = &boundary->cusp_region_end;
    boundary->cusp_region_begin.prev    = NULL;
    boundary->cusp_region_end.next      = NULL;
    boundary->cusp_region_end.prev      = &boundary->cusp_region_begin;

    for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
        // Intersection vertex doesn't have a center
        if (tri->tetIndex == boundary->intersectTetIndex && tri->tetVertex == boundary->intersectTetVertex) {
            index = init_intersect_cusp_region(boundary, tri, index);
            continue;
        }

        index = init_normal_cusp_region(boundary, tri, index);
    }

    update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
    boundary->numCuspRegions = index;
}

int init_intersect_cusp_region(ManifoldBoundary *boundary, CuspTriangle *tri, int index) {
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

                set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
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

            set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
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

        set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
        index++;

        // Find vertex with non-zero flow
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

        set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
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

            set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }
    }

    return index;
}

int init_normal_cusp_region(ManifoldBoundary *boundary, CuspTriangle *tri, int index) {
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

            set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
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

    set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which
 * allocates the attributes of the cusp region
 */

void set_cusp_region_data(CuspRegion *cusp_region_end, CuspTriangle *tri, int dist[4], int adjTri[4], int index) {
    int i, j, v1, v2, v3;
    CuspRegion *pRegion = NEW_STRUCT( CuspRegion );
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

void update_adj_region_data(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end) {
    int j;
    CuspRegion *pRegion;

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

CuspRegion *find_adj_region(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end,
        CuspRegion *xRegion, int face) {
    int vertex1, vertex2, yVertex1, yVertex2, yFace, distV1, distV2, distV3, tetIndex, tetVertex;
    CuspTriangle *tri = xRegion->tri;
    CuspRegion *yRegion;

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

OscillatingCurves *init_oscillating_curves(Triangulation *manifold, int *edgeClasses) {
    int i, j;
    OscillatingCurves *oscCurv = NEW_STRUCT( OscillatingCurves );

    oscCurv->numCurves = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++)
        if (edgeClasses[i])
            oscCurv->numCurves++;

    oscCurv->dual_curve_begin          = NEW_ARRAY(oscCurv->numCurves, DualCurves );
    oscCurv->dual_curve_end            = NEW_ARRAY(oscCurv->numCurves, DualCurves );
    oscCurv->edgeClass                 = NEW_ARRAY(oscCurv->numCurves, int);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edgeClasses[i])
            continue;

        oscCurv->edgeClass[j] = i;
        j++;
    }

    // which curve
    for (i = 0; i < oscCurv->numCurves; i++) {
        oscCurv->dual_curve_begin[i].next    = &oscCurv->dual_curve_end[i];
        oscCurv->dual_curve_begin[i].prev    = NULL;
        oscCurv->dual_curve_end[i].next      = NULL;
        oscCurv->dual_curve_end[i].prev      = &oscCurv->dual_curve_begin[i];
    }

    return oscCurv;
}

void free_oscillating_curves(OscillatingCurves *oscCurv) {
    int i;
    DualCurves *path;
    EdgeNode *node;

    for (i = 0; i < oscCurv->numCurves; i++) {
        while (oscCurv->dual_curve_begin[i].next != &oscCurv->dual_curve_end[i]) {
            path = oscCurv->dual_curve_begin[i].next;
            REMOVE_NODE(path);

            while (path->curves_begin.next != &path->curves_end) {
                node = path->curves_begin.next;
                REMOVE_NODE(node);
                my_free(node);
            }

            my_free(path);
        }
    }

    my_free(oscCurv->dual_curve_begin);
    my_free(oscCurv->dual_curve_end);
    my_free(oscCurv->edgeClass);
    my_free(oscCurv);
}

/*
 * Construct the graph dual to the cusp regions,
 * using region->index to label each vertex, and
 * adding edges using region->adjRegions[].
 */

Graph *construct_cusp_region_dual_graph(ManifoldBoundary *boundary) {
    int i, face;
    CuspRegion *region;

    Graph *graph1 = init_graph(boundary->numCuspRegions, FALSE);

    int *visited = NEW_ARRAY(graph1->nVertices, int);

    for (i = 0; i < graph1->nVertices; i++)
        visited[i] = FALSE;

    // Walk around the cusp triangulation inserting edges
    for (region = boundary->cusp_region_begin.next; region != &boundary->cusp_region_end; region = region->next) {
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

/*
 * flag
 *  - 0: Gluing info
 *  - 1: Dual Curves
 *  - 2: Cusp Regions
 *  - 3: Homology
 *  - 4: Edge Indices
 *  - 5: Dual Curve Paths
 *  - 6: Inside Edge
 *  - 7: Graph
 *  - 8: End points
 */

void print_debug_info(Triangulation *manifold, ManifoldBoundary **cusps, OscillatingCurves *oscCurv, int flag) {
    int i, j, k, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

    CuspTriangle *tri;
    EdgeNode *edge ;
    CuspRegion *region;
    DualCurves *path;
    Graph *g;
    ManifoldBoundary *boundary;

    if (!debug)
        return;

    if (!flag) {
        // Gluing Info
        printf("Triangle gluing info\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);
            boundary = cusps[i];

            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                for (j = 0; j < 4; j++) {
                    if (j == tri->tetVertex)
                        continue;

                    x_vertex1 = (int) remaining_face[tri->tetVertex][j];
                    x_vertex2 = (int) remaining_face[j][tri->tetVertex];
                    y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                    y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                    printf("    (Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
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
        }
    } else if (flag == 1) {
        // Dual Curve endpoints
        printf("Dual Curve Endpoint info\n");
//        for (path = boundary->dual_curve_begin->next; path != boundary->dual_curve_end; path = path->next) {
//            if (path->edgeClass == e0)
//                continue;
//
//            for (j = 0; j < 2; j++) {
//                edge = path->curves[j][START];
//                printf("Tet Vertex: %d (%d %d %d) Tet Vertex: %d (%d %d %d)\n",
//                       path->endpoints[j][START]->region->tetVertex,
//                       path->endpoints[j][START]->vertex,
//                       path->endpoints[j][START]->face,
//                       path->curves[j][START]->next->nextFace,
//                       path->endpoints[j][FINISH]->region->tetVertex,
//                       path->endpoints[j][FINISH]->vertex,
//                       path->endpoints[j][FINISH]->face,
//                       path->curves[j][FINISH]->prev->prevFace);
//            }
//        }
    } else if (flag == 2) {
        // Region Info
        printf("Cusp Region info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            boundary = cusps[i];
            for (region = boundary->cusp_region_begin.next; region != &boundary->cusp_region_end; region = region->next) {
                v1 = edgesThreeToFour[region->tetVertex][0];
                v2 = edgesThreeToFour[region->tetVertex][1];
                v3 = edgesThreeToFour[region->tetVertex][2];

                printf("    Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
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
        }

    } else if (flag == 3) {
        // Homology Info
        printf("Homology info\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            boundary = cusps[i];

            printf("Boundary %d\n", i);
            printf("Intersect Tet Index %d, Intersect Tet Vertex %d\n", boundary->intersectTetIndex, boundary->intersectTetVertex);
            printf("    Meridian\n");

            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tetIndex,
                       tri->tetVertex,
                       tri->tet->curve[M][right_handed][tri->tetVertex][0],
                       tri->tet->curve[M][right_handed][tri->tetVertex][1],
                       tri->tet->curve[M][right_handed][tri->tetVertex][2],
                       tri->tet->curve[M][right_handed][tri->tetVertex][3]
                );
            }
            printf("    Longitude\n");
            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tetIndex,
                       tri->tetVertex,
                       tri->tet->curve[L][right_handed][tri->tetVertex][0],
                       tri->tet->curve[L][right_handed][tri->tetVertex][1],
                       tri->tet->curve[L][right_handed][tri->tetVertex][2],
                       tri->tet->curve[L][right_handed][tri->tetVertex][3]
                );
            }
        }

    } else if (flag == 4) {
        // Edge indices
        printf("Edge classes\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            boundary = cusps[i];
            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                v1 = edgesThreeToFour[tri->tetVertex][0];
                v2 = edgesThreeToFour[tri->tetVertex][1];
                v3 = edgesThreeToFour[tri->tetVertex][2];

                printf("    (Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
                       "Vertex %d: (%d %d), Vertex %d: (%d %d)\n",
                       tri->tetIndex, tri->tetVertex,
                       v1, tri->vertices[v1].edgeClass, tri->vertices[v1].edgeIndex,
                       v2, tri->vertices[v2].edgeClass, tri->vertices[v2].edgeIndex,
                       v3, tri->vertices[v3].edgeClass, tri->vertices[v3].edgeIndex
                );
            }
        }

    } else if (flag == 5) {
        // Dual Curve Paths
        printf("Oscillating curve paths\n");

        // which dual curve
        for (i = 0; i < oscCurv->numCurves; i++) {
            j = 0;

            printf("Dual Curve %d\n", i);
            // which curve component
            for (path = oscCurv->dual_curve_begin[i].next; path != &oscCurv->dual_curve_end[i]; path = path->next) {
                printf("    Part %d: ", j);

                edge = &path->curves_begin;

                while ((edge = edge->next)->next != NULL) {
                    printf("%d ", edge->y);
                }

                printf("\n");
                j++;
            }
        }
    } else if (flag == 6) {
        // Inside Edge Info
        printf("Inside edge info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            boundary = cusps[i];
            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                printf("    (Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
                       tri->tetIndex,               // Tet Index
                       tri->tetVertex,                // Tet Vertex
                       edge3_between_faces[edgesThreeToFour[tri->tetVertex][1]][edgesThreeToFour[tri->tetVertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tetVertex][0]][edgesThreeToFour[tri->tetVertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tetVertex][0]][edgesThreeToFour[tri->tetVertex][1]]
                );
            }
        }
    } else if (flag == 7) {
        printf("Graph info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            boundary = cusps[i];

            printf("Boundary %d\n", i);
            g = boundary->dual_graph;
            for (j = 0; j < g->nVertices; j++) {
                if (g->pRegion[j] == NULL)
                        continue;

                printf("    Vertex %d (Tet Index: %d, Tet Vertex: %d): ",
                       j,
                       g->pRegion[j]->tetIndex,
                       g->pRegion[j]->tetVertex
                );
                edge = &g->edge_list_begin[j];
                while ((edge = edge->next)->next != NULL) {
                    printf("%d ", edge->y);
                }
                printf("\n");
            }
        }
    } else if (flag == 8) {
        // end point info
        printf("EndPoint Info\n");

        // which curve
        for (i = 0; i < oscCurv->numCurves; i++) {
            printf("Dual Curve %d\n", i);

            j = 0;
            // which component
            for (path = oscCurv->dual_curve_begin[i].next; path != &oscCurv->dual_curve_end[i]; path = path->next) {
                printf("    Part %d Cusp %d\n", j, path->endpoints[0].region->tri->tet->cusp[path->endpoints[0].region->tetVertex]->index);
                for (k = 0; k < 2; k++) {
                    if (k == 0)
                        printf("        Start: ");
                    else
                        printf("        End:   ");

                    printf("Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d (%d, %d)\n",
                           path->endpoints[k].regionIndex, path->endpoints[k].region->tetIndex,
                           path->endpoints[k].region->tetVertex, path->endpoints[k].face, path->endpoints[k].vertex,
                           path->endpoints[k].region->tri->vertices[path->endpoints[k].vertex].edgeClass,
                           path->endpoints[k].region->tri->vertices[path->endpoints[k].vertex].edgeIndex);
                }

                j++;
            }
        }
    }
    printf("-------------------------------\n");
}

/*
 * Calculate the number of curves passing around a vertex in the cusp triangulation.
 */

int flow(CuspTriangle *tri, int vertex) {
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

CuspTriangle *find_cusp_triangle(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end,
        CuspTriangle *tri, int face) {
    int tetIndex, tetVertex;
    CuspTriangle *pTri;

    tetIndex = tri->tet->neighbor[face]->index;
    tetVertex = EVALUATE(tri->tet->gluing[face], tri->tetVertex);

    for (pTri = cusp_triangle_begin->next; pTri != cusp_triangle_end; pTri = pTri->next) {
        if (pTri->tetIndex == tetIndex && pTri->tetVertex == tetVertex)
            return pTri;
    }

    // Didn't find a neighbour
    return NULL;
}

void do_oscillating_curves(ManifoldBoundary **cusps, OscillatingCurves *oscCurv, EndMultiGraph *multiGraph) {
    int i;

    for (i = 0; i < oscCurv->numCurves; i++) {
        do_one_dual_curve(cusps, &oscCurv->dual_curve_begin[i], &oscCurv->dual_curve_end[i],
                          multiGraph, oscCurv->edgeClass[i]);

        print_debug_info(cusps[0]->manifold, cusps, oscCurv, 2);
        print_debug_info(cusps[0]->manifold, cusps, oscCurv, 7);
        print_debug_info(cusps[0]->manifold, cusps, oscCurv, 5);
        print_debug_info(cusps[0]->manifold, cusps, oscCurv, 8);
    }

}

void do_one_dual_curve(ManifoldBoundary **cusps, DualCurves *dual_curve_begin,
        DualCurves *dual_curve_end, EndMultiGraph *multiGraph, int edgeClass) {
    int i, pathLen;
    CuspEndPoint *pCuspEndPoint;
    DualCurves *newPath;
    int orientation = START;

    // find paths through cusps
    pCuspEndPoint = find_multi_graph_path(multiGraph->multi_graph, cusps[0]->manifold, edgeClass, multiGraph->e0, &pathLen);

    // find paths inside each cusp
    dual_curve_begin->edgeClass[FINISH] = edgeClass;
    dual_curve_end->edgeClass[START]    = edgeClass;

    for (i = 0; i < pathLen; i++) {
//        print_debug_info(cusps[0]->manifold, cusps, NULL, 2);
//        print_debug_info(cusps[0]->manifold, cusps, NULL, 7);

        newPath = init_dual_curve(dual_curve_end, pCuspEndPoint[i].edgeClass1, pCuspEndPoint[i].edgeClass2);
        find_path_endpoints(cusps[pCuspEndPoint[i].cuspIndex]->dual_graph, dual_curve_begin, newPath, pCuspEndPoint[i].pos, orientation);
        do_one_cusp(cusps[pCuspEndPoint[i].cuspIndex], newPath, edgeClass);
        orientation = (orientation == START ? FINISH : START);
    }

    my_free(pCuspEndPoint);
}

/*
 * Construct oscillating curves on the boundary components
 */

void do_one_cusp(ManifoldBoundary *boundary, DualCurves *path, int edgeClass) {
    int *parent;
    bool *processed, *discovered;

    processed = NEW_ARRAY(boundary->dual_graph->nVertices, bool);
    discovered = NEW_ARRAY(boundary->dual_graph->nVertices, bool);
    parent = NEW_ARRAY(boundary->dual_graph->nVertices, int);

    // Find path using bfs
    init_search(boundary->dual_graph, processed, discovered, parent);
    bfs(boundary->dual_graph, path->endpoints[START].regionIndex, processed, discovered, parent);
    find_path(path->endpoints[START].regionIndex,path->endpoints[FINISH].regionIndex,
              parent, &path->curves_begin);
    update_path_info(boundary->dual_graph, path, edgeClass);

    // Reallocate memory
    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // Split the regions along the path
    split_cusp_regions_along_path(boundary, path);

    free_graph(boundary->dual_graph);
    boundary->dual_graph = construct_cusp_region_dual_graph(boundary);
}

DualCurves *init_dual_curve(DualCurves *curve_end, int edgeClassStart, int edgeClassFinish) {
    DualCurves *newPath = NEW_STRUCT( DualCurves );
    INSERT_BEFORE(newPath, curve_end);

    newPath->curves_begin.next = &newPath->curves_end;
    newPath->curves_begin.prev = NULL;
    newPath->curves_end.next = NULL;
    newPath->curves_end.prev = &newPath->curves_begin;

    newPath->edgeClass[START] = edgeClassStart;
    newPath->edgeClass[FINISH] = edgeClassFinish;

    return newPath;
}

/*
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass. If copy is true, find the path end point
 * corresponding to path1 and store it in path2. Else find a valid path end
 * point and store in path2
 */

void find_path_endpoints(Graph *g, DualCurves *pathStart, DualCurves *path, int pos, int orientation) {
    int endpoint1_edge_class, endpoint1_edge_index, endpoint1_copy, endpoint2_edge_class, endpoint2_edge_index, endpoint2_copy;
    PathEndPoint *endpoint1_path1, *endpoint1_path2, *endpoint2_path1, *endpoint2_path2;
    if (pos == FIRST) {
        endpoint1_path1 = NULL;
        endpoint1_path2 = &path->endpoints[START];
        endpoint1_edge_class = path->edgeClass[START];
        endpoint1_edge_index = START;
        endpoint1_copy = FALSE;

        endpoint2_path1 = NULL;
        endpoint2_path2 = &path->endpoints[FINISH];
        endpoint2_edge_class = path->edgeClass[FINISH];
        endpoint2_edge_index = START;
        endpoint2_copy = FALSE;
    } else if (pos == MIDDLE && orientation == START) {
        endpoint1_path1 = &path->prev->endpoints[START];
        endpoint1_path2 = &path->endpoints[START];
        endpoint1_edge_class = path->edgeClass[START];
        endpoint1_edge_index = FINISH;
        endpoint1_copy = TRUE;

        endpoint2_path1 = &path->prev->endpoints[FINISH];
        endpoint2_path2 = &path->endpoints[FINISH];
        endpoint2_edge_class = path->edgeClass[FINISH];
        endpoint2_edge_index = START;
        endpoint2_copy = FALSE;
    } else if (pos == MIDDLE && orientation == FINISH) {
        endpoint1_path1 = &path->prev->endpoints[FINISH];
        endpoint1_path2 = &path->endpoints[FINISH];
        endpoint1_edge_class = path->edgeClass[START];
        endpoint1_edge_index = FINISH;
        endpoint1_copy = TRUE;

        endpoint2_path1 = &path->prev->endpoints[START];
        endpoint2_path2 = &path->endpoints[START];
        endpoint2_edge_class = path->edgeClass[FINISH];
        endpoint2_edge_index = START;
        endpoint2_copy = FALSE;
    } else if (pos == LAST) {
        endpoint1_path1 = &pathStart->next->endpoints[START];
        endpoint1_path2 = &path->endpoints[START];
        endpoint1_edge_class = pathStart->next->edgeClass[START];
        endpoint1_edge_index = FINISH;
        endpoint1_copy = TRUE;

        endpoint2_path1 = &path->prev->endpoints[FINISH];
        endpoint2_path2 = &path->endpoints[FINISH];
        endpoint2_edge_class = path->prev->edgeClass[FINISH];
        endpoint2_edge_index = FINISH;
        endpoint2_copy = TRUE;
    } else {
        uFatalError("do_one_cusp", "symplectic_basis");
        return;
    }

    find_single_endpoints(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index, endpoint1_copy);
    find_single_endpoints(g, endpoint2_path1, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index, endpoint2_copy);
}

void find_single_endpoints(Graph *g, PathEndPoint *path1, PathEndPoint *path2,
        int edgeClass, int edgeIndex, bool copy) {
    int i, vertex, face1, face2, face;
    bool edge_lies_in_one_cusp;
    CuspRegion *pRegion;

    // which cusp region
    for (i = 0; i < g->nVertices; i++) {
        if (g->pRegion[i] == NULL)
            continue;

        pRegion = g->pRegion[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == pRegion->tetVertex)
                continue;

            if (pRegion->tri->vertices[vertex].edgeClass != edgeClass)
                continue;

            edge_lies_in_one_cusp = (pRegion->tri->tet->cusp[pRegion->tetVertex]->index == pRegion->tri->tet->cusp[vertex]->index);

            if (edge_lies_in_one_cusp && pRegion->tri->vertices[vertex].edgeIndex != edgeIndex) {
                continue;
            }

            face1 = (int) remaining_face[pRegion->tetVertex][vertex];
            face2 = (int) remaining_face[vertex][pRegion->tetVertex];

            // are we in the correct region for copy
            if (copy && (pRegion->tetIndex != path1->region->tetIndex ||
                pRegion->tetVertex != path1->vertex ||
                vertex != path1->region->tetVertex ||
                !pRegion->adjTri[path1->face] ||
                pRegion->curve[path1->face][vertex]))
                continue;

            if (copy) {
                face = path1->face;
            } else if (pRegion->adjTri[face1] && !pRegion->curve[face1][vertex]) {
                face = face1;
            } else if (pRegion->adjTri[face2] && !pRegion->curve[face2][vertex]) {
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

void update_path_info(Graph *g, DualCurves *path, int edgeClass) {
    int face;
    EdgeNode *node = &path->curves_begin;

    // path len 0
    if (node->next->next == NULL)
        return;

    node = node->next;
    // path len 1
    if (node->next->next == NULL) {
        for (face = 0; face < 4; face++)
            if (g->pRegion[node->y]->tetVertex != face &&
                path->endpoints[START].vertex != face &&
                path->endpoints[FINISH].vertex != face)
                break;

        g->pRegion[node->y]->tri->tet->extra[edgeClass].curve[g->pRegion[node->y]->tetVertex]
        [path->endpoints[FINISH].face]++;

        g->pRegion[node->y]->tri->tet->extra[edgeClass].curve[g->pRegion[node->y]->tetVertex]
        [path->endpoints[START].face]--;

        node->insideVertex = face;
        return;
    }


    // Set Header node
    update_path_endpoint_info(g->pRegion[node->y], node, &path->endpoints[START], edgeClass, START);

    for (node = node->next; node->next->next != NULL; node = node->next) {
        inside_vertex(g->pRegion[node->y], node);

        g->pRegion[node->y]->tri->tet->extra[edgeClass].curve[g->pRegion[node->y]->tetVertex][node->nextFace]++;
        g->pRegion[node->y]->tri->tet->extra[edgeClass].curve[g->pRegion[node->y]->tetVertex][node->prevFace]--;
    }

    // Set Tail node
    update_path_endpoint_info(g->pRegion[node->y], node, &path->endpoints[FINISH], edgeClass, FINISH);
}

void update_path_endpoint_info(CuspRegion *pRegion, EdgeNode *node, PathEndPoint *endPoint,
        int edgeClass, int pos) {
    int i, vertex1, vertex2;

    EdgeNode *nextNode;

    if (pos == START)
        nextNode = node->next;
    else
        nextNode = node->prev;

    node->nextFace = -1;
    node->prevFace = -1;

    for (i = 0; i < 4; i++) {
        if (i == pRegion->tetVertex || !pRegion->adjTri[i])
            continue;

        if (pRegion->adjRegions[i]->index == nextNode->y) {
            node->prevFace = i;
            node->nextFace = i;
        }
    }

    // next node isn't in an adjacent region
    if (node->nextFace == -1 || node->prevFace == -1)
        uFatalError("update_path_endpoint_info", "symplectic_basis");

    if (pos == START)
        node->prevFace = endPoint->face;
    else
        node->nextFace = endPoint->face;

    vertex1 = (int) remaining_face[pRegion->tetVertex][endPoint->vertex];
    vertex2 = (int) remaining_face[endPoint->vertex][pRegion->tetVertex];

    if (node->prevFace == endPoint->vertex) {
        node->insideVertex = endPoint->face == vertex1 ? vertex2 : vertex1;
    } else if (node->prevFace == endPoint->face) {
        node->insideVertex = -1;
    } else {
        node->insideVertex = endPoint->vertex;
    }

    if (node->nextFace == node->prevFace)
        return;

    pRegion->tri->tet->extra[edgeClass].curve[pRegion->tetVertex][node->nextFace]++;
    pRegion->tri->tet->extra[edgeClass].curve[pRegion->tetVertex][node->prevFace]--;

}

/*
 * The oscillating curve splits the region it passes through
 * into two regions. Split each region in two and update
 * attributes
 */

void split_cusp_regions_along_path(ManifoldBoundary *boundary, DualCurves *path) {
    int face, index = boundary->numCuspRegions;
    EdgeNode *node = path->curves_begin.next;
    CuspRegion *newRegion;
    Graph *g = boundary->dual_graph;

    // empty path
    if (node->next == NULL)
        return ;

    // path of len 1
    if (node->next->next == NULL) {
        newRegion = NEW_STRUCT(CuspRegion);
        INSERT_BEFORE(newRegion, &boundary->cusp_region_end)
        copy_region(g->pRegion[node->y], newRegion);

        face = node->insideVertex;

        newRegion->index = index;
        newRegion->adjTri[path->endpoints[START].vertex]  = 0;
        newRegion->adjTri[path->endpoints[FINISH].vertex] = 0;
        newRegion->dive[path->endpoints[START].vertex][path->endpoints[FINISH].vertex] = (face != path->endpoints[FINISH].face);
        newRegion->dive[path->endpoints[FINISH].vertex][path->endpoints[START].vertex] = (face != path->endpoints[START].face);

        g->pRegion[node->y]->adjTri[face] = 0;
        g->pRegion[node->y]->dive[face][path->endpoints[START].vertex]  = (face == path->endpoints[START].face);
        g->pRegion[node->y]->dive[face][path->endpoints[FINISH].vertex] = (face == path->endpoints[FINISH].face);

        update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
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
    update_cusp_triangle_endpoints(&boundary->cusp_region_begin, &boundary->cusp_region_end,
                                   g->pRegion[node->y], &path->endpoints[START], node, START);
    newRegion = update_cusp_region(g->pRegion[node->y], node,&path->endpoints[START], index, -1);
    INSERT_BEFORE(newRegion, &boundary->cusp_region_end);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        update_cusp_triangle(&boundary->cusp_region_begin, &boundary->cusp_region_end,g->pRegion[node->y], node);
        newRegion = update_cusp_region(g->pRegion[node->y], node, &path->endpoints[START], index, 0);
        INSERT_BEFORE(newRegion, &boundary->cusp_region_end);
        index++;
    }

    // update last region
    update_cusp_triangle_endpoints(&boundary->cusp_region_begin, &boundary->cusp_region_end,
                                   g->pRegion[node->y], &path->endpoints[FINISH], node, FINISH);
    newRegion = update_cusp_region(g->pRegion[node->y], node, &path->endpoints[FINISH], index, 1);
    INSERT_BEFORE(newRegion, &boundary->cusp_region_end);
    index++;

    update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
    boundary->numCuspRegions = index;
}

/*
 * Set the new and old region data. Draw a picture to see how
 * the attributes change in each case
 *   - flag = -1 : starting node
 *   - flag = 0  : interior node
 *   - flag = 1  : end node
 */

CuspRegion *update_cusp_region(CuspRegion *region, EdgeNode *node,
                                      PathEndPoint *endPoint, int index, int flag) {
    int face, vertex1, vertex2;
    CuspRegion *newRegion = NEW_STRUCT(CuspRegion);

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

void copy_region(CuspRegion *region1, CuspRegion *region2) {
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

void update_cusp_triangle(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end,
        CuspRegion *region, EdgeNode *node) {
    int face1, face2;
    CuspRegion *pRegion;

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

void update_cusp_triangle_endpoints(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end, CuspRegion *region,
        PathEndPoint *endPoint, EdgeNode *node, int pos) {
    int face, face1, face2;
    CuspRegion *pRegion;

    face1 = (int) remaining_face[region->tetVertex][endPoint->vertex];
    face2 = (int) remaining_face[endPoint->vertex][region->tetVertex];

    if (pos == START) {
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

void calculate_holonomy(Triangulation *manifold, int **symp_eqns, int numCurves) {
    int curve, v, f, ff;
    Tetrahedron *tet;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {

        // which curve
        for (curve = 0; curve < numCurves; curve++) {
            // which tet v
            for (v = 0; v < 4; v++) {

                for (f = 0; f < 4; f++) {
                    if (f == v)
                        continue;

                    ff = (int) remaining_face[v][f];

                    symp_eqns[curve][3 * tet->index + edge3_between_faces[f][ff]]
                                += FLOW(tet->extra[curve].curve[v][f], tet->extra[curve].curve[v][ff]);
                }
            }
        }
    }
}

/*
 * node lies in region midNode, find the vertex which the subpath
 * node->prev->y --> node->y --> node->next->y cuts off of the
 * cusp triangle midNode->tri.
 */

void inside_vertex(CuspRegion *midNode, EdgeNode *node) {
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

void init_search(Graph *g, bool *processed, bool *discovered, int *parent) {
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

void bfs(Graph *g, int start, bool *processed, bool *discovered, int *parent) {
    Queue *q = init_queue(10);
    int v, y;
    EdgeNode *p;

    enqueue(q, start);
    discovered[start] = true;

    while (!empty_queue(q)) {
        v = dequeue(q);
        processed[v] = true;
        p = &g->edge_list_begin[v];

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = true;
                parent[y] = v;
            }
        }
    }

    free_queue(q);
}

/*
 * Recover the path through the graph from the parents array
 */

void find_path(int start, int end, int *parents, EdgeNode *node) {
    EdgeNode *new_node = NEW_STRUCT(EdgeNode);
    new_node->edgeClass = -1;

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

Queue *init_queue(int size) {
    Queue *q = NEW_STRUCT( Queue );

    q->front = 0;
    q->rear = -1;
    q->len = 0;
    q->size = MIN(size, 256);
    q->array = NEW_ARRAY(q->size, int);

    return q;
}

Queue *enqueue(Queue *q, int i) {
    // Queue is full
    if ( q->size == q->len ) {
        q = resize_queue(q);
        q = enqueue(q, i);
    } else {
        q->rear = (q->rear + 1) % q->size;
        q->array[q->rear] = i;
        q->len++;
    }

    return q;
}

int dequeue(Queue *q) {
    // User to verify queue is not empty
    int i = q->array[q->front];

    q->front = (q->front + 1) % q->size;
    q->len--;

    return i;
}

int empty_queue(Queue *q) {
    return (!q->len);
}

Queue *resize_queue(Queue *q) {
    int i;
    Queue *p = init_queue(2 * q->size);

    // Copy elements to new array
    while (!empty_queue(q)) {
        i = dequeue(q);
        enqueue(p, i);
    }

    free_queue(q);
    return p;
}

void free_queue(Queue *q) {
    my_free(q->array);
    my_free(q);
}

// Graph Data Structure

/*
 * Initialise Graph
 *
 * Initialise the arrays of the graph 'g' to their default values
 */

Graph *init_graph(int maxVertices, bool directed) {
    int i;
    Graph *g = NEW_STRUCT(Graph);

    g->nVertices = maxVertices;
    g->directed = directed;

    g->edge_list_begin      = NEW_ARRAY(maxVertices, EdgeNode);
    g->edge_list_end        = NEW_ARRAY(maxVertices, EdgeNode);
    g->degree               = NEW_ARRAY(maxVertices, int);
    g->color                = NEW_ARRAY(maxVertices, int);
    g->pRegion              = NEW_ARRAY(maxVertices, CuspRegion *);

    for (i = 0; i < maxVertices; i++) {
        g->degree[i] = 0;
        g->color[i] = -1;

        g->edge_list_begin[i].next     = &g->edge_list_end[i];
        g->edge_list_begin[i].prev     = NULL;
        g->edge_list_end[i].next       = NULL;
        g->edge_list_end[i].prev       = &g->edge_list_begin[i];
        g->pRegion[i]                  = NULL;
    }

    return g;
}

/*
 * Free Graph
 *
 * Free the memory of the graph data structure
 */

void free_graph(Graph *g) {
    int i;
    EdgeNode *edge;

    for (i = 0; i < g->nVertices; i++) {
        while (g->edge_list_begin[i].next != &g->edge_list_end[i]) {
            edge = g->edge_list_begin[i].next;
            REMOVE_NODE(edge);
            my_free(edge);
        }
    }

    my_free(g->edge_list_begin);
    my_free(g->edge_list_end);
    my_free(g->degree);
    my_free(g->color);
    my_free(g->pRegion);
    my_free(g);
}

/*
 * Insert Edge
 *
 * Insert an edge into the graph 'g' from vertex x to y.
 */

int insert_edge(Graph *g, int x, int y, bool directed) {
    // Ignore edge if it already exists
    if (edge_exists(g, x, y))
        return x;

    EdgeNode *p = NEW_STRUCT( EdgeNode);
    INSERT_AFTER(p, &g->edge_list_begin[x]);
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

void delete_edge(Graph *g, int vertex_x, int vertex_y, bool directed) {
    EdgeNode *node, *deleted_node;

    node = &g->edge_list_begin[vertex_x];

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

int edge_exists(Graph *g, int v1, int v2) {
    EdgeNode *node = &g->edge_list_begin[v1];

    while ((node = node->next)->next != NULL) {
        if (node->y == v2) {
            return TRUE;
        }
    }

    return FALSE;
}

// -------------------------------------------------------

// End Multi Graph

EndMultiGraph *init_end_multi_graph(Triangulation *manifold, int *edgeClasses) {
    int i = 0;

    EndMultiGraph *multiGraph = NEW_STRUCT( EndMultiGraph );
    Graph *g = init_graph(manifold->num_cusps, FALSE);
    cusp_graph(manifold, g);
    print_graph(g, 2);

    multiGraph->multi_graph = init_graph(g->nVertices, g->directed);
    int *parent = NEW_ARRAY(multiGraph->multi_graph->nVertices, int);
    spanning_tree(g, multiGraph->multi_graph, 0, parent);
    color_graph(multiGraph->multi_graph);
    int *edges = find_tree_edges(multiGraph->multi_graph, manifold->num_tetrahedra);
    multiGraph->e0 = find_same_color_edge(multiGraph->multi_graph, g, edges);

    for (i = 0; i < manifold->num_tetrahedra; i++) {
        edgeClasses[i] = !edges[i];
    }
    edgeClasses[multiGraph->e0] = 0;

    my_free(edges);
    free_graph(g);
    my_free(parent);
    return multiGraph;
}

void free_end_multi_graph(EndMultiGraph *multiGraph) {
    free_graph(multiGraph->multi_graph);

    my_free(multiGraph);
}

int insert_edge_end_multi_graph(Graph *g, int x, int y, int edgeClass, bool directed) {
    // Ignore edge if it already exists
    if (edge_exists(g, x, y))
        return x;

    EdgeNode *p = NEW_STRUCT( EdgeNode);
    INSERT_AFTER(p, &g->edge_list_begin[x]);
    p->y = y;
    p->edgeClass = edgeClass;
    g->degree[x]++;

    if (!directed) {
        insert_edge_end_multi_graph(g, y, x, edgeClass, TRUE);
    }

    return x;
}


void cusp_graph(Triangulation *manifold, Graph *g) {
    int vertex1, vertex2;
    Tetrahedron *tet;
    EdgeClass *edge;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (vertex1 = 0; vertex1 < 4; vertex1++) {
            // which vertex of the cusp triangle at vertex1
            for (vertex2 = 0; vertex2 < 4; vertex2++) {
                if (vertex1 == vertex2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[vertex1][vertex2]];
                insert_edge_end_multi_graph(g, tet->cusp[vertex1]->index, tet->cusp[vertex2]->index, edge->index, g->directed);
            }
        }
    }
}

/*
 * Find a spanning tree of graph1
 */

void spanning_tree(Graph *graph1, Graph *graph2, int start, int *parent) {
    int i, edgeClass;
    EdgeNode *node;

    bool *processed = NEW_ARRAY(graph1->nVertices, bool);
    bool *discovered = NEW_ARRAY(graph1->nVertices, bool);

    // Find path using bfs
    init_search(graph1, processed, discovered, parent);
    bfs(graph1, start, processed, discovered, parent);

    for (i = 0; i < graph1->nVertices; i++) {
        if (parent[i] == -1)
            continue;

        edgeClass = -1;
        for (node = graph1->edge_list_begin[i].next; node != &graph1->edge_list_end[i]; node = node->next)
            if (node->y == parent[i]) {
                edgeClass = node->edgeClass;
                break;
            }

        if (edgeClass == -1)
            continue;

        insert_edge_end_multi_graph(graph2, i, parent[i], edgeClass, graph2->directed);
    }

    my_free(processed);
    my_free(discovered);
}

void color_graph(Graph *g) {
    int currentColor = 0, v;
    Queue *q = init_queue(g->nVertices);
    EdgeNode *node;

    g->color[0] = currentColor;
    q = enqueue(q, 0);

    while (!empty_queue(q)) {
        v = dequeue(q);
        currentColor = g->color[v];

        for (node = g->edge_list_begin[v].next; node != &g->edge_list_end[v]; node = node->next) {
            // graph is not bipartite
            if (g->color[node->y] == currentColor)
                uFatalError("color_graph", "symplectic_basis");

            if (g->color[node->y] != -1)
                continue;

            g->color[node->y] = !currentColor;
            q = enqueue(q, node->y);
        }
    }

    free_queue(q);
}

int *find_tree_edges(Graph *g, int numEdges) {
    int i, vertex;
    EdgeNode *node;
    int *edges = NEW_ARRAY(numEdges, int);

    for (i = 0; i < numEdges; i++) {
        edges[i] = 0;
    }

    // which vertex
    for (vertex = 0; vertex < g->nVertices; vertex++) {
        // which edge
        for (node = g->edge_list_begin[vertex].next; node != &g->edge_list_end[vertex]; node = node->next) {
            edges[node->edgeClass] = 1;
        }
    }

    return edges;
}

/*
 * g1 is the colored spanning tree of g2, return the
 * edge class of the edge in g2 which connects
 * vertices in g1 of the same color
 */

int find_same_color_edge(Graph *g1, Graph *g2, int *edgeClasses) {
    int v;
    EdgeNode *node;

    for (v = 0; v < g1->nVertices; v++) {
        for (node = g2->edge_list_begin[v].next; node != &g2->edge_list_end[v]; node = node->next) {
            if (g1->color[v] == g1->color[node->y] && !edgeClasses[node->edgeClass])
                // we found an edge
                return node->edgeClass;
        }
    }

    // we didn't find an edge connecting vertices of the same color
    uFatalError("find_same_color_edge", "symplectic_basis");
    return -1;
}

/*
 * Find the length of a path between start and end
 */

int find_path_len(int start, int end, int *parents, int pathLen) {
    if ((start == end) || (end == -1)) {
        return pathLen;
    } else {
        return find_path_len(start, parents[end], parents, pathLen + 1);
    }
}

CuspEndPoint *find_multi_graph_path(Graph *g, Triangulation *manifold, int edgeClass, int e0, int *pathLen) {
    bool *processed     = NEW_ARRAY(g->nVertices, bool);
    bool *discovered    = NEW_ARRAY(g->nVertices, bool);
    int *parent         = NEW_ARRAY(g->nVertices, int);
    int start, end, startE0, endE0;
    EdgeNode *node_begin = NEW_STRUCT( EdgeNode ), *node_end = NEW_STRUCT( EdgeNode), *tempNode;

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    find_edge_ends(g, manifold, edgeClass, &start, &end);
    find_edge_ends(g, manifold, e0, &startE0, &endE0);

    init_search(g, processed, discovered, parent);
    bfs(g, start, processed, discovered, parent);

    *pathLen = 0;
    *pathLen = find_path_len(start, end, parent, *pathLen);
    tempNode = NEW_STRUCT( EdgeNode );

    if (*pathLen % 2 == 1) {
        find_path(start, end, parent, node_begin);
    } else {
        init_search(g, processed, discovered, parent);
        bfs(g, start, processed, discovered, parent);

        find_path(start, startE0, parent, node_begin);

        my_free(tempNode);
        tempNode = node_end->prev;
//        tempNode->y = startE0;
//        tempNode->edgeClass = e0;
//        INSERT_BEFORE(tempNode, node_end);

        init_search(g, processed, discovered, parent);
        bfs(g, endE0, processed, discovered, parent);

        find_path(endE0, end, parent, node_end->prev);
        tempNode->next->edgeClass = e0;
    }

    CuspEndPoint *pCuspEndPoint = find_cusp_endpoint_edge_classes(g, node_begin, node_end, edgeClass, pathLen);

    // remove header and tailer nodes
    while (node_begin->next != node_end) {
        tempNode = node_begin->next;
        REMOVE_NODE(tempNode);
        my_free(tempNode);
    }

    my_free(node_begin);
    my_free(node_end);

    my_free(parent);
    my_free(discovered);
    my_free(processed);

    return pCuspEndPoint;
}

CuspEndPoint *find_cusp_endpoint_edge_classes(Graph *g, EdgeNode *node_begin,
        EdgeNode *node_end, int edgeClass, int *pathLen) {
    int i, cusp;
    EdgeNode *edge, *node;

    *pathLen = 0;
    for (node = node_begin->next; node != node_end; node = node->next)
        (*pathLen)++;

    CuspEndPoint *pCuspEndPoint = NEW_ARRAY(*pathLen, CuspEndPoint);

    pCuspEndPoint[0].edgeClass1 = edgeClass;

    i = 0;
    for (node = node_begin->next; node->next != node_end; node = node->next) {
        cusp = node->y;

        pCuspEndPoint[i].cuspIndex = node->y;
        pCuspEndPoint[i].pos = MIDDLE;

        if (i != 0)
            pCuspEndPoint[i].edgeClass1 = pCuspEndPoint[i - 1].edgeClass2;

        if (node->next->edgeClass != -1) {
            pCuspEndPoint[i].edgeClass2 = node->next->edgeClass;
            i++;
            continue;
        }

        for (edge = g->edge_list_begin[cusp].next; edge != &g->edge_list_end[cusp]; edge = edge->next) {
            if (node->next->y == edge->y) {
                pCuspEndPoint[i].edgeClass2 = edge->edgeClass;
                break;
            }
        }

        i++;
    }

    pCuspEndPoint[0].pos = FIRST;
    pCuspEndPoint[*pathLen - 1].cuspIndex = node->y;
    pCuspEndPoint[*pathLen - 1].edgeClass1 = pCuspEndPoint[*pathLen - 2].edgeClass2;
    pCuspEndPoint[*pathLen - 1].edgeClass2 = edgeClass;
    pCuspEndPoint[*pathLen - 1].pos = LAST;

    return pCuspEndPoint;
}

void find_edge_ends(Graph *g, Triangulation *manifold, int edgeClass, int *start, int *end) {
    int vertex1, vertex2;
    Tetrahedron *tet;
    EdgeClass *edge;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (vertex1 = 0; vertex1 < 4; vertex1++) {
            // which vertex of the cusp triangle at vertex1
            for (vertex2 = 0; vertex2 < 4; vertex2++) {
                if (vertex1 == vertex2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[vertex1][vertex2]];
                if (edge->index != edgeClass)
                    continue;

                *start = tet->cusp[vertex1]->index;
                *end   = tet->cusp[vertex2]->index;
                return;
            }
        }
    }


    // didn't find the edge class in the graph
    uFatalError("find_edge_ends", "symplectic_basis");
}

void print_graph(Graph *g, int flag) {
    int i;
    EdgeNode *edge;

    if (!debug)
        return;

    if (flag == 0) {
        printf("Multi Graph Info\n");
    } else if (flag == 1) {
        printf("Double Graph Info\n");
    } else if (flag == 2) {
        printf("Cusp Graph Info\n");
    }

    for (i = 0; i < g->nVertices; i++) {
        printf("(Cusp Index: %d) ", i);
        edge = &g->edge_list_begin[i];
        while ((edge = edge->next)->next != NULL) {
            printf("%d ", edge->y);
        }
        printf("\n");
    }

    printf("--------------------------\n");
}
