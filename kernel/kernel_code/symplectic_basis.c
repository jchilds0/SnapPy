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
    int                     edge_class;              /** edge_class of the edge for edges in the end multi graph */
    int                     next_face;               /** face the path crosses to the next node */
    int                     prev_face;               /** face the path crosses to the prev node */
    int                     inside_vertex;           /** inside vertex of the path */
    struct EdgeNode         *next;                  /** next node in doubly linked list */
    struct EdgeNode         *prev;                  /** prev node in doubly linked list */
} EdgeNode;

typedef struct Graph {
    EdgeNode                *edge_list_begin;       /** header node of doubly linked list */
    EdgeNode                *edge_list_end;         /** tail node of doubly linked list */
    struct CuspRegion       **regions;              /** list of regions in the graph */
    int                     *degree;                /** degree of each vertex */
    int                     *color;                 /** color a tree bipartite */
    int                     num_vertices;              /** number of vertices in the graph */
    int                     directed;               /** is the graph directed */
} Graph;

typedef struct CuspEndPoint {
    int                     cusp_index;
    int                     edge_class1;
    int                     edge_class2;
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
    int                     region_index;            /** index of the region the endpoint lies in */
    int                     num_adj_curves[4][4];         /** where the curve dives into the manifold */
    struct CuspRegion       *region;                /** pointer to the region the endpoint lies in */
} PathEndPoint;

typedef struct DualCurves {
    int                     edge_class[2];
    int                     cusp_index;             /** which cusp does the curve lie in */
    EdgeNode                curves_begin;           /** header node of doubbly linked list */
    EdgeNode                curves_end;             /** tailer node of doubbly linked list */
    PathEndPoint            endpoints[2];           /** path end points */
    struct DualCurves       *next;                  /** next dual curve in doubly linked list */
    struct DualCurves       *prev;                  /** prev dual curve in doubly linked list */
} DualCurves;

typedef struct OscillatingCurves {
    int                     num_curves;
    int                     *edge_class;
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
    int                     edge_class;
    int                     edge_index;
    EdgeClass               *edge;
    int                     v1;
    int                     v2;
} CuspVertex;

typedef struct CuspTriangle {
    Tetrahedron             *tet;                   /** tetrahedron the triangle comes from */
    int                     tet_index;               /** tet->index */
    int                     tet_vertex;              /** vertex the triangle comes from */
    int                     num_curves;              /** number of curves on the triangle */
    CuspVertex              vertices[4];            /** information about each vertex */
    struct CuspTriangle     *neighbours[4];         /** triangle neighbouring a face */
    struct CuspTriangle     *next;                  /** next cusp triangle on doubly linked list */
    struct CuspTriangle     *prev;                  /** prev cusp triangle on doubly linkled list */
} CuspTriangle;

typedef struct CuspRegion {
    CuspTriangle            *tri;                   /** cusp triangle the region lies on */
    int                     tet_index;               /** tri->tetIndex */
    int                     tet_vertex;              /** tri->tet_vertex */
    int                     index;                  /** index of the cusp region */
    int                     curve[4][4];            /** looking at face, number of curves between the region and vertex */
    int                     adj_cusp_triangle[4];              /** does the region meet this edge of the cusp triangle */
    int                     dive[4][4];             /** can we dive along the face into this vertex */
    int                     num_adj_curves[4][4];            /** stores the number of curves between a region and a face */
    int                     temp_adj_curves[4][4];        /** store the adj curve until pathfinding is complete */
    struct CuspRegion       *adj_cusp_regions[4];         /** index of the adjacent regions */
    struct CuspRegion       *adj_dive_regions[4];        /** regions which are adjacent by diving through the manifold */
    struct CuspRegion       *next;                  /** next cusp region in doubly linked list */
    struct CuspRegion       *prev;                  /** prev cusp region in doubly linked list */
} CuspRegion;

typedef struct ManifoldBoundary {
    int                     intersect_tet_index;      /** index of the intersection triangle */
    int                     intersect_tet_vertex;     /** vertex of the intersection triangle */
    int                     num_edge_classes;         /** number of edge classes in the boundary */
    int                     num_cusp_triangles;       /** number of cusp triangle in the boundary */
    int                     num_cusp_regions;         /** number of cusp regions in the boundary */
    Triangulation           *manifold;              /** manifold */
    Cusp                    *cusp;                  /** which cusp is the boundary in */
    Graph                   *dual_graph;            /** dual graph of the cusp region */
    CuspTriangle            cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    CuspTriangle            cusp_triangle_end;      /** tail node of doubly linked list of cusp triangles */
    CuspRegion              cusp_region_begin;      /** header node of doubly linked list of cusp regions */
    CuspRegion              cusp_region_end;        /** tail node of doubly linked list of cusp regions */
} ManifoldBoundary;

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
 * Graph
 */

Graph                   *init_graph(int, bool);
void                    free_graph(Graph *);
int                     insert_edge(Graph *, int, int, bool);
void                    delete_edge(Graph *, int, int, bool);
int                     edge_exists(Graph *, int, int);

/**
 * Breadth First Search
 */

void                    init_search(Graph *, bool *, bool *, int *);
void                    bfs(Graph *, int, bool *, bool *, int *);
void                    find_path(int, int, int *, EdgeNode *);

/**
 * Symplectic Basis
 */
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
void                    do_one_dual_curve(ManifoldBoundary **, OscillatingCurves *, DualCurves *, DualCurves *, EndMultiGraph *, int, int);
void                    do_one_cusp(ManifoldBoundary *, DualCurves *, int);
Graph *                 construct_cusp_region_dual_graph(ManifoldBoundary *);
void                    print_debug_info(Triangulation *, ManifoldBoundary **, OscillatingCurves *, int);
DualCurves              *init_dual_curve(DualCurves *, int, int);
void                    find_path_endpoints(Graph *, DualCurves *, DualCurves *, int, int);
PathEndPoint            *find_single_endpoint(Graph *, PathEndPoint *, int, int);
PathEndPoint            *find_single_matching_endpoint(Graph *, PathEndPoint *, PathEndPoint *, int, int);
void                    update_path_info(Graph *g, DualCurves *, int);
void                    update_path_endpoint_info(CuspRegion *, EdgeNode *, PathEndPoint *, int, int);
void                    split_cusp_regions_along_path(ManifoldBoundary *, DualCurves *);
CuspRegion              *split_cusp_region_path_interior(CuspRegion *, EdgeNode *, PathEndPoint *, int);
CuspRegion              *split_cusp_region_path_endpoint(CuspRegion *, EdgeNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle_path_interior(CuspRegion *, CuspRegion *, CuspRegion *, EdgeNode *);
void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, EdgeNode *, int);
void                    copy_region(CuspRegion *, CuspRegion *);
void                    update_adj_curve_along_path(ManifoldBoundary **, OscillatingCurves *, DualCurves *, DualCurves *, int);
void                    update_adj_curve_at_endpoint(PathEndPoint *, DualCurves *, DualCurves *);
void                    calculate_holonomy(Triangulation *, int **, int);
void                    find_interior_vertex(CuspRegion *, EdgeNode *);


/**
 * End Multi Graph
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

Graph *init_graph(int max_vertices, bool directed) {
    int i;
    Graph *g = NEW_STRUCT(Graph);

    g->num_vertices = max_vertices;
    g->directed = directed;

    g->edge_list_begin      = NEW_ARRAY(max_vertices, EdgeNode);
    g->edge_list_end        = NEW_ARRAY(max_vertices, EdgeNode);
    g->degree               = NEW_ARRAY(max_vertices, int);
    g->color                = NEW_ARRAY(max_vertices, int);
    g->regions              = NEW_ARRAY(max_vertices, CuspRegion *);

    for (i = 0; i < max_vertices; i++) {
        g->degree[i] = 0;
        g->color[i] = -1;

        g->edge_list_begin[i].next     = &g->edge_list_end[i];
        g->edge_list_begin[i].prev     = NULL;
        g->edge_list_end[i].next       = NULL;
        g->edge_list_end[i].prev       = &g->edge_list_begin[i];
        g->regions[i]                  = NULL;
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

    for (i = 0; i < g->num_vertices; i++) {
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
    my_free(g->regions);
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

// ---------------------------------------------------------------

// Breadth First Search from Skiena Algorithm Design Manual

/*
 * Initialise default values for bfs arrays
 */

void init_search(Graph *g, bool *processed, bool *discovered, int *parent) {
    int i;

    for (i = 0; i < g->num_vertices; i ++) {
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
    new_node->edge_class = -1;

    INSERT_AFTER(new_node, node);

    if ((start == end) || (end == -1)) {
        new_node->y = start;
    } else {
        find_path(start, parents[end], parents, node);
        new_node->y = end;
    }
}

/*
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the get_gluing_equations and get_symplectic_equations functions
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i, j;
    int *edge_classes = NEW_ARRAY(manifold->num_tetrahedra, int);
    peripheral_curves(manifold);

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, edge_classes, &symp_num_rows, *num_cols);

    // Construct return array
    *num_rows = 2 * (manifold->num_tetrahedra - manifold->num_cusps);
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i]) {
            my_free(symp_eqns[i]);
            continue;
        }

        eqns[2 * j] = gluing_equations_for_edge_class(manifold, i);
        eqns[2 * j + 1] = symp_eqns[i];
        j++;
    }

    my_free(edge_classes);
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

    ManifoldBoundary **cusps                = NEW_ARRAY(manifold->num_cusps, ManifoldBoundary *);
    Cusp *cusp                              = manifold->cusp_list_begin.next;
    EndMultiGraph *end_multi_graph          = init_end_multi_graph(manifold, edgeClasses);
    OscillatingCurves *oscillating_curves   = init_oscillating_curves(manifold, edgeClasses);

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

    print_debug_info(manifold, cusps, oscillating_curves, 0);       // Gluing
    print_debug_info(manifold, cusps, oscillating_curves, 3);       // Homology
    print_debug_info(manifold, cusps, oscillating_curves, 4);       // Edge classes
    print_debug_info(manifold, cusps, oscillating_curves, 6);       // Inside Edge
    print_debug_info(manifold, cusps, oscillating_curves, 2);       // Regions

    // Allocate Symplectic Equations Array
    *num_rows = oscillating_curves->num_curves;
    int **symp_eqns = NEW_ARRAY(manifold->num_tetrahedra, int *);

    for (i = 0; i < manifold->num_tetrahedra; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    do_oscillating_curves(cusps, oscillating_curves, end_multi_graph);
    calculate_holonomy(manifold, symp_eqns, manifold->num_tetrahedra);

    free_end_multi_graph(end_multi_graph);
    free_oscillating_curves(oscillating_curves);
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
    boundary->num_edge_classes = manifold->num_tetrahedra;
    boundary->num_cusp_triangles = 0;
    boundary->num_cusp_regions = 0;

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

                        boundary->intersect_tet_index = cusp->basepoint_tet->index;
                        boundary->intersect_tet_vertex = (int) cusp->basepoint_vertex;
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
            tri->tet_index = tri->tet->index;
            tri->tet_vertex = vertex;

            tri->num_curves = flow(tri, edgesThreeToFour[tri->tet_vertex][0])
                              + flow(tri, edgesThreeToFour[tri->tet_vertex][1])
                              + flow(tri, edgesThreeToFour[tri->tet_vertex][2]);

            for (face = 0; face < 4; face ++) {
                if (tri->tet_vertex == face)
                    continue;

                tri->vertices[face].v1 = tri->tet_vertex;
                tri->vertices[face].v2 = face;
                tri->vertices[face].edge = tri->tet->edge_class[edge_between_vertices[tri->vertices[face].v1][tri->vertices[face].v2]];
                tri->vertices[face].edge_class = tri->vertices[face].edge->index;
                tri->vertices[face].edge_index = -1;
            }
        }
    }

    // which cusp triangle
    for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
        // which vertex
        for (face = 0; face < 4; face++) {
            if (face == tri->tet_vertex)
                continue;

            tri->neighbours[face] = find_cusp_triangle(&boundary->cusp_triangle_begin, &boundary->cusp_triangle_end, tri, face);
        }
    }

    label_cusp_vertex_indices(&boundary->cusp_triangle_begin, &boundary->cusp_triangle_end, boundary->num_edge_classes);
    boundary->num_cusp_triangles = index;
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

    int *current_index = NEW_ARRAY(numEdgeClasses, int);

    for (i = 0; i < numEdgeClasses; i++)
        current_index[i] = 0;

    for (tri = cusp_triangle_begin->next; tri != cusp_triangle_end; tri = tri->next) {
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tet_vertex || tri->vertices[vertex].edge_index != -1)
                continue;

            walk_around_cusp_vertex(tri, vertex, current_index[tri->vertices[vertex].edge_class]);
            current_index[tri->vertices[vertex].edge_class]++;
        }
    }

    my_free(current_index);
}

/*
 * Walk around vertex cusp_vertex of triangle *tri
 * and set edge_index to index.
 */

void walk_around_cusp_vertex(CuspTriangle *tri, int cusp_vertex, int index) {
    int gluing_vertex, outside_vertex, old_gluing_vertex, old_cusp_vertex, old_outside_vertex;
    gluing_vertex = (int) remaining_face[cusp_vertex][tri->tet_vertex];
    outside_vertex = (int) remaining_face[tri->tet_vertex][cusp_vertex];

    while (tri->vertices[cusp_vertex].edge_index == -1) {
        tri->vertices[cusp_vertex].edge_index = index;

        // Move to the next cusp triangle
        old_cusp_vertex         = cusp_vertex;
        old_gluing_vertex       = gluing_vertex;
        old_outside_vertex      = outside_vertex;

        cusp_vertex          = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
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
        if (tri->tet_index == boundary->intersect_tet_index && tri->tet_vertex == boundary->intersect_tet_vertex) {
            index = init_intersect_cusp_region(boundary, tri, index);
            continue;
        }

        index = init_normal_cusp_region(boundary, tri, index);
    }

    update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
    boundary->num_cusp_regions = index;
}

int init_intersect_cusp_region(ManifoldBoundary *boundary, CuspTriangle *tri, int index) {
    int i, curve_index, vertex, v1, v2, v3;
    int adjTri[4], dist[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 1; i < flow(tri, vertex); i++) {
            for (curve_index = 0; curve_index < 2; curve_index++) {
                dist[v1]           = i;
                dist[v2]           = MIN(dist[v1], 2 * flow(tri, vertex) - dist[v1])
                                     + flow(tri, v2) + flow(tri, v1);
                dist[vertex]            = flow(tri, vertex) - dist[v1] + flow(tri, v1);
                dist[tri->tet_vertex]    = -1;

                adjTri[v1]         = 1;
                adjTri[v2]         = 0;
                adjTri[vertex]          = 0;
                adjTri[tri->tet_vertex]  = -1;

                set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
                index++;

                // Swap vertices
                v1 = (int) remaining_face[vertex][tri->tet_vertex];
                v2 = (int) remaining_face[tri->tet_vertex][vertex];
            }
        }

        // Region in the middle of face vertex
        if (flow(tri, v1) && flow(tri, v2)) {
            dist[v1]           = flow(tri, v2);
            dist[v2]           = flow(tri, v1);
            dist[vertex]            = MIN(flow(tri, v1) + dist[v1],
                                          flow(tri, v2) + dist[v2]) + flow(tri, vertex);
            dist[tri->tet_vertex]    = -1;

            adjTri[v1]         = 0;
            adjTri[v2]         = 0;
            adjTri[vertex]          = 1;
            adjTri[tri->tet_vertex]  = -1;

            set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }
    }

    // Region of dist 0 to vertex
    v1 = edgesThreeToFour[tri->tet_vertex][0];
    v2 = edgesThreeToFour[tri->tet_vertex][1];
    v3 = edgesThreeToFour[tri->tet_vertex][2];

    // Case 1; two regions
    if (atleast_two(!flow(tri, v1), !flow(tri, v2), !flow(tri, v3))) {
        dist[v1]           = flow(tri, v1);
        dist[v2]           = flow(tri, v2);
        dist[v3]           = flow(tri, v3);
        dist[tri->tet_vertex]    = -1;

        adjTri[v1]         = 1;
        adjTri[v2]         = 1;
        adjTri[v3]         = 1;
        adjTri[tri->tet_vertex]  = -1;

        set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
        index++;

        // Find vertex with non-zero flow
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tet_vertex)
                continue;

            if (flow(tri, vertex)) {
                v1 = vertex;
                v2 = (int) remaining_face[tri->tet_vertex][v1];
                v3 = (int) remaining_face[v1][tri->tet_vertex];
                break;
            }
        }
        dist[v1]           = 0;
        dist[v2]           = flow(tri, v1);
        dist[v3]           = flow(tri, v1);
        dist[tri->tet_vertex]    = -1;

        adjTri[v1]         = 0;
        adjTri[v2]         = 1;
        adjTri[v3]         = 1;
        adjTri[tri->tet_vertex]  = 0;

        set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
        index++;
    } else {
        // Case 2: three regions
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tet_vertex)
                continue;

            v1 = (int) remaining_face[tri->tet_vertex][vertex];
            v2 = (int) remaining_face[vertex][tri->tet_vertex];

            dist[vertex]            = 0;
            dist[v1]           = flow(tri, vertex) + flow(tri, v1);
            dist[v2]           = flow(tri, vertex) + flow(tri, v2);
            dist[tri->tet_vertex]    = -1;

            adjTri[vertex]          = 0;
            adjTri[v1]         = 1;
            adjTri[v2]         = 1;
            adjTri[tri->tet_vertex]  = 0;

            set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }
    }

    return index;
}

int init_normal_cusp_region(ManifoldBoundary *boundary, CuspTriangle *tri, int index) {
    int i, vertex, v1, v2;
    int dist[4], adjTri[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 0; i < flow(tri, vertex); i++) {
            dist[vertex] = i;
            dist[v1] = flow(tri, v1) + (flow(tri, vertex) - dist[vertex]);
            dist[v2] = flow(tri, v2) + (flow(tri, vertex) - dist[vertex]);
            dist[tri->tet_vertex] = -1;

            adjTri[vertex] = 0;
            adjTri[v1] = 1;
            adjTri[v2] = 1;
            adjTri[tri->tet_vertex] = 0;

            set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
            index++;
        }

    }

    // Center vertex
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        dist[vertex] = flow(tri, vertex);
        adjTri[vertex] = 1;
    }

    dist[tri->tet_vertex] = -1;
    adjTri[tri->tet_vertex] = 0;

    set_cusp_region_data(&boundary->cusp_region_end, tri, dist, adjTri, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which
 * allocates the attributes of the cusp region
 */

void set_cusp_region_data(CuspRegion *cusp_region_end, CuspTriangle *tri, int distance[4], int adj_cusp_triangle[4], int index) {
    int i, j, v1, v2, v3;
    CuspRegion *pRegion = NEW_STRUCT( CuspRegion );
    INSERT_BEFORE(pRegion, cusp_region_end);

    pRegion->tri                 = tri;
    pRegion->tet_index            = pRegion->tri->tet_index;
    pRegion->tet_vertex           = pRegion->tri->tet_vertex;
    pRegion->index               = index;

    // default values
    for (i = 0; i < 4; i++) {
        pRegion->adj_cusp_triangle[i] = 0;

        for (j = 0; j < 4; j++) {
            pRegion->curve[i][j] = -1;
            pRegion->dive[i][j] = 0;
            pRegion->num_adj_curves[i][j] = 0;
            pRegion->temp_adj_curves[i][j] = 0;
        }
    }

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tet_vertex][i];
        v2 = edgesThreeToFour[tri->tet_vertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tet_vertex][(i + 2) % 3];

        pRegion->curve[v2][v1] = distance[v1];
        pRegion->curve[v3][v1] = distance[v1];
        pRegion->dive[v2][v1] = (distance[v1] ? 0 : 1);
        pRegion->dive[v3][v1] = (distance[v1] ? 0 : 1);

        pRegion->adj_cusp_triangle[v1] = adj_cusp_triangle[v1];
    }
}

/*
 * Calculate which regions are located across cusp edges
 * and store the result in the adj_cusp_regions attribute
 */

void update_adj_region_data(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end) {
    int j;
    CuspRegion *pRegion;

    // Add adjacent region info
    for (pRegion = cusp_region_begin->next; pRegion != cusp_region_end; pRegion = pRegion->next) {
        for (j = 0; j < 4; j++) {
            if (!pRegion->adj_cusp_triangle[j] || pRegion->tet_vertex == j) {
                pRegion->adj_cusp_regions[j] = NULL;
                continue;
            }

            pRegion->adj_cusp_regions[j] = find_adj_region(cusp_region_begin, cusp_region_end, pRegion, j);
        }
    }
}

/*
 * Find the cusp region which is adjacent
 * to x across face.
 */

CuspRegion *find_adj_region(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end,
                            CuspRegion *x, int face) {
    int v1, v2, y_vertex1, y_vertex2, y_face, distance_v1, distance_v2, adj_face, tet_index, tet_vertex;
    CuspTriangle *tri = x->tri;
    CuspRegion *y;

    v1 = (int) remaining_face[tri->tet_vertex][face];
    v2 = (int) remaining_face[face][tri->tet_vertex];

    for (y = cusp_region_begin->next; y != cusp_region_end; y = y->next) {
        y_vertex1    = EVALUATE(tri->tet->gluing[face], v1);
        y_vertex2    = EVALUATE(tri->tet->gluing[face], v2);
        y_face       = EVALUATE(tri->tet->gluing[face], face);

        tet_index    = (tri->neighbours[face]->tet_index == y->tet_index);
        tet_vertex   = (tri->neighbours[face]->tet_vertex == y->tet_vertex);

        if (!tet_index || !tet_vertex)
            continue;

        distance_v1      = (x->curve[face][v1] == y->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == y->curve[y_face][y_vertex2]);
        adj_face         = y->adj_cusp_triangle[y_face];

        // missing distance
        if (y->curve[y_face][y_vertex1] == -1 || y->curve[y_face][y_vertex2] == -1)
            uFatalError("find_adj_region", "symplectic_basis");

        if (distance_v1 && distance_v2 && adj_face)
            return y;
    }

    // We didn't find a cusp y
    //uFatalError("find_cusp_region", "symplectic_basis");
    return NULL;
}

/*
 * Initialise dual curve doubly linked list which
 * stores the oscillating curves on the cusp
 */

OscillatingCurves *init_oscillating_curves(Triangulation *manifold, int *edge_classes) {
    int i, j;
    OscillatingCurves *curves = NEW_STRUCT(OscillatingCurves );

    curves->num_curves = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++)
        if (edge_classes[i])
            curves->num_curves++;

    curves->dual_curve_begin          = NEW_ARRAY(curves->num_curves, DualCurves );
    curves->dual_curve_end            = NEW_ARRAY(curves->num_curves, DualCurves );
    curves->edge_class                 = NEW_ARRAY(curves->num_curves, int);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i])
            continue;

        curves->edge_class[j] = i;
        j++;
    }

    // which curve
    for (i = 0; i < curves->num_curves; i++) {
        curves->dual_curve_begin[i].next    = &curves->dual_curve_end[i];
        curves->dual_curve_begin[i].prev    = NULL;
        curves->dual_curve_end[i].next      = NULL;
        curves->dual_curve_end[i].prev      = &curves->dual_curve_begin[i];
    }

    return curves;
}

void free_oscillating_curves(OscillatingCurves *curves) {
    int i;
    DualCurves *path;
    EdgeNode *node;

    for (i = 0; i < curves->num_curves; i++) {
        while (curves->dual_curve_begin[i].next != &curves->dual_curve_end[i]) {
            path = curves->dual_curve_begin[i].next;
            REMOVE_NODE(path);

            while (path->curves_begin.next != &path->curves_end) {
                node = path->curves_begin.next;
                REMOVE_NODE(node);
                my_free(node);
            }

            my_free(path);
        }
    }

    my_free(curves->dual_curve_begin);
    my_free(curves->dual_curve_end);
    my_free(curves->edge_class);
    my_free(curves);
}

/*
 * Construct the graph dual to the cusp regions,
 * using region->index to label each vertex, and
 * adding edges using region->adj_cusp_regions[].
 */

Graph *construct_cusp_region_dual_graph(ManifoldBoundary *boundary) {
    int i, face;
    CuspRegion *region;

    Graph *graph1 = init_graph(boundary->num_cusp_regions, FALSE);

    int *visited = NEW_ARRAY(graph1->num_vertices, int);

    for (i = 0; i < graph1->num_vertices; i++)
        visited[i] = FALSE;

    // Walk around the cusp triangulation inserting edges
    for (region = boundary->cusp_region_begin.next; region != &boundary->cusp_region_end; region = region->next) {
        if (visited[region->index])
            continue;

        for (face = 0; face < 4; face++) {
            if (!region->adj_cusp_triangle[face])
                continue;

            // Missing adj region data
            if (region->adj_cusp_regions[face] == NULL)
                uFatalError("construct_cusp_region_dual_graph", "symplectic_basis");

            insert_edge(graph1, region->index, region->adj_cusp_regions[face]->index, graph1->directed);
            graph1->regions[region->index] = region;
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

void print_debug_info(Triangulation *manifold, ManifoldBoundary **cusps, OscillatingCurves *curves, int flag) {
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
                    if (j == tri->tet_vertex)
                        continue;

                    x_vertex1 = (int) remaining_face[tri->tet_vertex][j];
                    x_vertex2 = (int) remaining_face[j][tri->tet_vertex];
                    y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                    y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                    printf("    (Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                           "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                           tri->tet_index,               // Tet Index
                           tri->tet_vertex,                // Tet Vertex
                           j,      // Cusp Edge
                           tri->tet->neighbor[j]->index,                              // Tet Index
                           EVALUATE(tri->tet->gluing[j], tri->tet_vertex),             // Tet Vertex
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
//            if (path->edge_class == e0)
//                continue;
//
//            for (j = 0; j < 2; j++) {
//                edge = path->curves[j][START];
//                printf("Tet Vertex: %d (%d %d %d) Tet Vertex: %d (%d %d %d)\n",
//                       path->endpoints[j][START]->region->tet_vertex,
//                       path->endpoints[j][START]->vertex,
//                       path->endpoints[j][START]->face,
//                       path->curves[j][START]->next->next_face,
//                       path->endpoints[j][FINISH]->region->tet_vertex,
//                       path->endpoints[j][FINISH]->vertex,
//                       path->endpoints[j][FINISH]->face,
//                       path->curves[j][FINISH]->prev->prev_face);
//            }
//        }
    } else if (flag == 2) {
        // Region Info
        printf("Cusp Region info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            boundary = cusps[i];
            for (region = boundary->cusp_region_begin.next; region != &boundary->cusp_region_end; region = region->next) {
                v1 = edgesThreeToFour[region->tet_vertex][0];
                v2 = edgesThreeToFour[region->tet_vertex][1];
                v3 = edgesThreeToFour[region->tet_vertex][2];

                printf("    Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
                       " (Adj Curves: [%d %d] [%d %d] [%d %d]) (Curves: [[%d %d] [%d %d] [%d %d]]) (Dive: [[%d %d] [%d %d] [%d %d]])\n",
                       region->index, region->tet_index, region->tet_vertex,
                       region->adj_cusp_triangle[v1], region->adj_cusp_triangle[v2], region->adj_cusp_triangle[v3],
                       region->adj_cusp_regions[v1] == NULL ? -1 : region->adj_cusp_regions[v1]->index,
                       region->adj_cusp_regions[v2] == NULL ? -1 : region->adj_cusp_regions[v2]->index,
                       region->adj_cusp_regions[v3] == NULL ? -1 : region->adj_cusp_regions[v3]->index,
                       region->num_adj_curves[v2][v1], region->num_adj_curves[v3][v1],
                       region->num_adj_curves[v1][v2], region->num_adj_curves[v3][v2],
                       region->num_adj_curves[v1][v3], region->num_adj_curves[v2][v3],
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
            printf("Intersect Tet Index %d, Intersect Tet Vertex %d\n", boundary->intersect_tet_index, boundary->intersect_tet_vertex);
            printf("    Meridian\n");

            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tet_index,
                       tri->tet_vertex,
                       tri->tet->curve[M][right_handed][tri->tet_vertex][0],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][1],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][2],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][3]
                );
            }
            printf("    Longitude\n");
            for (tri = boundary->cusp_triangle_begin.next; tri != &boundary->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tet_index,
                       tri->tet_vertex,
                       tri->tet->curve[L][right_handed][tri->tet_vertex][0],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][1],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][2],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][3]
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
                v1 = edgesThreeToFour[tri->tet_vertex][0];
                v2 = edgesThreeToFour[tri->tet_vertex][1];
                v3 = edgesThreeToFour[tri->tet_vertex][2];

                printf("    (Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
                       "Vertex %d: (%d %d), Vertex %d: (%d %d)\n",
                       tri->tet_index, tri->tet_vertex,
                       v1, tri->vertices[v1].edge_class, tri->vertices[v1].edge_index,
                       v2, tri->vertices[v2].edge_class, tri->vertices[v2].edge_index,
                       v3, tri->vertices[v3].edge_class, tri->vertices[v3].edge_index
                );
            }
        }

    } else if (flag == 5) {
        // Dual Curve Paths
        printf("Oscillating curve paths\n");

        // which dual curve
        for (i = 0; i < curves->num_curves; i++) {
            j = 0;

            printf("Dual Curve %d\n", i);
            // which curve component
            for (path = curves->dual_curve_begin[i].next; path != &curves->dual_curve_end[i]; path = path->next) {
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
                       tri->tet_index,               // Tet Index
                       tri->tet_vertex,                // Tet Vertex
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][1]][edgesThreeToFour[tri->tet_vertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][1]]
                );
            }
        }
    } else if (flag == 7) {
        printf("Graph info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            boundary = cusps[i];

            printf("Boundary %d\n", i);
            g = boundary->dual_graph;
            for (j = 0; j < g->num_vertices; j++) {
                if (g->regions[j] == NULL)
                        continue;

                printf("    Vertex %d (Tet Index: %d, Tet Vertex: %d): ",
                       j,
                       g->regions[j]->tet_index,
                       g->regions[j]->tet_vertex
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
        for (i = 0; i < curves->num_curves; i++) {
            printf("Dual Curve %d\n", i);

            j = 0;
            // which component
            for (path = curves->dual_curve_begin[i].next; path != &curves->dual_curve_end[i]; path = path->next) {
                printf("    Part %d Cusp %d\n", j, path->endpoints[0].region->tri->tet->cusp[path->endpoints[0].region->tet_vertex]->index);
                for (k = 0; k < 2; k++) {
                    if (k == 0)
                        printf("        Start: ");
                    else
                        printf("        End:   ");

                    x_vertex1 = (int) remaining_face[path->endpoints[k].region->tet_vertex][path->endpoints[k].vertex];
                    x_vertex2 = (int) remaining_face[path->endpoints[k].vertex][path->endpoints[k].region->tet_vertex];

                    printf("Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Num Curves (%d, %d) Edge Class (%d, %d)\n",
                           path->endpoints[k].region_index, path->endpoints[k].region->tet_index,
                           path->endpoints[k].region->tet_vertex, path->endpoints[k].face, path->endpoints[k].vertex,
                           path->endpoints[k].num_adj_curves[x_vertex1][path->endpoints[k].vertex],
                           path->endpoints[k].num_adj_curves[x_vertex2][path->endpoints[k].vertex],
                           path->endpoints[k].region->tri->vertices[path->endpoints[k].vertex].edge_class,
                           path->endpoints[k].region->tri->vertices[path->endpoints[k].vertex].edge_index);
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
    mflow = FLOW(tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    // Contribution from longitudinal curves
    lflow = FLOW(tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    retval = ABS(mflow) + ABS(lflow);
    return retval;
}

/*
 * Returns a pointer to the triangle with a given tetIndex and tet_vertex.
 */

CuspTriangle *find_cusp_triangle(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end,
        CuspTriangle *tri, int face) {
    int tet_index, tet_vertex;
    CuspTriangle *pTri;

    tet_index = tri->tet->neighbor[face]->index;
    tet_vertex = EVALUATE(tri->tet->gluing[face], tri->tet_vertex);

    for (pTri = cusp_triangle_begin->next; pTri != cusp_triangle_end; pTri = pTri->next) {
        if (pTri->tet_index == tet_index && pTri->tet_vertex == tet_vertex)
            return pTri;
    }

    // Didn't find a neighbour
    return NULL;
}

void do_oscillating_curves(ManifoldBoundary **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph) {
    int i;

    for (i = 0; i < curves->num_curves; i++) {
        do_one_dual_curve(cusps, curves, &curves->dual_curve_begin[i], &curves->dual_curve_end[i],
                          multi_graph, curves->edge_class[i], i);

        print_debug_info(cusps[0]->manifold, cusps, curves, 2);
        print_debug_info(cusps[0]->manifold, cusps, curves, 7);
        print_debug_info(cusps[0]->manifold, cusps, curves, 5);
        print_debug_info(cusps[0]->manifold, cusps, curves, 8);
    }
}

void do_one_dual_curve(ManifoldBoundary **cusps, OscillatingCurves *curves, DualCurves *dual_curve_begin,
                       DualCurves *dual_curve_end, EndMultiGraph *multi_graph, int edge_class, int curve_index) {
    int i, path_length;
    CuspEndPoint *cusp_end_point;
    DualCurves *path;
    int orientation = START;

    // find paths through cusps
    cusp_end_point = find_multi_graph_path(multi_graph->multi_graph, cusps[0]->manifold, edge_class, multi_graph->e0, &path_length);

    // find paths inside each cusp
    dual_curve_begin->edge_class[FINISH] = edge_class;
    dual_curve_end->edge_class[START]    = edge_class;

    for (i = 0; i < path_length; i++) {
//        print_debug_info(cusps[0]->manifold, cusps, NULL, 2);
//        print_debug_info(cusps[0]->manifold, cusps, NULL, 7);

        path = init_dual_curve(dual_curve_end, cusp_end_point[i].edge_class1, cusp_end_point[i].edge_class2);
        path->cusp_index = cusp_end_point[i].cusp_index;
        find_path_endpoints(cusps[path->cusp_index]->dual_graph, dual_curve_begin, path, cusp_end_point[i].pos, orientation);
        do_one_cusp(cusps[path->cusp_index], path, edge_class);
        orientation = (orientation == START ? FINISH : START);
    }

    update_adj_curve_along_path(cusps, curves, dual_curve_begin, dual_curve_end, curve_index);
    my_free(cusp_end_point);
}

/*
 * Construct oscillating curves on the boundary components
 */

void do_one_cusp(ManifoldBoundary *boundary, DualCurves *path, int edge_class) {
    int *parent;
    bool *processed, *discovered;

    processed = NEW_ARRAY(boundary->dual_graph->num_vertices, bool);
    discovered = NEW_ARRAY(boundary->dual_graph->num_vertices, bool);
    parent = NEW_ARRAY(boundary->dual_graph->num_vertices, int);

    // Find path using bfs
    init_search(boundary->dual_graph, processed, discovered, parent);
    bfs(boundary->dual_graph, path->endpoints[START].region_index, processed, discovered, parent);
    find_path(path->endpoints[START].region_index, path->endpoints[FINISH].region_index,
              parent, &path->curves_begin);
    update_path_info(boundary->dual_graph, path, edge_class);

    // Reallocate memory
    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // Split the regions along the path
    split_cusp_regions_along_path(boundary, path);

    free_graph(boundary->dual_graph);
    boundary->dual_graph = construct_cusp_region_dual_graph(boundary);
}

DualCurves *init_dual_curve(DualCurves *curve_end, int edge_class_start, int edge_class_finish) {
    DualCurves *path = NEW_STRUCT(DualCurves );
    INSERT_BEFORE(path, curve_end);

    path->curves_begin.next = &path->curves_end;
    path->curves_begin.prev = NULL;
    path->curves_end.next = NULL;
    path->curves_end.prev = &path->curves_begin;

    path->edge_class[START] = edge_class_start;
    path->edge_class[FINISH] = edge_class_finish;

    return path;
}

/*
 * Find the indicies of the cusp triangles which dive through the manifold
 * along the given edgeclass. If copy is true, find the path end point
 * corresponding to path1 and store it in path2. Else find a valid path end
 * point and store in path2
 */

void find_path_endpoints(Graph *g, DualCurves *path_start, DualCurves *path, int pos, int orientation) {
    int endpoint1_edge_class, endpoint1_edge_index, endpoint2_edge_class, endpoint2_edge_index;
    PathEndPoint *endpoint1_path1, *endpoint1_path2, *endpoint2_path1, *endpoint2_path2;
    if (pos == FIRST) {
        endpoint1_path2 = &path->endpoints[START];
        endpoint1_edge_class = path->edge_class[START];
        endpoint1_edge_index = START;
        find_single_endpoint(g, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path2 = &path->endpoints[FINISH];
        endpoint2_edge_class = path->edge_class[FINISH];
        endpoint2_edge_index = START;
        find_single_endpoint(g, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else if (pos == MIDDLE && orientation == START) {
        endpoint1_path1 = &path->prev->endpoints[START];
        endpoint1_path2 = &path->endpoints[START];
        endpoint1_edge_class = path->edge_class[START];
        endpoint1_edge_index = FINISH;
        find_single_matching_endpoint(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path2 = &path->endpoints[FINISH];
        endpoint2_edge_class = path->edge_class[FINISH];
        endpoint2_edge_index = START;
        find_single_endpoint(g, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else if (pos == MIDDLE && orientation == FINISH) {
        endpoint1_path1 = &path->prev->endpoints[FINISH];
        endpoint1_path2 = &path->endpoints[FINISH];
        endpoint1_edge_class = path->edge_class[START];
        endpoint1_edge_index = FINISH;
        find_single_matching_endpoint(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path2 = &path->endpoints[START];
        endpoint2_edge_class = path->edge_class[FINISH];
        endpoint2_edge_index = START;
        find_single_endpoint(g, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else if (pos == LAST) {
        endpoint1_path1 = &path_start->next->endpoints[START];
        endpoint1_path2 = &path->endpoints[START];
        endpoint1_edge_class = path_start->next->edge_class[START];
        endpoint1_edge_index = FINISH;
        find_single_matching_endpoint(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path1 = &path->prev->endpoints[FINISH];
        endpoint2_path2 = &path->endpoints[FINISH];
        endpoint2_edge_class = path->prev->edge_class[FINISH];
        endpoint2_edge_index = FINISH;
        find_single_matching_endpoint(g, endpoint2_path1, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else {
        uFatalError("do_one_cusp", "symplectic_basis");
        return;
    }
}

PathEndPoint *find_single_endpoint(Graph *g, PathEndPoint *path_endpoint, int edge_class, int edge_index) {
    int i, j, k, vertex, face1, face2, face;
    bool edge_lies_in_one_cusp;
    CuspRegion *pRegion;

    // which cusp region
    for (i = 0; i < g->num_vertices; i++) {
        if (g->regions[i] == NULL)
            continue;

        pRegion = g->regions[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == pRegion->tet_vertex)
                continue;

            if (pRegion->tri->vertices[vertex].edge_class != edge_class)
                continue;

            edge_lies_in_one_cusp = (pRegion->tri->tet->cusp[pRegion->tet_vertex]->index == pRegion->tri->tet->cusp[vertex]->index);

            if (edge_lies_in_one_cusp && pRegion->tri->vertices[vertex].edge_index != edge_index) {
                continue;
            }

            face1 = (int) remaining_face[pRegion->tet_vertex][vertex];
            face2 = (int) remaining_face[vertex][pRegion->tet_vertex];

            if (pRegion->dive[face1][vertex]) {
                face = face1;
            } else if (pRegion->dive[face2][vertex]) {
                face = face2;
            } else
                continue;

            path_endpoint->region = pRegion;
            path_endpoint->vertex = vertex;
            path_endpoint->face = face;
            path_endpoint->region_index = i;

            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    path_endpoint->num_adj_curves[j][k] = pRegion->num_adj_curves[j][k];
            return path_endpoint;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_single_endpoints", "symplectic_basis");
    return NULL;
}

PathEndPoint *find_single_matching_endpoint(Graph *g, PathEndPoint *path_endpoint1, PathEndPoint *path_endpoint2, int edge_class, int edge_index) {
    int i, j, k, vertex;
    bool edge_lies_in_one_cusp, region_index, region_vertex, dive_vertex, region_dive, region_curve;
    CuspRegion *pRegion;

    // which cusp region
    for (i = 0; i < g->num_vertices; i++) {
        if (g->regions[i] == NULL)
            continue;

        pRegion = g->regions[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == pRegion->tet_vertex)
                continue;

            if (pRegion->tri->vertices[vertex].edge_class != edge_class)
                continue;

            edge_lies_in_one_cusp = (pRegion->tri->tet->cusp[pRegion->tet_vertex]->index == pRegion->tri->tet->cusp[vertex]->index);

            if (edge_lies_in_one_cusp && pRegion->tri->vertices[vertex].edge_index != edge_index) {
                continue;
            }

            // are we in the correct region for copy
            region_index = (pRegion->tet_index != path_endpoint1->region->tet_index);
            region_vertex = (pRegion->tet_vertex != path_endpoint1->vertex);
            dive_vertex = (vertex != path_endpoint1->region->tet_vertex);
            region_dive = !pRegion->dive[path_endpoint1->face][vertex];
            region_curve = (pRegion->num_adj_curves[path_endpoint1->face][vertex] != path_endpoint1->num_adj_curves[path_endpoint1->face][path_endpoint1->vertex]);

            if (region_index || region_vertex || dive_vertex || region_dive || region_curve)
                continue;

            path_endpoint2->region = pRegion;
            path_endpoint2->vertex = vertex;
            path_endpoint2->face = path_endpoint1->face;
            path_endpoint2->region_index = i;

            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    path_endpoint2->num_adj_curves[j][k] = pRegion->num_adj_curves[j][k];

            return path_endpoint2;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_single_endpoints", "symplectic_basis");
    return NULL;
}

/*
 * After finding a path, each node contains the index of the
 * region it lies in. Update path info calculates the face
 * the path crosses to get to the next node and the vertex
 * it cuts off to simplify combinatorial holonomy calculation.
 */

void update_path_info(Graph *g, DualCurves *path, int edge_class) {
    int face;
    EdgeNode *node = &path->curves_begin;

    // path len 0
    if (node->next->next == NULL)
        return;

    node = node->next;
    // path len 1
    if (node->next->next == NULL) {
        for (face = 0; face < 4; face++)
            if (g->regions[node->y]->tet_vertex != face &&
                path->endpoints[START].vertex != face &&
                path->endpoints[FINISH].vertex != face)
                break;

        g->regions[node->y]->tri->tet->extra[edge_class].curve[g->regions[node->y]->tet_vertex]
        [path->endpoints[FINISH].face]++;

        g->regions[node->y]->tri->tet->extra[edge_class].curve[g->regions[node->y]->tet_vertex]
        [path->endpoints[START].face]--;

        node->inside_vertex = face;
        return;
    }


    // Set Header node
    update_path_endpoint_info(g->regions[node->y], node, &path->endpoints[START], edge_class, START);

    for (node = node->next; node->next->next != NULL; node = node->next) {
        find_interior_vertex(g->regions[node->y], node);

        g->regions[node->y]->tri->tet->extra[edge_class].curve[g->regions[node->y]->tet_vertex][node->next_face]++;
        g->regions[node->y]->tri->tet->extra[edge_class].curve[g->regions[node->y]->tet_vertex][node->prev_face]--;
    }

    // Set Tail node
    update_path_endpoint_info(g->regions[node->y], node, &path->endpoints[FINISH], edge_class, FINISH);
}

void update_path_endpoint_info(CuspRegion *pRegion, EdgeNode *node, PathEndPoint *path_endpoint,
                               int edge_class, int pos) {
    int i, vertex1, vertex2;

    EdgeNode *next_node;

    if (pos == START)
        next_node = node->next;
    else
        next_node = node->prev;

    node->next_face = -1;
    node->prev_face = -1;

    for (i = 0; i < 4; i++) {
        if (i == pRegion->tet_vertex || !pRegion->adj_cusp_triangle[i])
            continue;

        if (pRegion->adj_cusp_regions[i]->index == next_node->y) {
            node->prev_face = i;
            node->next_face = i;
        }
    }

    // next node isn't in an adjacent region
    if (node->next_face == -1 || node->prev_face == -1)
        uFatalError("update_path_endpoint_info", "symplectic_basis");

    if (pos == START)
        node->prev_face = path_endpoint->face;
    else
        node->next_face = path_endpoint->face;

    vertex1 = (int) remaining_face[pRegion->tet_vertex][path_endpoint->vertex];
    vertex2 = (int) remaining_face[path_endpoint->vertex][pRegion->tet_vertex];

    if (node->prev_face == path_endpoint->vertex) {
        node->inside_vertex = path_endpoint->face == vertex1 ? vertex2 : vertex1;
    } else if (node->prev_face == path_endpoint->face) {
        node->inside_vertex = -1;
    } else {
        node->inside_vertex = path_endpoint->vertex;
    }

    if (node->next_face == node->prev_face)
        return;

    pRegion->tri->tet->extra[edge_class].curve[pRegion->tet_vertex][node->next_face]++;
    pRegion->tri->tet->extra[edge_class].curve[pRegion->tet_vertex][node->prev_face]--;

}

/*
 * The oscillating curve splits the region it passes through
 * into two regions. Split each region in two and update
 * attributes
 */

void split_cusp_regions_along_path(ManifoldBoundary *boundary, DualCurves *path) {
    int face, index = boundary->num_cusp_regions;
    EdgeNode *node = path->curves_begin.next;
    CuspRegion *region;
    Graph *g = boundary->dual_graph;

    // empty path
    if (node->next == NULL)
        return ;

    // path of len 1
    if (node->next->next == NULL) {
        region = NEW_STRUCT(CuspRegion);
        INSERT_BEFORE(region, &boundary->cusp_region_end)
        copy_region(g->regions[node->y], region);

        face = node->inside_vertex;

        region->index = index;
        region->adj_cusp_triangle[path->endpoints[START].vertex]  = 0;
        region->adj_cusp_triangle[path->endpoints[FINISH].vertex] = 0;
        region->temp_adj_curves[path->endpoints[START].vertex][path->endpoints[FINISH].vertex]++;
        region->temp_adj_curves[path->endpoints[FINISH].vertex][path->endpoints[START].vertex]++;
        region->dive[path->endpoints[START].vertex][path->endpoints[FINISH].vertex] = (face != path->endpoints[FINISH].face);
        region->dive[path->endpoints[FINISH].vertex][path->endpoints[START].vertex] = (face != path->endpoints[START].face);

        g->regions[node->y]->adj_cusp_triangle[face] = 0;
        g->regions[node->y]->temp_adj_curves[face][path->endpoints[START].vertex]++;
        g->regions[node->y]->temp_adj_curves[face][path->endpoints[FINISH].vertex]++;
        g->regions[node->y]->dive[face][path->endpoints[START].vertex]  = (face == path->endpoints[START].face);
        g->regions[node->y]->dive[face][path->endpoints[FINISH].vertex] = (face == path->endpoints[FINISH].face);

        update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
        boundary->num_cusp_regions++;
        return;
    }

    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and region to the left of the curve.
     */
    update_cusp_triangle_endpoints(&boundary->cusp_region_begin, &boundary->cusp_region_end,
                                   g->regions[node->y], &path->endpoints[START], node, START);
    region = split_cusp_region_path_endpoint(g->regions[node->y], node, &path->endpoints[START], index, START);
    INSERT_BEFORE(region, &boundary->cusp_region_end);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        update_cusp_triangle_path_interior(&boundary->cusp_region_begin, &boundary->cusp_region_end,
                                           g->regions[node->y], node);
        region = split_cusp_region_path_interior(g->regions[node->y], node, &path->endpoints[START], index);
        INSERT_BEFORE(region, &boundary->cusp_region_end);
        index++;
    }

    // update last region
    update_cusp_triangle_endpoints(&boundary->cusp_region_begin, &boundary->cusp_region_end,
                                   g->regions[node->y], &path->endpoints[FINISH], node, FINISH);
    region = split_cusp_region_path_endpoint(g->regions[node->y], node, &path->endpoints[FINISH], index, FINISH);
    INSERT_BEFORE(region, &boundary->cusp_region_end);
    index++;

    update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
    boundary->num_cusp_regions = index;
}

/*
 * Set the new and old region data. Draw a picture to see how
 * the attributes change in each case
 *   - flag = -1 : starting node
 *   - flag = 0  : interior node
 *   - flag = 1  : end node
 */

CuspRegion *split_cusp_region_path_interior(CuspRegion *region, EdgeNode *node, PathEndPoint *path_endpoint, int index) {
    int v1, v2;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    v1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    v2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;

    // Update new region
    new_region->curve[v1][node->inside_vertex]++;
    new_region->curve[v2][node->inside_vertex]++;
    new_region->dive[node->inside_vertex][v1]    = region->dive[node->inside_vertex][v1];
    new_region->dive[v2][v1]               = region->dive[v2][v1];
    new_region->dive[node->inside_vertex][v2]    = region->dive[node->inside_vertex][v2];
    new_region->dive[v1][v2]               = region->dive[v1][v2];

    // Update region
    region->curve[v2][v1]++;
    region->curve[v1][v2]++;
    region->dive[v2][v1]                  = 0;
    region->dive[node->inside_vertex][v1]       = 0;
    region->dive[v1][v2]                  = 0;
    region->dive[node->inside_vertex][v2]       = 0;
    region->adj_cusp_triangle[node->inside_vertex]              = 0;

    return new_region;
}

CuspRegion *split_cusp_region_path_endpoint(CuspRegion *region, EdgeNode *node, PathEndPoint *path_endpoint, int index, int pos) {
    int face, vertex1, vertex2;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    vertex1 = (int) remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = (int) remaining_face[path_endpoint->vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;

    if (pos == START) {
        face = node->next_face;
    } else {
        face = node->prev_face;
    }

    if (face == path_endpoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        new_region->curve[path_endpoint->vertex][vertex2]++;
        new_region->temp_adj_curves[vertex1][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]      = (path_endpoint->face == vertex1);
        new_region->dive[vertex2][path_endpoint->vertex]      = region->dive[vertex2][path_endpoint->vertex];
        new_region->dive[vertex2][vertex1]               = region->dive[vertex2][vertex1];
        new_region->dive[path_endpoint->vertex][vertex1]      = region->dive[path_endpoint->vertex][vertex1];
        new_region->adj_cusp_triangle[vertex1]                      = 0;

        region->curve[path_endpoint->vertex][vertex1]++;
        region->temp_adj_curves[vertex2][path_endpoint->vertex]++;
        region->dive[vertex2][path_endpoint->vertex]         = (path_endpoint->face == vertex2);
        region->dive[vertex2][vertex1]                  = 0;
        region->dive[path_endpoint->vertex][vertex1]         = 0;
        region->adj_cusp_triangle[vertex2]                         = 0;
    } else if (face == path_endpoint->face) {
        // curve passes through the face that carries it
        new_region->curve[path_endpoint->face][path_endpoint->face == vertex1 ? vertex2 : vertex1]++;
        new_region->temp_adj_curves[face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex]++;
        new_region->dive[path_endpoint->face][path_endpoint->vertex]                   = region->dive[path_endpoint->face][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->vertex]                                 = 0;
        new_region->adj_cusp_triangle[path_endpoint->face == vertex1 ? vertex2 : vertex1]    = 0;

        region->curve[path_endpoint->face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face][path_endpoint->vertex]                      = 0;
    } else {
        // Curve goes around the vertex
        new_region->curve[face][path_endpoint->face]++;
        new_region->temp_adj_curves[path_endpoint->face][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]              = region->dive[vertex1][path_endpoint->vertex];
        new_region->dive[vertex2][path_endpoint->vertex]              = region->dive[vertex2][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->face]                       = 0;
        new_region->adj_cusp_triangle[path_endpoint->vertex]                     = 0;

        region->curve[face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex] = 0;
    }

    return new_region;
}

/*
 * region1 splits into region1 and region2, set them up
 * to be split
 */

void copy_region(CuspRegion *region1, CuspRegion *region2) {
    int i, j;

    if (region1 == NULL || region2 == NULL || region1->tri == NULL)
        uFatalError("copy_region", "symplectic_basis");

    region2->tri                            = region1->tri;
    region2->tet_index                       = region1->tet_index;
    region2->tet_vertex                      = region1->tet_vertex;

    for (i = 0; i < 4; i++) {
        region2->adj_cusp_triangle[i]                  = region1->adj_cusp_triangle[i];
        region2->adj_cusp_regions[i]              = NULL;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]            = region1->curve[i][j];
            region2->dive[i][j]             = 0;
            region2->num_adj_curves[i][j]         = region1->num_adj_curves[i][j];
            region2->temp_adj_curves[i][j]     = region1->temp_adj_curves[i][j];
        }
    }
}

/*
 * After splitting each region the path travels through,
 * the attributes for other regions in the same cusp
 * triangle is now out of date. Update cusp triangles
 * for nodes in the interior of the path.
 */

void update_cusp_triangle_path_interior(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end,
                                        CuspRegion *region, EdgeNode *node) {
    int face1, face2;
    CuspRegion *pRegion;

    face1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    face2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    for (pRegion = cusp_region_start->next; pRegion != cusp_region_end; pRegion = pRegion->next) {
        // is the region initialised?
        if (pRegion == NULL || pRegion->tet_index == -1)
            continue;

        // which triangle are we in?
        if (pRegion->tet_index != region->tet_index || pRegion->tet_vertex != region->tet_vertex)
            continue;

        if (pRegion->curve[face1][node->inside_vertex] > region->curve[face1][node->inside_vertex]) {
            pRegion->curve[face1][node->inside_vertex]++;
            pRegion->dive[face1][node->inside_vertex] = 0;
        }
        else if (pRegion->curve[face1][node->inside_vertex] < region->curve[face1][node->inside_vertex]) {
            pRegion->curve[face1][face2]++;
        }

        if (pRegion->curve[face2][node->inside_vertex] > region->curve[face2][node->inside_vertex]) {
            pRegion->curve[face2][node->inside_vertex]++;
            pRegion->dive[face2][node->inside_vertex] = 0;
        }
        else if (pRegion->curve[face2][node->inside_vertex] < region->curve[face2][node->inside_vertex]) {
            pRegion->curve[face2][face1]++;
        }
    }
}

/*
 * After splitting each curveRegion the path travels through,
 * the attributes for other regions in the same cusp
 * triangle is now out of date. Update cusp triangles
 * for nodes at the end of the path.
 */

void update_cusp_triangle_endpoints(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end, CuspRegion *region,
                                    PathEndPoint *path_endpoint, EdgeNode *node, int pos) {
    int face, face1, face2;
    CuspRegion *pRegion;

    face1 = (int) remaining_face[region->tet_vertex][path_endpoint->vertex];
    face2 = (int) remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (pos == START) {
        face = node->next_face;
    } else {
        face = node->prev_face;
    }

    for (pRegion = cusp_region_start->next; pRegion != cusp_region_end; pRegion = pRegion->next) {
        if (pRegion == NULL || pRegion->tet_index == -1)
            continue;

        // which triangle are we in?
        if (pRegion->tet_index != region->tet_index || pRegion->tet_vertex != region->tet_vertex)
            continue;

        if (face == path_endpoint->vertex) {
            // curve passes through the face opposite the vertex it dives through
            if (pRegion->curve[path_endpoint->vertex][face1] > region->curve[path_endpoint->vertex][face1]) {
                pRegion->curve[face][face1]++;
                pRegion->dive[face][face1] = 0;
                pRegion->temp_adj_curves[face2][path_endpoint->vertex]++;

            } else if (pRegion->curve[path_endpoint->vertex][face1] < region->curve[path_endpoint->vertex][face1]) {
                pRegion->curve[face][face2]++;
                pRegion->dive[face][face2] = 0;
                pRegion->temp_adj_curves[face1][path_endpoint->vertex]++;

            }

            continue;
        }

        // Curve goes around the vertex or passes through the face that carries it
        if (pRegion->curve[face][path_endpoint->vertex] > region->curve[face][path_endpoint->vertex]) {
            pRegion->curve[face][path_endpoint->vertex]++;
            pRegion->temp_adj_curves[face][path_endpoint->vertex]++;
            pRegion->dive[face][path_endpoint->vertex] = 0;

        } else if (pRegion->curve[face][path_endpoint->vertex] < region->curve[face][path_endpoint->vertex]) {
            pRegion->curve[face][face == face1 ? face2 : face1]++;
            pRegion->temp_adj_curves[face == face1 ? face2 : face1][path_endpoint->vertex]++;
            pRegion->dive[face][face == face1 ? face2 : face1] = 0;
        }
    }
}

void update_adj_curve_along_path(ManifoldBoundary **cusps, OscillatingCurves *curves, DualCurves *dual_curve_begin, DualCurves *dual_curve_end, int curve_index) {
    int i, j;
    DualCurves *curve;
    PathEndPoint *path_endpoint;
    CuspRegion *pRegion;
    ManifoldBoundary *cusp;

    // Update regions curve data
    for (curve = dual_curve_begin->next; curve != dual_curve_end; curve = curve->next) {
        // which cusp
        cusp = cusps[curve->cusp_index];

        for (pRegion = cusp->cusp_region_begin.next; pRegion != &cusp->cusp_region_end; pRegion = pRegion->next) {
            // which cusp region
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    pRegion->num_adj_curves[i][j] += pRegion->temp_adj_curves[i][j];
                    pRegion->temp_adj_curves[i][j] = 0;
                }
            }
        }
    }

    // update endpoint curve data
    for (i = 0; i < curve_index; i++) {
        // which oscillating curve

        for (curve = curves->dual_curve_begin[i].next; curve != &curves->dual_curve_end[i]; curve = curve->next) {
            // which component of the curve

            for (j = 0; j < 2; j++) {
                // which end point
                path_endpoint = &curve->endpoints[j];

                update_adj_curve_at_endpoint(path_endpoint, dual_curve_begin, dual_curve_end);
            }
        }
    }
}

void update_adj_curve_at_endpoint(PathEndPoint *path_endpoint, DualCurves *dual_curve_begin, DualCurves *dual_curve_end) {
    int i, face1, face2;
    DualCurves *path;
    PathEndPoint *curve_end_point;

    for (path = dual_curve_begin->next; path != dual_curve_end; path = path->next) {
        for (i = 0; i < 2; i++) {
            curve_end_point = &path->endpoints[i];

            // Cusp Triangle
            if (curve_end_point->region->tet_index != path_endpoint->region->tet_index ||
                curve_end_point->region->tet_vertex != path_endpoint->region->tet_vertex)
                continue;

            // Dive vertex
            if (curve_end_point->vertex != path_endpoint->vertex)
                continue;

            // update path data
            face1 = (int) remaining_face[path_endpoint->region->tet_vertex][path_endpoint->vertex];
            face2 = (int) remaining_face[path_endpoint->vertex][path_endpoint->region->tet_vertex];

            if (curve_end_point->face == face1 && path_endpoint->face == face2)
                path_endpoint->num_adj_curves[face1][path_endpoint->vertex]++;
            else if (curve_end_point->face == face2 && path_endpoint->face == face1)
                path_endpoint->num_adj_curves[face2][path_endpoint->vertex]++;
            else if (curve_end_point->face == face1 && path_endpoint->face == face1) {
                if (path_endpoint->num_adj_curves[face1][path_endpoint->vertex] >= curve_end_point->num_adj_curves[face1][path_endpoint->vertex])
                    path_endpoint->num_adj_curves[face1][path_endpoint->vertex]++;
                else
                    path_endpoint->num_adj_curves[face2][path_endpoint->vertex]++;
            } else if (curve_end_point->face == face2 && path_endpoint->face == face2) {
                if (path_endpoint->num_adj_curves[face2][path_endpoint->vertex] >= curve_end_point->num_adj_curves[face2][path_endpoint->vertex])
                    path_endpoint->num_adj_curves[face2][path_endpoint->vertex]++;
                else
                    path_endpoint->num_adj_curves[face1][path_endpoint->vertex]++;
            } else
                uFatalError("update_adj_curve_at_endpoint", "symplectic_basis");
        }
    }
}

// ----------------------------------

/*
 * Construct the symplectic equations from the dual curves
 */

void calculate_holonomy(Triangulation *manifold, int **symp_eqns, int num_curves) {
    int curve, v, f, ff;
    Tetrahedron *tet;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {

        // which curve
        for (curve = 0; curve < num_curves; curve++) {
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
 * node lies in region region, find the vertex which the subpath
 * node->prev->y --> node->y --> node->next->y cuts off of the
 * cusp triangle region->tri.
 */

void find_interior_vertex(CuspRegion *region, EdgeNode *node) {
    int i, vertex1, vertex2;

    for (i = 0; i < 4; i++) {
        if (i == region->tet_vertex)
            continue;

        vertex1 = (int) remaining_face[region->tet_vertex][i];
        vertex2 = (int) remaining_face[i][region->tet_vertex];

        if (!region->adj_cusp_triangle[vertex1] || !region->adj_cusp_triangle[vertex2])
            continue;

        if (region->adj_cusp_regions[vertex1]->index == node->next->y
            && region->adj_cusp_regions[vertex2]->index == node->prev->y) {
            node->inside_vertex = i;
            node->next_face = vertex1;
            node->prev_face = vertex2;
            return;
        } else if (region->adj_cusp_regions[vertex2]->index == node->next->y && region->adj_cusp_regions[vertex1]->index == node->prev->y) {
            node->inside_vertex = i;
            node->next_face = vertex2;
            node->prev_face = vertex1;
            return;
        }
    }

    // where does the next node go?
    uFatalError("find_interior_vertex", "symplectic_basis");
}

// -------------------------------------------------------

// End Multi Graph

EndMultiGraph *init_end_multi_graph(Triangulation *manifold, int *edge_classes) {
    int i = 0;
    int *parent, *edges;
    EndMultiGraph *multi_graph = NEW_STRUCT(EndMultiGraph );

    Graph *g = init_graph(manifold->num_cusps, FALSE);
    cusp_graph(manifold, g);
    print_graph(g, 2);

    multi_graph->multi_graph = init_graph(g->num_vertices, g->directed);
    parent = NEW_ARRAY(multi_graph->multi_graph->num_vertices, int);
    spanning_tree(g, multi_graph->multi_graph, 0, parent);
    color_graph(multi_graph->multi_graph);
    edges = find_tree_edges(multi_graph->multi_graph, manifold->num_tetrahedra);
    multi_graph->e0 = find_same_color_edge(multi_graph->multi_graph, g, edges);

    for (i = 0; i < manifold->num_tetrahedra; i++) {
        edge_classes[i] = !edges[i];
    }
    edge_classes[multi_graph->e0] = 0;

    my_free(edges);
    free_graph(g);
    my_free(parent);
    return multi_graph;
}

void free_end_multi_graph(EndMultiGraph *multiGraph) {
    free_graph(multiGraph->multi_graph);

    my_free(multiGraph);
}

int insert_edge_end_multi_graph(Graph *g, int x, int y, int edge_class, bool directed) {
    // Ignore edge if it already exists
    if (edge_exists(g, x, y))
        return x;

    EdgeNode *p = NEW_STRUCT( EdgeNode);
    INSERT_AFTER(p, &g->edge_list_begin[x]);
    p->y = y;
    p->edge_class = edge_class;
    g->degree[x]++;

    if (!directed) {
        insert_edge_end_multi_graph(g, y, x, edge_class, TRUE);
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
    int i, edge_class;
    EdgeNode *node;

    bool *processed = NEW_ARRAY(graph1->num_vertices, bool);
    bool *discovered = NEW_ARRAY(graph1->num_vertices, bool);

    // Find path using bfs
    init_search(graph1, processed, discovered, parent);
    bfs(graph1, start, processed, discovered, parent);

    for (i = 0; i < graph1->num_vertices; i++) {
        if (parent[i] == -1)
            continue;

        edge_class = -1;
        for (node = graph1->edge_list_begin[i].next; node != &graph1->edge_list_end[i]; node = node->next)
            if (node->y == parent[i]) {
                edge_class = node->edge_class;
                break;
            }

        if (edge_class == -1)
            continue;

        insert_edge_end_multi_graph(graph2, i, parent[i], edge_class, graph2->directed);
    }

    my_free(processed);
    my_free(discovered);
}

void color_graph(Graph *g) {
    int color = 0, v;
    Queue *q = init_queue(g->num_vertices);
    EdgeNode *node;

    g->color[0] = color;
    q = enqueue(q, 0);

    while (!empty_queue(q)) {
        v = dequeue(q);
        color = g->color[v];

        for (node = g->edge_list_begin[v].next; node != &g->edge_list_end[v]; node = node->next) {
            // graph is not bipartite
            if (g->color[node->y] == color)
                uFatalError("color_graph", "symplectic_basis");

            if (g->color[node->y] != -1)
                continue;

            g->color[node->y] = !color;
            q = enqueue(q, node->y);
        }
    }

    free_queue(q);
}

int *find_tree_edges(Graph *g, int num_edges) {
    int i, vertex;
    EdgeNode *node;
    int *edges = NEW_ARRAY(num_edges, int);

    for (i = 0; i < num_edges; i++) {
        edges[i] = 0;
    }

    // which vertex
    for (vertex = 0; vertex < g->num_vertices; vertex++) {
        // which edge
        for (node = g->edge_list_begin[vertex].next; node != &g->edge_list_end[vertex]; node = node->next) {
            edges[node->edge_class] = 1;
        }
    }

    return edges;
}

/*
 * g1 is the colored spanning tree of g2, return the
 * edge class of the edge in g2 which connects
 * vertices in g1 of the same color
 */

int find_same_color_edge(Graph *g1, Graph *g2, int *edge_classes) {
    int v;
    EdgeNode *node;

    for (v = 0; v < g1->num_vertices; v++) {
        for (node = g2->edge_list_begin[v].next; node != &g2->edge_list_end[v]; node = node->next) {
            if (g1->color[v] == g1->color[node->y] && !edge_classes[node->edge_class])
                // we found an edge
                return node->edge_class;
        }
    }

    // we didn't find an edge connecting vertices of the same color
    uFatalError("find_same_color_edge", "symplectic_basis");
    return -1;
}

/*
 * Find the length of a path between start and end
 */

int find_path_len(int start, int end, int *parents, int path_length) {
    if ((start == end) || (end == -1)) {
        return path_length;
    } else {
        return find_path_len(start, parents[end], parents, path_length + 1);
    }
}

CuspEndPoint *find_multi_graph_path(Graph *g, Triangulation *manifold, int edge_class, int e0, int *path_length) {
    bool *processed     = NEW_ARRAY(g->num_vertices, bool);
    bool *discovered    = NEW_ARRAY(g->num_vertices, bool);
    int *parent         = NEW_ARRAY(g->num_vertices, int);
    int start, end, startE0, endE0;
    EdgeNode *node_begin = NEW_STRUCT( EdgeNode ), *node_end = NEW_STRUCT( EdgeNode), *tempNode;
    CuspEndPoint *cusp_end_point;

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    find_edge_ends(g, manifold, edge_class, &start, &end);
    find_edge_ends(g, manifold, e0, &startE0, &endE0);

    init_search(g, processed, discovered, parent);
    bfs(g, start, processed, discovered, parent);

    *path_length = 0;
    *path_length = find_path_len(start, end, parent, *path_length);
    tempNode = NEW_STRUCT( EdgeNode );

    if (*path_length % 2 == 1) {
        find_path(start, end, parent, node_begin);
    } else {
        init_search(g, processed, discovered, parent);
        bfs(g, start, processed, discovered, parent);

        find_path(start, startE0, parent, node_begin);

        my_free(tempNode);
        tempNode = node_end->prev;

        init_search(g, processed, discovered, parent);
        bfs(g, endE0, processed, discovered, parent);

        find_path(endE0, end, parent, node_end->prev);
        tempNode->next->edge_class = e0;
    }

    cusp_end_point = find_cusp_endpoint_edge_classes(g, node_begin, node_end, edge_class, path_length);

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

    return cusp_end_point;
}

CuspEndPoint *find_cusp_endpoint_edge_classes(Graph *g, EdgeNode *node_begin,
                                              EdgeNode *node_end, int edge_classes, int *path_length) {
    int i, cusp;
    EdgeNode *edge, *node;

    *path_length = 0;
    for (node = node_begin->next; node != node_end; node = node->next)
        (*path_length)++;

    CuspEndPoint *pCuspEndPoint = NEW_ARRAY(*path_length, CuspEndPoint);

    pCuspEndPoint[0].edge_class1 = edge_classes;

    i = 0;
    for (node = node_begin->next; node->next != node_end; node = node->next) {
        cusp = node->y;

        pCuspEndPoint[i].cusp_index = node->y;
        pCuspEndPoint[i].pos = MIDDLE;

        if (i != 0)
            pCuspEndPoint[i].edge_class1 = pCuspEndPoint[i - 1].edge_class2;

        if (node->next->edge_class != -1) {
            pCuspEndPoint[i].edge_class2 = node->next->edge_class;
            i++;
            continue;
        }

        for (edge = g->edge_list_begin[cusp].next; edge != &g->edge_list_end[cusp]; edge = edge->next) {
            if (node->next->y == edge->y) {
                pCuspEndPoint[i].edge_class2 = edge->edge_class;
                break;
            }
        }

        i++;
    }

    pCuspEndPoint[0].pos = FIRST;
    pCuspEndPoint[*path_length - 1].cusp_index = node->y;
    pCuspEndPoint[*path_length - 1].edge_class1 = pCuspEndPoint[*path_length - 2].edge_class2;
    pCuspEndPoint[*path_length - 1].edge_class2 = edge_classes;
    pCuspEndPoint[*path_length - 1].pos = LAST;

    return pCuspEndPoint;
}

void find_edge_ends(Graph *g, Triangulation *manifold, int edge_class, int *start, int *end) {
    int v1, v2;
    Tetrahedron *tet;
    EdgeClass *edge;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (v1 = 0; v1 < 4; v1++) {
            // which vertex of the cusp triangle at v1
            for (v2 = 0; v2 < 4; v2++) {
                if (v1 == v2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[v1][v2]];
                if (edge->index != edge_class)
                    continue;

                *start = tet->cusp[v1]->index;
                *end   = tet->cusp[v2]->index;
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

    for (i = 0; i < g->num_vertices; i++) {
        printf("(Cusp Index: %d) ", i);
        edge = &g->edge_list_begin[i];
        while ((edge = edge->next)->next != NULL) {
            printf("%d ", edge->y);
        }
        printf("\n");
    }

    printf("--------------------------\n");
}
