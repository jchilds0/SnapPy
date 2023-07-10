/*
 *  Symplectic Basis
 *
 *  Computes the symplectic basis of a cusped 3-manifold
 */

#include <stdio.h>
#include "SnapPea.h"
#include "kernel.h"
#include "addl_code.h"

#define at_least_two(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))

#define FIRST                   0
#define MIDDLE                  1
#define LAST                    2
#define START                   0
#define FINISH                  1

/*
 * Queue
 */

typedef struct Queue {
    int                         front;      // First element of queue
    int                         rear;       // Last element of queue
    int                         len;        // num of elements
    int                         size;       // array size
    int                         *array;
} Queue ;

/**
 * Graph
 */

typedef struct EdgeNode {
    int                         y;                       /** cusp region index */
    struct EdgeNode             *next;                   /** next node in doubly linked list */
    struct EdgeNode             *prev;                   /** prev node in doubly linked list */
} EdgeNode;

typedef struct Graph {
    EdgeNode                    *edge_list_begin;        /** header node of doubly linked list */
    EdgeNode                    *edge_list_end;          /** tail node of doubly linked list */
    struct CuspRegion           **regions;               /** list of regions in the graph */
    int                         *degree;                 /** degree of each vertex */
    int                         *color;                  /** color a tree bipartite */
    int                         num_vertices;            /** number of vertices in the graph */
    Boolean                     directed;                /** is the graph directed */
} Graph;

typedef struct CuspEndPoint {
    int                         cusp_index;
    int                         edge_class[2];
    struct CuspEndPoint         *next;
} CuspEndPoint;

typedef struct EndMultiGraph {
    int                         e0;                      /** edge connecting vertices of the same color */
    int                         num_edge_classes;
    int                         num_cusps;
    int                         **edges;                 /** edge_class[u][v] is the edge class of the edge u->v */
    Boolean                     *edge_classes;           /** which edge classes are in the multigraph */
    Graph                       *multi_graph;            /** tree with extra edge of cusps */
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
    int                         curve[4][4];            /** oscillating curve holonomy for a cusp triangle */
};

typedef struct PathEndPoint {
    FaceIndex                   face;                   /** face containg the short rectangle carrying the curve */
    VertexIndex                 vertex;                 /** vertex we dive through the manifold along */
    int                         region_index;           /** index of the region the endpoint lies in */
    int                         num_adj_curves[4][4];   /** where the curve dives into the manifold */
    struct CuspRegion           *region;                /** pointer to the region the endpoint lies in */
} PathEndPoint;

typedef struct PathNode { 
    int                         cusp_region_index;
    FaceIndex                   next_face;               /** face the path crosses to the next node */
    FaceIndex                   prev_face;               /** face the path crosses to the prev node */
    VertexIndex                 inside_vertex;           /** inside vertex of the path */
    struct CuspTriangle         *tri;                    /** cusp triangle the node lies in */
    struct PathNode             *next;                   /** next node in doubly linked list */
    struct PathNode             *prev;
} PathNode;

typedef struct DualCurves {
    int                         edge_class[2];
    int                         cusp_index;             /** which cusp does the curve lie in */
    PathNode                    curves_begin;           /** header node of doubbly linked list */
    PathNode                    curves_end;             /** tailer node of doubbly linked list */
    PathEndPoint                endpoints[2];           /** path end points */
    struct DualCurves           *next;                  /** next dual curve in doubly linked list */
    struct DualCurves           *prev;                  /** prev dual curve in doubly linked list */
} DualCurves;

typedef struct OscillatingCurves {
    int                         num_curves;
    int                         *edge_class;
    DualCurves                  *dual_curve_begin;      /** array of doubly linked lists of dual curves */
    DualCurves                  *dual_curve_end;        /** array of doubly linkek lists of dual curves */
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
    int                         edge_class;
    int                         edge_index;
    EdgeClass                   *edge;
    VertexIndex                 v1;
    VertexIndex                 v2;
} CuspVertex;

typedef struct CuspTriangle {
    Tetrahedron                 *tet;                   /** tetrahedron the triangle comes from */
    Cusp                        *cusp;                  /** cusp the triangle lies in */
    int                         tet_index;              /** tet->index */
    VertexIndex                 tet_vertex;             /** vertex the triangle comes from */
    int                         num_curves;             /** number of curves on the triangle */
    CuspVertex                  vertices[4];            /** information about each vertex */
    struct CuspTriangle         *neighbours[4];         /** triangle neighbouring a face */
    struct CuspTriangle         *next;                  /** next cusp triangle on doubly linked list */
    struct CuspTriangle         *prev;                  /** prev cusp triangle on doubly linkled list */
} CuspTriangle;

typedef struct CuspRegion {
    CuspTriangle                *tri;                   /** cusp triangle the region lies on */
    int                         tet_index;              /** tri->tetIndex */
    VertexIndex                 tet_vertex;             /** tri->tet_vertex */
    int                         index;                  /** index of the cusp region */
    int                         curve[4][4];            /** looking at face, number of curves between the region and vertex */
    Boolean                     adj_cusp_triangle[4];   /** does the region meet this edge of the cusp triangle */
    Boolean                     dive[4][4];             /** can we dive along the face into this vertex */
    int                         num_adj_curves[4][4];   /** stores the number of curves between a region and a face */
    int                         temp_adj_curves[4][4];  /** store the adj curve until pathfinding is complete */
    struct CuspRegion           *adj_cusp_regions[4];   /** index of the adjacent regions */
    struct CuspRegion           *adj_dive_regions[4];   /** regions which are adjacent by diving through the manifold */
    struct CuspRegion           *next;                  /** next cusp region in doubly linked list */
    struct CuspRegion           *prev;                  /** prev cusp region in doubly linked list */
} CuspRegion;

typedef struct ManifoldBoundary {
    int                         intersect_tet_index;    /** index of the intersection triangle */
    VertexIndex                 intersect_tet_vertex;   /** vertex of the intersection triangle */
    int                         num_edge_classes;       /** number of edge classes in the boundary */
    int                         num_cusp_triangles;     /** number of cusp triangle in the boundary */
    int                         num_cusp_regions;       /** number of cusp regions in the boundary */
    Triangulation               *manifold;              /** manifold */
    Cusp                        *cusp;                  /** which cusp is the boundary in */
    Graph                       *dual_graph;            /** dual graph of the cusp region */
    CuspTriangle                cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    CuspTriangle                cusp_triangle_end;      /** tail node of doubly linked list of cusp triangles */
    CuspRegion                  cusp_region_begin;      /** header node of doubly linked list of cusp regions */
    CuspRegion                  cusp_region_end;        /** tail node of doubly linked list of cusp regions */
    int                         train_line_len;
    int                         *train_line_path;       /** train line */
    PathEndPoint                **train_line_endpoint;
} ManifoldBoundary;

/**
 * Queue Data Structure
 */

Queue                   *init_queue(int);
Queue                   *enqueue(Queue *, int);
int                     dequeue(Queue *);
Queue                   *resize_queue(Queue *);
Boolean                 empty_queue(Queue *);
void                    free_queue(Queue *);

/**
 * Graph
 */

Graph                   *init_graph(int, Boolean);
void                    free_graph(Graph *);
int                     insert_edge(Graph *, int, int, Boolean);
void                    delete_edge(Graph *, int, int, Boolean);
Boolean                 edge_exists(Graph *, int, int);

/**
 * Breadth First Search
 */

void                    init_search(Graph *, Boolean *, Boolean *, int *);
void                    bfs(Graph *, int, Boolean *, Boolean *, int *);
void                    find_path(int, int, int *, EdgeNode *);
Graph                   **find_connected_graph_components(Graph *);
int                     *dfs_target_edges(Graph *, int *, int *, int, int *);
Boolean                 check_neighbours(Graph *, int, int *, int *, int, int *, int *, const Boolean *);
Boolean                 is_path_an_edge_cover(const int *, const int *, int, const int *, int);

/**
 * Symplectic Basis
 */

int                     *gluing_equations_for_edge_class(Triangulation *, int);
int                     **get_symplectic_equations(Triangulation *, Boolean *, int *, int);

ManifoldBoundary        *init_boundary(Triangulation *, Cusp *);
void                    free_boundary(ManifoldBoundary **, int);
void                    init_cusp_triangulation(Triangulation *, ManifoldBoundary *);
int                     net_flow_around_vertex(CuspTriangle *, int);
void                    label_triangulation_edges(Triangulation *);
CuspTriangle            *find_cusp_triangle(CuspTriangle *, CuspTriangle *, CuspTriangle *, int);
void                    label_cusp_vertex_indices(CuspTriangle *, CuspTriangle *, int);
void                    walk_around_cusp_vertex(CuspTriangle *, int, int);
void                    init_cusp_region(ManifoldBoundary *);
int                     init_intersect_cusp_region(ManifoldBoundary *, CuspTriangle *, int);
int                     init_intersect_vertex_two_zero_flows(ManifoldBoundary *, CuspTriangle *, int);
int                     init_normal_cusp_region(ManifoldBoundary *, CuspTriangle *, int);
void                    set_cusp_region_data(CuspRegion *, CuspTriangle *, const int [4], const Boolean [4], int);
void                    update_adj_region_data(CuspRegion *, CuspRegion *);
CuspRegion              *find_adj_region(CuspRegion *, CuspRegion *, CuspRegion *, int);
void                    init_train_line(ManifoldBoundary *);
DualCurves              *init_dual_curve(int, int);
OscillatingCurves       *init_oscillating_curves(Triangulation *, const Boolean *);
void                    free_oscillating_curves(OscillatingCurves *);
void                    find_intersection_triangle(Triangulation *, ManifoldBoundary *);

/**
 * Train lines
 */

void                    do_manifold_train_lines(ManifoldBoundary **, EndMultiGraph *);
void                    find_edge_class_edges(ManifoldBoundary **, EndMultiGraph *, int **, int **);
void                    find_edge_class_edges_on_cusp(ManifoldBoundary *, EndMultiGraph *, Graph **, Graph *, int *, int *);
void                    update_edge_classes_on_cusp(ManifoldBoundary **, int, Boolean **, int, int, int **, int **);
Graph                   *construct_edge_class_graph(CuspTriangle *, CuspTriangle *, const Boolean *, int);
int                     *find_primary_train_line(ManifoldBoundary *, EndMultiGraph *, int *, int *, int *);
void                    find_train_line_endpoints(ManifoldBoundary **, EndMultiGraph *, PathEndPoint ***, int **, int **);
void                    find_train_line_endpoints_on_cusp(ManifoldBoundary *, EndMultiGraph *, PathEndPoint **, int *, int *, const Boolean *);
void                    update_cusp_regions_along_train_line(ManifoldBoundary *);
void                    convert_path_array_to_path_node(CuspRegion *, CuspRegion *, PathNode *, PathNode *, const int *, int);
PathNode                *endpoint_path_point_to_path_node(CuspRegion *, PathEndPoint *, const int *, int);
PathNode                *interior_path_point_to_path_node(CuspRegion *, const int *, int);
CuspRegion              *get_region_from_index(CuspRegion *, CuspRegion *, int);
CuspRegion              *find_adj_dive_region(CuspRegion *, CuspRegion *, CuspRegion *, int, int);


/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    do_oscillating_curves(ManifoldBoundary **, OscillatingCurves *, EndMultiGraph *);
void                    do_one_dual_curve(ManifoldBoundary **, OscillatingCurves *, DualCurves *, DualCurves *, EndMultiGraph *, int, int);
void                    do_one_cusp(ManifoldBoundary *, DualCurves *);
Graph *                 construct_cusp_region_dual_graph(ManifoldBoundary *);
void                    print_debug_info(Triangulation *, ManifoldBoundary **, OscillatingCurves *, int);
void                    find_path_endpoints(Graph *, DualCurves *, DualCurves *, int, int);
Boolean                 is_valid_endpoint(CuspRegion *, int, int, int, int);
PathEndPoint            *find_single_endpoint(Graph *, PathEndPoint *, int, int);
PathEndPoint            *find_single_matching_endpoint(Graph *, PathEndPoint *, PathEndPoint *, int, int);
void                    graph_path_to_dual_curve(Graph *g, EdgeNode *, EdgeNode *, DualCurves *);
PathNode                *endpoint_edge_node_to_path_node(CuspRegion *, EdgeNode *, PathEndPoint *, int);
PathNode                *interior_edge_node_to_path_node(CuspRegion *, EdgeNode *);
void                    split_cusp_regions_along_path(ManifoldBoundary *, DualCurves *);
CuspRegion              *split_cusp_region_path_interior(CuspRegion *, PathNode *, int);
CuspRegion              *split_cusp_region_path_endpoint(CuspRegion *, PathNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle_path_interior(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *);
void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, PathNode *, int);
void                    copy_region(CuspRegion *, CuspRegion *);
void                    update_adj_curve_along_path(ManifoldBoundary **, OscillatingCurves *, DualCurves *, DualCurves *, int);
void                    update_adj_curve_at_endpoint(PathEndPoint *, DualCurves *, DualCurves *);
CuspRegion              **find_endpoint_region(ManifoldBoundary *, int, int, int *);
CuspRegion              *find_endpoint_adj_region(ManifoldBoundary *, CuspRegion *, int);
void                    update_adj_dive(ManifoldBoundary **, int);
void                    update_path_holonomy(DualCurves *, int);
void                    calculate_holonomy(Triangulation *, int **, int);


/**
 * End Multi Graph
 */

EndMultiGraph           *init_end_multi_graph(Triangulation *);
void                    free_end_multi_graph(EndMultiGraph *);
Graph                   *spanning_tree(Graph *, int, int *);
int                     **find_end_multi_graph_edge_classes(EndMultiGraph *, Triangulation *);
int                     find_edge_class(Triangulation *, int, int);
void                    cusp_graph(Triangulation *, Graph *);
void                    color_graph(Graph *);
int                     find_same_color_edge(Triangulation *, EndMultiGraph *, Graph *);
int                     find_path_len(int, int, int *, int);
CuspEndPoint            *find_multi_graph_path(EndMultiGraph *, Triangulation *, int, int *);
CuspEndPoint            *graph_path_to_cusp_path(EndMultiGraph *, EdgeNode *, EdgeNode *, int);
void                    find_edge_ends(Graph *, Triangulation *, int, int *, int *);
Boolean                 *edge_classes_on_cusp(EndMultiGraph *, int);
void                    print_graph(Graph *);

int edgesThreeToFour[4][3] = {{1, 2, 3},
                              {0, 2, 3},
                              {0, 1, 3},
                              {0, 1, 2}};

// -------------------------------------------------

// Data Structures

// Queue Data Structure

Queue *init_queue(int size) {
    Queue *q    = NEW_STRUCT( Queue );

    q->front    = 0;
    q->rear     = -1;
    q->len      = 0;
    q->size     = MIN(size, 256);
    q->array    = NEW_ARRAY(q->size, int);

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

Boolean empty_queue(Queue *q) {
    if (q->len > 0)
        return FALSE;

    return TRUE;
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

Graph *init_graph(int max_vertices, Boolean directed) {
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

int insert_edge(Graph *g, int x, int y, Boolean directed) {
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

void delete_edge(Graph *g, int vertex_x, int vertex_y, Boolean directed) {
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

Boolean edge_exists(Graph *g, int v1, int v2) {
    EdgeNode *node = &g->edge_list_begin[v1];

    while ((node = node->next)->next != NULL) {
        if (node->y == v2) {
            return TRUE;
        }
    }

    return FALSE;
}

// ---------------------------------------------------------------

// Breadth First Search

/*
 * Initialise default values for bfs arrays
 */

void init_search(Graph *g, Boolean *processed, Boolean *discovered, int *parent) {
    int i;

    for (i = 0; i < g->num_vertices; i ++) {
        processed[i] = FALSE;
        discovered[i] = FALSE;
        parent[i] = -1;
    }
}

/*
 * Graph search algorithm starting at vertex 'start'.
 */

void bfs(Graph *g, int start, Boolean *processed, Boolean *discovered, int *parent) {
    Queue *q = init_queue(10);
    int v, y;
    EdgeNode *p;

    enqueue(q, start);
    discovered[start] = TRUE;

    while (!empty_queue(q)) {
        v = dequeue(q);
        processed[v] = TRUE;
        p = &g->edge_list_begin[v];

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = TRUE;
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
    //new_node->edge_class = -1;

    INSERT_AFTER(new_node, node);

    if ((start == end) || (end == -1)) {
        new_node->y = start;
    } else {
        find_path(start, parents[end], parents, node);
        new_node->y = end;
    }
}

/* 
 * Given a graph g, find the spanning tree of each connected component
 */

Graph **find_connected_graph_components(Graph *g) {
    int i, j, index = 0;
    EdgeNode *edge;
    int *parent             = NEW_ARRAY(g->num_vertices, int);
    Boolean *edge_classes   = NEW_ARRAY(g->num_vertices, Boolean); 
    Boolean *processed      = NEW_ARRAY(g->num_vertices, Boolean);
    Boolean *discovered     = NEW_ARRAY(g->num_vertices, Boolean);
    Graph *temp_graph;
    Graph **graphs          = NEW_ARRAY(g->num_vertices, Graph *);

    for (i = 0; i < g->num_vertices; i++)
        graphs[i] = NULL;

    for (i = 0; i < g->num_vertices; i++) {
        if (edge_classes[i])
            continue;

        init_search(g, processed, discovered, parent);
        bfs(g, i, processed, discovered, parent);

        temp_graph = init_graph(g->num_vertices, g->directed);

        for (j = 0; j < temp_graph->num_vertices; j++) {
            edge_classes[j] = edge_classes[j] || discovered[j];

            if (!discovered[j])
                continue;

            for (edge = g->edge_list_begin[j].next; edge != &g->edge_list_end[j]; edge = edge->next) {
                insert_edge(temp_graph, j, edge->y, temp_graph->directed);
            }
        }

        graphs[index++] = spanning_tree(temp_graph, i, parent);
    }

    return graphs;
}

/*
 * Find a path which covers the set of edges 'edges' using depth first search
 */

int *dfs_target_edges(Graph *g, int *edges_start, int *edges_end, int num_edges, int *path_len) {
    int i;
    int *path = NEW_ARRAY(g->num_vertices + 1, int);
    Boolean found;
    Boolean *visited = NEW_ARRAY(g->num_vertices, Boolean);

    for (i = 0; i < g->num_vertices; i++) {
        visited[i] = FALSE;
    }

    *path_len = 0;
    found = FALSE;

    for (i = 0; i < num_edges; i++) {
        if (edges_start[i] != -1 || edges_end[i] != -1) {
            found = TRUE;
        }
    }

    if (!found) {
        my_free(visited);
        return NULL;
    }

    if (check_neighbours(g, edges_start[0], edges_start, edges_end, num_edges, path, path_len, visited)) {
        my_free(visited);
        return path;
    }

    uFatalError("dfs_target_edges", "symplectic_basis");
    return NULL;
}

Boolean check_neighbours(Graph *g, int start, int *edges_start, int *edges_end, int num_edges, 
                         int *path, int *path_len, const Boolean *visited) {
    int i;
    EdgeNode *node;
    Boolean result, *current_visited = NEW_ARRAY(g->num_vertices, Boolean);

    path[(*path_len)++] = start;
    for (i = 0; i < g->num_vertices; i++) {
        current_visited[i] = visited[i];
    }

    if (is_path_an_edge_cover(edges_start, edges_end, num_edges, path, *path_len)) {
        my_free(current_visited);
        return TRUE;
    }

    for (node = g->edge_list_begin[start].next; node != &g->edge_list_end[start]; node = node->next) {
        if (!current_visited[node->y] || node->y == path[0]) {
            result = check_neighbours(g, node->y, edges_start, edges_end, num_edges, path, path_len, current_visited);

            if (result == FALSE)
                continue;

            my_free(current_visited);
            return TRUE;
        }
    }

    my_free(current_visited);
    (*path_len)--;
    return FALSE;
}

Boolean is_path_an_edge_cover(const int *edges_begin, const int *edges_end, int num_edges, const int *path, int path_len) {
    int i, j;
    Boolean edge_cover = TRUE;
    Boolean *edge_covered = NEW_ARRAY(num_edges, Boolean);

    for (i = 0; i < path_len - 1; i++) {
        for (j = 0; j < num_edges; j++) {
            if ((edges_begin[j] == path[i] && edges_end[j] == path[i + 1])
            || (edges_begin[j] == path[i + 1] && edges_end[j] == path[i]))
                edge_covered[j] = TRUE;
        }
    }

    for (j = 0; j < num_edges; j++) {
        if (!edge_covered[j])
            edge_cover = FALSE;
    }

    my_free(edge_covered);
    return edge_cover;
}

// ---------------------------------------------------

// Symplectic Basis

/*
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the get_gluing_equations and get_symplectic_equations functions
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols) {
    int i, j;
    Boolean *edge_classes = NEW_ARRAY(manifold->num_tetrahedra, Boolean);
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
 * and reduce the chance of memory leak since we don't need all
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

int **get_symplectic_equations(Triangulation *manifold, Boolean *edge_classes, int *num_rows, int num_cols) {
    int i, j, k;
    label_triangulation_edges(manifold);

    ManifoldBoundary **cusps                = NEW_ARRAY(manifold->num_cusps, ManifoldBoundary *);
    Cusp *cusp                              = manifold->cusp_list_begin.next;
    EndMultiGraph *end_multi_graph          = init_end_multi_graph(manifold);

    for (i = 0; i < end_multi_graph->num_edge_classes; i++)
        edge_classes[i] = end_multi_graph->edge_classes[i] == TRUE ? FALSE : TRUE;

    edge_classes[end_multi_graph->e0] = FALSE;

    OscillatingCurves *oscillating_curves   = init_oscillating_curves(manifold, edge_classes);

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

//    do_manifold_train_lines(cusps, end_multi_graph);
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

    boundary->manifold              = manifold;
    boundary->cusp                  = cusp;
    boundary->num_edge_classes      = manifold->num_tetrahedra;
    boundary->num_cusp_triangles    = 0;
    boundary->num_cusp_regions      = 0;

    find_intersection_triangle(manifold, boundary);
    init_cusp_triangulation(manifold, boundary);
    init_cusp_region(boundary);
    init_train_line(boundary);

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

                        boundary->intersect_tet_index  = cusp->basepoint_tet->index;
                        boundary->intersect_tet_vertex = cusp->basepoint_vertex;
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
    int index = 0;
    VertexIndex vertex;
    FaceIndex face;
    Tetrahedron *tet;
    CuspTriangle *tri;

    // Allocate Cusp Triangulation Header and Tail Null nodes
    boundary->cusp_triangle_begin.next      = &boundary->cusp_triangle_end;
    boundary->cusp_triangle_begin.prev      = NULL;
    boundary->cusp_triangle_end.next        = NULL;
    boundary->cusp_triangle_end.prev        = &boundary->cusp_triangle_begin;

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
            tri->cusp = tet->cusp[vertex];
            tri->tet_index = tri->tet->index;
            tri->tet_vertex = vertex;

            tri->num_curves = net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][0])
                              + net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][1])
                              + net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][2]);

            for (face = 0; face < 4; face ++) {
                if (tri->tet_vertex == face)
                    continue;

                tri->vertices[face].v1              = tri->tet_vertex;
                tri->vertices[face].v2              = face;
                tri->vertices[face].edge            = tri->tet->edge_class[
                        edge_between_vertices[tri->vertices[face].v1][tri->vertices[face].v2]];
                tri->vertices[face].edge_class      = tri->vertices[face].edge->index;
                tri->vertices[face].edge_index      = -1;
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

        cusp_vertex             = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
        gluing_vertex           = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_outside_vertex);
        outside_vertex          = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_gluing_vertex);
        tri                     = tri->neighbours[old_gluing_vertex];
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
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 1; i < net_flow_around_vertex(tri, vertex); i++) {
            for (curve_index = 0; curve_index < 2; curve_index++) {
                distance[v1]                    = i;
                distance[v2]                    = MIN(distance[v1], 2 * net_flow_around_vertex(tri, vertex) - distance[v1])
                                                + net_flow_around_vertex(tri, v2) + net_flow_around_vertex(tri, v1);
                distance[vertex]                = net_flow_around_vertex(tri, vertex)
                                                - distance[v1] + net_flow_around_vertex(tri, v1);
                distance[tri->tet_vertex]       = -1;

                adj_triangle[v1]                = 1;
                adj_triangle[v2]                = 0;
                adj_triangle[vertex]            = 0;
                adj_triangle[tri->tet_vertex]   = -1;

                set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);
                index++;

                // Swap vertices
                v1 = (int) remaining_face[vertex][tri->tet_vertex];
                v2 = (int) remaining_face[tri->tet_vertex][vertex];
            }
        }

        // Region in the middle of face vertex
        if (net_flow_around_vertex(tri, v1) && net_flow_around_vertex(tri, v2)) {
            distance[v1]                    = net_flow_around_vertex(tri, v2);
            distance[v2]                    = net_flow_around_vertex(tri, v1);
            distance[vertex]                = MIN(net_flow_around_vertex(tri, v1) + distance[v1],
                                              net_flow_around_vertex(tri, v2) + distance[v2]) + net_flow_around_vertex(tri, vertex);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[v1]                = 0;
            adj_triangle[v2]                = 0;
            adj_triangle[vertex]            = 1;
            adj_triangle[tri->tet_vertex]   = -1;

            set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);
            index++;
        }
    }

    // Region of distance 0 to vertex
    v1 = edgesThreeToFour[tri->tet_vertex][0];
    v2 = edgesThreeToFour[tri->tet_vertex][1];
    v3 = edgesThreeToFour[tri->tet_vertex][2];

    // Edge Case: Two vertices with 0 flow
    if (at_least_two(!net_flow_around_vertex(tri, v1), !net_flow_around_vertex(tri, v2), !net_flow_around_vertex(tri, v3)))
        return init_intersect_vertex_two_zero_flows(boundary, tri, index);

    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        distance[vertex]               = 0;
        distance[v1]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v1);
        distance[v2]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v2);
        distance[tri->tet_vertex]      = -1;

        adj_triangle[vertex]           = 0;
        adj_triangle[v1]               = 1;
        adj_triangle[v2]               = 1;
        adj_triangle[tri->tet_vertex]  = 0;

        set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);
        index++;
    }

    return index;
}

int init_intersect_vertex_two_zero_flows(ManifoldBoundary *boundary, CuspTriangle *tri, int index) {
    int vertex, v1, v2, v3, distance[4];
    Boolean adj_triangle[4];

    v1 = (int) edgesThreeToFour[tri->tet_vertex][0];
    v2 = (int) edgesThreeToFour[tri->tet_vertex][1];
    v3 = (int) edgesThreeToFour[tri->tet_vertex][2];

    distance[v1]                   = net_flow_around_vertex(tri, v1);
    distance[v2]                   = net_flow_around_vertex(tri, v2);
    distance[v3]                   = net_flow_around_vertex(tri, v3);
    distance[tri->tet_vertex]      = -1;

    adj_triangle[v1]               = 1;
    adj_triangle[v2]               = 1;
    adj_triangle[v3]               = 1;
    adj_triangle[tri->tet_vertex]  = -1;

    set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);
    index++;

    // Find vertex with non-zero flow
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        if (net_flow_around_vertex(tri, vertex)) {
            v1 = vertex;
            v2 = (int) remaining_face[tri->tet_vertex][v1];
            v3 = (int) remaining_face[v1][tri->tet_vertex];
            break;
        }
    }
    distance[v1]                    = 0;
    distance[v2]                    = net_flow_around_vertex(tri, v1);
    distance[v3]                    = net_flow_around_vertex(tri, v1);
    distance[tri->tet_vertex]       = -1;

    adj_triangle[v1]                = 0;
    adj_triangle[v2]                = 1;
    adj_triangle[v3]                = 1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);

    return index + 1;
}

int init_normal_cusp_region(ManifoldBoundary *boundary, CuspTriangle *tri, int index) {
    int i, vertex, v1, v2;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 0; i < net_flow_around_vertex(tri, vertex); i++) {
            distance[vertex]                = i;
            distance[v1]                    = net_flow_around_vertex(tri, v1)
                                            + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[v2]                    = net_flow_around_vertex(tri, v2)
                                            + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[vertex]            = 0;
            adj_triangle[v1]                = 1;
            adj_triangle[v2]                = 1;
            adj_triangle[tri->tet_vertex]   = 0;

            set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);
            index++;
        }

    }

    // center region
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        distance[vertex]        = net_flow_around_vertex(tri, vertex);
        adj_triangle[vertex]    = 1;
    }

    distance[tri->tet_vertex]       = -1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(&boundary->cusp_region_end, tri, distance, adj_triangle, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which
 * allocates the attributes of the cusp region
 */

void set_cusp_region_data(CuspRegion *cusp_region_end, CuspTriangle *tri, const int distance[4],
                          const Boolean adj_cusp_triangle[4], int index) {
    int i, j, v1, v2, v3;
    CuspRegion *region = NEW_STRUCT( CuspRegion );
    INSERT_BEFORE(region, cusp_region_end);

    region->tri             = tri;
    region->tet_index       = region->tri->tet_index;
    region->tet_vertex      = region->tri->tet_vertex;
    region->index           = index;

    // default values
    for (i = 0; i < 4; i++) {
        region->adj_cusp_triangle[i] = 0;
        region->adj_dive_regions[i]  = NULL;

        for (j = 0; j < 4; j++) {
            region->curve[i][j]             = -1;
            region->dive[i][j]              = 0;
            region->num_adj_curves[i][j]    = 0;
            region->temp_adj_curves[i][j]   = 0;
        }
    }

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tet_vertex][i];
        v2 = edgesThreeToFour[tri->tet_vertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tet_vertex][(i + 2) % 3];

        region->curve[v2][v1]   = distance[v1];
        region->curve[v3][v1]   = distance[v1];
        region->dive[v2][v1]    = distance[v1] ? FALSE : TRUE;
        region->dive[v3][v1]    = distance[v1] ? FALSE : TRUE;

        region->adj_cusp_triangle[v1] = adj_cusp_triangle[v1];
    }
}

/*
 * Calculate which regions are located across cusp edges
 * and store the result in the adj_cusp_regions attribute
 */

void update_adj_region_data(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end) {
    int j;
    CuspRegion *region;

    // Add adjacent region info
    for (region = cusp_region_begin->next; region != cusp_region_end; region = region->next) {
        for (j = 0; j < 4; j++) {
            if (!region->adj_cusp_triangle[j] || region->tet_vertex == j) {
                region->adj_cusp_regions[j] = NULL;
                continue;
            }

            region->adj_cusp_regions[j] = find_adj_region(cusp_region_begin, cusp_region_end, region, j);
        }
    }
}

/*
 * Find the cusp region which is adjacent
 * to x across face.
 */

CuspRegion *find_adj_region(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end,
                            CuspRegion *x, int face) {
    int v1, v2, y_vertex1, y_vertex2, y_face, distance_v1, distance_v2, tet_index, tet_vertex;
    Boolean adj_face;
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

void init_train_line(ManifoldBoundary *cusp) {

}

DualCurves *init_dual_curve(int edge_class_start, int edge_class_finish) {
    DualCurves *path = NEW_STRUCT( DualCurves );

    path->curves_begin.next     = &path->curves_end;
    path->curves_begin.prev     = NULL;
    path->curves_end.next       = NULL;
    path->curves_end.prev       = &path->curves_begin;

    path->edge_class[START]     = edge_class_start;
    path->edge_class[FINISH]    = edge_class_finish;

    return path;
}

/*
 * Initialise dual curve doubly linked list which
 * stores the oscillating curves on the cusp
 */

OscillatingCurves *init_oscillating_curves(Triangulation *manifold, const Boolean *edge_classes) {
    int i, j;
    OscillatingCurves *curves = NEW_STRUCT(OscillatingCurves );

    curves->num_curves = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++)
        if (edge_classes[i])
            curves->num_curves++;

    curves->dual_curve_begin          = NEW_ARRAY(curves->num_curves, DualCurves );
    curves->dual_curve_end            = NEW_ARRAY(curves->num_curves, DualCurves );
    curves->edge_class                = NEW_ARRAY(curves->num_curves, int);

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
    PathNode *path_node;

    for (i = 0; i < curves->num_curves; i++) {
        while (curves->dual_curve_begin[i].next != &curves->dual_curve_end[i]) {
            path = curves->dual_curve_begin[i].next;
            REMOVE_NODE(path);

            while (path->curves_begin.next != &path->curves_end) {
                path_node = path->curves_begin.next;
                REMOVE_NODE(path_node);
                my_free(path_node);
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
 *  - 1: Train Lines
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
    CuspRegion *region;
    EdgeNode *edge_node;
    PathNode *path_node;
    DualCurves *path;
    Graph *g;
    ManifoldBoundary *boundary;

    if (!debug)
        return;

    if (flag != 1)
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
        // Cusp Train Lines
        printf("Train Lines\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            boundary = cusps[i];
            printf("    Train Line Regions:");
            for (j = 0; j < boundary->train_line_len; j++) {
                printf(" %d", boundary->train_line_path[j]);
            }

            printf("\n");
        }
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

                for (path_node = path->curves_begin.next; 
                    path_node != &path->curves_end; 
                    path_node = path_node->next)
                    printf("%d ", path_node->cusp_region_index);

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
                for (edge_node = g->edge_list_begin[j].next;
                    edge_node != &g->edge_list_end[j];
                    edge_node = edge_node->next)
                    printf("%d ", edge_node->y);

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

int net_flow_around_vertex(CuspTriangle *tri, int vertex) {
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

// ------------------------------------

/*
 * Train lines 
 *
 * Each cusp contains a number of edge classes which are split into two partitions, 
 * those in the end multi graph and those not in the end multi graph.
 *
 * We want a collection of curves which connect any edge class in the end multi graph
 * to any other on a single cusp for use when finding oscillating curves. 
 *
 * Firstly we construct a graph consisting of vertices which have an edge class in 
 * the end multi graph and edges consisting of cusp edges which connect vertices 
 * of the graph. 
 *
 * Then for each connected component we find a spanning tree. The edges of each 
 * spanning tree represent the faces we dive along to dive into the manifold. For 
 * singletons we add any edge around the vertex. 
 *
 * This choice on one cusp induces choices on the other cusps. 
 *
 * For each tree on a cusp, we construct a curve which passses through all edges 
 * in the tree. Then join each of these curves together to form the final train tracks.
 * The train traks give a path to connect any edge class to any other, ensuring curves 
 * along gamma have intersection number 0. 
 *
 * After finding the train tracks, we need to split each cusp region near a vertex 
 * corresponding to an edge class since we don't do this when we find train tracks, 
 * nor when we find a path along the train tracks.
 */

void do_manifold_train_lines(ManifoldBoundary **cusps, EndMultiGraph *multi_graph) {
    int cusp_index, edge_class;
    int **edges_begin = NEW_ARRAY(multi_graph->num_cusps, int *),
        **edges_end = NEW_ARRAY(multi_graph->num_cusps, int *);
    PathEndPoint ***endpoints = NEW_ARRAY(multi_graph->num_cusps, PathEndPoint **);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        endpoints[cusp_index] = NEW_ARRAY(multi_graph->num_edge_classes, PathEndPoint *);
        edges_begin[cusp_index] = NEW_ARRAY(multi_graph->num_edge_classes, int);
        edges_end[cusp_index] = NEW_ARRAY(multi_graph->num_edge_classes, int);

        for (edge_class = 0; edge_class < multi_graph->num_edge_classes; edge_class++) {
            endpoints[cusp_index][edge_class] = NULL;
        }
    }

    find_edge_class_edges(cusps, multi_graph, edges_begin, edges_end);
    find_train_line_endpoints(cusps, multi_graph, endpoints, edges_begin, edges_end);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        cusps[cusp_index]->train_line_path = find_primary_train_line(cusps[cusp_index],
                                                                     multi_graph,
                                                                     edges_begin[cusp_index],
                                                                     edges_end[cusp_index],
                                                                     &cusps[cusp_index]->train_line_len);
        cusps[cusp_index]->train_line_endpoint = endpoints[cusp_index];

        update_cusp_regions_along_train_line(cusps[cusp_index]);
    }

    print_debug_info(cusps[0]->manifold, cusps, NULL, 1);

    my_free(endpoints);
}

void find_edge_class_edges(ManifoldBoundary **cusps, EndMultiGraph *multi_graph, int **edges_begin, int **edges_end) {
    int edge_class, cusp_index;
    Boolean found_edge_class;
    Boolean **edge_classes = NEW_ARRAY(multi_graph->num_cusps, Boolean *);
    Graph *cusp_graph, **cusp_connected_graphs;
    Queue *queue = init_queue(multi_graph->num_cusps);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        edge_classes[cusp_index] = edge_classes_on_cusp(multi_graph, cusp_index);
        edge_classes[cusp_index][multi_graph->e0] = FALSE;
    }

    enqueue(queue, 0);

    while (!empty_queue(queue)) {
        cusp_index = dequeue(queue);

        found_edge_class = FALSE;
        for (edge_class = 0; edge_class < multi_graph->num_edge_classes; edge_class++) {
            if (!edge_classes[cusp_index][edge_class])
                found_edge_class = TRUE;
        }

        if (!found_edge_class)
            continue;

        // edge_classes[cusp_index] contains some true entries
        cusp_graph = construct_edge_class_graph(&cusps[cusp_index]->cusp_triangle_begin,
                                                &cusps[cusp_index]->cusp_triangle_end,
                                                edge_classes[cusp_index],
                                                multi_graph->num_edge_classes);

        cusp_connected_graphs = find_connected_graph_components(cusp_graph);

        // assign edges to edge classes
        find_edge_class_edges_on_cusp(cusps[cusp_index],
                                      multi_graph,
                                      cusp_connected_graphs,
                                      cusps[cusp_index]->dual_graph,
                                      edges_begin[cusp_index],
                                      edges_end[cusp_index]);

        // update dive edges classes
        update_edge_classes_on_cusp(cusps, cusp_index, edge_classes, multi_graph->num_cusps, multi_graph->num_edge_classes, edges_begin, edges_end);
    }
}

Graph *construct_edge_class_graph(CuspTriangle *triangle_begin, CuspTriangle *triangle_end,
                                  const Boolean *edge_classes, int num_edge_classes) {
    Graph *g = init_graph(num_edge_classes, FALSE);
    CuspTriangle *tri;
    VertexIndex v1, v2;
    CuspVertex *vertex1, *vertex2;

    for (tri = triangle_begin->next; tri != triangle_end; tri = tri->next) {
        for (v1 = 0; v1 < 4; v1++) {
            if (v1 == tri->tet_vertex)
                continue;

            for (v2 = 0; v2 < 4; v2++) {
                if (v2 == v1 || v2 == tri->tet_vertex)
                    continue;

                vertex1 = &tri->vertices[v1];
                vertex2 = &tri->vertices[v2];

                if (edge_classes[vertex1->edge_class] && edge_classes[vertex2->edge_class]) {
                    insert_edge(g, vertex1->edge_class, vertex2->edge_class, g->directed);
                }
            }
        }
    }

    return g;
}

void find_edge_class_edges_on_cusp(ManifoldBoundary *cusp, EndMultiGraph *multi_graph, 
                                   Graph **cusp_graphs, Graph *dual_graph,
                                   int *edges_begin, int *edges_end) {
    int edge_class, cusp_graph_index;
    VertexIndex v1, v2;
    FaceIndex face;
    Boolean found;
    Boolean *edge_classes = edge_classes_on_cusp(multi_graph, cusp->cusp->index);
    CuspRegion *region;
    CuspVertex *vertex1, *vertex2;

    for (edge_class = 0; edge_class < multi_graph->num_edge_classes; edge_class++) {
        edges_begin[edge_class] = -1;
        edges_end[edge_class] = -1;

        if (!edge_classes[edge_class])
            continue;

        // find edge in cusp_graphs incident to edge_class
        for (region = cusp->cusp_region_begin.next; region != &cusp->cusp_region_begin; region = region->next) {
            for (face = 0; face < 4; face++) {
                if (face == region->tet_vertex)
                    continue;

                if (!region->adj_cusp_triangle[face])
                    continue;

                v1 = remaining_face[region->tet_vertex][face];
                v2 = remaining_face[face][region->tet_vertex];

                vertex1 = &region->tri->vertices[v1];
                vertex2 = &region->tri->vertices[v2];

                if (vertex1->edge_class != edge_class && vertex2->edge_class != edge_class)
                    continue;

                if ((vertex1->edge_class == edge_class && !region->dive[face][v1]) ||
                (vertex2->edge_class == edge_class && !region->dive[face][v2]))
                    continue;

                found = FALSE;
                // is 'face' an edge of the cusp graphs?
                for (cusp_graph_index = 0; cusp_graph_index < multi_graph->num_edge_classes; cusp_graph_index++) {
                    if (cusp_graphs[cusp_graph_index] == NULL)
                        continue;

                    if (edge_exists(cusp_graphs[cusp_graph_index], vertex1->edge_class, vertex2->edge_class))
                        found = TRUE;
                }

                if (!found)
                    continue;

                // region lies adjacent to an edge in the cusp graphs and can dive into the target edge class along face
                edges_begin[edge_class] = region->index;
                edges_end[edge_class]   = region->adj_cusp_regions[face]->index;
            }
        }
    }
}

void update_edge_classes_on_cusp(ManifoldBoundary **cusps, int current_cusp_index, Boolean **edge_classes, 
                                 int num_cusps, int num_edge_classes, int **edges_begin, int **edges_end) {
    int cusp_index, other_cusp_index, edge_class;
    VertexIndex v1, v2, vertex;
    FaceIndex face;
    CuspRegion *region1, *region2, *adj_region1, *adj_region2;
    ManifoldBoundary *cusp = cusps[current_cusp_index];

    for (edge_class = 0; edge_class < num_edge_classes; edge_class++) {
        if (!edge_classes[current_cusp_index][edge_class])
            continue;

        region1 = get_region_from_index(&cusp->cusp_region_begin, &cusp->cusp_region_end, edges_begin[cusp->cusp->index][edge_class]);
        region2 = get_region_from_index(&cusp->cusp_region_begin, &cusp->cusp_region_end, edges_end[cusp->cusp->index][edge_class]);

        for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
            if (edge_classes[cusp_index][edge_class] && cusp_index != cusp->cusp->index)
                other_cusp_index = cusp_index;
        }

        for (face = 0; face < 4; face++) {
            if (face == region1->tet_vertex)
                continue;

            if (region1->adj_cusp_regions[face]->index == region2->index)
                break;
        }

        v1 = remaining_face[region1->tet_vertex][face];
        v2 = remaining_face[face][region1->tet_vertex];

        if (region1->tri->vertices[v1].edge_class == edge_class)
            vertex = v1;
        else
            vertex = v2;

        edge_classes[other_cusp_index][edge_class] = FALSE;

        adj_region1 = find_adj_dive_region(&cusps[other_cusp_index]->cusp_region_begin,
                                           &cusps[other_cusp_index]->cusp_region_end,
                                           region1, face, vertex);
        adj_region2 = find_adj_dive_region(&cusps[other_cusp_index]->cusp_region_begin,
                                           &cusps[other_cusp_index]->cusp_region_end,
                                           region2,
                                           EVALUATE(region1->tri->tet->gluing[face], face),
                                           EVALUATE(region1->tri->tet->gluing[face], vertex));

        edges_begin[other_cusp_index][edge_class] = adj_region1->index;
        edges_end[other_cusp_index][edge_class] = adj_region2->index;
    }
}

int *find_primary_train_line(ManifoldBoundary *cusp, EndMultiGraph *multi_graph, int *edges_start, int *edges_end, int *path_len) {
    int *paths;

    paths = dfs_target_edges(cusp->dual_graph, edges_start, edges_end, multi_graph->num_edge_classes, path_len);

    return paths;
}

void find_train_line_endpoints(ManifoldBoundary **cusps, EndMultiGraph *multi_graph, PathEndPoint ***endpoints, 
                               int **edges_begin, int **edges_end) {
    int cusp_index;
    Queue *queue = init_queue(multi_graph->num_cusps);
    Boolean **edge_classes = NEW_ARRAY(multi_graph->num_cusps, Boolean *);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        edge_classes[cusp_index] = edge_classes_on_cusp(multi_graph, cusp_index);
        edge_classes[cusp_index][multi_graph->e0] = FALSE;
    }

    enqueue(queue, 0);

    while (!empty_queue(queue)) {
        cusp_index = dequeue(queue);

        find_train_line_endpoints_on_cusp(cusps[cusp_index], 
                                          multi_graph, 
                                          endpoints[cusp_index], 
                                          edges_begin[cusp_index], 
                                          edges_end[cusp_index], 
                                          edge_classes[cusp_index]);
    }
}

void find_train_line_endpoints_on_cusp(ManifoldBoundary *cusp, EndMultiGraph *multi_graph, PathEndPoint **endpoints,
                                       int *edges_begin, int *edges_end, const Boolean *edge_classes) {
    int edge_class;
    CuspRegion *region1, *region2;
    FaceIndex face, f;
    VertexIndex vertex, v, v1, v2;

    for (edge_class = 0; edge_class < multi_graph->num_edge_classes; edge_class++) {
        if (!edge_classes[edge_class]) {
            endpoints[edge_class] = NULL;
            continue;
        }

        if (edges_begin[edge_class] == -1)
            uFatalError("find_train_line_endpoints_on_cusp", "symplectic_basis");

        region1 = get_region_from_index(&cusp->cusp_region_begin, &cusp->cusp_region_end, edges_begin[edge_class]);
        region2 = get_region_from_index(&cusp->cusp_region_begin, &cusp->cusp_region_end, edges_end[edge_class]);

        face = -1;
        for (f = 0; f < 4; f++) {
            if (f == region1->tet_vertex)
                continue;

            if (region1->adj_cusp_regions[f]->index == region2->index)
                face = f;
        }

        v1 = remaining_face[region1->tet_vertex][face];
        v2 = remaining_face[face][region1->tet_vertex];

        if (region1->tri->vertices[v1].edge_class == edge_class)
            vertex = v1;
        else 
            vertex = v2;

        endpoints[edge_class] = NEW_STRUCT( PathEndPoint );
        endpoints[edge_class]->vertex = vertex;
        endpoints[edge_class]->face = face;
        endpoints[edge_class]->region = region1;
        endpoints[edge_class]->region_index = region1->index;

        for (f = 0; f < 4; f++)
            for (v = 0; v < 4; v++) 
                endpoints[edge_class]->num_adj_curves[f][v] = 0;
    }
}

void update_cusp_regions_along_train_line(ManifoldBoundary *cusp) {
    int index = 0;
    CuspRegion *region, *p_region;
    PathNode *path_begin = NEW_STRUCT( PathNode ),
             *path_end   = NEW_STRUCT( PathNode );
    PathNode *node;

    convert_path_array_to_path_node(&cusp->cusp_region_begin, &cusp->cusp_region_end,
                                    path_begin, path_end,cusp->train_line_path,cusp->train_line_len);

    path_begin->next = path_end;
    path_begin->prev = NULL;
    path_end->next = NULL;
    path_end->prev = path_begin;

    // interior edges
    for (node = path_begin->next; node != path_end; node = node->next) {
        p_region = get_region_from_index(&cusp->cusp_region_begin, &cusp->cusp_region_end, node->cusp_region_index);
        update_cusp_triangle_path_interior(&cusp->cusp_region_begin, &cusp->cusp_region_end, p_region, node);
        region = split_cusp_region_path_interior(p_region, node, index);
        INSERT_BEFORE(region, &cusp->cusp_region_end);
        index++;
    }
}

void convert_path_array_to_path_node(CuspRegion *region_begin, CuspRegion *region_end, PathNode *path_begin,
                                     PathNode *path_end, const int *path, int path_len) {
    int i;
    PathNode *node;
    CuspRegion *region;

    for (i = 1; i < path_len - 1; i++) {
        region = get_region_from_index(region_begin, region_end, path[i]);
        node = interior_path_point_to_path_node(region, path, i);
        INSERT_BEFORE(node, path_end);
    }
}

PathNode *endpoint_path_point_to_path_node(CuspRegion *region, PathEndPoint *path_endpoint, const int *path, int pos) {
    int next_node;
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = path[pos];
    path_node->tri = region->tri;

    if (pos == 0)
        next_node = path[1];
    else
        next_node = path[pos - 1];

    path_node->next_face = -1;
    path_node->prev_face = -1;

    for (face = 0; face < 4; face++) {
        if (face == region->tet_vertex || !region->adj_cusp_triangle[face])
            continue;

        if (region->adj_cusp_regions[face]->index == next_node) {
            path_node->prev_face = face;
            path_node->next_face = face;
        }
    }

    // next node isn't in an adjacent region
    if (path_node->next_face == -1 || path_node->prev_face == -1)
        uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

    if (pos == 0)
        path_node->prev_face = path_endpoint->face;
    else
        path_node->next_face = path_endpoint->face;

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];


    if (path_node->prev_face == path_endpoint->vertex) {
        if (path_endpoint->face == vertex1)
            path_node->inside_vertex = vertex2;
        else
            path_node->inside_vertex = vertex1;
    } else if (path_node->prev_face == path_endpoint->face) {
        path_node->inside_vertex = -1;
    } else {
        path_node->inside_vertex = path_endpoint->vertex;
    }

    if (path_node->next_face == path_node->prev_face)
        return path_node;
        //uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

    return path_node;
}

PathNode *interior_path_point_to_path_node(CuspRegion *region, const int *path, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = path[pos];
    path_node->tri = region->tri;

    for (face = 0; face < 4; face++) {
        if (face == region->tet_vertex)
            continue;

        vertex1 = remaining_face[region->tet_vertex][face];
        vertex2 = remaining_face[face][region->tet_vertex];

        if (!region->adj_cusp_triangle[vertex1] || !region->adj_cusp_triangle[vertex2])
            continue;

        if (region->adj_cusp_regions[vertex1]->index == path[pos + 1] 
         && region->adj_cusp_regions[vertex2]->index == path[pos - 1]) {
            path_node->inside_vertex = face;
            path_node->next_face = vertex1;
            path_node->prev_face = vertex2;
            return path_node;

        } else if (region->adj_cusp_regions[vertex2]->index == path[pos + 1] 
                && region->adj_cusp_regions[vertex1]->index == path[pos - 1]) {
            path_node->inside_vertex = face;
            path_node->next_face = vertex2;
            path_node->prev_face = vertex1;
            return path_node;

        }
    }

    uFatalError("interior_path_point_to_path_node", "symplectic_basis");
    return NULL;
}

CuspRegion *get_region_from_index(CuspRegion *region_begin, CuspRegion *region_end, int index) {
    CuspRegion *region, *retval = NULL;

    for (region = region_begin->next; region != region_end; region = region->next) {
        if (region->index != index)
            continue;

        retval = region;
    }

    if (retval == NULL) {
        uFatalError("get_region_from_index", "symplectic_basis");
        return NULL;
    }

    return retval;
}

CuspRegion *find_adj_dive_region(CuspRegion *region_begin, CuspRegion *region_end, CuspRegion *start_region, int face, int vertex) {
    CuspRegion *region;

    for (region = region_begin->next; region != region_end; region = region->next) {
        if (region->tet_index != start_region->tet_index)
            continue;

        if (region->tet_vertex != vertex)
            continue;

        if (!region->adj_cusp_triangle[face])
            continue;

        if (!region->dive[face][start_region->tet_vertex])
            continue;

        return region;
    }

    uFatalError("find_adj_dive_region", "symplectic_basis");
    return NULL;
}

// ------------------------------------

void do_oscillating_curves(ManifoldBoundary **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph) {
    int i;

    for (i = 0; i < curves->num_curves; i++) {
        do_one_dual_curve(cusps, curves,
                          &curves->dual_curve_begin[i],
                          &curves->dual_curve_end[i],
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
    CuspEndPoint *cusp_end_point, *endpoint;
    DualCurves *path;
    int orientation = START;

    // find paths through cusps
    cusp_end_point = find_multi_graph_path(multi_graph, cusps[0]->manifold, edge_class, &path_length);

    // find paths inside each cusp
    dual_curve_begin->edge_class[FINISH] = edge_class;
    dual_curve_end->edge_class[START]    = edge_class;
    
    i = FIRST;
    for (endpoint = cusp_end_point; endpoint != NULL; endpoint = endpoint->next) {
        print_debug_info(cusps[0]->manifold, cusps, NULL, 2);
        print_debug_info(cusps[0]->manifold, cusps, NULL, 7);
        
        if (endpoint->next == NULL)
            i = LAST;

        path = init_dual_curve(endpoint->edge_class[START], endpoint->edge_class[FINISH]);
        INSERT_BEFORE(path, dual_curve_end);
        path->cusp_index = endpoint->cusp_index;
        find_path_endpoints(cusps[path->cusp_index]->dual_graph, dual_curve_begin, path, i, orientation);
        do_one_cusp(cusps[path->cusp_index], path);
        update_path_holonomy(path, edge_class);
        orientation = (orientation == START ? FINISH : START);
        i = MIDDLE;
    }

    update_adj_curve_along_path(cusps, curves, dual_curve_begin, dual_curve_end, curve_index);

    endpoint = cusp_end_point;
    while (endpoint != NULL) {
        cusp_end_point = endpoint;
        endpoint = endpoint->next;
        my_free(cusp_end_point);
    }
}

/*
 * Construct oscillating curves on the boundary components
 */

void do_one_cusp(ManifoldBoundary *boundary, DualCurves *path) {
    int *parent;
    Boolean *processed, *discovered;
    EdgeNode *node_begin, *node_end, *edge_node;

    processed   = NEW_ARRAY(boundary->dual_graph->num_vertices, Boolean);
    discovered  = NEW_ARRAY(boundary->dual_graph->num_vertices, Boolean);
    parent      = NEW_ARRAY(boundary->dual_graph->num_vertices, int);

    node_begin  = NEW_STRUCT( EdgeNode );
    node_end    = NEW_STRUCT( EdgeNode );

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    // Find path using bfs
    init_search(boundary->dual_graph, processed, discovered, parent);
    bfs(boundary->dual_graph, path->endpoints[START].region_index, processed, discovered, parent);
    find_path(path->endpoints[START].region_index, path->endpoints[FINISH].region_index,
              parent, node_begin);
    graph_path_to_dual_curve(boundary->dual_graph, node_begin, node_end, path);

    // Reallocate memory
    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // Split the regions along the path
    split_cusp_regions_along_path(boundary, path);

    free_graph(boundary->dual_graph);
    boundary->dual_graph = construct_cusp_region_dual_graph(boundary);

    while (node_begin->next != node_end) {
        edge_node = node_begin->next;
        REMOVE_NODE(edge_node);
        my_free(edge_node);
    }

    my_free(node_begin);
    my_free(node_end);
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
        endpoint1_path2         = &path->endpoints[START];
        endpoint1_edge_class    = path->edge_class[START];
        endpoint1_edge_index    = START;
        find_single_endpoint(g, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path2         = &path->endpoints[FINISH];
        endpoint2_edge_class    = path->edge_class[FINISH];
        endpoint2_edge_index    = START;
        find_single_endpoint(g, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else if (pos == MIDDLE && orientation == START) {
        endpoint1_path1         = &path->prev->endpoints[START];
        endpoint1_path2         = &path->endpoints[START];
        endpoint1_edge_class    = path->edge_class[START];
        endpoint1_edge_index    = FINISH;
        find_single_matching_endpoint(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path2         = &path->endpoints[FINISH];
        endpoint2_edge_class    = path->edge_class[FINISH];
        endpoint2_edge_index    = START;
        find_single_endpoint(g, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else if (pos == MIDDLE && orientation == FINISH) {
        endpoint1_path1         = &path->prev->endpoints[FINISH];
        endpoint1_path2         = &path->endpoints[FINISH];
        endpoint1_edge_class    = path->edge_class[START];
        endpoint1_edge_index    = FINISH;
        find_single_matching_endpoint(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path2         = &path->endpoints[START];
        endpoint2_edge_class    = path->edge_class[FINISH];
        endpoint2_edge_index    = START;
        find_single_endpoint(g, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else if (pos == LAST) {
        endpoint1_path1         = &path_start->next->endpoints[START];
        endpoint1_path2         = &path->endpoints[START];
        endpoint1_edge_class    = path_start->next->edge_class[START];
        endpoint1_edge_index    = FINISH;
        find_single_matching_endpoint(g, endpoint1_path1, endpoint1_path2, endpoint1_edge_class, endpoint1_edge_index);

        endpoint2_path1         = &path->prev->endpoints[FINISH];
        endpoint2_path2         = &path->endpoints[FINISH];
        endpoint2_edge_class    = path->prev->edge_class[FINISH];
        endpoint2_edge_index    = FINISH;
        find_single_matching_endpoint(g, endpoint2_path1, endpoint2_path2, endpoint2_edge_class, endpoint2_edge_index);
    } else {
        uFatalError("do_one_cusp", "symplectic_basis");
        return;
    }
}

Boolean is_valid_endpoint(CuspRegion *region, int edge_class, int edge_index, int face, int vertex) {
    Boolean edge_lies_in_one_cusp;

    if (vertex == region->tet_vertex)
        return FALSE;

    if (region->tri->vertices[vertex].edge_class != edge_class)
        return FALSE;

    if (region->tri->tet->cusp[region->tet_vertex]->index == region->tri->tet->cusp[vertex]->index)
        edge_lies_in_one_cusp = TRUE;
    else
        edge_lies_in_one_cusp = FALSE;

    if (edge_lies_in_one_cusp && region->tri->vertices[vertex].edge_index != edge_index) {
        return FALSE;
    }

    if (region->adj_cusp_triangle[face] == TRUE && region->curve[face][vertex] == 0)
        return TRUE;
    else
        return FALSE;

    return region->dive[face][vertex];
}

PathEndPoint *find_single_endpoint(Graph *g, PathEndPoint *path_endpoint, int edge_class, int edge_index) {
    int i, j, k;
    VertexIndex vertex;
    FaceIndex face1, face2, face;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < g->num_vertices; i++) {
        if (g->regions[i] == NULL) {
            continue;
        }

        region = g->regions[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            face1 = remaining_face[region->tet_vertex][vertex];
            face2 = remaining_face[vertex][region->tet_vertex];

            if (is_valid_endpoint(region, edge_class, edge_index, face1, vertex))
                face = face1;
            else if (is_valid_endpoint(region, edge_class, edge_index, face2, vertex))
                face = face2;
            else
                continue;

            path_endpoint->region           = region;
            path_endpoint->vertex           = vertex;
            path_endpoint->face             = face;
            path_endpoint->region_index     = i;

            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    path_endpoint->num_adj_curves[j][k] = region->num_adj_curves[j][k];
            return path_endpoint;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_single_endpoints", "symplectic_basis");
    return NULL;
}

PathEndPoint *find_single_matching_endpoint(Graph *g, PathEndPoint *path_endpoint1, PathEndPoint *path_endpoint2, int edge_class, int edge_index) {
    int i, j, k;
    Boolean region_index, region_vertex, dive_vertex, region_dive, region_curve;
    VertexIndex vertex;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < g->num_vertices; i++) {
        if (g->regions[i] == NULL)
            continue;

        region = g->regions[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            if (!is_valid_endpoint(region, edge_class, edge_index, path_endpoint1->face, vertex))
                continue;

            // are we in the correct region for copy
            region_index    = (Boolean) (region->tet_index != path_endpoint1->region->tet_index);
            region_vertex   = (Boolean) (region->tet_vertex != path_endpoint1->vertex);
            dive_vertex     = (Boolean) (vertex != path_endpoint1->region->tet_vertex);
            region_dive     = (Boolean) !region->adj_cusp_triangle[path_endpoint1->face]; //!region->dive[path_endpoint1->face][vertex];
            region_curve    = (Boolean) region->curve[path_endpoint1->face][vertex]; //(region->num_adj_curves[path_endpoint1->face][vertex] != path_endpoint1->num_adj_curves[path_endpoint1->face][path_endpoint1->vertex]);

            if (region_index || region_vertex || dive_vertex || region_dive || region_curve)
                continue;

            path_endpoint2->region          = region;
            path_endpoint2->vertex          = vertex;
            path_endpoint2->face            = path_endpoint1->face;
            path_endpoint2->region_index    = i;

            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    path_endpoint2->num_adj_curves[j][k] = region->num_adj_curves[j][k];

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

void graph_path_to_dual_curve(Graph *g, EdgeNode *node_begin, EdgeNode *node_end, DualCurves *path) {
    FaceIndex face;
    EdgeNode *edge_node;
    PathNode *path_node;
    CuspRegion *region;

    // path len 0
    if (node_begin->next == node_end)
        return;

    edge_node = node_begin->next;
    // path len 1
    if (edge_node->next == node_end) {
        for (face = 0; face < 4; face++)
            if (g->regions[edge_node->y]->tet_vertex != face &&
                path->endpoints[START].vertex != face &&
                path->endpoints[FINISH].vertex != face)
                break;

        region = g->regions[edge_node->y];

        path_node = NEW_STRUCT( PathNode );
        INSERT_AFTER(path_node, &path->curves_begin);
        path_node->next_face = path->endpoints[FINISH].face;
        path_node->prev_face = path->endpoints[START].face;
        path_node->cusp_region_index = edge_node->y;
        path_node->tri = region->tri;
        path_node->inside_vertex = face;
        return;
    }

    // Set Header node
    path_node = endpoint_edge_node_to_path_node(g->regions[edge_node->y], edge_node, &path->endpoints[START], START);
    INSERT_BEFORE(path_node, &path->curves_end);

    for (edge_node = node_begin->next->next; edge_node->next != node_end; edge_node = edge_node->next) {
        region = g->regions[edge_node->y];
        path_node = interior_edge_node_to_path_node(region, edge_node);
        INSERT_BEFORE(path_node, &path->curves_end);
    }

    // Set Tail node
    path_node = endpoint_edge_node_to_path_node(g->regions[edge_node->y], edge_node, &path->endpoints[FINISH],FINISH);
    INSERT_BEFORE(path_node, &path->curves_end);
}

PathNode *endpoint_edge_node_to_path_node(CuspRegion *region, EdgeNode *edge_node, PathEndPoint *path_endpoint, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    EdgeNode *next_node;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    if (pos == START)
        next_node = edge_node->next;
    else
        next_node = edge_node->prev;

    path_node->next_face = -1;
    path_node->prev_face = -1;

    for (face = 0; face < 4; face++) {
        if (face == region->tet_vertex || !region->adj_cusp_triangle[face])
            continue;

        if (region->adj_cusp_regions[face]->index == next_node->y) {
            path_node->prev_face = face;
            path_node->next_face = face;
        }
    }

    // next node isn't in an adjacent region
    if (path_node->next_face == -1 || path_node->prev_face == -1)
        uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

    if (pos == START)
        path_node->prev_face = path_endpoint->face;
    else
        path_node->next_face = path_endpoint->face;

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (path_node->prev_face == path_endpoint->vertex) {
        if (path_endpoint->face == vertex1)
            path_node->inside_vertex = vertex2;
        else
            path_node->inside_vertex = vertex1;
    } else if (path_node->prev_face == path_endpoint->face) {
        path_node->inside_vertex = -1;
    } else {
        path_node->inside_vertex = path_endpoint->vertex;
    }

    if (path_node->next_face == path_node->prev_face)
        return path_node;
        //uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

    return path_node;
}

/*
 * node lies in 'region', find the vertex which the subpath
 * node->prev->y --> node->y --> node->next->y cuts off of the
 * cusp triangle region->tri.
 */

PathNode *interior_edge_node_to_path_node(CuspRegion *region, EdgeNode *edge_node) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    for (face = 0; face < 4; face++) {
        if (face == region->tet_vertex)
            continue;

        vertex1 = remaining_face[region->tet_vertex][face];
        vertex2 = remaining_face[face][region->tet_vertex];

        if (!region->adj_cusp_triangle[vertex1] || !region->adj_cusp_triangle[vertex2])
            continue;

        if (region->adj_cusp_regions[vertex1]->index == edge_node->next->y 
         && region->adj_cusp_regions[vertex2]->index == edge_node->prev->y) {
            path_node->inside_vertex = face;
            path_node->next_face = vertex1;
            path_node->prev_face = vertex2;
            return path_node;

        } else if (region->adj_cusp_regions[vertex2]->index == edge_node->next->y 
                && region->adj_cusp_regions[vertex1]->index == edge_node->prev->y) {
            path_node->inside_vertex = face;
            path_node->next_face = vertex2;
            path_node->prev_face = vertex1;
            return path_node;

        }
    }

    // where does the next node go?
    uFatalError("interior_edge_node_to_path_node", "symplectic_basis");
    return NULL;
}

/*
 * The oscillating curve splits the region it passes through
 * into two regions. Split each region in two and update
 * attributes
 */

void split_cusp_regions_along_path(ManifoldBoundary *boundary, DualCurves *path) {
    int index = boundary->num_cusp_regions;
    FaceIndex face;
    PathNode *node = path->curves_begin.next;
    CuspRegion *region, *p_region;
    Graph *g = boundary->dual_graph;

    // empty path
    if (node->next == NULL)
        return ;

    // path of len 1
    if (node->next->next == NULL) {
        region = NEW_STRUCT(CuspRegion);
        p_region = g->regions[node->cusp_region_index];
        INSERT_BEFORE(region, &boundary->cusp_region_end)
        copy_region(p_region, region);

        face = node->inside_vertex;

        region->index = index;
        region->adj_cusp_triangle[path->endpoints[START].vertex]  = 0;
        region->adj_cusp_triangle[path->endpoints[FINISH].vertex] = 0;
        region->temp_adj_curves[path->endpoints[START].vertex][path->endpoints[FINISH].vertex]++;
        region->temp_adj_curves[path->endpoints[FINISH].vertex][path->endpoints[START].vertex]++;
        region->dive[path->endpoints[START].vertex][path->endpoints[FINISH].vertex] = (Boolean) (face != path->endpoints[FINISH].face);
        region->dive[path->endpoints[FINISH].vertex][path->endpoints[START].vertex] = (Boolean) (face != path->endpoints[START].face);

        p_region->adj_cusp_triangle[face] = 0;
        p_region->temp_adj_curves[face][path->endpoints[START].vertex]++;
        p_region->temp_adj_curves[face][path->endpoints[FINISH].vertex]++;
        p_region->dive[face][path->endpoints[START].vertex]  = (Boolean) (face == path->endpoints[START].face);
        p_region->dive[face][path->endpoints[FINISH].vertex] = (Boolean) (face == path->endpoints[FINISH].face);

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
    p_region = g->regions[node->cusp_region_index];
    update_cusp_triangle_endpoints(&boundary->cusp_region_begin, &boundary->cusp_region_end,
                                   p_region, &path->endpoints[START], node, START);
    region = split_cusp_region_path_endpoint(p_region, node, &path->endpoints[START], index, START);
    INSERT_BEFORE(region, &boundary->cusp_region_end);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        p_region = g->regions[node->cusp_region_index];
        update_cusp_triangle_path_interior(&boundary->cusp_region_begin, &boundary->cusp_region_end, p_region, node);
        region = split_cusp_region_path_interior(p_region, node, index);
        INSERT_BEFORE(region, &boundary->cusp_region_end);
        index++;
    }

    // update last region
    p_region = g->regions[node->cusp_region_index];
    update_cusp_triangle_endpoints(&boundary->cusp_region_begin, &boundary->cusp_region_end, p_region, &path->endpoints[FINISH], node, FINISH);
    region = split_cusp_region_path_endpoint(p_region, node, &path->endpoints[FINISH], index, FINISH);
    INSERT_BEFORE(region, &boundary->cusp_region_end);
    index++;

    update_adj_region_data(&boundary->cusp_region_begin, &boundary->cusp_region_end);
    boundary->num_cusp_regions = index;
}

/*
 * Set the new and old region data. Draw a picture to see how
 * the attributes change in each case
 */

CuspRegion *split_cusp_region_path_interior(CuspRegion *region, PathNode *node, int index) {
    int v1, v2;
    CuspRegion *new_region = NEW_STRUCT( CuspRegion );

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
    new_region->dive[node->inside_vertex][v1]       = region->dive[node->inside_vertex][v1];
    new_region->dive[v2][v1]                        = region->dive[v2][v1];
    new_region->dive[node->inside_vertex][v2]       = region->dive[node->inside_vertex][v2];
    new_region->dive[v1][v2]                        = region->dive[v1][v2];

    // Update region
    region->curve[v2][v1]++;
    region->curve[v1][v2]++;
    region->dive[v2][v1]                            = 0;
    region->dive[node->inside_vertex][v1]           = 0;
    region->dive[v1][v2]                            = 0;
    region->dive[node->inside_vertex][v2]           = 0;
    region->adj_cusp_triangle[node->inside_vertex]  = 0;

    return new_region;
}

CuspRegion *split_cusp_region_path_endpoint(CuspRegion *region, PathNode *path_node, PathEndPoint *path_endpoint, int index, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;

    if (pos == START) {
        face = path_node->next_face;
    } else {
        face = path_node->prev_face;
    }

    if (face == path_endpoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        new_region->curve[path_endpoint->vertex][vertex2]++;
        new_region->temp_adj_curves[vertex1][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]      = (Boolean) (path_endpoint->face == vertex1);
        new_region->dive[vertex2][path_endpoint->vertex]      = region->dive[vertex2][path_endpoint->vertex];
        new_region->dive[vertex2][vertex1]                    = region->dive[vertex2][vertex1];
        new_region->dive[path_endpoint->vertex][vertex1]      = region->dive[path_endpoint->vertex][vertex1];
        new_region->adj_cusp_triangle[vertex1]                = 0;

        region->curve[path_endpoint->vertex][vertex1]++;
        region->temp_adj_curves[vertex2][path_endpoint->vertex]++;
        region->dive[vertex2][path_endpoint->vertex]         = (Boolean) (path_endpoint->face == vertex2);
        region->dive[vertex2][vertex1]                       = 0;
        region->dive[path_endpoint->vertex][vertex1]         = 0;
        region->adj_cusp_triangle[vertex2]                   = 0;
    } else if (face == path_endpoint->face) {
        // curve passes through the face that carries it
        new_region->curve[path_endpoint->face][path_endpoint->face == vertex1 ? vertex2 : vertex1]++;
        new_region->temp_adj_curves[face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex]++;
        new_region->dive[path_endpoint->face][path_endpoint->vertex]                         = region->dive[path_endpoint->face][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->vertex]                                 = 0;
        new_region->adj_cusp_triangle[path_endpoint->face == vertex1 ? vertex2 : vertex1]    = 0;

        region->curve[path_endpoint->face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face][path_endpoint->vertex] = 0;
    } else {
        // Curve goes around the vertex
        new_region->curve[face][path_endpoint->face]++;
        new_region->temp_adj_curves[path_endpoint->face][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]              = region->dive[vertex1][path_endpoint->vertex];
        new_region->dive[vertex2][path_endpoint->vertex]              = region->dive[vertex2][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->face]            = 0;
        new_region->adj_cusp_triangle[path_endpoint->vertex]          = 0;

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

    region2->tri            = region1->tri;
    region2->tet_index      = region1->tet_index;
    region2->tet_vertex     = region1->tet_vertex;

    for (i = 0; i < 4; i++) {
        region2->adj_cusp_triangle[i]   = region1->adj_cusp_triangle[i];
        region2->adj_cusp_regions[i]    = NULL;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]            = region1->curve[i][j];
            region2->dive[i][j]             = 0;
            region2->num_adj_curves[i][j]   = region1->num_adj_curves[i][j];
            region2->temp_adj_curves[i][j]  = region1->temp_adj_curves[i][j];
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
                                        CuspRegion *region, PathNode *node) {
    int face1, face2;
    CuspRegion *current_region;

    face1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    face2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    for (current_region = cusp_region_start->next; current_region != cusp_region_end; current_region = current_region->next) {
        // is the region initialised?
        if (current_region == NULL || current_region->tet_index == -1)
            continue;

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (current_region->curve[face1][node->inside_vertex] > region->curve[face1][node->inside_vertex]) {
            current_region->curve[face1][node->inside_vertex]++;
            current_region->dive[face1][node->inside_vertex] = 0;
        }
        else if (current_region->curve[face1][node->inside_vertex] < region->curve[face1][node->inside_vertex]) {
            current_region->curve[face1][face2]++;
        }

        if (current_region->curve[face2][node->inside_vertex] > region->curve[face2][node->inside_vertex]) {
            current_region->curve[face2][node->inside_vertex]++;
            current_region->dive[face2][node->inside_vertex] = 0;
        }
        else if (current_region->curve[face2][node->inside_vertex] < region->curve[face2][node->inside_vertex]) {
            current_region->curve[face2][face1]++;
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
                                    PathEndPoint *path_endpoint, PathNode *node, int pos) {
    FaceIndex face, face1, face2;
    CuspRegion *current_region;

    face1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    face2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (pos == START) {
        face = node->next_face;
    } else {
        face = node->prev_face;
    }

    for (current_region = cusp_region_start->next; current_region != cusp_region_end; current_region = current_region->next) {
        if (current_region == NULL || current_region->tet_index == -1)
            continue;

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (face == path_endpoint->vertex) {
            // curve passes through the face opposite the vertex it dives through
            if (current_region->curve[path_endpoint->vertex][face1] > region->curve[path_endpoint->vertex][face1]) {
                current_region->curve[face][face1]++;
                current_region->dive[face][face1] = 0;
                current_region->temp_adj_curves[face2][path_endpoint->vertex]++;

            } else if (current_region->curve[path_endpoint->vertex][face1] < region->curve[path_endpoint->vertex][face1]) {
                current_region->curve[face][face2]++;
                current_region->dive[face][face2] = 0;
                current_region->temp_adj_curves[face1][path_endpoint->vertex]++;

            }

            continue;
        }

        // Curve goes around the vertex or passes through the face that carries it
        if (current_region->curve[face][path_endpoint->vertex] > region->curve[face][path_endpoint->vertex]) {
            current_region->curve[face][path_endpoint->vertex]++;
            current_region->temp_adj_curves[face][path_endpoint->vertex]++;
            current_region->dive[face][path_endpoint->vertex] = 0;

        } else if (current_region->curve[face][path_endpoint->vertex] < region->curve[face][path_endpoint->vertex]) {
            current_region->curve[face][face == face1 ? face2 : face1]++;
            current_region->temp_adj_curves[face == face1 ? face2 : face1][path_endpoint->vertex]++;
            current_region->dive[face][face == face1 ? face2 : face1] = 0;
        }
    }
}

void update_adj_curve_along_path(ManifoldBoundary **cusps, OscillatingCurves *curves, DualCurves *dual_curve_begin, DualCurves *dual_curve_end, int curve_index) {
    int i, j;
    DualCurves *curve;
    PathEndPoint *path_endpoint;
    CuspRegion *region;
    ManifoldBoundary *cusp;

    // Update regions curve data
    for (curve = dual_curve_begin->next; curve != dual_curve_end; curve = curve->next) {
        // which cusp
        cusp = cusps[curve->cusp_index];

        for (region = cusp->cusp_region_begin.next; region != &cusp->cusp_region_end; region = region->next) {
            // which cusp region
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    region->num_adj_curves[i][j] += region->temp_adj_curves[i][j];
                    region->temp_adj_curves[i][j] = 0;
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

/*
 * Find a list of cusp regions which can dive through the manifold along the
 * edge with edge class 'edge_class'.
 */

CuspRegion **find_endpoint_regions(ManifoldBoundary *cusp, int edge_class, int *num_regions) {
    int i, v, face1, face2;
    CuspVertex *vertex;
    CuspRegion *region;
    CuspRegion **retval, **regions = NEW_ARRAY(cusp->num_cusp_regions, CuspRegion *);

    *num_regions = 0;

    for (region = cusp->cusp_region_begin.next; region != &cusp->cusp_region_end; region = region->next) {
        for (v = 0; v < 4; v++) {
            vertex = &region->tri->vertices[v];

            if (vertex->edge_class != edge_class)
                continue;

            face1 = (int) remaining_face[region->tet_vertex][v];
            face2 = (int) remaining_face[v][region->tet_vertex];

            if (!region->dive[face1][v] && !region->dive[face2][v])
                continue;

            regions[*num_regions] = region;
            (*num_regions)++;
        }
    }

    // reduce array size
    retval = NEW_ARRAY(*num_regions, CuspRegion *);
    for (i = 0; i < *num_regions; i++)
        retval[i] = regions[i];

    my_free(regions);
    return retval;
}

/*
 * Find the cusp region adjacent to 'region' by diving through the manifold
 * into vertex 'vertex' of the cusp region.
 */

CuspRegion *find_endpoint_adj_region(ManifoldBoundary *cusp, CuspRegion *region, int vertex) {
    int face1, face2, edge_class, edge_index;
    Boolean dive_along_face1, dive_along_face2, curve_face1, curve_face2, endpoint_face1, endpoint_face2;
    CuspRegion *current_region;

    face1 = (int) remaining_face[region->tet_vertex][vertex];
    face2 = (int) remaining_face[vertex][region->tet_vertex];

    edge_class = region->tri->vertices[vertex].edge_class;
    edge_index = region->tri->vertices[vertex].edge_index;

    dive_along_face1 = region->dive[face1][vertex];
    dive_along_face2 = region->dive[face2][vertex];

    // which cusp region
    for (current_region = cusp->cusp_region_begin.next;
         current_region != &cusp->cusp_region_end;
         current_region = current_region->next) {

        if (current_region->tet_index != region->tet_index ||
            current_region->tet_vertex != vertex)
            continue;

        curve_face1 = (Boolean) (region->curve[face1][vertex] == current_region->curve[face1][region->tet_vertex]);
        curve_face2 = (Boolean) (region->curve[face2][vertex] == current_region->curve[face2][region->tet_vertex]);

        if (!curve_face1 || !curve_face2)
            continue;

        endpoint_face1 = is_valid_endpoint(current_region, edge_class, edge_index, face1, region->tet_vertex);
        endpoint_face2 = is_valid_endpoint(current_region, edge_class, edge_index, face2, region->tet_vertex);

        if ( ( dive_along_face1 && endpoint_face1 ) || ( dive_along_face2 && endpoint_face2 ) )
            return region;
    }

    // didn't find valid path endpoints
    return NULL;
}

void update_adj_dive(ManifoldBoundary **cusps, int num_cusps) {
    int i, v;
    ManifoldBoundary *cusp;
    CuspRegion *region;

    for (i = 0; i < num_cusps; i++) {
        cusp = cusps[i];

        for (region = cusp->cusp_region_begin.next; region != &cusp->cusp_region_end; region = region->next) {
            for (v = 0; v < 4; v++) {
                if (v == region->tet_vertex)
                    continue;

                region->adj_dive_regions[v] = find_endpoint_adj_region(cusps[region->tri->tet->cusp[v]->index], region, v);
            }
        }
    }
}

void update_path_holonomy(DualCurves *path, int edge_class) {
    PathNode *path_node;
    CuspTriangle *tri;

    for (path_node = path->curves_begin.next; path_node != &path->curves_end; path_node = path_node->next) {
        tri = path_node->tri;
        tri->tet->extra[edge_class].curve[tri->tet_vertex][path_node->next_face]++;
        tri->tet->extra[edge_class].curve[tri->tet_vertex][path_node->prev_face]--;
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

// -------------------------------------------------------

// End Multi Graph

EndMultiGraph *init_end_multi_graph(Triangulation *manifold) {
    int i, j;
    int *parent;
    EndMultiGraph *multi_graph = NEW_STRUCT( EndMultiGraph );
    
    multi_graph->num_cusps = manifold->num_cusps;
    multi_graph->num_edge_classes = manifold->num_tetrahedra;

    Graph *g = init_graph(multi_graph->num_cusps, FALSE);
    cusp_graph(manifold, g);
//    printf("Cusp Graph\n");
//    print_graph(g);

    parent = NEW_ARRAY(g->num_vertices, int);

    multi_graph->multi_graph = spanning_tree(g, 0, parent);
    color_graph(multi_graph->multi_graph);

    multi_graph->edges = find_end_multi_graph_edge_classes(multi_graph, manifold);
    multi_graph->e0 = find_same_color_edge(manifold, multi_graph, g);

    multi_graph->edge_classes = NEW_ARRAY(multi_graph->num_edge_classes, Boolean);
    for (i = 0; i < multi_graph->num_edge_classes; i++) {
        multi_graph->edge_classes[i] = FALSE;
    }

    for (i = 0; i < multi_graph->num_cusps; i++) { 
        for (j = 0; j < multi_graph->num_cusps; j++) { 
            if (multi_graph->edges[i][j] == -1)
                continue;

            multi_graph->edge_classes[multi_graph->edges[i][j]] = TRUE;
        }
    }

    free_graph(g);
    my_free(parent);
    return multi_graph;
}

void free_end_multi_graph(EndMultiGraph *multi_graph) {
    int i;

    free_graph(multi_graph->multi_graph);
    
    for (i = 0; i < multi_graph->num_cusps; i++)
        my_free(multi_graph->edges[i]);
    
    my_free(multi_graph->edges);
    my_free(multi_graph);
}

void cusp_graph(Triangulation *manifold, Graph *g) {
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

Graph *spanning_tree(Graph *graph1, int start, int *parent) {
    int i;

    Boolean *processed = NEW_ARRAY(graph1->num_vertices, Boolean);
    Boolean *discovered = NEW_ARRAY(graph1->num_vertices, Boolean);

    Graph *graph2 = init_graph(graph1->num_vertices, graph1->directed);

    // Find path using bfs
    init_search(graph1, processed, discovered, parent);
    bfs(graph1, start, processed, discovered, parent);

    for (i = 0; i < graph1->num_vertices; i++) {
        if (parent[i] == -1)
            continue;

        insert_edge(graph2, i, parent[i], graph2->directed);
    }

    my_free(processed);
    my_free(discovered);

    return graph2;
}

/*
 * Assign an edge class to each edge of the graph g and return an array of 
 * Booleans indicating if an edge class is in the graph.
 */

int **find_end_multi_graph_edge_classes(EndMultiGraph *multi_graph, Triangulation *manifold) {
    int i, j, edge_class, **cusps;
    EdgeNode *edge_node;
    Graph *g = multi_graph->multi_graph;

    cusps = NEW_ARRAY(multi_graph->num_cusps, int *);

    for (i = 0; i < multi_graph->num_cusps; i++) {
        cusps[i] = NEW_ARRAY(multi_graph->num_cusps, int);

        for (j = 0; j < multi_graph->num_cusps; j++)
            cusps[i][j] = -1;
    }

    for (i = 0; i < g->num_vertices; i++) {
        for (edge_node = g->edge_list_begin[i].next; edge_node != &g->edge_list_end[i]; edge_node = edge_node->next) {
            edge_class = find_edge_class(manifold, i, edge_node->y);
            cusps[i][edge_node->y] = edge_class; 
            cusps[edge_node->y][i] = edge_class; 
        }
    }
    
    return cusps;
}

/*
 * Find an edge class whose edge connects cusp1 and cusp2
 */

int find_edge_class(Triangulation *manifold, int cusp1, int cusp2) {
    int v1, v2;
    EdgeClass *edge;
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        for (v1 = 0; v1 < 4; v1++) {
            for (v2 = 0; v2 < 4; v2++) {
                if (v1 == v2)
                    continue;

                if (tet->cusp[v1]->index != cusp1 || tet->cusp[v2]->index != cusp2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[v1][v2]];
                return edge->index;
            }
        }
    }

    uFatalError("find_edge_class", "symplectic_basis");
    return 0;
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

/*
 * g1 is the colored spanning tree of g2, return the
 * edge class of the edge in g2 which connects
 * vertices in g1 of the same color
 */

int find_same_color_edge(Triangulation  *manifold, EndMultiGraph *multi_graph, Graph *g2) {
    int cusp;
    EdgeNode *node;
    Graph *g1 = multi_graph->multi_graph;

    for (cusp = 0; cusp < g2->num_vertices; cusp++) {
        for (node = g2->edge_list_begin[cusp].next; node != &g2->edge_list_end[cusp]; node = node->next) {
            if (g1->color[cusp] == g1->color[node->y] && multi_graph->edges[cusp][node->y] == -1) 
                // we found an edge
                return find_edge_class(manifold, cusp, node->y);
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

CuspEndPoint *find_multi_graph_path(EndMultiGraph *multi_graph, Triangulation *manifold, int edge_class, int *path_length) {
    Graph *g = multi_graph->multi_graph;
    Boolean *processed     = NEW_ARRAY(g->num_vertices, Boolean);
    Boolean *discovered    = NEW_ARRAY(g->num_vertices, Boolean);
    int *parent         = NEW_ARRAY(g->num_vertices, int);
    int start, end, startE0, endE0;
    EdgeNode *node_begin = NEW_STRUCT( EdgeNode ), *node_end = NEW_STRUCT( EdgeNode), *tempNode;
    CuspEndPoint *cusp_end_point;

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    find_edge_ends(g, manifold, edge_class, &start, &end);
    find_edge_ends(g, manifold, multi_graph->e0, &startE0, &endE0);

    init_search(g, processed, discovered, parent);
    bfs(g, start, processed, discovered, parent);

    *path_length = 0;
    *path_length = find_path_len(start, end, parent, *path_length);

    /*
     * Construct a path of even legnth through the end multi graph. 
     * As a convention the path stored in node_begin -> node_end is 
     * the cusps the path visits, the edge class of the edge is 
     * determined by the edge classes in the end multi graph and e0.
     */
    if (*path_length % 2 == 1) {
        find_path(start, end, parent, node_begin);
    } else {
        init_search(g, processed, discovered, parent);
        bfs(g, start, processed, discovered, parent);

        find_path(start, startE0, parent, node_begin);

        init_search(g, processed, discovered, parent);
        bfs(g, endE0, processed, discovered, parent);

        find_path(endE0, end, parent, node_end->prev);
    }

    cusp_end_point = graph_path_to_cusp_path(multi_graph, node_begin, node_end, edge_class);

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

/*
 * Converts the EdgeNode path through the cusps in the end multigraph to 
 * a CuspEndPoint path which contains the edge classes on each cusp. 
 * A CuspEndPoint corresponding to one section of an oscillating curve, and
 * constructing such a section for all CuspEndPoints gives the whole curve
 *
 * node_begin is a doubly linked list with NULL header and tailer nodes and
 * the return CuspEndPoint is a pointer to a singly linked list with a NULL 
 * next value for the last node. 
 */

CuspEndPoint *graph_path_to_cusp_path(EndMultiGraph *multi_graph, EdgeNode *node_begin, EdgeNode *node_end, int edge_class) {
    int cusp, prev_edge_class;
    EdgeNode *node;
    CuspEndPoint *start_endpoint, *endpoint;

    start_endpoint = NEW_STRUCT( CuspEndPoint );
    start_endpoint->next = NULL;
    endpoint = start_endpoint;

    prev_edge_class = edge_class;
    for (node = node_begin->next; node->next != node_end; node = node->next) {
        cusp = node->y;

        endpoint->next = NEW_STRUCT( CuspEndPoint );
        endpoint->next->next = NULL;
        endpoint->cusp_index = cusp;
        endpoint->edge_class[START] = prev_edge_class;
        endpoint->edge_class[FINISH] = multi_graph->edges[node->y][node->next->y];

        if (endpoint->edge_class[FINISH] == -1)
            endpoint->edge_class[FINISH] = multi_graph->e0;
        
        prev_edge_class = endpoint->edge_class[FINISH];
        endpoint = endpoint->next;
    }

    endpoint->cusp_index = node->y;
    endpoint->edge_class[START] = prev_edge_class;
    endpoint->edge_class[FINISH] = edge_class;
    return start_endpoint;
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

Boolean *edge_classes_on_cusp(EndMultiGraph *multi_graph, int cusp_index) {
    int i;
    EdgeNode *edge_node;
    Boolean *edge_classes = NEW_ARRAY(multi_graph->num_edge_classes, Boolean);

    for (i = 0; i < multi_graph->num_edge_classes; i++)
        edge_classes[i] = FALSE;

    for (edge_node = multi_graph->multi_graph->edge_list_begin[cusp_index].next;
         edge_node != &multi_graph->multi_graph->edge_list_end[cusp_index];
         edge_node = edge_node->next) {
        edge_classes[multi_graph->edges[cusp_index][edge_node->y]] = TRUE;
    }

    return edge_classes;
}

void print_graph(Graph *g) {
    int i;
    EdgeNode *edge;

    if (!debug)
        return;

    for (i = 0; i < g->num_vertices; i++) {
        printf("(Vertex: %d) ", i);
        edge = &g->edge_list_begin[i];
        while ((edge = edge->next)->next != NULL) {
            printf("%d ", edge->y);
        }
        printf("\n");
    }

    printf("--------------------------\n");
}
