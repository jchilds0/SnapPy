//
// Created by joshu on 16/05/2023.
//

#ifndef SNAPPEA_SYMPLECTIC_BASIS_H
#define SNAPPEA_SYMPLECTIC_BASIS_H

#endif //SNAPPEA_SYMPLECTIC_BASIS_H

#include "triangulation.h"
#include <stdbool.h>


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
 * Graph
 */

struct EdgeNode {
    int                     y;                      /** cusp region index */
    int                     edgeClass;              /** edgeClass of the edge for edges in the end multi graph */
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
    int                     *color;                 /** color a tree bipartite */
    int                     nVertices;              /** number of vertices in the graph */
    int                     directed;               /** is the graph directed */
};

struct EndMultiGraph {
    int                     e0;                     /** edge connecting vertices of the same color */
    struct Graph            *multi_graph;           /** tree with extra edge of cusps */
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

struct extra {
    int                     curve[4][4];            /** oscillating curve holonomy for a cusp triangle */
};

struct PathEndPoint {
    int                     face;                   /** face containg the short rectangle carrying the curve */
    int                     vertex;                 /** vertex we dive through the manifold along */
    int                     regionIndex;            /** index of the region the endpoint lies in */
    struct CuspRegion       *region;                /** pointer to the region the endpoint lies in */
};

struct DualCurves {
    int                     edgeClass[2];
    struct EdgeNode         curves_begin;           /** header node of doubbly linked list */
    struct EdgeNode         curves_end;             /** tailer node of doubbly linked list */
    struct PathEndPoint     endpoints[2];           /** path end points */
    struct DualCurves       *next;                  /** next dual curve in doubly linked list */
    struct DualCurves       *prev;                  /** prev dual curve in doubly linked list */
};

struct OscillatingCurves {
    int                     numCurves;
    int                     *edgeClass;
    struct DualCurves       *dual_curve_begin;      /** array of doubly linked lists of dual curves */
    struct DualCurves       *dual_curve_end;        /** array of doubly linkek lists of dual curves */
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
    Triangulation           *manifold;              /** manifold */
    Cusp                    *cusp;                  /** which cusp is the boundary in */
    struct Graph            *dual_graph;            /** dual graph of the cusp region */
    struct CuspTriangle     cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    struct CuspTriangle     cusp_triangle_end;      /** tail node of doubly linked list of cusp triangles */
    struct CuspRegion       cusp_region_begin;      /** header node of doubly linked list of cusp regions */
    struct CuspRegion       cusp_region_end;        /** tail node of doubly linked list of cusp regions */
};

// Graph
struct Graph            *init_graph(int, bool);
void                    free_graph(struct Graph *);
int                     insert_edge(struct Graph *, int, int, bool);
void                    delete_edge(struct Graph *, int, int, bool);
int                     edge_exists(struct Graph *, int, int);

// Symplectic Basis
int                     **get_symplectic_equations(Triangulation *, int *, int *, int);

struct ManifoldBoundary *init_boundary(Triangulation *, Cusp *);
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
struct OscillatingCurves*init_oscillating_curves(Triangulation *, int *);
void                    free_oscillating_curves(struct OscillatingCurves *);
void                    find_intersection_triangle(Triangulation *, struct ManifoldBoundary *);

/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    do_oscillating_curves(struct ManifoldBoundary **, struct OscillatingCurves *, struct EndMultiGraph *);
void                    do_one_dual_curve(struct ManifoldBoundary **, struct DualCurves *, struct DualCurves *, struct EndMultiGraph *, int);
void                    do_one_cusp(struct ManifoldBoundary *, struct DualCurves *, struct DualCurves *, int);
struct Graph *          construct_cusp_region_dual_graph(struct ManifoldBoundary *);
void                    print_debug_info(Triangulation *, struct ManifoldBoundary **, struct OscillatingCurves *, int);
void                    find_path_endpoints(struct Graph *, struct PathEndPoint *, struct PathEndPoint *, int, int, bool);
void                    update_path_info(struct Graph *g, struct DualCurves *, int);
void                    update_path_endpoint_info(struct CuspRegion *, struct EdgeNode *, struct PathEndPoint *, int, int);
void                    split_cusp_regions_along_path(struct ManifoldBoundary *, struct DualCurves *);
struct CuspRegion       *update_cusp_region(struct CuspRegion *region, struct EdgeNode *, struct PathEndPoint *, int, int);
void                    update_cusp_triangle(struct CuspRegion *, struct CuspRegion *, struct CuspRegion *, struct EdgeNode *);
void                    update_cusp_triangle_endpoints(struct CuspRegion *, struct CuspRegion *, struct CuspRegion *, struct PathEndPoint *, struct EdgeNode *, int);
void                    copy_region(struct CuspRegion *, struct CuspRegion *);
void                    calculate_holonomy(Triangulation *, int **, int);
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

void                    init_search(struct Graph *, bool *, bool *, int *, int *);
void                    bfs(struct Graph *, int, bool *, bool *, int *, int *);
void                    find_path(int, int, int *, struct EdgeNode *);


/**
 * Spanning Tree for End Multi Graph
 */

struct EndMultiGraph    *init_end_multi_graph(Triangulation *, int *);
void                    free_end_multi_graph(struct EndMultiGraph *);
int                     insert_edge_end_multi_graph(struct Graph *, int, int, int, bool);
void                    spanning_tree(struct Graph *, struct Graph *, int, int *);
void                    cusp_graph(Triangulation *, struct Graph *);
void                    color_graph(struct Graph *);
int                     *find_tree_edges(struct Graph *, int);
int                     find_same_color_edge(struct Graph *, struct Graph *, int *);
int                     find_path_len(int, int, int *, int);
void                    find_multi_graph_path(struct Graph *, int, int, struct EdgeNode *);
void                    print_graph(struct Graph *, int);
