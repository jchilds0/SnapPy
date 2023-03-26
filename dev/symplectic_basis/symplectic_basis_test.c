//
// Created by joshu on 19/03/2023.
//

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "kernel.h"
#include "kernel_namespace.h"
#include "SnapPea.h"
#include "symplectic_basis.h"
#include "unix_cusped_census.h"
#include "../addl_code/addl_code.h"

void testQueue(void);
void testBreadthFirstSearch(void);
void testDual(void);

int main() {
    printf("Testing queue\n");
    testQueue();

    printf("Testing breadth first search: \n");
    testBreadthFirstSearch();

    printf("Testing Symplectic Basis: \n");
    testDual();
}

void testQueue(void) {
    int i, j;
    char *tests[] = {"Enqueue", "Dequeue", "Circular", "Resize", "Empty"};

    struct queue q, p;

    for (i = 0; i < 5; i++) {
        printf("    %s: ", tests[i]);

        initialise_queue(&q, 10);
        initialise_queue(&p, 3);

        int index[] = {0, 0, 0, 0};

        if (i == 0) {
            // Enqueue
            enqueue(&q, 0);
            enqueue(&q, 0);
            enqueue(&q, 1);
            enqueue(&q, 0);

            index[0] = dequeue(&q);
            index[1] = dequeue(&q);
            index[2] = dequeue(&q);
            index[3] = dequeue(&q);

            if (index[0] == 0 && index[1] == 0 && index[2] == 1 && index[3] == 0) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 1) {
            // Dequeue
            enqueue(&q, 0);
            enqueue(&q, 0);
            enqueue(&q, 1);
            enqueue(&q, 0);

            index[0] = dequeue(&q);
            index[1] = dequeue(&q);
            index[2] = dequeue(&q);
            index[3] = dequeue(&q);

            if (index[0] == 0 && index[1] == 0 && index[2] == 1 && index[3] == 0) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 2) {
            // Circular
            for (j = 0; j < 9; j++)
                enqueue(&q, 0);

            for (j = 0; j < 8; j++)
                dequeue(&q);

            enqueue(&q, 1);
            enqueue(&q, 0);
            enqueue(&q, 1);

            index[0] = dequeue(&q);
            index[1] = dequeue(&q);
            index[2] = dequeue(&q);
            index[3] = dequeue(&q);

            if (index[0] == 0 && index[1] == 1 && index[2] == 0 && index[3] == 1) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 3) {
            // Resize
            p = *enqueue(&p, 1);
            p = *enqueue(&p, 0);
            p = *enqueue(&p, 1);
            p = *enqueue(&p, 0);

            index[0] = dequeue(&p);
            index[1] = dequeue(&p);
            index[2] = dequeue(&p);
            index[3] = dequeue(&p);

            if (index[0] == 1 && index[1] == 0 && index[2] == 1 && index[3] == 0) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 4) {
            enqueue(&p, 0);
            enqueue(&p, 0);
            dequeue(&p);
            enqueue(&p, 0);
            enqueue(&p, 0);
            dequeue(&p);
            dequeue(&p);
            dequeue(&p);

            if (empty_queue(&q) && empty_queue(&p)) {
                printf("Passed");
            } else if (empty_queue(&q)){
                printf("Failed - Queue p not empty");
            } else if (empty_queue(&p)){
                printf("Failed - Queue q not empty");
            } else {
                printf("Failed");
            }
        }
        printf("\n");

        free_queue(&p);
        free_queue(&q);
    }
}

void testBreadthFirstSearch(void) {
    int i, j, nvertices = 4;
    graph g;

    initialise_graph(&g, maxVertices, 10, TRUE);

    bool *processed = malloc(sizeof(bool *) * maxVertices);
    bool *discovered = malloc(sizeof(bool *) * maxVertices);
    int *parent = malloc(sizeof(int *) * maxVertices);

    int *path = malloc(sizeof(int *) * maxVertices);

    int matrix[4][4] = {
            {0, 1, 0, 1},
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {1, 0, 0, 0}
    };

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            if (matrix[i][j] == 1)
                insert_edge(&g, i, j, TRUE);

    initialise_search(&g, processed, discovered, parent);
    bfs(&g, 1, processed, discovered, parent);

    find_path(1, 3, parent, path, 0);
    if (path[0] == 3 && path[1] == 0 && path[2] == 1) {
        printf("    Passed\n");
    } else {
        printf("    Failed - Path returned [%d, %d, %d, %d]\n", path[0], path[1], path[2], path[3]);
    }

    free_graph(&g);
    free(processed);
    free(discovered);
    free(parent);
    free(path);
}

void testDual(void) {
    int eqns[2][6] = {
            {0, 1, 2, 3, 4, 5},
            {1, 2, 3, 4, 5, 6}
    };
    int **basis, dual_rows, i, j;

    Triangulation *theTriangulation;

    theTriangulation = GetCuspedCensusManifold("", 5, oriented_manifold, 4);

    basis = get_symplectic_basis(theTriangulation, &dual_rows);

    for (i = 0; i < dual_rows; i ++) {
        for (j = 0; j < 3 * (get_num_tetrahedra(theTriangulation)); j ++) {
            if (basis[i][j] != eqns[i][j]) {

            }
        }
    }
}