//
// Created by joshu on 19/03/2023.
//

#include <stdio.h>
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
    int i;
    char *tests[] = {"Enqueue", "Dequeue", "Circular", "Resize", "Empty"};

    Triangulation *theTriangulation;
    Tetrahedron *tet0, *tet1;

    theTriangulation = GetCuspedCensusManifold("", 5, oriented_manifold, 4);
    tet0 = theTriangulation->tet_list_begin.next;
    tet1 = tet0->next;

    for (i = 0; i < 5; i++) {
        printf("    %s: ", tests[i]);

        if (i == 0) {
            int index[] = {0, 0, 0, 0};
            struct queue q;
            initialise_queue(&q, 10);

            // Enqueue
            enqueue(&q, tet0);
            enqueue(&q, tet0);
            enqueue(&q, tet1);
            enqueue(&q, tet0);

            index[0] = dequeue(&q)->index;
            index[1] = dequeue(&q)->index;
            index[2] = dequeue(&q)->index;
            index[3] = dequeue(&q)->index;

            if (index[0] == 0 && index[1] == 0 && index[2] == 1 && index[3] == 0) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 1) {
            int index[] = {0, 0, 0, 0};
            struct queue q;
            initialise_queue(&q, 10);

            // Dequeue
            enqueue(&q, tet0);
            enqueue(&q, tet0);
            enqueue(&q, tet1);
            enqueue(&q, tet0);

            index[0] = dequeue(&q)->index;
            index[1] = dequeue(&q)->index;
            index[2] = dequeue(&q)->index;
            index[3] = dequeue(&q)->index;

            if (index[0] == 0 && index[1] == 0 && index[2] == 1 && index[3] == 0) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 2) {
            int j, index[] = {0, 0, 0, 0};
            struct queue q;
            Tetrahedron *tet;
            initialise_queue(&q, 10);

            // Circular
            for (j = 0; j < 9; j++)
                enqueue(&q, tet0);

            for (j = 0; j < 8; j++)
                dequeue(&q);

            enqueue(&q, tet1);
            enqueue(&q, tet0);
            enqueue(&q, tet1);

            index[0] = dequeue(&q)->index;
            index[1] = dequeue(&q)->index;
            index[2] = dequeue(&q)->index;
            index[3] = dequeue(&q)->index;

            if (index[0] == 0 && index[1] == 1 && index[2] == 0 && index[3] == 1) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 3) {
            int index[] = {0, 0, 0, 0};
            struct queue q;
            initialise_queue(&q, 2);

            // Resize
            q = *enqueue(&q, tet1);
            q = *enqueue(&q, tet0);
            q = *enqueue(&q, tet1);
            q = *enqueue(&q, tet0);

            index[0] = dequeue(&q)->index;
            index[1] = dequeue(&q)->index;
            index[2] = dequeue(&q)->index;
            index[3] = dequeue(&q)->index;

            if (index[0] == 1 && index[1] == 0 && index[2] == 1 && index[3] == 0) {
                printf("Passed");
            } else {
                printf("Failed - Queue returned [%d, %d, %d, %d]", index[0], index[1], index[2], index[3]);
            }
        } else if (i == 4) {
            struct queue q, p;
            initialise_queue(&q, 3);
            initialise_queue(&p, 3);

            enqueue(&q, tet0);
            enqueue(&q, tet0);
            dequeue(&q);
            enqueue(&q, tet0);
            enqueue(&q, tet0);
            dequeue(&q);
            dequeue(&q);
            dequeue(&q);

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
    }
}

void testBreadthFirstSearch(void) {

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