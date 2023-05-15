/*
 *
 */

#include "SnapPea.h"
#include "unix_cusped_census.h"
#include "../addl_code/addl_code.h"
#include <stdio.h>

void printMatrix(int**, int, int);

int main(void) {
    int i, j, **eqns, num_rows, num_cols;
    Triangulation *theTriangulation;

    int count = 1;
    int numTet[] = {7};
    int index[][2] = {{1, 10}};

    for (i = 0; i < count; i++) {
        for (j = index[i][0]; j < index[i][1]; j++) {
            theTriangulation = GetCuspedCensusManifold("", numTet[i], oriented_manifold, j);

            if (get_orientability(theTriangulation) == nonorientable_manifold)
                continue;

            if (get_num_cusps(theTriangulation) != 1)
                continue;

            printf("Num Tet: %d Index: %d\n", numTet[i], j);
            if (theTriangulation != NULL) {
                eqns = get_symplectic_basis(theTriangulation, &num_rows, &num_cols);
                printMatrix(eqns, num_cols, num_rows);
                free_symplectic_basis(eqns, num_rows);
            } else
                printf("Couldn't read census manifold.\n");
        }
    }

    return 0;
}


void printMatrix(int** M, int numCols, int numRows) {
    int i, j;

    for (i = 0; i < numRows; i ++) {
        printf("[");

        for (j = 0; j < numCols - 1; j ++) {
            printf(" %d,", M[i][j]);
        }

        printf(" %d]\n", M[i][numCols - 1]);
    }
}