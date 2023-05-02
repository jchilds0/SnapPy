/*
 *	unix_cusped_census_main.c
 *
 *	The main() in this file illustrates the use of the unix-style
 *	GetCuspedCensusManifold().  It reads all the manifolds of
 *	5 or fewer tetrahedra and prints their volumes.
 */

#include "SnapPea.h"
#include "unix_cusped_census.h"
#include "../addl_code/addl_code.h"
#include <stdio.h>

void printMatrix(int**, int, int);

int main(void) {
    int **eqns, num_rows, num_cols;
    Triangulation *theTriangulation;

    theTriangulation = GetCuspedCensusManifold("", 5, oriented_manifold, 4);

    if (theTriangulation != NULL) {
        eqns = get_symplectic_basis(theTriangulation, &num_rows, &num_cols);
        printMatrix(eqns, num_cols, num_rows);
        free_symplectic_basis(eqns, num_rows);
    } else
        printf("Couldn't read census manifold.\n");

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