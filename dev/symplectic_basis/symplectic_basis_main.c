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
void mergeRows(Triangulation *);

int main(void) {

    Triangulation *theTriangulation;

    theTriangulation = GetCuspedCensusManifold("", 5, oriented_manifold, 4);

    if (theTriangulation != NULL) {
        mergeRows(theTriangulation);
    } else
        printf("Couldn't read census manifold.\n");

    return 0;
}


void mergeRows(Triangulation *manifold) {
    int numRows, numCols;
    int** eqns;

    // Edge Equations
    eqns = get_gluing_equations(manifold, &numRows, &numCols);
    printMatrix(eqns, numCols, numRows);
    free_gluing_equations(eqns, numRows);

    // Dual Equations
    eqns = get_symplectic_basis(manifold);
    printMatrix(eqns, numCols, numRows);
    free_symplectic_basis(eqns, numRows);

    free_triangulation(manifold);
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