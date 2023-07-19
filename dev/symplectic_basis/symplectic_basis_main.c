/*
 *
 */

#include "SnapPea.h"
#include "unix_cusped_census.h"
#include "unix_file_io.h"
#include "addl_code.h"
#include <stdio.h>

void printMatrix(int**, int, int);

int main(void) {
    int i, **eqns, num_rows, num_cols;
    Triangulation *theTriangulation;

    int fromFile = 0;

    int count = 1;
    int numTet[] = {6, 7, 7, 7, 7, 7};
    int index[] = {443, 2208, 2652, 2942, 3140, 3507};

    char *error[] = {"CuspedCensusData/1.tri",
                     "CuspedCensusData/4.tri",
                     "CuspedCensusData/35.tri",
                     "CuspedCensusData/76.tri",
                     "CuspedCensusData/7703.tri",
                     "CuspedCensusData/7855.tri",
                     "CuspedCensusData/12442.tri"};

    char *link[] = {"CuspedCensusData/link-0.tri", 
                    "CuspedCensusData/link-1.tri",
                    "CuspedCensusData/link-2.tri",
                    "CuspedCensusData/link-3.tri",
                    "CuspedCensusData/link-4.tri"
    };


    for (i = 0; i < count; i++) {
        if (fromFile == 1) {
            printf("Triangulation: %s\n", error[i]);
            theTriangulation = read_triangulation(error[i]);
        } else if (fromFile == 2) {
            printf("Triangulation: %s\n", link[i]);
            theTriangulation = read_triangulation(link[i]);
        } else {
            theTriangulation = GetCuspedCensusManifold("", numTet[i], oriented_manifold, index[i]);
            printf("Num Tet: %d Index: %d\n", numTet[i], index[i]);
        }

        if (theTriangulation != NULL) {
            if (get_orientability(theTriangulation) == nonorientable_manifold)
                continue;

            eqns = get_symplectic_basis(theTriangulation, &num_rows, &num_cols);
            printMatrix(eqns, num_cols, num_rows);
            printf("---------------------------\n");
            free_symplectic_basis(eqns, num_rows);
            free_triangulation(theTriangulation);
        } else
            printf("Couldn't read census manifold.\n");
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
