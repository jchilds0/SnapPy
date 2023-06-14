//
// Created by joshu on 19/03/2023.
//

#include <stdio.h>
#include "kernel.h"
#include "kernel_namespace.h"
#include "unix_cusped_census.h"
#include "addl_code.h"
#include "symplectic_basis_old.h"

void testDual(void);
int omega(int *, int *, int);
void printMatrix(int **, int, int);

int main() {
    printf("Testing Symplectic Basis: \n");
    testDual();
}

void testDual(void) {
    int **basis, dual_rows, dual_cols, i, j, k;
    Triangulation *theTriangulation;

    int failed[] = {0, 0, 0};
    int count[] = {0, 0, 0};
    int index[][3] = {{5, 0, 110},
                      {6, 0, 950},
                      {7, 0, 3550}
    };
    int **passed = NEW_ARRAY(3, int *);
    passed[0] = NEW_ARRAY(index[0][2] - index[0][1], int);
    passed[1] = NEW_ARRAY(index[1][2] - index[1][1], int);
    passed[2] = NEW_ARRAY(index[2][2] - index[2][1], int);

    for (i = 0; i < 3; i++)
        for (j = 0; j < index[i][2] - index[i][1]; j++)
            passed[i][j] = 1;

    for (i = 0; i < 3; i++) {
        for (j = index[i][1]; j < index[i][2]; j++) {
            theTriangulation = GetCuspedCensusManifold("", index[i][0], oriented_manifold, j);

            if (get_orientability(theTriangulation) == nonorientable_manifold)
                continue;

            if (j == 2208 || j == 2652 || j == 2942 || j == 3140 || j == 3507)
                continue;

            printf("Num Tet: %d Index: %d \n", index[i][0], j);

            basis = get_symplectic_basis(theTriangulation, &dual_rows, &dual_cols);

            for (k = 0; k < dual_rows / 2; k ++) {
                if (ABS(omega(basis[2 * k], basis[2 * k + 1], dual_cols)) == 2) {
                    printf("Passed\n");
                    continue;
                }

                printf("Failed Num Tet: %d Index: %d \n", index[i][0], j);
                printMatrix(basis, dual_cols, dual_rows);
                passed[i][j - index[i][1]] = 0;
                failed[i]++;
                break;
            }
            count[i]++;

            free_triangulation(theTriangulation);
            free_symplectic_basis(basis, dual_rows);
        }
    }

    for (i = 0; i < 3; i++) {
        printf("    (Num. of Tet %d) Failed: %d out of %d tests\n", index[i][0], failed[i], count[i]);
    }
}

int omega(int *v1, int *v2, int numCols) {
    int i, yyval = 0;

    for (i = 0; i < numCols / 3; i++) {
        yyval += ((v1[3 * i] - v1[3 * i + 2]) * (v2[3 * i + 1] - v2[3 * i + 2])
                  - (v1[3 * i + 1] - v1[3 * i + 2]) * (v2[3 * i] - v2[3 * i + 2]));
    }

    return yyval;
}

void printMatrix(int** M, int numCols, int numRows) {
    int i, j;

    for (i = 0; i < numRows; i++) {
        printf("[");

        for (j = 0; j < numCols - 1; j++) {
            printf(" %d,", M[i][j]);
        }

        printf(" %d]\n", M[i][numCols - 1]);
    }
}
