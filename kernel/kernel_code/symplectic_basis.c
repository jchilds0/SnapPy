/*
 *  Symplectic Basis
 *
 *  Computes the symplectic basis of a cusped 3-manifold
*/

#include "kernel.h"
#include "kernel_namespace.h"

int** get_symplectic_basis(Triangulation *manifold, int* dual_rows)
{
    int i, j, T;
    int** eqns;

    T = manifold -> num_tetrahedra;

    *dual_rows = 2;

    eqns = NEW_ARRAY(*dual_rows, int*);

    for (i = 0; i < *dual_rows; i ++) {
        eqns[i] = NEW_ARRAY(3*T, int);

        for (j = 0; j < 3 * T; j ++)
            eqns[i][j] = i + j;
    }

    return eqns;
}