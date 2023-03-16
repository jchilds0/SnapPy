/*
 *	unix_cusped_census_main.c
 *
 *	The main() in this file illustrates the use of the unix-style
 *	GetCuspedCensusManifold().  It reads all the manifolds of
 *	5 or fewer tetrahedra and prints their volumes.
 */

#include "SnapPea.h"
#include "unix_cusped_census.h"
#include <stdio.h>


int main(void) {
    Triangulation *theTriangulation;

    theTriangulation = GetCuspedCensusManifold("", 5, oriented_manifold, 125);

    if (theTriangulation != NULL) {
        printf("%s  %d\n", get_triangulation_name(theTriangulation), get_symplectic_basis(theTriangulation));
        free_triangulation(theTriangulation);
    } else
        printf("Couldn't read census manifold.\n");

    return 0;
}
