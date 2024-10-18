/***********************************************************************
/
/  PROTOSUBGRID CLASS (Force grids to be less than the maximum side)
/
/  written by: dcollins
/  date:       fall 2024.
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int ProtoSubgrid::ForceSubgridEdge(int dim, int &NumberOfNewGrids,
        int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2])
{
    /* Error check */
    fprintf(stderr,"CLOWN forcer\n");

    if (dim >= GridRank) {
        ENZO_VFAIL("Passed dim(%"ISYM") > GridRank(%"ISYM")\n", dim, GridRank)
    }

    if (Signature[dim] == NULL) {
        ENZO_VFAIL("Signature %"ISYM" not yet computed.\n", dim)
    }

    /* Initialize */

    int i = 0;
    NumberOfNewGrids = 0;

    for ( i=0; i < GridDimension[dim]; i++){

        /* Look for the start of a new subgrid. */
        if ( i % ForceSubgridEdgeSize == 0) {
            if ( NumberOfNewGrids > 0 ){
                GridEnds[NumberOfNewGrids-1][1] = StartIndex[dim] + i-1;
            }
            GridEnds[NumberOfNewGrids++][0] = StartIndex[dim] + i;
        }
        if ( i == GridDimension[dim]-1){
                GridEnds[NumberOfNewGrids-1][1] = StartIndex[dim] + i;
        }
    }

    return SUCCESS;
}
