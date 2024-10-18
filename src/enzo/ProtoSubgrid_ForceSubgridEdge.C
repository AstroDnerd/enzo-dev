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
 
  if (dim >= GridRank) {
    ENZO_VFAIL("Passed dim(%"ISYM") > GridRank(%"ISYM")\n", dim, GridRank)
  }
 
  if (Signature[dim] == NULL) {
    ENZO_VFAIL("Signature %"ISYM" not yet computed.\n", dim)
  }
 
  /* Initialize */
 
  int i = 0;
  NumberOfNewGrids = 0;
 
  /* Loop over signature. */
 
  //fprintf(stderr, "CLOWN grid dim %d\n", GridDimension[dim]);
  while (i < GridDimension[dim]) {
 
    /* Look for the start of a new subgrid. */
      //fprintf(stderr,"CLOWN wtf i mod %d mod %d %d\n", i, ForceSubgridEdgeSize, i%ForceSubgridEdgeSize);
 
    if ( i % ForceSubgridEdgeSize == 0) {
      GridEnds[NumberOfNewGrids][0] = StartIndex[dim] + i;
      i++;
 
      /* Now find the end of the subgrid. */
 
      while (i < GridDimension[dim] && i % ForceSubgridEdgeSize > 0 )
	i++;
      GridEnds[NumberOfNewGrids++][1] = StartIndex[dim] + i-1;
      //fprintf(stderr,"CLOWN Gridends[%d] = %d %d\n",NumberOfNewGrids, GridEnds[NumberOfNewGrids-1][0], GridEnds[NumberOfNewGrids-1][1]);

      if ( NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS ) {
        ENZO_VFAIL("PE %"ISYM" NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS in ProtoSubgrid_FindGridsByZeroSignature\n", MyProcessorNumber)

      }

    }
 
    /* Next zone in signature. */
 
    i++;
  }
 
  return SUCCESS;
}
