//
// ExtraFunction
// Generic grid member function for debugging.
//

#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

int grid::ExtraFunction(char * message){

    if ( ProcessorNumber == MyProcessorNumber )
        return SUCCESS;

//fprintf(stderr,"CLOWN %s %p\n",message,AccelerationField);
    fprintf(stderr,"CLOWN %s \n",message);

    return SUCCESS;
}
