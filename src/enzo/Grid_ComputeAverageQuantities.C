
//#include <stdio.h>
#include <iostream>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "AverageQuantities.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
/*
enum quan 
        self.stuff['mach'].append( np.sqrt(reg.quantities['WeightedAverageQuantity']('mean_square_velocity','cell_volume').in_units('code_velocity**2').v ) )
        self.stuff['vx'].append(reg.quantities['WeightedAverageQuantity']('velocity_x','cell_volume').in_units('code_velocity').v)
        self.stuff['vy'].append(reg.quantities['WeightedAverageQuantity']('velocity_y','cell_volume').in_units('code_velocity').v)
        self.stuff['vz'].append(reg.quantities['WeightedAverageQuantity']('velocity_z','cell_volume').in_units('code_velocity').v)
        self.stuff['px'].append(reg.quantities['WeightedAverageQuantity']('momentum_x','cell_volume').in_units('code_density*code_velocity').v)
        self.stuff['py'].append(reg.quantities['WeightedAverageQuantity']('momentum_y','cell_volume').in_units('code_density*code_velocity').v)
        self.stuff['pz'].append(reg.quantities['WeightedAverageQuantity']('momentum_z','cell_volume').in_units('code_density*code_velocity').v)
        self.stuff['ex'].append(reg.quantities['WeightedAverageQuantity']('eng_x','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['ey'].append(reg.quantities['WeightedAverageQuantity']('eng_y','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['ez'].append(reg.quantities['WeightedAverageQuantity']('eng_z','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['ke_tot'].append(reg.quantities['WeightedAverageQuantity']('kinetic_energy','cell_volume').in_units('code_density*code_velocity**2').v)
        reg.set_field_parameter('bulk_velocity',ds.arr([self.stuff['vx'][-1],self.stuff['vy'][-1], self.stuff['vz'][-1]],'code_velocity'))
        self.stuff['ke_rel'].append(reg.quantities['WeightedAverageQuantity']('rel_kinetic_energy','cell_volume').in_units('code_density*code_velocity**2').v)
        if self.potential_written:
            self.stuff['grav_pot'].append(reg.quantities['WeightedAverageQuantity']('grav_pot','cell_volume').in_units('code_density*code_velocity**2').v)
        self.stuff['gas_work'].append(reg.quantities['WeightedAverageQuantity']('gas_work','cell_volume').in_units('code_density*code_velocity**2').v)

        if car.ds['HydroMethod'] in [4,6]:
            self.stuff['Bx'].append( (reg['Bx']*volume).sum()/total_volume)
            self.stuff['By'].append( (reg['By']*volume).sum()/total_volume)
            self.stuff['Bz'].append( (reg['Bz']*volume).sum()/total_volume)

            self.stuff['bx'].append( ( (reg['Bx']-self.stuff['Bx'][-1])*volume).sum()/total_volume)
            self.stuff['by'].append( ( (reg['By']-self.stuff['By'][-1])*volume).sum()/total_volume)
            self.stuff['bz'].append( ( (reg['Bz']-self.stuff['Bz'][-1])*volume).sum()/total_volume)

            self.stuff['bx2'].append( np.sqrt(( (reg['Bx']-self.stuff['Bx'][-1])**2*volume).sum()/total_volume) )
            self.stuff['by2'].append( np.sqrt(( (reg['By']-self.stuff['By'][-1])**2*volume).sum()/total_volume) )
            self.stuff['bz2'].append( np.sqrt(( (reg['Bz']-self.stuff['Bz'][-1])**2*volume).sum()/total_volume) )

            self.stuff['Bfield_strength'].append( (reg['magnetic_field_strength']*volume).sum()/total_volume)
            self.stuff['AlfvenSpeed'].append( (volume*reg['magnetic_field_strength']/np.sqrt(np.pi*4*reg['density']) ).sum()/total_volume)

            self.stuff['AlfMach'].append( self.stuff['mach'][-1]/self.stuff['AlfvenSpeed'][-1])
            self.stuff['beta'].append( 2*(self.stuff['AlfMach'][-1]/self.stuff['mach'][-1])**2 )
*/


int grid::ComputeAverageQuantities(){
    if ( ProcessorNumber != MyProcessorNumber ){
        return SUCCESS;
    }
    BaryonFieldIndex BF;
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
            TENum, B1Num, B2Num, B3Num);
    BF.Dens= DensNum;
    BF.GEN= GENum;
    BF.Vel1= Vel1Num;
    BF.Vel2= Vel2Num;
    BF.Vel3= Vel3Num;
    BF.TE= TENum;
    BF.B1= B1Num;
    BF.B2= B2Num; 
    BF.B3= B3Num;
    BF.Volume = CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0];

    int i,j,k, index, size=GridDimension[0]*GridDimension[1]*GridDimension[2];
    for( k=GridStartIndex[2];k<=GridEndIndex[2]; k++)
        for( j=GridStartIndex[1];j<=GridEndIndex[1]; j++)
            for( i=GridStartIndex[0];i<=GridEndIndex[0]; i++){
                index = i+GridDimension[0]*(j+GridDimension[1]*k);
                for ( AvgQuanMapType::iterator it1=AverageList.begin(); it1!= AverageList.end(); it1++){
                    it1->second->addup(BaryonField,index,GridDimension,BF);
                }
            }

    return SUCCESS;
}
