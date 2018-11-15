
/***********************************************************************
/
/  COMPUTE Average Quantities
/
/  written by: collins.
/  date:       October 2018.
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <map>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "AverageQuantities.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
void SetupAverageQuantities(){
    if (  AvgMap.empty() ){
      AvgMap["root"]   =&avg_root;
      AvgMap["cycle"]  =&avg_cycle;
      AvgMap["time"]   =&avg_time;
      AvgMap["volume"] =&tot_vol;
      FirstMoments["density"]=&avg_density;
      SecondMoments["density"]=&std_density;
      for ( AvgQuanMapType::iterator it1=FirstMoments.begin(); it1!= FirstMoments.end(); it1++){
          char * label = new  char[MAX_LINE_LENGTH];
          sprintf( label, "%s_avg", it1->first);
          AvgMap[label] = it1->second;
      }
      for ( AvgQuanMapType::iterator it1=SecondMoments.begin(); it1!= SecondMoments.end(); it1++){
         char * label = new  char[MAX_LINE_LENGTH];
         sprintf( label, "%s_std", it1->first);
         AvgMap[label] = it1->second;
      }
    }
    for ( AvgQuanMapType::iterator it1=AvgMap.begin(); it1!= AvgMap.end(); it1++){
        it1->second->current_mean=0;
    }
}
void SerializeAverageQuantities( float * AllQuantities){
    int n=0;
    for ( AvgQuanMapType::iterator it1=AvgMap.begin(); it1!= AvgMap.end(); it1++, n++){
        AllQuantities[n] = it1->second->current_mean;
    }
}
void DeSerializeAverageQuantities( float * AllQuantities, int cycle){
    int n=0;
    for ( AvgQuanMapType::iterator it1=AvgMap.begin(); it1!= AvgMap.end(); it1++,n++){
        if ( MyProcessorNumber == ROOT_PROCESSOR ){
            it1->second->list[cycle] = AllQuantities[n];
        }

    }
}
float * DeSerializeQ( char * quan){
    float * out = new float [ AvgMap[quan]->list.size()];
    int n=0;
    for ( QuanMap::iterator it1=AvgMap[quan]->list.begin(); it1!= AvgMap[quan]->list.end(); it1++,n++){
        out[n] = it1->second;
    }
    return out;

}
float * DeSerialize( QuanMap * thisone){
    float * out = new float [ thisone->size()];
    int n=0;
    for ( QuanMap::iterator it1=thisone->begin(); it1!= thisone->end(); it1++,n++){
        out[n] = it1->second;
    }
    return out;

}
void AverageQuantityWrite( char * name ){
    if ( MyProcessorNumber != ROOT_PROCESSOR || AvgMap.empty() ){
        return;
    }

    hid_t       file_id, dset, dataspace_id;
    herr_t      status, h5_status, h5_error = -1;

    char filename[100];
    sprintf(filename, "%s.AverageQuantities.h5",name);

    //file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if( file_id == -1 ){
        fprintf(stderr,"ERROR IN ERROR: ignore previous warning.  Opening hdf5 file.\n");
        file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    hsize_t number[1];
    int n=0;
    float * output;
    char label[MAX_LINE_LENGTH];
    for ( AvgQuanMapType::iterator it1=AvgMap.begin(); it1!= AvgMap.end(); it1++,n++){
        //Create Dataspace 
        number[0] = it1->second->list.size();

        output = DeSerialize( &it1->second->list );
        dataspace_id=H5Screate_simple(1, number, NULL);
        sprintf(label,"%s",it1->first);
        dset = H5Dcreate(file_id, label, HDF5_PREC, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, output);


        delete [] output;
        status = H5Sclose(dataspace_id);
        status = H5Dclose(dset);
    }
    status = H5Fclose(file_id);
}
int CommunicationBroadcastValues(float *Values, int Number, int BroadcastProcessor);
int ComputeAverageQuantities( LevelHierarchyEntry *LevelArray[],
        TopGridData *MetaData,int level, int cycle,
        int * LevelCycleCount, int * LevelSubCycleCount){

    SetupAverageQuantities();
    if ( MyProcessorNumber == ROOT_PROCESSOR ){
        AvgMap["time"]->current_mean = MetaData->Time;
        AvgMap["cycle"]->current_mean = MetaData->CycleNumber;
    }
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
        Temp->GridData->ComputeAverageQuantities();
        Temp             = Temp->NextGridThisLevel;
    }

    float * AllQuantities = new float[ AvgMap.size()];

    SerializeAverageQuantities( AllQuantities);
    CommunicationSumValues(AllQuantities, AvgMap.size());
    int nCycle=MetaData->CycleNumber;
    DeSerializeAverageQuantities( AllQuantities, nCycle);
    //Store standard deviations, not variance
    float d2, db, drms;
    if ( MyProcessorNumber == ROOT_PROCESSOR ){
        for ( AvgQuanMapType::iterator it1=SecondMoments.begin(); it1!= SecondMoments.end(); it1++){
            d2 = SecondMoments[it1->first]->list[nCycle];
            db=  FirstMoments[it1->first]->list[nCycle];
            drms = POW(d2 - db*db,0.5);
            it1->second->list[nCycle] = drms;


        }
    }


    delete AllQuantities;


}
