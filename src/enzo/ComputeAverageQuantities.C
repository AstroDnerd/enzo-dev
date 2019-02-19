
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
#include <string.h>
#include "EnzoTiming.h"
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
    if (  FullList.empty() ){
     // Scalars["root"]   =&avg_root;
     Scalars["cycle"]  =&avg_cycle;
     Scalars["time"]   =&avg_time;
     Scalars["volume"] =&tot_vol;
     Scalars["dt"] =&avg_dt;
      FirstMoments["density"]=&avg_density;
      //StringMapTest["density"]=&avg_density;
      FirstMoments["vx"]=&avg_vx;
      FirstMoments["vy"]=&avg_vy;
      FirstMoments["vz"]=&avg_vz;
      FirstMoments["px"]=&avg_px;
      FirstMoments["py"]=&avg_py;
      FirstMoments["pz"]=&avg_pz;
      FirstMoments["bx"]=&avg_bx;
      FirstMoments["by"]=&avg_by;
      FirstMoments["bz"]=&avg_bz;
      FirstMoments["ex"]=&avg_ex;
      FirstMoments["ey"]=&avg_ey;
      FirstMoments["ez"]=&avg_ez;
      FirstMoments["sq_alf_x"]=&sq_alf_x;
      FirstMoments["sq_alf_y"]=&sq_alf_y;
      FirstMoments["sq_alf_z"]=&sq_alf_z;
      SecondMoments["density"]=&std_density;
      SecondMoments["vx"]=&std_vx;
      SecondMoments["vy"]=&std_vy;
      SecondMoments["vz"]=&std_vz;
      SecondMoments["bx"]=&std_bx;
      SecondMoments["by"]=&std_by;
      SecondMoments["bz"]=&std_bz;
      for ( AvgQuanMapType::iterator it1=FirstMoments.begin(); it1!= FirstMoments.end(); it1++){
          std::string label = it1->first + "_avg";
          AverageList[label] = it1->second;
      }
      for ( AvgQuanMapType::iterator it1=SecondMoments.begin(); it1!= SecondMoments.end(); it1++){
          std::string label = it1->first + "_std";
         AverageList[label] = it1->second;
      }
      for ( AvgQuanMapType::iterator it1=AverageList.begin(); it1!= AverageList.end(); it1++){
          std::string label = it1->first;
         FullList[label] = it1->second;
      }
      for ( AvgQuanMapType::iterator it1=Scalars.begin(); it1!= Scalars.end(); it1++){
          std::string label = it1->first;
         FullList[label] = it1->second;
      }
    }
    for ( AvgQuanMapType::iterator it1=FullList.begin(); it1!= FullList.end(); it1++){
        it1->second->current_mean=0;
    }
}
void SerializeAverageQuantities( float * AllQuantities){
    int n=0;
    for ( AvgQuanMapType::iterator it1=AverageList.begin(); it1!= AverageList.end(); it1++, n++){
        AllQuantities[n] = (float) it1->second->current_mean;
    }
}
void DeSerializeAverageQuantities( float * AllQuantities, int cycle){
    int n=0;
    for ( AvgQuanMapType::iterator it1=AverageList.begin(); it1!= AverageList.end(); it1++,n++){
        if ( MyProcessorNumber == ROOT_PROCESSOR ){
            it1->second->list[cycle] = AllQuantities[n];
        }

    }
}
float * DeSerializeQ( std::string quan){
    float * out = new float [ FullList[quan]->list.size()];
    int n=0;
    for ( QuanMap::iterator it1=FullList[quan]->list.begin(); it1!= FullList[quan]->list.end(); it1++,n++){
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
    if ( MyProcessorNumber != ROOT_PROCESSOR || FullList.empty() ){
        return;
    }

    hid_t       file_id, dset, dataspace_id;
    herr_t      status, h5_status, h5_error = -1;

    char filename[100];
    sprintf(filename, "%s.AverageQuantities.h5",name);
    fprintf(stderr,"Write Average %s\n",filename);

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
    for ( AvgQuanMapType::iterator it1=FullList.begin(); it1!= FullList.end(); it1++,n++){
        //Create Dataspace 
        number[0] = it1->second->list.size();

        output = DeSerialize( &it1->second->list );
        dataspace_id=H5Screate_simple(1, number, NULL);
        //sprintf(label,"%s",it1->first);
        dset = H5Dcreate(file_id, it1->first.c_str(), HDF5_PREC, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dset, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, output);


        delete [] output;
        it1->second->list.clear(); //EMPTY THE RUNNING LIST.
        status = H5Sclose(dataspace_id);
        status = H5Dclose(dset);
    }
    status = H5Fclose(file_id);
}
int CommunicationBroadcastValues(float *Values, int Number, int BroadcastProcessor);
int ComputeAverageQuantities( LevelHierarchyEntry *LevelArray[],
        TopGridData *MetaData,int level, float dt){
    TIMER_START("AverageQuantities");
    SetupAverageQuantities();
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
        Temp->GridData->ComputeAverageQuantities();
        Temp             = Temp->NextGridThisLevel;
    }
    
    float * AllQuantities = new float[ AverageList.size()];
  
    SerializeAverageQuantities( AllQuantities);
    CommunicationSumValues(AllQuantities, AverageList.size());
    int nCycle=MetaData->CycleNumber;
    DeSerializeAverageQuantities( AllQuantities, nCycle);
    if ( MyProcessorNumber == ROOT_PROCESSOR ){
        Scalars["time"]->list[nCycle] = MetaData->Time;
        Scalars["cycle"]->list[nCycle] = MetaData->CycleNumber;
        Scalars["dt"]->list[nCycle] = dt;
    }
    float d2, db, drms;
    if ( MyProcessorNumber == ROOT_PROCESSOR ){
        for ( AvgQuanMapType::iterator it1=SecondMoments.begin(); it1!= SecondMoments.end(); it1++){
            //d2 = SecondMoments[it1->first]->list[nCycle];
            //db=  FirstMoments[it1->first]->list[nCycle];
            //drms = POW(d2 - db*db,0.5);
            //it1->second->list[nCycle] = drms;
            it1->second->square_to_std(FirstMoments[it1->first]->list[nCycle], nCycle);
        }
    }


    delete AllQuantities;
    TIMER_STOP("AverageQuantities");


}
