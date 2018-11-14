#ifndef AVG_QUANTITY_
#define AVG_QUANTITY_
#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif


struct BaryonFieldIndex{
    int Dens, GEN, Vel1, Vel2, Vel3, TE, B1, B2, B3;
    FLOAT Volume;
};

typedef std::map<int , float> QuanMap;

class average_root{
    public:
    QuanMap list;
    float current_mean;
    int id;
    char * label;
    virtual void addup(float ** BaryonField, int index, int * GridDimenision,
                        BaryonFieldIndex Ind){
        this->current_mean += 1;
    }
};
typedef std::map<char *,average_root *> AvgQuanMapType;
EXTERN AvgQuanMapType AvgMap;
EXTERN average_root avg_root;

class average_density: public average_root{
    public:
    virtual void addup(float ** BaryonField, int index, int * GridDimenision,
                        BaryonFieldIndex Ind){
        this->current_mean += BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_density avg_density;
class square_density: public average_root{
    public:
    virtual void addup(float ** BaryonField, int index, int * GridDimenision,
                        BaryonFieldIndex Ind){
        this->current_mean += BaryonField[Ind.Dens][index]*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN square_density std_density;

class total_volume: public average_root{
    public:
    virtual void addup(float ** BaryonField, int index, int * GridDimenision, BaryonFieldIndex Ind){
        this->current_mean += Ind.Volume;
    }
};
EXTERN total_volume tot_vol;

class avg_scalar: public average_root{
    //this is for scalar quantities like cycle and time.
    public:
    virtual void addup(float ** BaryonField, int index, int * GridDimenision, BaryonFieldIndex Ind){
        //doesn't do anything.
    }
};
EXTERN avg_scalar avg_time;
EXTERN avg_scalar avg_cycle;



void SetupAverageQuantities();
void SerializeAverageQuantities( float * AllQuantities);
void DeSerializeAverageQuantities( float * AllQuantities);
  // thistest
#endif //AVG_QUANTITY_
