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
#define AVERAGE_ARGS float ** BaryonField, int index, int * GridDimenision, BaryonFieldIndex Ind

class average_root{
    public:
    QuanMap list;
    float current_mean;
    int id;
    char * label;
    virtual void addup(AVERAGE_ARGS){
        this->current_mean += 1;
    }
};
typedef std::map<char *,average_root *> AvgQuanMapType;
EXTERN AvgQuanMapType AvgMap;
EXTERN AvgQuanMapType FirstMoments;
EXTERN AvgQuanMapType SecondMoments;
EXTERN average_root avg_root;
class average_density: public average_root{
    public:
    virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_density avg_density;

class square_density: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Dens][index]*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN square_density std_density;

class average_vx: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel1][index]*Ind.Volume;
    }
};
EXTERN average_vx avg_vx;

class average_vy: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel2][index]*Ind.Volume;
    }
};
EXTERN average_vy avg_vy;

class average_vz: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel3][index]*Ind.Volume;
    }
};
EXTERN average_vz avg_vz;


class average_px: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel1][index]*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_px avg_px;

class average_py: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel2][index]*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_py avg_py;

class average_pz: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel3][index]*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_pz avg_pz;

class average_ex: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean +=0.5* POW(BaryonField[Ind.Vel1][index],2)*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_ex avg_ex;

class average_ey: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean +=0.5* POW(BaryonField[Ind.Vel2][index],2)*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_ey avg_ey;

class average_ez: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean +=0.5* POW(BaryonField[Ind.Vel3][index],2)*BaryonField[Ind.Dens][index]*Ind.Volume;
    }
};
EXTERN average_ez avg_ez;


class square_vx: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel1][index]*BaryonField[Ind.Vel1][index]*Ind.Volume;
    }
};
EXTERN square_vx std_vx;


class square_vy: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel2][index]*BaryonField[Ind.Vel2][index]*Ind.Volume;
    }
};
EXTERN square_vy std_vy;


class square_vz: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.Vel3][index]*BaryonField[Ind.Vel3][index]*Ind.Volume;
    }
};
EXTERN square_vz std_vz;

class average_bx: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.B1][index]*Ind.Volume;
    }
};
EXTERN average_bx avg_bx;

class average_by: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.B2][index]*Ind.Volume;
    }
};
EXTERN average_by avg_by;

class average_bz: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.B3][index]*Ind.Volume;
    }
};
EXTERN average_bz avg_bz;

class square_bx: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.B1][index]*BaryonField[Ind.B1][index]*Ind.Volume;
    }
};
EXTERN square_bx std_bx;


class square_by: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.B2][index]*BaryonField[Ind.B2][index]*Ind.Volume;
    }
};
EXTERN square_by std_by;


class square_bz: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += BaryonField[Ind.B3][index]*BaryonField[Ind.B3][index]*Ind.Volume;
    }
};
EXTERN square_bz std_bz;



class total_volume: public average_root{
    public: virtual void addup(AVERAGE_ARGS){
        this->current_mean += Ind.Volume;
    }
};
EXTERN total_volume tot_vol;

class avg_scalar: public average_root{
    //this is for scalar quantities like cycle and time.
    public: virtual void addup(AVERAGE_ARGS){
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
