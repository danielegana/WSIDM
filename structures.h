// ****************** STRUCTURE DEFINITIONS *********

#define npointsv 30
#define nummom 30
#define nummomnu 30


//Structure of perturbations
struct pertstruct{
    
    double phi,deltarad,viscrad,vbrad,deltaDM,vDM,viscDM,deltab,viscnu,psi,vrad,vb,rhoPvisctot; //visc here is sigma in Bertschinger's paper. the velocities are not theta, are theta/k^2.
    double Psi[nummom][npointsv]; // // The last moment, Psi[nummom-1] is defined as Psi*qoverepsilon
    double PsideltaDM[npointsv]; //Define also special vectors for the first moments. Useful for integral interpolation/
    double PsivDM[npointsv];
    double Psivisc[npointsv];
    double Fnu[nummomnu];

    
};


//Structure of derivatives of perturbations
struct dstruct{
    
    double phi,deltarad,vbrad,viscrad,deltaDM,vDM,deltab,vb,vrad;
    double Psi[nummom-1][npointsv];
    double Fnu[nummomnu];

    
};

// Super-structure given as an input to RKF routine.
struct inputRKF{
    
    struct pertstruct pert;
    double s;
    struct pertstruct y;
    struct pertstruct yp1;
    struct pertstruct zp1;
    struct dstruct k1;
    struct dstruct k2;
    struct dstruct k3;
    struct dstruct k4;
    struct dstruct k5;
    struct dstruct k6;
    struct dstruct ksum;
    int rkfcounter;
    
    double RKFalg[7][6];

};


// Structure for background and other parameters

struct paramsstepstruct{
    
    int decoupl;
    int wdmflag;
    int recomb;
    int irec;
    double k;
    double deta;
    double eta;
    double a;
    double Heta;
    double tol;

    
    double rhoDMa;
    double rhoba;
    double rhorada;
    double rhonua;
    
    double TM;
    double dlogTM;
    double csbsq;
    
    double usDMsq;
    double wDM;
    double wDMp;
    double usbradsq;
    double RB;
    
    double neSigmaTcomov;
    double neSigmaTa;

    
    int n;
    int nmom;
    int nmomnu;
        
};

//Function to sum derivatives for RKF algorithm
void sumstruct( struct inputRKF *inputRKF, struct paramsstepstruct *paramsstep){
        
        (*inputRKF).ksum.phi=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.phi+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.phi+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.phi+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.phi+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.phi+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.phi;

        (*inputRKF).ksum.deltarad=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.deltarad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.deltarad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.deltarad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.deltarad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.deltarad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.deltarad;

        (*inputRKF).ksum.deltab=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.deltab+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.deltab+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.deltab+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.deltab+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.deltab+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.deltab;
    

     //Neutrinos. Sum only to nummom-1 because derivatives are only calculated up to that moment.
    for ((*paramsstep).nmomnu=0; (*paramsstep).nmomnu<nummomnu-1; (*paramsstep).nmomnu++) {
        
        (*inputRKF).ksum.Fnu[(*paramsstep).nmomnu]=
          (*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.Fnu[(*paramsstep).nmomnu]+
          (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.Fnu[(*paramsstep).nmomnu]+
          (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.Fnu[(*paramsstep).nmomnu]+
          (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.Fnu[(*paramsstep).nmomnu]+
          (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.Fnu[(*paramsstep).nmomnu]+
          (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.Fnu[(*paramsstep).nmomnu];
        
    }
    
    
    
    if ((*paramsstep).recomb==0) {
        
            (*inputRKF).ksum.vbrad=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.vbrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.vbrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.vbrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.vbrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.vbrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.vbrad;
        
        (*inputRKF).ksum.viscrad=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.viscrad+
                                  (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.viscrad+
                                  (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.viscrad+
                                  (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.viscrad+
                                  (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.viscrad+
                                  (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.viscrad;


        
    } else{
        
          (*inputRKF).ksum.vb=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.vb+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.vb+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.vb+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.vb+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.vb+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.vb;

                
            (*inputRKF).ksum.vrad=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.vrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.vrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.vrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.vrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.vrad+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.vrad;

        
    }
    
    if ((*paramsstep).decoupl==0) {
        
        
        (*inputRKF).ksum.deltaDM=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.deltaDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.deltaDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.deltaDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.deltaDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.deltaDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.deltaDM;


        (*inputRKF).ksum.vDM=(*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.vDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.vDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.vDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.vDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.vDM+
                            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.vDM;
        
    } else{

    for ((*paramsstep).nmom=0; (*paramsstep).nmom<nummom-1; (*paramsstep).nmom++) {
        for ((*paramsstep).n=0; (*paramsstep).n<npointsv; (*paramsstep).n++) {
        
            (*inputRKF).ksum.Psi[(*paramsstep).nmom][(*paramsstep).n]=
            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][0]*(*inputRKF).k1.Psi[(*paramsstep).nmom][(*paramsstep).n]+
            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][1]*(*inputRKF).k2.Psi[(*paramsstep).nmom][(*paramsstep).n]+
            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][2]*(*inputRKF).k3.Psi[(*paramsstep).nmom][(*paramsstep).n]+
            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][3]*(*inputRKF).k4.Psi[(*paramsstep).nmom][(*paramsstep).n]+
            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][4]*(*inputRKF).k5.Psi[(*paramsstep).nmom][(*paramsstep).n]+
            (*inputRKF).RKFalg[(*inputRKF).rkfcounter][5]*(*inputRKF).k6.Psi[(*paramsstep).nmom][(*paramsstep).n];

        }
    }
    }

};
