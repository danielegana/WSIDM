
// ************************* INTEGRALS AND INTERPOLATIONS **********************

extern double Vmax;
extern double Varr[npointsv];
extern double Vint;
extern double deltaDMpref[npointsv];
extern double vDMpref[npointsv];
extern double viscpref[npointsv];
extern double mrel;




double deltaDMintegral2(struct inputRKF *inputRKF){
    
    double deltaDM=0;
    int n;
    
    for (n=0; n<npointsv; n++) {
        
        deltaDM=deltaDM+(*inputRKF).y.PsideltaDM[n];
        
    }
    
    return deltaDM;
}


double vDMintegral2(struct inputRKF *inputRKF){
    
    double vDM=0;;
    
    int n;
    
    for (n=0; n<npointsv; n++) {
        
        vDM=vDM+(*inputRKF).y.PsivDM[n];
        
    }
    
    return vDM;
}


double viscintegral2(struct inputRKF *inputRKF){
    
    double visc=0;
        
    int n;
    
    for (n=0; n<npointsv; n++) {
        
        visc=visc+(*inputRKF).y.Psivisc[n];
        
    }
    
    return visc;
}

double rhointegral(double fm4arr[npointsv],struct paramsstepstruct *paramsstep){
    
    double rho=0;
        
    int n;
    
    for (n=0; n<npointsv; n++) {
        
        if((*paramsstep).a<anr){
            rho=rho+4*M_PI/mrel*(fm4arr[n])*square(Varr[n])*Varr[n]*mrel/pow((*paramsstep).a,4.);}
        else{
            rho=rho+4*M_PI/mrel*(fm4arr[n])*square(Varr[n])*(*paramsstep).a*mrel/pow((*paramsstep).a,4.);
            
        }

        }
        
    
    return rho;
}



