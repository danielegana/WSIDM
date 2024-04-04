extern double Vint;
extern double Vmax;
extern double anr;

extern double Varr[npointsv];
extern double stepfactorlate;

extern double deltaDMprefNR[npointsv];
extern double vDMprefNR[npointsv];
extern double viscprefNR[npointsv];

extern double deltaDMprefwdm[npointsv];
extern double vDMprefwdm[npointsv];
extern double viscprefwdm[npointsv];

extern double dlogfNRarr[npointsv];
extern double dlogfoverVNRarr[npointsv];
extern double dlogfrelarr[npointsv];
extern double dlogfoverVrelarr[npointsv];

extern double mrel;
extern double T0r;


// CHECK DM VISCOSITY FACTORS

// ************* FUNCTION TO UPDATE THE VALUES OF THE PERTURBATIONS GIVEN DERIVATIVES ****************
void pertstep(struct paramsstepstruct *paramsstepnext, struct inputRKF *inputRKF){
        
    // ************* Phi, photon and baryon density perturbations ****************
    
    (*inputRKF).y.phi=(*inputRKF).pert.phi+(*paramsstepnext).deta*(*inputRKF).ksum.phi;
    (*inputRKF).y.deltarad=(*inputRKF).pert.deltarad+(*paramsstepnext).deta*(*inputRKF).ksum.deltarad;
    (*inputRKF).y.deltab=(*inputRKF).pert.deltab+(*paramsstepnext).deta*(*inputRKF).ksum.deltab;

    // ************* Neutrino moments ****************
    for ((*paramsstepnext).nmomnu=0; (*paramsstepnext).nmomnu<nummomnu-1; (*paramsstepnext).nmomnu++){
        (*inputRKF).y.Fnu[(*paramsstepnext).nmomnu]=(*inputRKF).pert.Fnu[(*paramsstepnext).nmomnu]+(*paramsstepnext).deta*(*inputRKF).ksum.Fnu[(*paramsstepnext).nmomnu];
       
              }
    // Truncation
    (*inputRKF).y.Fnu[nummomnu-1]=Fnunummom((*paramsstepnext).k, (*paramsstepnext).eta, (*inputRKF).y.Fnu[nummomnu-2], (*inputRKF).y.Fnu[nummomnu-3]);
    
    (*inputRKF).y.viscnu=1./2.*(*inputRKF).y.Fnu[2];
    
    // ************* Photon and baryon velocity perturbations and photon viscosity ****************
    if ((*paramsstepnext).recomb==0) { //Pre-recombination
        
        //Photon-baryon plasma velocity perturbation and photon viscosity
        (*inputRKF).y.vbrad=(*inputRKF).pert.vbrad+(*paramsstepnext).deta*(*inputRKF).ksum.vbrad;
        (*inputRKF).y.viscrad=viscrad(paramsstepnext, &((*inputRKF).y));

        
    } else {  //Post-recombination
         
        //Separate photon and baryon velocity peturbations
        (*inputRKF).y.vb=(*inputRKF).pert.vb+(*paramsstepnext).deta*(*inputRKF).ksum.vb;
        (*inputRKF).y.vrad=(*inputRKF).pert.vrad+(*paramsstepnext).deta*(*inputRKF).ksum.vrad;
    }
    
    // *************             ****************
    // ************* Dark Matter ****************
    // *************             ****************

    if ((*paramsstepnext).decoupl==0) {
        
            // ***************************************************************************************************************
            // *************  Pre-decoupling evolution
            // ***************************************************************************************************************

            (*inputRKF).y.deltaDM=(*inputRKF).pert.deltaDM+(*paramsstepnext).deta*(*inputRKF).ksum.deltaDM;
        
       // printf("%d",(*paramsstepnext).decoupl);
            (*inputRKF).y.vDM=(*inputRKF).pert.vDM+(*paramsstepnext).deta*(*inputRKF).ksum.vDM;
        
        if ((*paramsstepnext).a<anr) {
            (*inputRKF).y.viscDM=-square((*paramsstepnext).k)*4./15.*(*inputRKF).pert.vDM*taucomov((*paramsstepnext).a);

        } else {
            (*inputRKF).y.viscDM=-square((*paramsstepnext).k)*4./15.*(*inputRKF).pert.vDM*taucomov(anr)*((*paramsstepnext).a)/anr;
        }

                   
            // ***************************************************************************************************************
            // *************  Post-decoupling evolution.
            // ***************************************************************************************************************
        
            // *******
            // Decoupled evolution,  both with NR and R decoupling ***.
            // *******
        
            } else {
                            
            //q-loop (q=m*v, q as in Ma-Bertschinger's notation)
            for ((*paramsstepnext).n=0; (*paramsstepnext).n<npointsv; (*paramsstepnext).n++) {
                
                //Moment-loop
                for ((*paramsstepnext).nmom=0; (*paramsstepnext).nmom<nummom-1; (*paramsstepnext).nmom++) {
                    
                    (*inputRKF).y.Psi[(*paramsstepnext).nmom][(*paramsstepnext).n]=
                                    (*inputRKF).pert.Psi[(*paramsstepnext).nmom][(*paramsstepnext).n]+
                                    (*paramsstepnext).deta*(*inputRKF).ksum.Psi[(*paramsstepnext).nmom][(*paramsstepnext).n];
                }
                
               //Last moment: truncation. Last moment in my array is defined multiplied by q/epsilon.
               (*inputRKF).y.Psi[nummom-1][(*paramsstepnext).n]=Psinummomtimesqovere(paramsstepnext,&((*inputRKF).y));
                
            }

            //Vectors containing first three moments and some prefactors. Only used to perform the integrals below.
            
            for ((*paramsstepnext).n=0; (*paramsstepnext).n<npointsv; (*paramsstepnext).n++){
                
               (*inputRKF).y.PsideltaDM[(*paramsstepnext).n]=deltaDMpreffun(paramsstepnext)*(*inputRKF).y.Psi[0][(*paramsstepnext).n];
               (*inputRKF).y.PsivDM[(*paramsstepnext).n]=vDMpreffun(paramsstepnext)*(*inputRKF).y.Psi[1][(*paramsstepnext).n];
               (*inputRKF).y.Psivisc[(*paramsstepnext).n]=viscpreffun(paramsstepnext)*
                                                        (*inputRKF).y.Psi[2][(*paramsstepnext).n];
                
            }

            // Calculate density, velocity and viscosity integrals.
            (*inputRKF).y.deltaDM= deltaDMintegral2(inputRKF);
            (*inputRKF).y.vDM= vDMintegral2(inputRKF);
            (*inputRKF).y.viscDM=viscintegral2(inputRKF);
                
            }
        
    
    // ************* After all viscosities are obtained, we can calculate the psi potential ****************
    
    
    //Total viscosity times density and pressure. Used to calculate the psi potential below.
    (*inputRKF).y.rhoPvisctot=(*paramsstepnext).rhoDMa*(1+(*paramsstepnext).wDM)*(*inputRKF).y.viscDM+4./3.*(*paramsstepnext).rhonua* (*inputRKF).y.viscnu+4./3.*(*paramsstepnext).rhorada*(*inputRKF).y.viscrad;

    (*inputRKF).y.psi= psi(paramsstepnext,&((*inputRKF).y));
 
};


// ************* FUNCTION TO CALCULATE THE  PERTURBATION DERIVATIVES****************
void dstep(struct paramsstepstruct *paramsstepnext, struct pertstruct *pert, struct dstruct *de){
                            
    if ((*paramsstepnext).recomb==0) {
        
        (*de).phi=dphi(paramsstepnext, pert);
        
        //Baryons and photons
        
        (*de).deltarad=dradtight( (*de).phi, paramsstepnext, pert);

        (*de).vbrad=dvbrad(paramsstepnext, pert);
        
        (*de).deltab=dbtight( (*de).phi, paramsstepnext, pert);
        
        //Neutrinos
        
        (*de).Fnu[0]= dFnu0((*de).phi,paramsstepnext, pert);
        
        (*de).Fnu[1]= dFnu1(paramsstepnext, pert);

        
           for ((*paramsstepnext).nmomnu=2; (*paramsstepnext).nmomnu<nummomnu-1; (*paramsstepnext).nmomnu++){
               
                (*de).Fnu[(*paramsstepnext).nmomnu]=dFnul(paramsstepnext, pert);
               
                     }


    } else {
        
        (*de).phi=dphilate(paramsstepnext, pert);
        
        //Baryons and photons
        
        (*de).deltab=db((*de).phi, paramsstepnext, pert);
        
        (*de).vb=dvb( paramsstepnext, pert);  
        
        (*de).deltarad=drad((*de).phi,  paramsstepnext, pert);
        
        (*de).vrad=dvrad( paramsstepnext, pert);
        
         //Neutrinos (same as above)
               
        (*de).Fnu[0]= dFnu0((*de).phi,paramsstepnext, pert);

        (*de).Fnu[1]= dFnu1(paramsstepnext, pert);

                    
        for ((*paramsstepnext).nmomnu=2; (*paramsstepnext).nmomnu<nummomnu-1; (*paramsstepnext).nmomnu++){
           
            (*de).Fnu[(*paramsstepnext).nmomnu]=dFnul(paramsstepnext, pert);
           
        }
        
    }
    
    // *************
    // ************* Dark Matter ****************
    // *************
    
    if ((*paramsstepnext).decoupl==0) {
    
    // ************* Pre-decoupling evolution, both for relativistic and non-relativistic DM *************

    (*de).deltaDM=dDM((*de).phi, paramsstepnext, pert);
    (*de).vDM=dvDM(  paramsstepnext, pert);
                   
    //************* Post-decoupling evolution *************
     } else {
                  
         for ((*paramsstepnext).n=0; (*paramsstepnext).n<npointsv; (*paramsstepnext).n++) {
                          
             (*de).Psi[0][(*paramsstepnext).n]= dPsi0((*de).phi,paramsstepnext,pert);
             (*de).Psi[1][(*paramsstepnext).n]= dPsi1(paramsstepnext,pert);
             (*de).Psi[nummom-2][(*paramsstepnext).n]= dPsillast(paramsstepnext,pert);

             for ((*paramsstepnext).nmom=2; (*paramsstepnext).nmom<nummom-2; (*paramsstepnext).nmom++) {
                 
             (*de).Psi[(*paramsstepnext).nmom][(*paramsstepnext).n]= dPsil(paramsstepnext,pert);
                 
             }
             
         }
    
         }

};


// This function calculates the Runge-Kutta-Fehlberg derivatives given a perturbation
void dstepRKF(struct paramsstepstruct *paramsstep, struct paramsstepstruct *paramsstepnext, struct inputRKF *inputRKF){
    
    (*inputRKF).rkfcounter=0;
    
    dstep(paramsstep, &(*inputRKF).pert, &(*inputRKF).k1); // Calculates derivative
    
    updateparams(paramsstep, paramsstepnext, 1./4.); // Updates time parameter for next derivative
    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstepnext, inputRKF); // Calculates next perturbation for next derivative
    dstep(paramsstepnext, &(*inputRKF).y, &(*inputRKF).k2); // Calculates next derivative
    (*inputRKF).rkfcounter++; // rkfcounter=1

    updateparams(paramsstep, paramsstepnext, 3./8.);
    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstepnext, inputRKF);
    dstep(paramsstepnext, &(*inputRKF).y , &(*inputRKF).k3);
    (*inputRKF).rkfcounter++; //rkfcounter=2


    updateparams(paramsstep, paramsstepnext, 12./13.);
    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstepnext, inputRKF);
    dstep(paramsstepnext, &(*inputRKF).y , &(*inputRKF).k4);
    (*inputRKF).rkfcounter++; //rkfcounter=3

    
    updateparams(paramsstep, paramsstepnext, 1.);
    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstepnext, inputRKF);
    dstep(paramsstepnext,  &(*inputRKF).y , &(*inputRKF).k5);
    (*inputRKF).rkfcounter++; //rkfcounter=4

    
    updateparams(paramsstep, paramsstepnext, 1./2.);
    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstepnext, inputRKF);
    dstep(paramsstepnext,  &(*inputRKF).y ,&(*inputRKF).k6);
    (*inputRKF).rkfcounter++; //rkfcounter=5


    updateparams(paramsstep, paramsstep, 1.); // Take advantage of last RKF step to update paramsstep for next iteration
    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstep, inputRKF);
    (*inputRKF).rkfcounter++; //rkfcounter=6

    (*inputRKF).yp1=(*inputRKF).y;

    sumstruct(inputRKF, paramsstep);
    pertstep(paramsstepnext, inputRKF); //
    
    (*inputRKF).zp1=(*inputRKF).y;

     (*inputRKF).s=
    pow((*paramsstep).tol/(2*pow(
                                 pow((((*inputRKF).zp1).phi-((*inputRKF).yp1).phi),2)+
                                 pow((((*inputRKF).zp1).psi-((*inputRKF).yp1).psi),2)+
                                 pow((((*inputRKF).zp1).deltaDM-((*inputRKF).yp1).deltaDM),2)+
                                 pow((((*inputRKF).zp1).Fnu[0]-((*inputRKF).yp1).Fnu[0]),2)+
                                 pow((((*inputRKF).zp1).deltarad-((*inputRKF).yp1).deltarad),2)+
                                 pow((((*inputRKF).zp1).viscrad-((*inputRKF).yp1).viscrad),2)+
                                 pow((((*inputRKF).zp1).deltab-((*inputRKF).yp1).deltab),2)+
                                 pow((((*inputRKF).zp1).viscDM-((*inputRKF).yp1).viscDM),2)
                                 ,1./2.)
                           ),1./4.);
    
 //    (*inputRKF).s=1;
    
    (*inputRKF).pert=(*inputRKF).y;

};
 
