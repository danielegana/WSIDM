#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include "structures.h"
#include "background.h"
#include "thermodynamics.h"
#include "derivatives.h"
#include "integrals.h"
#include "RKF.h"


// ****************************************************************
// ****************************************************************
// ****************************************************************
// ****************************************************************
// ********************** BEGINNING OF CODE  **********************
// ****************************************************************
// ****************************************************************
// ****************************************************************
// ****************************************************************
// ****************************************************************


/*
 Boltzmann code for the evolution of self-interacting warm dark-matter.
 We work in conformal Newtonian gauge.
 
 Daniel Egana-Ugrinovic
 */

/*
 
 ATTENTION:
 TO CHANGE FERMIONS INTO BOSONS
 Modify mrel in this file.
 Modify xirho and xin in background.h
 Modify f0m4rel, dlogfrel and dlogfoverVrel in thermodynamics.h.
 Modify background.h line 78

 */

double Vmax;
double Varr[npointsv];
double stepfactorlate;
double Vint;

double deltaDMprefNR[npointsv];
double vDMprefNR[npointsv];
double viscprefNR[npointsv];

double deltaDMprefwdm[npointsv];
double vDMprefwdm[npointsv];
double viscprefwdm[npointsv];

double dlogfNRarr[npointsv];
double dlogfoverVNRarr[npointsv];

double dlogfrelarr[npointsv];
double dlogfoverVrelarr[npointsv];

double frelm4arr[npointsv];
double fNRm4arr[npointsv];

double anr;
double mrel;
double T0r;

double sigmaoverm;
double sigmaovermcm2g;


int main(int argc, char** argv) {
    
    // ************* DEFINE FLOATS AND STRUCTURES ****************

    double etamin,etamax,amin=1.e-10;
    int n=0,i=0, ik=0;
    int irecit=0;


    //inputRKF is the super-structure containing all perturbations and derivatives.
    struct inputRKF inputRKF;
    //paramsstep and paramsstepnext contain the background variables and other parameters. paramsstep keeps track of each time interation. paramsstepnext is used for the RKF algorithm, which within each time iteration takes fractional time iterations.
    struct paramsstepstruct paramsstep={0.};
    struct paramsstepstruct paramsstepnext={0.};

    double khinvMPc;
    
    // Photon free-streaming scale. Used to suppress photon perturbations after recombination.
    double kphotonfs=1E-38;

    //Initial and final conformal time
    etamin=pow(8*M_PI*G*(omegar0+omeganu0)*rhoc/3.,-1./2.)*amin;
    etamax=2215114565788903483279384025350285037666304.;

    //Constants for Runge-Kutta-Fehlberg algorithm
    const double RKFalg[7][6]= {{1./4., 0., 0., 0. , 0., 0.},
          {3./32., 9./32., 0., 0. , 0., 0.},
          {1932./2197., -7200./2197., 7296./2197., 0 , 0. ,0.},
          {439./216., -8., 3680./513., -845./4104. , 0., 0.},
           {-8./27., 2., -3544/2565, 1859./4104.,-11./40., 0.  },
          {25./216., 0., 1408./2565., 2197./4101., -1./5. ,  0.},
          {16./135., 0.,  6656./12825., 28561./56430., -9./50. , 2./55.}
          };


    
    // ************* NAMES OF INPUT DATA TABLES ****************
    
    char directory[100]="/home/degana/wsidm/ssixteen_clean/";
    char temp[100];
        
    char kfilename[100]="ktable.txt";
    strcpy(temp,directory);
    strcat(temp,kfilename);
    strcpy(kfilename,temp);
    

    
    char nesigmaTRECFASTfilename[200]="neSigmaTRECFAST.txt"; //The table neSigmaTRECFAST has the comoving nesigmaT. it needs to be multiplied by 1/a^2 to get the physical nesigmaTa. 
         strcpy(temp,directory);
         strcat(temp,nesigmaTRECFASTfilename);
         strcpy(nesigmaTRECFASTfilename,temp);
    
    char atableRECFASTfilename[200]="atableRECFAST.txt";
       strcpy(temp,directory);
       strcat(temp,atableRECFASTfilename);
       strcpy(atableRECFASTfilename,temp);
    
    char TMtableRECFASTfilename[200]="TMRECFAST.txt";
       strcpy(temp,directory);
       strcat(temp,TMtableRECFASTfilename);
       strcpy(TMtableRECFASTfilename,temp);
    
    
    char dlogTMtableRECFASTfilename[200]="dlogTMdlogaRECFAST.txt";
       strcpy(temp,directory);
       strcat(temp,dlogTMtableRECFASTfilename);
       strcpy(dlogTMtableRECFASTfilename,temp);

    
    char paramfile[100]="params.txt";
    strcpy(temp,directory);
    strcat(temp,paramfile);
    strcpy(paramfile,temp);
    
    
    // ************* POINTERS OF I/0 DATA TABLES ****************
    
    FILE *paramfilepointer = fopen(paramfile, "r");
    
    FILE *kdatapointer;
    kdatapointer = fopen(kfilename, "r");
    
    FILE *adatapointerRECFAST;
    adatapointerRECFAST = fopen(atableRECFASTfilename, "r");
    
    FILE *nesigmaTdatapointerRECFAST;
    nesigmaTdatapointerRECFAST = fopen(nesigmaTRECFASTfilename, "r");
    
    FILE *TMtabledatapointerRECFAST;
    TMtabledatapointerRECFAST = fopen(TMtableRECFASTfilename, "r");
    
    FILE *dlogTMtabledatapointerRECFAST;
    dlogTMtabledatapointerRECFAST = fopen(dlogTMtableRECFASTfilename, "r");
    
    // ************* COUNT LINES IN MOMENTUM TABLES ****************

      int npointsk=0;
      char ch;

      while(!feof(kdatapointer))
      {
        ch = fgetc(kdatapointer);
        if(ch == '\n')
        {
          npointsk++;
        }
      }
      rewind(kdatapointer);
 
      int skipLines;
      if (argc == 1)
         skipLines = 0;
      else
         skipLines = 2*atoi(argv[1]);
    
    // ************* READ TABLE OF MOMENTA AND PARAMETERS. ****************
    
    double k[npointsk];
    double temp2;

    for (ik = 0; ik < npointsk; ik++){
    fscanf(kdatapointer, "%lf", &k[ik] );
    }
    fclose(kdatapointer);

    for (int ip = 0; ip < skipLines; ip++){
    fscanf(paramfilepointer, "%lf", &temp2 );
    }
       
    fscanf(paramfilepointer, "%lf", &anr );
    fscanf(paramfilepointer, "%lf", &sigmaovermcm2g );
    
    sigmaoverm=sigmaovermcm2g*cm2overgtoinvGeV2;
    
    fclose(paramfilepointer);
    
    // ************* COUNT LINES IN RECFAST TABLES ****************

        int npointsrecfast=1;

        while(!feof(nesigmaTdatapointerRECFAST))
        {
          ch = fgetc(nesigmaTdatapointerRECFAST);
          if(ch == '\n')
          {
            npointsrecfast++;
          }
        }
        rewind(nesigmaTdatapointerRECFAST);
    
    
    // ************* READ RECFAST tables. ****************
    
    double *nesigmaT = malloc(sizeof(double)*npointsrecfast);
    double *arecfast = malloc(sizeof(double)*npointsrecfast);
    double *TMrecfast = malloc(sizeof(double)*npointsrecfast);
    double *dlogTMrecfast = malloc(sizeof(double)*npointsrecfast);

    
    for (irecit = 0; irecit < npointsrecfast; irecit++){
   fscanf(nesigmaTdatapointerRECFAST, "%lf",nesigmaT+irecit);
    }
    fclose(nesigmaTdatapointerRECFAST);
    
    for (irecit = 0; irecit < npointsrecfast; irecit++){
    fscanf((adatapointerRECFAST), "%lf", arecfast+irecit);
    }
    fclose((adatapointerRECFAST));
    
    for (irecit = 0; irecit < npointsrecfast; irecit++){
    fscanf((TMtabledatapointerRECFAST), "%lf", TMrecfast+irecit);
    }
    fclose((TMtabledatapointerRECFAST));
    
    for (irecit = 0; irecit < npointsrecfast; irecit++){
    fscanf((dlogTMtabledatapointerRECFAST), "%lf", dlogTMrecfast+irecit);
    }
    fclose((dlogTMtabledatapointerRECFAST));
    
    // ************* DM Mass ****************

            
    /* Fermions. Comment line below for bosons */
    
    mrel=pow(2.,7./8.)*sqrt(M_PI)*pow(rhoDM,1./4.)*pow(zeta5*xirho*gamma,3./8.)/(pow(gf,1./4.)*pow(xin,3./8.)*pow(3.,1./4.)*pow(anr,3./4.)*pow(zeta3,5./8.));  // DM mass given a relativistic Fermi distro
        
    /* Uncomment for bosons
     
     mrel=sqrt(M_PI)*pow(rhoDM,1./4.)*pow(2*zeta5*xirho*gamma/xin,3./8.)
        /(pow(gf,1./4.)*pow(xin,3./8.)*pow(anr,3./4.)*pow(zeta3,5./8.)); // DM mass given a relativistic bose distro
     
     */
    
    
    //************* NAME OF OUTPUT DATA TABLES
    
    char phifilename[100]="phidata.txt";
    
    strcpy(temp,directory);
    strcat(temp,phifilename);
    strcpy(phifilename,temp);
    
    char deltaradfilename[100]="deltaraddata.txt";
    strcpy(temp,directory);
    strcat(temp,deltaradfilename);
    strcpy(deltaradfilename,temp);
       
    char vbradfilename[100]="vbraddata.txt";
    strcpy(temp,directory);
    strcat(temp,vbradfilename);
    strcpy(vbradfilename,temp);
    
    char deltaDMfilename[100];
    strcpy(deltaDMfilename,"deltaDM anr=");
    gcvt(anr,10,temp);
    strcat(deltaDMfilename,temp);
    strcat(deltaDMfilename," sm=");
    gcvt(sigmaovermcm2g,10,temp);
    strcat(deltaDMfilename,temp);
    strcat(deltaDMfilename,".txt");
     
    strcpy(temp,directory);
    strcat(temp,deltaDMfilename);
    strcpy(deltaDMfilename,temp);

    
    char deltabfilename[100];
    strcpy(deltabfilename,"deltab anr=");
    gcvt(anr,10,temp);
    strcat(deltabfilename,temp);
    strcat(deltabfilename," sm=");
    gcvt(sigmaovermcm2g,10,temp);
    strcat(deltabfilename,temp);
    strcat(deltabfilename,".txt");
    strcpy(temp,directory);
    strcat(temp,deltabfilename);
    strcpy(deltabfilename,temp);
    
    //  FILE *phidatapointer = fopen(phifilename, "w");
    FILE *deltaDMdatapointer = fopen(deltaDMfilename, "w");
    FILE *deltabdatapointer = fopen(deltabfilename, "w");
    
    
    //Calculate velocity-dependent arrays.

    Vmax=6*anr;
    Vint=Vmax/(npointsv-1);
    
    for (i=0; i<npointsv; i++) {
        
        Varr[i]=Vint*i;
        
        deltaDMprefNR[i]=Vint*4*M_PI/(rhoDM)*square(Varr[i])*f0m4(Varr[i]);
        vDMprefNR[i]=Vint*(-4*M_PI/(rhoDM))*Varr[i]*Varr[i]*Varr[i]*f0m4(Varr[i]);
        viscprefNR[i]=Vint*(8*M_PI/(3*rhoDM))*square(Varr[i])*square(Varr[i])*f0m4(Varr[i]);
        
        dlogfNRarr[i]=dlogfNR(Varr[i]);
        dlogfoverVNRarr[i]=dlogfoverVNR(Varr[i]);
        
        deltaDMprefwdm[i]=Vint*4*M_PI*square(Varr[i])*f0m4rel(Varr[i])/mrel; // epsilon/(rho a^4) in thermodynamics.h)
        vDMprefwdm[i]=Vint*(-4*M_PI)*Varr[i]*Varr[i]*Varr[i]*f0m4rel(Varr[i]); // 1/(a^4*rho(1+wDM)*k) in thermodynamics.h
        viscprefwdm[i]=Vint*(8*M_PI/3.)*square(Varr[i])*square(Varr[i])*f0m4rel(Varr[i])*mrel; // 1/(epsilon*a^4*rho(1+wDM)) in thermodynamics.h
        
        dlogfrelarr[i]=dlogfrel(Varr[i]);
        dlogfoverVrelarr[i]=dlogfoverVrel(Varr[i]);
        
        frelm4arr[i]=Vint*f0m4rel(Varr[i]);
        fNRm4arr[i]=Vint*f0m4(Varr[i]);
        
    }
    
    // ****************************************************************
    // ****************************************************************
    // ****************************************************************
    // ******************* INITIALIZE MOMENTUM LOOP *******************
    // ****************************************************************
    // ****************************************************************
    // ****************************************************************

    
    #pragma omp parallel private(i,n, khinvMPc, paramsstep, paramsstepnext, inputRKF)
    {
    #pragma omp for
    for(ik=0;ik<npointsk;ik++){

        khinvMPc=k[ik]*GeVtoinvMPc/hcosmo;
        
        //Zero-out structures.
        paramsstep=(const struct paramsstepstruct){ 0. };
        paramsstepnext=(const struct paramsstepstruct){ 0. };
        inputRKF.pert=(const struct pertstruct){ 0. };
        inputRKF.s=1.;
        inputRKF.y=(const struct pertstruct){ 0. };
        inputRKF.yp1=(const struct pertstruct){ 0. };
        inputRKF.zp1=(const struct pertstruct){ 0. };
        inputRKF.k1=(const struct dstruct){ 0. };
        inputRKF.k2=(const struct dstruct){ 0. };
        inputRKF.k3=(const struct dstruct){ 0. };
        inputRKF.k4=(const struct dstruct){ 0. };
        inputRKF.k5=(const struct dstruct){ 0. };
        inputRKF.k6=(const struct dstruct){ 0. };
        inputRKF.ksum=(const struct dstruct){ 0. };
        inputRKF.rkfcounter=0;
                
        //Set tolerances. Tolerance for khinvMPc>10 hMPc^-1 sometimes needs to be smaller.
        if (khinvMPc<10.) {
            paramsstep.tol=1E-5;//.25E-5;
        }
        
        if (khinvMPc<=30.&&khinvMPc>=10.) {
            paramsstep.tol=1E-6;//5E-7;
                }
        
        printf("Loop %d, calculating k=%.5e, deta=%.5e, etamax=%.5e \n ", ik , khinvMPc,paramsstep.deta,etamax);
        
        //Copy Runge-Kutta-Fehlberg algorithm constants to the super-structure inputRKF.
        memcpy(inputRKF.RKFalg, RKFalg, sizeof(RKFalg));
        
        //Set initial background variables.
        paramsstep.k=k[ik];
        paramsstep.eta=etamin;
        paramsstep.deta=(etamin)/20.;
        paramsstep.a=1.e-10;
        paramsstep.Heta=calcH(paramsstep.a);
        paramsstep.rhoDMa=rhoDMafun(paramsstep.a);
        paramsstep.rhoba=rhomafun(rhob ,paramsstep.a);
        paramsstep.rhorada=rhoradafun(rhorad ,paramsstep.a);
        paramsstep.rhonua=rhoradafun(rhonu ,paramsstep.a);
        paramsstep.usDMsq=ussqfun( paramsstep.decoupl, paramsstep.a);
        paramsstep.wDM=paramsstep.usDMsq;
        paramsstep.wDMp=0.;
        paramsstep.usbradsq=usbradsqfun( paramsstep.rhoba, paramsstep.rhorada);
        paramsstep.RB=RBfun(paramsstep.rhoba, paramsstep.rhorada);
        paramsstepnext.k=k[ik];
        paramsstep.neSigmaTcomov=*nesigmaT;
        paramsstep.neSigmaTa=0;
        paramsstep.TM=*TMrecfast;
        paramsstep.dlogTM=*dlogTMrecfast;
        paramsstep.csbsq=cbsfun(paramsstep.TM,paramsstep.dlogTM);
        paramsstepnext.TM=*TMrecfast;
        paramsstepnext.dlogTM=*dlogTMrecfast;
        paramsstepnext.csbsq=cbsfun(paramsstep.TM,paramsstep.dlogTM);

        
        
        //If the decoupling happens before the particle becomes NR, we're in the usual warm dark matter regime
         
         if(3.16033*sigmaovermcm2g<anr){
             paramsstep.wdmflag=1;
             paramsstepnext.wdmflag=1;
             printf("\n Warm Dark Matter case \n");
         }
        
        printf("\n WDM flag=%d (1=WDM, 0=NRDM)\n",paramsstep.wdmflag);

        // ****************************************************************
        // ******************* INITIAL CONDITIONS FOR PERTURBATIONS *******
        // ****************************************************************

        
        inputRKF.pert.psi=-1/((4*Rnu+15.)/10.); //Normalization as in CAMB and CLASS see page 12 in https://cosmologist.info/notes/CAMB.pdf
        inputRKF.pert.phi=(1+2./5*Rnu)*inputRKF.pert.psi;
        inputRKF.pert.deltarad=-2.*inputRKF.pert.psi;
        inputRKF.pert.Fnu[0]=-2.*inputRKF.pert.psi;
        inputRKF.pert.deltaDM=-3./2.*inputRKF.pert.psi;
        inputRKF.pert.deltab=-3./2.*inputRKF.pert.psi;
        inputRKF.pert.vbrad=-1./2.*paramsstep.eta*inputRKF.pert.psi;
        inputRKF.pert.Fnu[1]=4./3.*k[ik]*1./2.*paramsstep.eta*inputRKF.pert.psi;
        inputRKF.pert.Fnu[2]=2./15.*square(k[ik]*paramsstep.eta)*inputRKF.pert.psi;
        inputRKF.pert.vDM=-1./2.*paramsstep.eta*inputRKF.pert.psi;
        inputRKF.pert.viscDM=0.;
        inputRKF.pert.viscnu=1./15.*square(k[ik]*paramsstep.eta)*inputRKF.pert.psi;
        inputRKF.pert.viscrad=0.;
        inputRKF.pert.rhoPvisctot=4./3.*(paramsstep).rhonua*inputRKF.pert.viscnu;
        
        // ****************************************************************
        // ****************************************************************
        // ****************************************************************
        // ******************* INITIALIZE TIME LOOP ***********************
        // ****************************************************************
        // ****************************************************************
        // ****************************************************************
        
        
        for(i=0;paramsstep.eta<=etamax;i++){
            
            // Find closest element in recombination array
          
            for (irecit=paramsstep.irec; irecit<npointsrecfast; irecit++) {
                if((*(arecfast+paramsstep.irec))<paramsstep.a){
                        paramsstep.irec++;
                        
                        paramsstep.neSigmaTcomov=(*(nesigmaT+paramsstep.irec));
                        paramsstep.TM=(*(TMrecfast+paramsstep.irec));
                        paramsstep.dlogTM=(*(dlogTMrecfast+paramsstep.irec));
                        paramsstep.csbsq=cbsfun(paramsstep.TM,paramsstep.dlogTM);
                        
                        paramsstepnext.TM=(*(TMrecfast+paramsstep.irec));
                        paramsstepnext.dlogTM=(*(dlogTMrecfast+paramsstep.irec));
                        paramsstepnext.csbsq=cbsfun(paramsstep.TM,paramsstep.dlogTM);

                    }
                
                if((*(arecfast+paramsstep.irec))>paramsstep.a){
                    break;
                    }
            }
            
            // Check DM decoupling conditions. paramsstep.decoupl=0 if DM is decoupled, >0 otherwise.
            if (oneovertaucomov(paramsstep.a) < paramsstep.Heta) {
                (paramsstep.decoupl)++;
                (paramsstepnext.decoupl)++;
            }
            
   /*         double rhotest;  TEST if rho integral coincides with rhoDM
            rhotest=rhointegral(frelm4arr,&paramsstep);
            
          printf("\n a value= %.3e, rhotest=%.3e,rho=%.3e", paramsstep.a, rhotest,paramsstep.rhoDMa);*/

            
            if (paramsstepnext.decoupl==1){

                printf("\n DM decouples at a=%.5e, oneovertau=%.5e, Hcomov=%.5e, sigmaoverm=%.5e, i=%d \n", paramsstep.a,oneovertaucomov(paramsstep.a),paramsstep.Heta,sigmaoverm/cm2overgtoinvGeV2,i);
                

                if(paramsstep.wdmflag==0){
                       
                       for (n=0; n<npointsv; n++) {
                        inputRKF.pert.Psi[0][n]=-1./3.*inputRKF.pert.deltaDM*dlogfNRarr[n];
                        inputRKF.pert.Psi[1][n]=k[ik]*paramsstep.a/3.*inputRKF.pert.vDM*dlogfoverVNRarr[n];
                        inputRKF.pert.Psi[2][n]=square(k[ik])*2./(15*paramsstep.Heta)*inputRKF.pert.vDM*dlogfNRarr[n];

                                          }
                   } else if(paramsstep.wdmflag==1){
                       
                       for (n=0; n<npointsv; n++) {
                           
                           if (paramsstep.a>anr) {
                           inputRKF.pert.Psi[0][n]=-1./3.*inputRKF.pert.deltaDM*dlogfrelarr[n]; // This is here because the model may be WDM, but if by the time the simulation starts we're already in the NR regime, then epsilon=a*m, so the matching condition gets the factor of 1/3 in front, not 1/4 like when epsilon=q.
                           }
                           else{
                               inputRKF.pert.Psi[0][n]=-1./3.*inputRKF.pert.deltaDM*dlogfrelarr[n];
                           }
                         inputRKF.pert.Psi[1][n]=k[ik]*paramsstep.a/3.*inputRKF.pert.vDM*dlogfoverVrelarr[n];
        
                       //    inputRKF.pert.Psi[1][n]=k[ik]/3.*inputRKF.pert.vDM*dlogfrelarr[n];
                           inputRKF.pert.Psi[2][n]=square(k[ik])*2./(15*paramsstep.Heta)*inputRKF.pert.vDM*dlogfrelarr[n];
                                                    
                       }
                   }
                 
            }
                        
            // ****************************************************************
            // ******************* PRE-RECOMBINATION EVOLUTION ****************
            // ****************************************************************
            
            if(paramsstep.eta<etarec+paramsstep.deta){
                
                //Call RKF algorithm to perform one time-step, and adapt step-size.
                dstepRKF(&paramsstep, &paramsstepnext, &inputRKF);
                paramsstep.deta=inputRKF.s*paramsstep.deta;
                
             // ******************* TRANSITION TO POST-RECOMBINATION LOOP *******************
                
                if(etarec-paramsstep.deta<=paramsstep.eta && paramsstep.eta<=etarec && paramsstep.recomb==0){
                                            
                    printf("Viscosity at recombination: nesigmaTcomov=%.5e\n",paramsstep.neSigmaTcomov);
                    paramsstep.recomb=1;
                    paramsstepnext.recomb=1;
                    inputRKF.pert.vb=inputRKF.pert.vbrad;
                    inputRKF.pert.deltarad=inputRKF.pert.deltarad;
                    inputRKF.pert.vrad=inputRKF.pert.vbrad;
                }
                
            } // End of pre-recombination if
                           
            // ****************************************************************
            // ******************* POST-RECOMBINATION EVOLUTION ***************
            // ****************************************************************

            if(paramsstep.recomb>0){
                                
                //Call RKF algorithm to perform one time-step,  and adapt step-size.
                dstepRKF( &paramsstep, &paramsstepnext, &inputRKF);
                paramsstep.deta=inputRKF.s*paramsstep.deta;
                
                //Suppress photon perturbations due to free-streaming.
                inputRKF.pert.deltarad=inputRKF.pert.deltarad*exp(-square(k[ik]/kphotonfs));
                inputRKF.pert.vrad=inputRKF.pert.vrad*exp(-square(k[ik]/kphotonfs));

                     
            } // End of if post-recombination if
                    
            
             // ************* WRITE DM AND BARYON PERTURBATIONS TO FILES ****************
            if(paramsstep.eta>etamax){
               printf("\n deltaDM=%.5e \n",inputRKF.pert.deltaDM);
               printf("\n deltab=%.5e \n",inputRKF.pert.deltab);
               fprintf(deltaDMdatapointer, "%.5e   %.5e \n", khinvMPc,inputRKF.pert.deltaDM);
               fprintf(deltabdatapointer, "%.5e   %.5e \n", khinvMPc,inputRKF.pert.deltab);

            }
            
        }
        
        // End of time loop
        
    } // End of momentum loop
    } // End of parallelization routine

    free(nesigmaT);
    free(arecfast);
    free(TMrecfast);
    free(dlogTMrecfast);



    // ************* CLOSE FILE POINTERS ****************
    
//    fclose(phidatapointer);
    fclose(deltaDMdatapointer);
    fclose(deltabdatapointer);

    return 0;
    
    
}




    




