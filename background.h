/*    *************** DEFINITION OF BACKGROUND CONSTANTS AND FUNCTIONS ************ */

extern double anr;

#define G 1./148406323681458953382492023727134198357.
#define wrad 1./3.
#define usradsq 1./3.
#define hcosmo 0.6732117
#define H0GeV 1.436041e-42 //Number derived from hcosmo with unit conversions
#define rhoc 3.*H0GeV*H0GeV/(8*M_PI*G)
#define omegar0 0.0000545651 // Here the omegas are capital omegas
#define omeganu0 3.046*pow(4./11.,4./3.)*7/8*omegar0
#define Rnu (omeganu0/(omeganu0+omegar0))
#define omegaDM0 0.2650128
#define omegab0 0.0493868
#define omegam0 (omegab0+omegaDM0)
#define omegaLambda0 (1.-omegam0-omegar0-omeganu0)
#define rhoDM omegaDM0*rhoc
#define rhoM omegam0*rhoc
#define rhob omegab0*rhoc
#define rhorad omegar0*rhoc
#define rhonu omeganu0*rhoc
#define etarec 4.38005e40 //***
#define arec 1/(1089.914+1)
#define muGeV 1.22


#define gamma 5./3.
#define gf 2.

#define zeta3 1.20206
#define zeta5 1.03693

// For fermions (comment for bosons)

#define xirho 45./8.
#define xin 3./4.

/* For bosons (comment for fermions)
 
#define xirho 6.
#define xin 1.
 
*/
 
#define GeVtoinvMPc 1.566497462E+38
#define cm2overgtoinvGeV2 4586.56

#define square(x) (x*x)
#define cube(x) (x*x*x)
#define etastep( detavar,  etavar) (etavar+detavar)
#define astep(detavar, avar, Hetavar) (avar*(1+detavar*Hetavar))

// Comoving Hubble scale (H(a)*a=da/deta/a=Heta)
#define calcH(avar) (H0GeV*sqrt(square(avar)*omegaLambda0+omegam0/avar+(omegar0+omeganu0)/square(avar)))

#define rhomafun(rhom, a) (rhom/(a*a*a))

#define rhoradafun(rhor, a) (rhor/(a*a*a*a))

#define RBfun(rhoba, rhorada) (3.*rhoba/(4.*rhorada))

#define usbradsqfun(rhoba, rhorada) (1./(3.*(1+(3.*rhoba/(4.*rhorada)))))

#define knufs0 sqrt(4*M_PI*G*rhonu)

#define knufsfun(a) knufs0/a
    
double astepfun(double detavar,double etavar,double avar, double Hetavar){ //This function is not being used.
    if(avar<1e-8){
        return etavar*pow(8*M_PI*G*rhoc*(omeganu0+omegar0)/3,1./2.);
    }
    else{
        return avar*(1+detavar*Hetavar);
    }
}

double rhoDMafun(double a){
    
    double rhoDMa;
    
    if (a>0.876078*anr/(pow(gamma,1./2.))) {  // FERMIONS
        //       if (a>1.11941*anr/(pow(gamma,1./2.))) {  // BOSONS
        rhoDMa=rhoDM/(a*a*a);
    } else {
        rhoDMa=rhoDM*0.876078*anr/(pow(gamma,1./2.))/(square(a)*square(a)); // Assumes Fermi Distro zero chem pot
 //      rhoDMa=rhoDM*1.11941*anr/(pow(gamma,1./2.))/(square(a)*square(a)); // BOSONS. Assumes Bose distro zero chem pot

    }
    
    return rhoDMa;
}






    




