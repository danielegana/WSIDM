// ############ THERMODYNAMICS ##############

extern double anr;
extern double sigmaoverm;
extern double Varr[npointsv];

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



#define csovar 1/sqrt(3.)*anr


#define tau(a) (1/(sigmaoverm*rhoDMafun(a)*velfun(a)))

#define taucomov(a) (1/(a*sigmaoverm*rhoDMafun(a)*velfun(a)))

#define oneovertaucomov(a) (velfun(a)*a*sigmaoverm*rhoDMafun(a))

#define f0m4(V) (rhoDM*pow(gamma,3./2.)/pow(2.*M_PI*square(anr)/3.,3./2.)*exp(-3*gamma*square(V)/(2*square(anr)))) //Note that it is mass independent

#define dlogfNR(V) (-3*gamma*square(V)/(square(anr)))   // This is dlogf/dlogV=dlogf/dV*V

#define dlogfoverVNR(V) (-3*gamma*V/(square(anr)))

// Old #define f0m4rel(V) ((1/cube(2*M_PI))*square(mrel)*square(mrel)*gf/(exp(3.*gamma*V/anr)+1.))


/* Relativistic fermi distro */

#define frelexppref sqrt(2.*gamma*xirho*zeta5/(xin*zeta3))

#define f0m4rel(V) ((1/cube(2*M_PI))*square(mrel)*square(mrel)*gf/(exp(frelexppref*V/anr)+1.))

#define dlogfrel(V) (-frelexppref*V/anr/(1.+exp(-frelexppref*V/anr)))

#define dlogfoverVrel(V) (-frelexppref/anr/(1.+exp(-frelexppref*V/anr)))

/* Relativistic bose distro. Uncomment for bosons
 
#define f0m4rel(V) ((1/cube(2*M_PI))*square(mrel)*square(mrel)*gf/(exp(frelexppref*V/anr)-1.))
 
#define dlogfrel(V) (-frelexppref*V/anr/(1.-exp(-frelexppref*V/anr)))
 
#define dlogfoverVrel(V) (-frelexppref/anr/(1.-exp(-frelexppref*V/anr)))
 

*/


double velfun(double a){
    
    double velvar;
    
    if (a>anr) {
        velvar=1./a*anr;
    } else {
        velvar=1.;
    }
    
    return velvar;
}

double cbsfun(double TM, double dlogTM){
    
    return TM/muGeV*(1.-1./3.*dlogTM);
    
}

double ussqfun(double decoupl, double a){
    
    // decoupl=0 when fluid is still coupled, taucomov(sigmaoverm, a , anr ) < 1./(Heta)
    
    double usvar;
    
    if (decoupl==0) {
        if (a>anr) {
            usvar=csovar/a;
        } else {
            usvar=csovar/anr;
        }
    } else {
        usvar=0.;
    }
    
    return square(usvar);
    
}

double epsilon(struct paramsstepstruct *paramsstep){
    return sqrt(square((*paramsstep).a*mrel)+square(Varr[(*paramsstep).n]*mrel));
}

double qoverepsilon(struct paramsstepstruct *paramsstep){
/*    if ((*paramsstep).a<anr){
        return 1.;
    } else {
        return Varr[(*paramsstep).n]/(*paramsstep).a;
    } */
    return Varr[(*paramsstep).n]*mrel/sqrt(square((*paramsstep).a*mrel)+square(Varr[(*paramsstep).n]*mrel));
}

double epsilonoverqtimesV(struct paramsstepstruct *paramsstep){
   /* if ((*paramsstep).a<anr){
        return Varr[(*paramsstep).n];
        
    } else {
        return (*paramsstep).a;
    }*/
    return sqrt(square((*paramsstep).a*mrel)+square(Varr[(*paramsstep).n]*mrel))/mrel;
}

double dlogf(struct paramsstepstruct *paramsstep){
    if ((*paramsstep).wdmflag==0){
        return dlogfNRarr[(*paramsstep).n];
    }
    else{
        return dlogfrelarr[(*paramsstep).n];
    }
 
}

double dlogfoverV(struct paramsstepstruct *paramsstep){
    if ((*paramsstep).wdmflag==0){
        return dlogfoverVNRarr[(*paramsstep).n];
    } else {
        return dlogfoverVrelarr[(*paramsstep).n];
    }
}

double deltaDMpreffun(struct paramsstepstruct *paramsstep){
    
    double ddmpref;
    
    if ((*paramsstep).wdmflag==0){
        ddmpref=deltaDMprefNR[(*paramsstep).n];
    } else {
        ddmpref= deltaDMprefwdm[(*paramsstep).n]*epsilon(paramsstep)/((*paramsstep).rhoDMa*square((*paramsstep).a)*square((*paramsstep).a));
        
    }
    return ddmpref;
}

double vDMpreffun(struct paramsstepstruct *paramsstep){
    if ((*paramsstep).wdmflag==0){
        return vDMprefNR[(*paramsstep).n]/
        ((*paramsstep).k*(*paramsstep).a);
    } else {
        return  vDMprefwdm[(*paramsstep).n]/
        ((*paramsstep).k*(*paramsstep).rhoDMa*(1+(*paramsstep).wDM)*square((*paramsstep).a)*square((*paramsstep).a));
    }
    
}

double viscpreffun
(struct paramsstepstruct *paramsstep){
    if ((*paramsstep).wdmflag==0){
        return 1/square((*paramsstep).a)*viscprefNR[(*paramsstep).n];
    } else {
        return viscprefwdm[(*paramsstep).n]/ (epsilon(paramsstep)*(*paramsstep).rhoDMa*(1+(*paramsstep).wDM)*square((*paramsstep).a)*square((*paramsstep).a));
    }
}




