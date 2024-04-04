/*    *************** DEFINITION OF PERTURBATION DERIVATIVES ************
 
Define all derivatives as macros (helps with code speed)
 
 */

#define psi(paramsptr, pertptr) (-1/square((*paramsptr).k)*(12*M_PI*G*square((*paramsptr).a)*(*pertptr).rhoPvisctot)+(*pertptr).phi)

#define dphi(paramsptr, pertptr) (1/(3.*(*paramsptr).Heta)*((-4.*M_PI*G*square((*paramsptr).a))*((*paramsptr).rhorada*(*pertptr).deltarad+(*paramsptr).rhonua*(*pertptr).Fnu[0]+(*paramsptr).rhoDMa*(*pertptr).deltaDM+(*paramsptr).rhoba*(*pertptr).deltab)-3.*square((*paramsptr).Heta)*(*pertptr).psi-square((*paramsptr).k)*(*pertptr).phi))

#define dphilate(paramsptr, pertptr) (-4.*M_PI*G*square((*paramsptr).a)*((*paramsptr).rhoDMa*(*pertptr).vDM + (*paramsptr).rhoba*(*pertptr).vb + 4./3.*(*paramsptr).rhorada*(*pertptr).vrad-1/(*paramsptr).k*(*paramsptr).rhonua*(*pertptr).Fnu[1])-(*paramsptr).Heta*(*pertptr).psi)

#define dradtight(dphiRHS, paramsptr, pertptr) (3.*(1+wrad)*dphiRHS-(3*(*paramsptr).Heta*(usradsq-wrad)*(*pertptr).deltarad-(1.+wrad)*square((*paramsptr).k)*(*pertptr).vbrad))

#define drad(dphiRHS, paramsptr, pertptr) (3.*(1+wrad)*dphiRHS-(3*(*paramsptr).Heta*(usradsq-wrad)*(*pertptr).deltarad-(1.+wrad)*square((*paramsptr).k)*(*pertptr).vrad))

#define viscrad(paramsptr, pertptr) (-8./27.*(*pertptr).vbrad*square((*paramsptr).k)/(*paramsptr).neSigmaTa)
//viscrad is sigma in Ma's notation. this is only the photon viscosity

#define dvbrad(paramsptr, pertptr) ( -(*pertptr).psi - ((*paramsptr).Heta*(*paramsptr).RB/(1.+(*paramsptr).RB)*(*pertptr).vbrad +3./4.*(*paramsptr).usbradsq*(*pertptr).deltarad)+(*pertptr).viscrad/(1+(*paramsptr).RB))

#define dvrad(paramsptr, pertptr) (-(*pertptr).psi - (1./4.*(*pertptr).deltarad)+(*pertptr).viscrad+(*paramsptr).neSigmaTa*((*pertptr).vb-(*pertptr).vrad))

#define dDM( dphiRHS, paramsptr, pertptr) ( 3*(1 + (*paramsptr).wDM)*dphiRHS -3*(*paramsptr).Heta*((*paramsptr).usDMsq - (*paramsptr).wDM)*(*pertptr).deltaDM + (1 + (*paramsptr).wDM)*square((*paramsptr).k)*(*pertptr).vDM)
            
#define dvDM( paramsptr, pertptr) (-(*paramsptr).Heta*(1-3*(*paramsptr).wDM)*(*pertptr).vDM-(*pertptr).psi+(*pertptr).viscDM+(-(*paramsptr).wDMp*(*pertptr).vDM -(*paramsptr).usDMsq*(*pertptr).deltaDM)/(1 + (*paramsptr).wDM))

#define dbtight( dphiRHS, paramsptr, pertptr) (3*dphiRHS + square((*paramsptr).k)*(*pertptr).vbrad )

#define db( dphiRHS,  paramsptr, pertptr) (3*dphiRHS + square((*paramsptr).k)*(*pertptr).vb )
    
#define dvb(paramsptr, pertptr) (-(*paramsptr).Heta*(*pertptr).vb -(*pertptr).psi+(*paramsptr).neSigmaTa/(*paramsptr).RB*((*pertptr).vrad-(*pertptr).vb)-(*paramsptr).csbsq*(*pertptr).deltab)
    
#define dPsi0(dphiRHS, paramsptr, pertptr) ( -qoverepsilon(paramsptr)*(*paramsptr).k*(*pertptr).Psi[1][(*paramsstepnext).n] - dphiRHS*dlogf(paramsptr) )
    
#define dPsi1(paramsptr, pertptr) ( qoverepsilon(paramsptr)*(*paramsptr).k/3*((*pertptr).Psi[0][(*paramsstepnext).n]-2*(*pertptr).Psi[2][(*paramsstepnext).n]) - epsilonoverqtimesV(paramsptr)*(*paramsptr).k/3.*(*pertptr).psi*dlogfoverV(paramsptr))
        
#define dPsil(paramsptr, pertptr) ( (*paramsptr).k/((2*(*paramsptr).nmom+1))*qoverepsilon(paramsptr)*((*paramsptr).nmom*(*pertptr).Psi[(*paramsptr).nmom-1][(*paramsptr).n]-((*paramsptr).nmom+1)*(*pertptr).Psi[(*paramsptr).nmom+1][(*paramsptr).n]))

#define dPsillast(paramsptr, pertptr) ( (*paramsptr).k/(2*(nummom-2)+1)*(qoverepsilon(paramsptr)*(nummom-2)*(*pertptr).Psi[(nummom-2)-1][(*paramsptr).n]-((nummom-2)+1)*(*pertptr).Psi[(nummom-2)+1][(*paramsptr).n]))

// The last moment, Psi[nummom-1] is defined as Psi*qoverepsilon

#define Psinummomtimesqovere(paramsptr, pertptr) ((2*(nummom-2)+1)/((*paramsptr).k*(*paramsptr).eta)*(*pertptr).Psi[nummom-2][(*paramsptr).n]-qoverepsilon(paramsptr)*(*pertptr).Psi[nummom-3][(*paramsptr).n])

#define dFnu0(dphi,paramsptr, pertptr) ( -(*paramsptr).k* (*pertptr).Fnu[1]+4.*dphi)

#define dFnu1(paramsptr, pertptr) ( 4./3.*(*paramsptr).k*(1./4.*(*pertptr).Fnu[0]-1./2.*(*pertptr).Fnu[2]+(*pertptr).psi))

#define dFnul(paramsptr, pertptr) ((*paramsstepnext).k/(2.*(*paramsstepnext).nmomnu+1.)*((*paramsstepnext).nmomnu*(*pertptr).Fnu[(*paramsstepnext).nmomnu-1]-((*paramsstepnext).nmomnu+1)*(*pertptr).Fnu[(*paramsstepnext).nmomnu+1]))

#define Fnunummom(k, eta,Fm1,Fm2) ((2*nummomnu-3)/(k*eta)*Fm1-Fm2)


//Function to update background values.

void updateparams(struct paramsstepstruct *paramsstep, struct paramsstepstruct *paramsstepnext, double fraction){
    
     (*paramsstepnext).deta=(*paramsstep).deta*fraction;
     (*paramsstepnext).eta=etastep((*paramsstepnext).deta,(*paramsstep).eta);
     (*paramsstepnext).a=astep((*paramsstepnext).deta,(*paramsstep).a,(*paramsstep).Heta);
     (*paramsstepnext).Heta=calcH((*paramsstepnext).a);
     
     (*paramsstepnext).rhoDMa=rhoDMafun((*paramsstepnext).a);
     (*paramsstepnext).rhoba=rhomafun(rhob ,(*paramsstepnext).a);
     (*paramsstepnext).rhorada=rhoradafun(rhorad ,(*paramsstepnext).a);
     (*paramsstepnext).rhonua=rhoradafun(rhonu ,(*paramsstepnext).a);

    
     (*paramsstepnext).usDMsq=ussqfun((*paramsstepnext).decoupl, (*paramsstepnext).a);
     (*paramsstepnext).wDM=(*paramsstepnext).usDMsq; // Add wDMp if needed
     (*paramsstepnext).usbradsq=usbradsqfun((*paramsstepnext).rhoba, (*paramsstepnext).rhorada);
     (*paramsstepnext).RB=RBfun((*paramsstepnext).rhoba, (*paramsstepnext).rhorada);
    
     (*paramsstepnext).deta=(*paramsstep).deta; // This is needed, because once the derivatives in the RKF algorithm are computed with fractional time steps, the variation of the function deltaf is computed with the full time step, h*k1, h*k2, etc
        
     (*paramsstepnext).neSigmaTa=((*paramsstep).neSigmaTcomov)/square((*paramsstepnext).a);

    
}


