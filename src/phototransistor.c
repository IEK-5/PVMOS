/*
    PhotoVoltaic Module Simulator (PVMOS), a finite difference ODE solver 
    for solar modules. 
    Copyright (C) 2014  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*****************************************************************              
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
 +  IEK-5 PHOTOVOLTAIK                                           *              
 *                                                               *              
 *        ########                _   _                          *              
 *     ##########                |_| |_|                         *              
 *    ##########     ##         _ _   _ _     ___ ____ _   _     *              
 *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *              
 *   #########     #####    _  | | | | | |    | | |   | |_| |    *              
 *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *              
 *    #          ######     \___/ \___/|_____|___\____|_| |_|    *              
 *     ##      #######      F o r s c h u n g s z e n t r u m    *              
 *       ##########                                              *              
 *                                                               *              
 *   http://www.fz-juelich.de/iek/iek-5/DE/Home/home_node.html   *              
 *****************************************************************
 *                                                               *
 *    Dr. Bart E. Pieters 2015                                   *
 *                                                               *             
 *****************************************************************/                                                                             

/*****************************************************************           
 * FUNCTION:                                                     *                
 * compute equivalent circuit for a solar cell                   *
 *                                                               *            
 *****************************************************************/     
/*****************************************************************
 * The following script implements the following circuit:	 *
 * Phototransistor Model with shunt and series resistance:	 *
 *			 					 *
 *                          2      Rs      Va                    *                     
 *                                ___                            *                  
 *                       +----+--|___|-----o                     *                  
 *                       |    |                                  *                  
 *                       |    |                                  *                  
 *              1        /C   |                                  *                  
 *                   B |/    .-.                                 *                  
 *               +-----|     | |  Rsh                            *                      
 *               |     |\E   | |                                 *                  
 *              .-.      \   '-'                                 *                  
 *         Iph (---)     |    |                                  *                  
 *              '-'      |    |                                  *                  
 *               |       |    |             0                    *                   
 *               +-------+----+------------o      		 *
 *								 *
 *****************************************************************/                                        
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phototransistor.h"


#define k_b 8.6173323478e-5

/* hard coded breakoff criterion for the kicrhoff current law error */
#define EPS 1e-6

/* Mostly simple and fast newton-raphson will suffice */
#define MAXITER 25

/* If it does not we go to slower fallback routines, i.e. make the number above to small and you end up in slow routines too often, make it too high and you iterate a lot before going into a robust routine */	
#define MAXITER_2 100

#define T0 300


#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))

void * InitPhotoTransistorStruct(int *size)
{
	PhotoTransistor *res;
	*size=sizeof(PhotoTransistor);
	res=malloc(sizeof(PhotoTransistor));
	res->Jsbc=5e-4;
	res->Jsbe=5e-15;
	res->Jse=0.0;
	res->Jsc=0.0;
	res->Nf=1.0;
	res->Nr=1.0;
	res->Ne=1.5;
	res->Nc=2.0;
	res->Vaf=INFINITY;
	res->Var=INFINITY;
	res->Bf=3.0;
	res->Jph=0.014;
	res->Rs=1e-7;
	res->Rsh=1e7;
	res->EgBE=1.2;
	res->PhiBC=0.3;
	res->XTIBE=3.0;
	res->XTIBC=3.0;
	res->XTB=0.0;
	return (void *) res;
}

static inline void DiodeEq(double V, double Is, double N, double T, double Eg, double XTI, double gmin, double *I, double *dI)
{
	/* 
	Diode current and derivative thereof.
	Input:
	V	- Applied voltage
	Is	- Saturation current at T0
	N	- ideality factor
	Eg	- Bandgap/barrier height
	XTI	- Temperature dependency if saturation current ((T/T0)^XTI)
	gmin	- in case of gmin stepping we can add a conductance in parallel
	*/
	 
	double nVt, fTs;
	
	nVt=N*k_b*T;
	fTs=exp((T/T0-1)*Eg/nVt)*pow((T/T0), XTI);
	
	if (I)
		(*I)=fTs*Is*(exp(V/nVt)-1)+V*gmin;
	if (dI)
		(*dI)=fTs*Is*exp(V/nVt)/nVt+gmin;
}

static inline void System(double Va, PhotoTransistor M, double V1, double V2, double gmin, double T, double *dV1, double *dV2, double *E, double *dJdVa, double *Gmax)
{
	/*
	 Create linearized system matrix and solve one NR step
	 Input:
	 	Va		- Applied voltage
		M		- Model parameters
		V1,V2		- node voltages
		gmin		- Gmin value
	 
	 Output:
		dV1, dV2	- Change in node voltages (NR step)
		E		- Error in Kirchhof law for current node voltages
		dJdVa		- Derivative of the current
		Gmax		- Maximum diagonal value in the system (used to initiate the gmin value when going to Gmin stepping to solve convergence issues)
	*/
	/* init Nodal Analysis Matrix */
	double S11=0, S12=0, S21=0, S22=0;
	/* RHS */
	double I1=0, I2=0;
	double Vbe, Vbc;
	double IF, gBEI, IBEN, gBEN, Bf, IBE, g_pi;
	double IBC, gBCI, IBCN, gBCN, g_mu;
	double gIF, Qb, IT, dQdVbe, dQdVbc, g_mf, g_mr;
	double IBEeq, IBCeq, ICEeq;
	double b1, b2, Dx, Dy, D, Gce;
	
	*E=1.0;
	
	I1=M.Jph;
	/* add NA-stamps for Rs and Rp */
	S22+=1/M.Rs;
	S22+=1/M.Rsh;
	I2+=Va/M.Rs;
	
	/* compute Vbe and Vbc */
	Vbe=V1;
	Vbc=V1-V2;
	
	/************ Gummel-Poon
	 * symbol names should correspond to common Gummel-Poon documentation */
	
	/* forward and reverse currents & conductivities, note we do not use gmin here, now I connect all nodes to ground with gmin */
	DiodeEq(Vbe,M.Jsbe,M.Nf, T, M.EgBE, M.XTIBE, 0, &IF, &gBEI);
	DiodeEq(Vbe,M.Jse*pow((T/T0),M.XTB),M.Ne, T, M.EgBE, M.XTIBE, 0, &IBEN, &gBEN);
	Bf=M.Bf*pow((T/T0), M.XTB);
	gBEI/=Bf;
	IBE=IF/Bf+IBEN;
	g_pi=gBEI+gBEN;
	
	
	/* Turns out that a Br unequal to 0 makes the transistor active, i.e. generate power
	 * Thus I want Br to be strictly 0, this model implementation makes Br 0
	 */
	 	
	DiodeEq(Vbc,M.Jsbc,M.Nr, T, M.PhiBC, M.XTIBC, 0, &IBC, &gBCI);
	DiodeEq(Vbc,M.Jsc*pow((T/T0),M.XTB),M.Nc, T, M.PhiBC, M.XTIBC, 0, &IBCN, &gBCN);
	g_mu=gBCI+gBCN;
		
	/* transfer current and conductivities */
	gIF=gBEI*M.Bf;
	
	Qb=1/(1-Vbc/M.Vaf-Vbe/M.Var);
	
	IT=IF/Qb;
	dQdVbe=Qb*Qb/M.Var;
	dQdVbc=Qb*Qb/M.Vaf;
	
	g_mf=(gIF-IT*dQdVbe)/Qb;
	g_mr=(-IT*dQdVbc)/Qb;
	
	/* Right Hand Side Terms */	
	IBEeq=IBE-g_pi*Vbe;
	IBCeq=IBC-g_mu*Vbc;
	ICEeq=IT-g_mf*Vbe+g_mr*Vbc;
	/* NA entries */
	S11+=(g_mu+g_pi);
	S12-=g_mu;
	I1-=(IBEeq+IBCeq);
	
	S21-=(g_mu-g_mf+g_mr);
	S22+=(g_mu+g_mr);
	I2+=(IBCeq-ICEeq);
	
	/* add gmin to diagonal */
	S11+=gmin;
	S22+=gmin;
	
	/* Comute KCL error and Newton step */
	/* KCL errors */
	b1=I1-S11*V1-S12*V2;
	b2=I2-S21*V1-S22*V2;
	/* Current MS KCL error (Root mean square minus the root) */
	*E=b1*b1+b2*b2;
	/* solve system */
	D=S11*S22-S12*S21;
	Dx=b1*S22-S12*b2;
	Dy=S11*b2-b1*S21;
	*dV1=Dx/D;
	*dV2=Dy/D;
	
	/* compute dJ/dVa */
	Gce=(g_mf*g_mu-g_mr*g_pi)/(g_pi+g_mu)+1/(1/g_pi+1/g_mu)+1/M.Rsh;
	*dJdVa=1/(M.Rs+1/Gce);
	
	/* maximum Conductivity */
	*Gmax=MAX(S11,S22)-gmin;
}



static inline int Solve_gmin(double Va, double *V1, double *V2, PhotoTransistor M, double gmin, double T, double *dJdVa, double *Gmax, double *err, double eps, int maxiter)
{
	/* Solve the circuit at applied bias of Va
	   Input:
	  	Va	- Applied bias voltage
	  	V1,V2	- initial node voltages
	  	M	- Model parameters
		gmin	- minimum conductance
		eps	- accuracy criterion (KCL error)
		maxiter - muximum number or iterations
	   Output:
	  	V1,V2	- updated node voltages
		dJdVa	- derivative current to applied voltage
		Gmax	- maximum diagonal value in the system
		err	- KCL error of the solution
	   return value:
	   	0 	- solution is not converges
		1	- solutiuon did converge
	
	*/
	int conv=1, iter=0;
	double dVmax, dV1, dV2, dv1, dv2, En, E;
	double dvm, a, v1, v2;
	
	/* maximum voltage step in primitive damping scheme, depends linearly on T*/
	dVmax=5*(T/300);
	
	/* Do newton step and determine KCL Error*/
	System(Va, M, *V1, *V2, gmin, T, &dv1, &dv2, &En, dJdVa, Gmax);
	eps*=eps;
	
	while ((En>eps) && (iter<maxiter))
	{
		/* Current KCL error*/
		E=En;
		
		/* Voltagge step*/
		dV1=dv1;
		dV2=dv2;
		
		/* simple damping */
		dvm=MAX(fabs(dV1), fabs(dV2));
		if (dvm>dVmax)
		{
			dV1*=dVmax/dvm;
			dV2*=dVmax/dvm;
		}
		
		/* do a linear search in the NR direction 
		  such that the new KCL error is smaller than the current Error
		  we start with the full NR step (a=1) */
		a=1.0;
		v1=(*V1)+dV1/a;
		v2=(*V2)+dV2/a;
		System(Va, M, v1, v2, gmin, T, &dv1, &dv2, &En, dJdVa, Gmax);
		while ((En>E)&&(a<9))
		{
			/* if the new error is larger than the current, reduce the step by a factor two
			  we stop at 4 iterations (a 0V step does not get us closer to a solution either)*/
			a*=2;
			v1=(*V1)+dV1/a;
			v2=(*V2)+dV2/a;
			System(Va, M, v1, v2, gmin, T, &dv1, &dv2, &En, dJdVa, Gmax);
		}
		*V1=v1;
		*V2=v2;
		iter++;
	}
	if (En>eps)
		conv=0;  /* not converged*/
	*err=sqrt(En);
	return conv;
}

#define NGMIN 10
#define NOGMIN 7
#define EPSGMIN 1e-9
static inline int Step_gmin(double Va, PhotoTransistor M, double T, double *V1, double *V2, int maxiter)
{
	/* Gmin stepping
	   Input:
	  	Va	- Applied bias voltage
	  	M	- Model parameters
		maxiter - muximum number or iterations
	   Output:
	  	V1,V2	- updated node voltages
	   return value:
	   	0 	- solution is not converges
		1	- solutiuon did converge
	
	*/
	
	/* Gmin stepping in cases convergence is poor*/
	double gmin, dj, gmax;
	int GminIter=0, conv=0;
	double a, E;
	double V1ref, V2ref, GminMin;
	
	/* initiate node voltages, discars current values as we only go here in case of problems */
	*V1=Va;
	*V2=Va;
	
	/* compute the system to get an appropriate value for gmin (gmax) */
	System(Va, M, *V1, *V2, 0, T, &V1ref, &V2ref, &E, &dj, &gmax);
	V1ref=*V1; /* reference voltages contain the last trusted solution */
	V2ref=*V2; /* reference voltages contain the last trusted solution */
	
	*V1=V1ref;
	*V2=V2ref;
	
	gmin=gmax; /* gmax is the maximum value on the diagonal, as gmin is placed from every node to ground, this gmin value is bound to make the system more stable */
	GminMin=2*gmin;  /* gmin value of the last trusted solution, only update voltages if the current gmin value is lower and the solution converged */
	
	a=pow(gmin/EPSGMIN, 1.0/NOGMIN); /* rate of stepping gmin */
	while((GminIter<NGMIN)&&((gmin>EPSGMIN)||(conv==0)))
	{
		conv=Solve_gmin(Va, V1, V2, M, gmin, T, &dj, &gmax, &E, EPS/10, maxiter);
		if (conv==0)
		{
			a=sqrt(a); /* gmin stepping not going so well, slow rates down */
			gmin*=a;   /* increase gmin */	
			*V1=V1ref; /* discard last solution */	
			*V2=V2ref; /* discard last solution */
		}
		else
		{
			if (GminMin>gmin)	/* only then we can have a better solution */
			{
				a=a*sqrt(a);  /* gmin stepping going well, try speeding up */
				V1ref=*V1;	/* update reference solution */
				V2ref=*V2;
				GminMin=gmin;	/* update the gmin value for the reference solution */
			}
			gmin/=a;	
		}
		GminIter++;			
	}
	/* Vref should now be close to the final solution */
	gmin=0;
	conv=Solve_gmin(Va, V1, V2, M, gmin, T, &dj, &gmax, &E, EPS, MAXITER);
	if (!conv)
	{
		/* return the best known solution */
		*V1=V1ref;
		*V2=V2ref;
	}
	return conv;
}


#define RSMIN 1e-9
static inline void BracketV2(double Va, double *V1, double *V2, PhotoTransistor M, double T)
{

	/* often convergence is poor when Rs is high In this routine we set Rs to a low value to solve the transitor part (using newton as before)*/
	/* we can then use simple bracketing to solve V2, slow but robust */
	/*
	   Input:
	  	Va	- Applied bias voltage
	  	M	- Model parameters
	  	V1,V2	- initial node voltages
	   Output:
	  	V1,V2	- updated node voltages
	*/
	
	int conv=1;
	double dVmax, dj, gmax;
	double v1, v2, v, vmax, emax, vmin, emin, E, v1min, v1max;
	PhotoTransistor M2;	
	
	if (M.Rs<RSMIN)
	{
		/* convergence problesm not due to series resistance, go for gmin stepping alone */
		conv=Step_gmin(Va, M, T, V1, V2, MAXITER_2);
		return;
	}
	/* Bracketing the root: */
	/* copy parameters, set low resistance value */
	M2=M;
	M2.Rs=RSMIN;
	
	dVmax=EPS; /* break off for bracketing */
	
	/* V2 must be between the Voc of the phototransistor and the shunt and the applied voltage */
	/* If Va is larger than 0V the maximum V2 is Va+Rs*Jph */
	/* if Va is negative the maximum V2 is the maximum Voc, We set this to Jph*Rsh */ 
	vmax=Va+M.Rs*M.Jph;
	if (Va<0)
		vmax=M.Rsh*M.Jph; 
	
	/* The smallest possible Voc is 0V. The minimum possible voltage thus lies either at 0V, in case Va is positive */
	/* or at Va */
	vmin=MIN(Va,0);
	
	v1=0;
	v2=vmax;	
	conv=Solve_gmin(vmax, &v1, &v2, M2, 0, T, &dj, &gmax, &E, EPS, MAXITER);
	if (conv==0)
		conv=Step_gmin(vmax, M2, T, &v1, &v2, MAXITER_2); 
	/* sunstract current through the series resistzance from the current through the transistor and shunt, should be 0 if we find the solution */
	emax=(vmax-v2)/M2.Rs-(Va-v2)/M.Rs;
	/* store v1 value for better initial guesses later */
	v1max=v1;
	
	
	/* minimum v2 is minimum of Va and 0 (no negative photocurrents allowed) */
		
	conv=Solve_gmin(vmin, &v1, &v2, M2, 0, T, &dj, &gmax, &E, EPS, MAXITER);
	if (conv==0)
		conv=Step_gmin(vmin, M2, T, &v1, &v2, MAXITER_2);
	emin=(vmin-v2)/M2.Rs-(Va-v2)/M.Rs;
	/* store v1 value for better initial guesses later */
	v1min=v1;
	
	if (emin*emax>0)
	{
		/* OK this really is not going well
		fprintf(stderr, "This should not happen, the bracketing interval does not contain root in BracketV2\n");
		fprintf(stderr, "If this happens please infrom the PVMOS developer(s) as it should not! \n");
		fprintf(stderr, "However, no worries, I'll just adapt the range for you\n"); */
		v=0;
		conv=Solve_gmin(v, &v1, &v2, M2, 0, T, &dj, &gmax, &E, EPS, MAXITER);
		E=(v-v2)/M2.Rs-(Va-v2)/M.Rs;
		if (E<0)
		{
			emax=E;
			vmax=v;
			while(emax<0)
			{
				/* vmax is too small */
				vmin=v;
				v1min=v1;
				emin=emax;
				
				vmax+=0.5;
				conv=Solve_gmin(vmax, &v1, &v2, M2, 0, T, &dj, &gmax, &E, EPS, MAXITER);
				if (conv==0)
					conv=Step_gmin(vmax, M2, T, &v1, &v2, MAXITER_2);
				emax=(vmax-v2)/M2.Rs-(Va-v2)/M.Rs;
				v1max=v1;
			}
		}
		else
		{
			emin=E;
			vmin=v;
			while(emin>0)
			{
				/* vmax is too small */
				vmax=v;
				v1max=v1;
				emax=emin;
				
				vmin-=0.5;
				conv=Solve_gmin(vmin, &v1, &v2, M2, 0, T, &dj, &gmax, &E, EPS, MAXITER);
				if (conv==0)
					conv=Step_gmin(vmin, M2, T, &v1, &v2, MAXITER_2);
				emin=(vmin-v2)/M2.Rs-(Va-v2)/M.Rs;
				v1min=v1;
			}
		}
	}
	/* simple bracketing the root for v2 */
	while (fabs(vmax-vmin)>dVmax)
	{
		v=(vmax+vmin)/2; /* new V2 estimate */
		v2=v;
		v1=(v1max+v1min)/2; /* new v1 estimate */
		conv=Solve_gmin(v, &v1, &v2, M2, 0, T, &dj, &gmax, &E, EPS, MAXITER);
		if (conv==0)
			conv=Step_gmin(v, M2, T, &v1, &v2, MAXITER_2);
		E=(v-v2)/M2.Rs-(Va-v2)/M.Rs;
		if (E*emax>0)
		{
			vmax=v;
			v1max=v1;
		}
		else
		{
			vmin=v;
			v1min=v1;
		}
	}
	*V1=v1;
	*V2=v;
}

void Phototransistor(double Va, PhotoTransistor M, double T, double *V1, double *V2, double *J, double *dJdVa)
{
	/* Solve the circuit at applied bias of Va
	 * Input:
	 *	Va	- Applied bias voltage
	 *	M	- Model parameters */
	/*   Output:
	  	V1,V2	- updated node voltages 
		J	- current density through the circuit
		dJdVa	- derivative of the above to the applied voltage
		
		*/
	double gmin, dj, gmax;
	int conv;
	double err;
	gmin=0;
	double v1, v2;
	/* fastest results will be obtained if you comment these two lines out, it will resuse the previous solution as a first guess to the next */
	v1=Va/2;
	v2=Va/2;
	
	/* this  routine accepts NULL vectors to the current or derivative values */
	if (!dJdVa)
		dJdVa=&dj;
		
	/* fisrt shot, does it converge? */
	conv=Solve_gmin(Va, &v1, &v2, M, gmin, T, dJdVa, &gmax, &err, EPS, MAXITER);
	if (conv==0)
	{
		/* newton raphson did not do it, try something more robust */
		double dv1, dv2;
		BracketV2(Va, &v1, &v2, M, T);
		System(Va, M, v1, v2, 0, T, &dv1, &dv2, &err, dJdVa, &gmax);
		if (err>EPS)	
			conv=Solve_gmin(Va, &v1, &v2, M, gmin, T, dJdVa, &gmax, &err, EPS, MAXITER_2);
		else
			conv=1;
			
	}
	if (!conv)
		fprintf(stderr, "Warning: solution did not converge. KCL error: %e\n", err);
	if (J)
		*J=(Va-(v2))/M.Rs;
	if (V1)
		*V1=v1;
	if (V2)
		*V2=v2;
}

