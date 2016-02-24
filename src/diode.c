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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "diode.h"
#define k_b 8.6173323478e-5
#define T0 300
#define EPS 1e-6
#define MAXITER 50


void * InitOneTwoDiodeStruct(int *size)
{
	OneTwoDiode *res;
	*size=sizeof(OneTwoDiode);
	res=malloc(sizeof(OneTwoDiode));
	res->J01=1e-12;
	res->J02=1e-8;
	res->Jph=0;
	res->nid1=1;
	res->nid2=2;
	res->Eg=1.12;
	res->Rs=1e-5;
	res->Rsh=1e4;
	return (void *) res;
}
/* generic diode equation */
double D(double V, double J0, double Eg, double T, double n, double Jph)
{
	return J0*(T*T*T/(T0*T0*T0))*exp(-(T0-T)*Eg/(n*k_b*T*T0))*(exp(V/(k_b*T*n))-1)-Jph;
}

/* derivative */
double dDdV(double V, double J0, double Eg, double T, double n)
{
	return J0*(T*T*T/(T0*T0*T0))*exp(-(T0-T)*Eg/(n*k_b*T*T0))*exp(V/(k_b*T*n))/(k_b*T*n);
}

/* Two Diode Model, note the fixed ideality factors */
/*                                                                                 
                                 +----+                                         
       +------+------+------+----|    |---o                                     
   .   |      |      |      |    +----+                                         
  /|\ _|_   __|__  __|__   _|_     Rp                                           
   | /___\   / \    / \   |   |                                                 
Jph| \_ _/  /_ _\  /_ _\  |_ _|                                                 
   |   |      |      |      |                                                   
       |      |      |      |                                                   
       +------+------+------+-------------o                                     
             n=1    n=2    Rsh                                                  
             J01    J02                                                                                                           
*/
void TwoDiode(double V, OneTwoDiode Pars, double T, double *I, double *dI, double *Vj)
{
	double V1, lI, E, dIdV;
	int iter=0;
	V1=V;
	do
	{
		lI=D(V1,Pars.J01,Pars.Eg,T,1.0,Pars.Jph)+D(V1,Pars.J02,Pars.Eg,T,2.0,0.0)+V1/Pars.Rsh;
		dIdV=dDdV(V1, Pars.J01, Pars.Eg, T, 1.0)+dDdV(V1, Pars.J02, Pars.Eg, T, 2.0)+1/Pars.Rsh;
		E=((V-V1)/Pars.Rs-lI);
		V1+=((V-V1)/Pars.Rs-lI)/(dIdV+1/Pars.Rs); /* newton-rapson */
		iter++;
	} while ((fabs(E/lI)>EPS)&&(iter<MAXITER));
	if (dI)
		(*dI)=1/((1/dIdV)+Pars.Rs);
	if (I)
		(*I)=lI;
	if (Vj)
		(*Vj)=V1;
}
/* One Diode Model, note the non-fixed ideality factor */
/*                                                                                 
                            +----+                                         
       +------+--------+----|    |---o                                     
   .   |      |        |    +----+                                         
  /|\ _|_   __|__     _|_     Rp                                           
   | /___\   / \     |   |                                                 
Jph| \_ _/  /_ _\    |_ _|                                                 
   |   |      |        |                                                   
       |      |        |                                                   
       +------+--------+-------------o                                     
             n        Rsh                                                  
             J0                                                                                                                
*/

void OneDiode(double V, OneTwoDiode Pars, double T, double *I, double *dI, double *Vj)
{
	double V1, lI, E, dIdV;
	int iter=0;
	V1=V;
	do
	{
		lI=D(V1,Pars.J01,Pars.Eg,T,Pars.nid1,Pars.Jph)+V1/Pars.Rsh;
		dIdV=dDdV(V1, Pars.J01, Pars.Eg, T, Pars.nid1)+1/Pars.Rsh;
		E=((V-V1)/Pars.Rs-lI);
		V1+=((V-V1)/Pars.Rs-lI)/(dIdV+1/Pars.Rs); /* newton-rapson */
		iter++;
	} while ((fabs(E/lI)>EPS)&&(iter<MAXITER));
	if (dI)
		(*dI)=1/((1/dIdV)+Pars.Rs);
	if (I)
		(*I)=lI;
	if (Vj)
		(*Vj)=V1;
}
