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
#define k_b 8.6173323478e-5
#define T0 300
#define EPS 1e-6
#define MAXITER 50
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
void TwoDiode(double V, double J01, double J02, double Eg, double T, double Jph, double Rs, double Rsh, double *I, double *dI, double *Vj)
{
	double V1, lI, E, dIdV;
	int iter=0;
	V1=V;
	do
	{
		lI=D(V1,J01,Eg,T,1.0,Jph)+D(V1,J02,Eg,T,2.0,0.0)+V1/Rsh;
		dIdV=dDdV(V1, J01, Eg, T, 1.0)+dDdV(V1, J02, Eg, T, 2.0)+1/Rsh;
		E=((V-V1)/Rs-lI);
		V1+=((V-V1)/Rs-lI)/(dIdV+1/Rs); /* newton-rapson */
		iter++;
	} while ((fabs(E/lI)>EPS)&&(iter<MAXITER));
	if (dI)
		(*dI)=1/((1/dIdV)+Rs);
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

void OneDiode(double V, double J0,double n, double Eg, double T, double Jph, double Rs, double Rsh, double *I, double *dI, double *Vj)
{
	double V1, lI, E, dIdV;
	int iter=0;
	V1=V;
	do
	{
		lI=D(V1,J0,Eg,T,n,Jph)+V1/Rsh;
		dIdV=dDdV(V1, J0, Eg, T, n)+1/Rsh;
		E=((V-V1)/Rs-lI);
		V1+=((V-V1)/Rs-lI)/(dIdV+1/Rs); /* newton-rapson */
		iter++;
	} while ((fabs(E/lI)>EPS)&&(iter<MAXITER));
	if (dI)
		(*dI)=1/((1/dIdV)+Rs);
	if (I)
		(*I)=lI;
	if (Vj)
		(*Vj)=V1;
}
