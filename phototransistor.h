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
typedef struct PhotoTransistor { 
	double Jsbc;	/*  A cm-2  	 	Saturation current density of the main BC junction */
	double Jsbe;	/*  A cm-2  	 	Saturation current density of the main BE junction */
	double Jse;	/*  A cm-2  	 	Saturation current density of the "leakage" BE junction */
	double Jsc;	/*  A cm-2  	 	Saturation current density of the "leakage" BC junction */
	double Nf;	/*  -     	 	Idiality factor forward operation (BE junction) */
	double Nr;	/*  -     	 	Idiality factor reverse operation (BC junction) */
	double Ne;	/*  -     	 	Idiality factor leakage diode BE junction */
	double Nc;	/*  -     	 	Idiality factor leakage diode BC junction */
	double Vaf;	/*  V     	 	forward operation Early Voltage */
	double Var;	/*  V     	 	reverse operation Early Voltage */
	double Bf;	/*  -     	 	forward beta */
	double Jph;	/*  A cm-2     	 	Photocurrent BE junction */
	double Rs;	/*  Ohm cm2    	 	series resistance */
	double Rsh;	/*  Ohm cm2    	 	shunt resistance */
	double EgBE;	/*  eV 	 	 	Bandgap for temperature dependencies BE junction */
	double PhiBC;	/*  eV 	 	 	Barrier height for temperature dependencies BC junction */
	double XTIBE;	/*  -          	 	Saturation Current temperature dependency  */
	double XTIBC;	/*  -          	 	Saturation Current temperature dependency  */
	double XTB;	/*  -          	 	Beta temperature dependency	 */
} PhotoTransistor;

void * InitPhotoTransistorStruct(int *size);
void Phototransistor(double Va, PhotoTransistor M, double T, double *V1, double *V2, double *J, double *dJdVa);
