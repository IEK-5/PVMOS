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
typedef struct NumSet {
	double tol_kcl_abs, tol_kcl_rel, tol_v_abs, tol_v_rel;
	int max_iter, N_lin_search; 
	int GminStep;
	int GminStart;
	double GminMax;
	double GminFac;
	double Gmin;
} NumSet; 

extern NumSet Numeric_Settings; /* initialized in solve.c */
                                                                          
void BubbleSortJV(int n, double *V, double *J);
void Resistance(mesh M, node N1, node N2, double *R);
void Diode(mesh M, node N, int inter_index, double V, double *I, double *dIdV, double *Vj);
int FindVa(double Va, double *list, int Nva);
void SolveVa(mesh *M, double Vstart, double Vend, int Nstep);
int SolPar(mesh *M, int *isc, int *imp_m, int *imp, int *imp_p, int *ioc_m, int *ioc_p);
void RefineOC(mesh *M, double tol_i, double tol_v, int Niter);
void RefineMPP(mesh *M, double tol_i, double tol_v, int Niter);
void AdaptiveSolveVa(mesh *M, double Va, double rel_threshold, int N);
double *LocalyCollectedCurrent(mesh *M, double Va, int diode_index, int *nodes, int diff, int NL);
