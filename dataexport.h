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

void PrintMesh(char *fn, mesh *M, int *selected);
void PrintConn(char *fn, mesh *M, int *selected);
void PrintSurfDef(char *fn, mesh *M, int *selected);
void PrintSurfV(char *fn, mesh *M, int *selected);
void PrintSolPar(char *fn, mesh *M);
void PrintIV(char *fn, mesh *M);
void PrintProbe(char *fn, mesh *M, double x, double y);
void PrintInIp(char *fn, mesh *M, int *selected);
void PrintPars(char *fn, mesh *M);
/*void SurfVPlotNearest(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny);*/
void SurfVPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny);
void SurfVjPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny);
void SurfPPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny);
void SurfJPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny);
void SurfEPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny);
void PrintLocallyCollectedCurrent(char *fn, mesh *M, double x1, double y1, double x2, double y2, int Nx, int Ny, double Va, int diode_index, int diff, int NL);
void PrintLocalJV(char *fn, mesh M, double x, double y, int inter_index, double Vstart, double Vend, int Nstep);
