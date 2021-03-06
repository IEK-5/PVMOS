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
 #include "polygon.h"                                                                           
int FindPos(mesh M, int id, double x, double y, int *NOTINMESH);
int * PolySelectNodes(polygon P, mesh M, int *sel_nodes);
int * PolyContourSelectNodes(double d, polygon P, int loop, mesh M, int *sel_nodes);
int * CircSelectNodes(double x, double y, double r, mesh M, int *sel_nodes);
int * CircContourSelectNodes(double x, double y, double r, double d, mesh M, int *sel_nodes);
int * RectSelectNodes(double x1, double y1, double x2, double y2, mesh M, int *sel_nodes);
int * RectContourSelectNodes(double x1, double y1, double x2, double y2, double d, mesh M, int *sel_nodes);
int * SelectArea(mesh M, int *sel_nodes, char *name);
void ResolvContour(polygon P, mesh *M, int loop, double D);
void ResolvCircle(double x, double y, double r, mesh *M, double D);
void ResolvRect(double x1, double y1, double x2, double y2, mesh *M, double D);
