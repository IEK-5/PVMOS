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
/* Data Structures */
typedef enum {JVD,ONED,TWOD, PHOTOT} diode_model;
/* struct describing the connection between two electrodes */

typedef struct ElConn {
	diode_model model;
	void * ParStruct;
	int ParSize;	/* using the type size_t will make the 32 bit version incompatible with the 64 bit version */
	double *V, *J; /* voltage versus current density A/cm^2*/
	int N; /* number of current-voltage pairs */
} ElConn;

/* struct with all local properties i.e. decribing the solar cell, sheet resistances, etc */
typedef struct local_prop {
	char * name;
	double *Rel; /* array of electrode sheet resistances */
	double *Rvp, *Rvn, *Rvg; /* resistance (in Ohm cm^2) of the electrodes to the positive (vp) /negative (vn) terminal  (negative values indicate no connection) */
	ElConn *conn;
	double T;
	double Vg;
	int SplitX, SplitY; /* whether the node can be split during adaptive meshing, in x and y direction  */
} local_prop ;

/* Result struct */
typedef struct results {
	double *Va;
	double *I;
	double ***Vn;
	int Nva;	
} results ;

/* the node struct */
typedef struct node {
	int *north;		/* Adjacent nodes north */
	int *south;		/* Adjacent nodes south */
	int *west;		/* Adjacent nodes west */
	int *east;		/* Adjacent nodes east */
	int id;			/* this node's id number */
	double x1,y1,x2,y2;
	int P;			/* local properties index */
} node;

/* mesh struct, i.e. a collection of nodes, local properties and results */
typedef struct mesh {
	node * nodes;			/* nodes */
	int Nn;				/* number of nodes */
	int Nel;			/* number of electrodes */
	local_prop *P;			/* local properties */
	int Na;				/* Number of area's with differing local properties */
	results res;			/* struct with calculated data for the mesh */
} mesh;

#ifndef _HAS_MESHHASH
extern int NMESHHASH;
extern unsigned const char MESHHASH[];
#endif
local_prop InitProperties (char *name, int Nel);
mesh InitMesh(char *name, double x1, double x2, double y1, double y2, int Nx, int Ny);
void FreeProperties(local_prop *P, int Nel);
void FreeMesh(mesh *M);
void DuplicateProperties(mesh *M, local_prop *dest, local_prop *source);
mesh JoinMeshes(mesh M1, mesh M2, double xoff, double yoff);
mesh JoinMeshes_H(mesh M1, mesh M2, double yoff);
mesh JoinMeshes_V(mesh M1, mesh M2, double xoff); 
void AddRowNorth(mesh *M, double dy);
void AddColEast(mesh * M, double dx);
void AddRowSouth(mesh *M, double dy);
void AddColWest(mesh * M, double dx);
mesh DuplicateMesh(mesh M);
void DuplicateNode(mesh M, node *d, int source_id);
node *SearchNode(mesh M, int id);
int FindProperties(mesh M, char *name);
void NewProperties(mesh * M, char *name);
int DeleteUnusedProperties(mesh *M);
void AssignProperties(mesh *M, int *select, int P);
void AssignPropertiesMesh(mesh *M, int P);
void AddElectrode(mesh *M);
void PurgeResults(mesh *M);
void PurgeResultAtIndex(mesh *M, int index);
void SplitNodeX(int id, mesh *M);
void SplitNodeY(int id, mesh *M);
void SplitNodeXY(int id, mesh *M);
void SplitMeshX(mesh *M);
void SplitMeshY(mesh *M);
void SplitMeshXY(mesh *M);
void SplitMeshLong(mesh *M);
void SplitMeshWhileCoarse(mesh *M, double dx, double dy);
void SplitListX(mesh *M, int *list);
void SplitListY(mesh *M, int *list);
void SplitListXY(mesh *M, int *list);
void SplitListLong(mesh *M, int *list);
void SplitListWhileCoarse(mesh *M, int *list, double dx, double dy);
void SplitXY_Ntimes(mesh *M, int id, int Nx, int Ny);
void Mesh_ScaleMove(mesh *M, double fx, double fy, double dx, double dy);
void RotateMesh(mesh *M, double x, double y, int d);
void GetMeshBB(mesh *M, double *x1, double *y1, double *x2, double *y2);
void SetMeshBB(mesh *M, double x1, double y1, double x2, double y2, int FixR, double *fx, double *fy, double *dx, double *dy);
void CleanUpMesh(mesh *M, int *merged);
void Chunkify(mesh *M);
void WriteMesh(char *fn, mesh *M);
void ReadMesh(char *fn, mesh *M);

