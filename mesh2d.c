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
 * routines to create and modify meshes a variable rectangular   *
 * mesh suitable for finite-differences                          *      
 * NOTE: The meshing routines are basically in 2D, however, the  *           
 *       PVMOS solver is 3D. The limitation is that in the       * 
 *       z-dimension  the mesh is regular (i.e. he 3rd dimension * 
 *       exists as a certain number of layers in the mesh)       *               
 *                                                               *       
 *****************************************************************/     
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "meshhash.h"
#include "mesh2d.h"
#include "main.h"
#include "list.h"
#include "utils.h"
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))
#define TINY 1e-10
#define TWOPI 6.28318530717959
#define MAXNAMELEN 255
#define T0 300

local_prop InitProperties (char *name, int Nel)
{
	local_prop P;
	int i;
	/* local properties struct */
	P.name=malloc(MAXNAMELEN*sizeof(char));
	if (strlen(name)+1<MAXNAMELEN)
		strncpy(P.name, name,strlen(name)+1);
	else
	{
		strncpy(P.name, name,MAXNAMELEN-1);
		P.name[MAXNAMELEN]='\0';
	}
	
	P.Rel=malloc((Nel+1)*sizeof(double));
	P.Rvp=malloc((Nel+1)*sizeof(double));
	P.Rvn=malloc((Nel+1)*sizeof(double));	
	P.conn=malloc((Nel)*sizeof(ElConn));
	
	for (i=0;i<Nel-1;i++)
	{
		P.Rel[i]=1;
		P.Rvp[i]=-1.0;
		P.Rvn[i]=-1.0;
		
		P.conn[i].model=JVD;
		P.conn[i].ParStruct=NULL;
		P.conn[i].ParSize=0;
		
		P.conn[i].V=malloc(2*sizeof(double));
		P.conn[i].J=malloc(2*sizeof(double));
		P.conn[i].V[0]=-1;
		P.conn[i].J[0]=-1;
		P.conn[i].V[1]=1;
		P.conn[i].J[1]=1;
		P.conn[i].N=2;		
	}
	
	P.Rel[i]=1;
	P.Rvp[i]=-1.0;
	P.Rvn[i]=-1.0;	
	
	P.T=300;
	P.SplitX=1;
	P.SplitY=1;
	return P;
}

mesh InitMesh(char *name, double x1, double x2, double y1, double y2, int Nx, int Ny)
/* returns a regularly spaced mesh with Nx columns and Ny rows*/
/* The coordinates x1,y1 represents the south-west corner (lower left) of the mesh and the coordinates x2,y2 the north-east corner (upper right)*/
{
	mesh M;
	int Nn, i, j;
	double x_step, y_step, xx1, xx2;
	
	/* results struct */
	M.res.Nva=0;
	M.res.Va=malloc(sizeof(double));
	M.res.I=malloc(sizeof(double));
	M.res.Vn=malloc(sizeof(double **));
	
	
	/* local properties struct */
	M.P=malloc((2)*sizeof(local_prop));
	M.Na=1;	
	M.P[0]=InitProperties (name,2);
	
	/* nodes */
	Nn=Nx*Ny;
	M.Nn=Nn;
	M.Nel=2;	/* two electrodes per default */
	M.nodes=malloc((Nn+1)*sizeof(node));
	x_step=(x2-x1)/Nx;
	y_step=(y2-y1)/Ny;
	Nn=0;
	for (i=0;i<Nx;i++)
	{
		xx1=x1+i*x_step;
		xx2=xx1+x_step;
		for (j=0;j<Ny;j++)
		{
			M.nodes[Nn].x1=xx1;
			M.nodes[Nn].x2=xx2;
			M.nodes[Nn].y1=y1+j*y_step;
			M.nodes[Nn].y2=M.nodes[Nn].y1+y_step;
			M.nodes[Nn].id=Nn;
			if (i!=0)
			{
				M.nodes[Nn].west=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].west[0]=1;
				M.nodes[Nn].west[1]=(i-1)*Ny+j;
			}
			else
			{
				M.nodes[Nn].west=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].west[0]=0;
			}
			if (i!=Nx-1)
			{
				M.nodes[Nn].east=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].east[0]=1;
				M.nodes[Nn].east[1]=(i+1)*Ny+j;
			}
			else
			{
				M.nodes[Nn].east=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].east[0]=0;
			}
			if (j!=0)
			{
				M.nodes[Nn].south=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].south[0]=1;
				M.nodes[Nn].south[1]=i*Ny+j-1;
			}
			else
			{
				M.nodes[Nn].south=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].south[0]=0;
			}
			if (j!=Ny-1)
			{
				M.nodes[Nn].north=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].north[0]=1;
				M.nodes[Nn].north[1]=i*Ny+j+1;
			}
			else
			{
				M.nodes[Nn].north=malloc(LISTBLOCK*sizeof(int));
				M.nodes[Nn].north[0]=0;
			}
			
			M.nodes[Nn].P=0; /* all nodes link to default proprties */
			Nn++;		
		}
	}
	return M;
}
void FreeProperties(local_prop *P, int Nel)
{
	int j;
	free(P->name);
	free(P->Rel);
	free(P->Rvp);
	free(P->Rvn);
	for (j=0;j<Nel-1;j++)
	{
		free(P->conn[j].V);
		free(P->conn[j].J);
		if (P->conn[j].ParStruct)
			free(P->conn[j].ParStruct);
	}
	free(P->conn);
}
void FreeMesh(mesh *M)
{
	int i,j;
	for (i=0;i<M->Nn;i++)
	{
		free(M->nodes[i].north);
		free(M->nodes[i].south);
		free(M->nodes[i].east);
		free(M->nodes[i].west);
	}
	free(M->nodes);
	for (i=0;i<M->Na;i++)
		FreeProperties(M->P+i, M->Nel);
	free(M->P);
	for (i=0;i<M->res.Nva;i++)
	{
		for (j=0;j<M->Nel;j++)
			free(M->res.Vn[i][j]);
		free(M->res.Vn[i]);
	}
	free(M->res.Vn);
	free(M->res.Va);
	free(M->res.I);
}

node *SearchNode(mesh M, int id)
/* finds a node id and returns a pointer to the node */
{
	int i=0, min, max;
	/* first shot, node-id corresponds to node position*/
	/* For robustness it is good to not require the node-id  is the same as the node position but for speed it is definately to be preferred */
	/* For sorting od the mesh is is absolutely required that this routine always works */
	if ((id<M.Nn)&&(id>0))
		if (M.nodes[id].id==id)
			return M.nodes+id;
	
	/* fallback, still assume the mesh is ordered (i.e. node id's increase with index) */
	min=0;
	max=M.Nn-1;
	while(max-min>1)
	{
		i=(min+max)/2;
		if (M.nodes[i].id==id)
			break;
		if (M.nodes[i].id<id)
			min=i;
		else
			max=i;
	}
	if (M.nodes[i].id==id)
		return M.nodes+i;	
	if (M.nodes[min].id==id)
		return M.nodes+min;
	if (M.nodes[max].id==id)
		return M.nodes+max;
	
	Print(DEBUG,"Chaos reighns the mesh!\nStarting a slow, exhaustive search for the requested node");	
	/* we can of course always search through the whole list  in case you make a mess of your mesh */
	while (i<M.Nn)
	{
		if (M.nodes[i].id==id)
			break;
		i++;
	}
	if (i==M.Nn)
	{
		Warning("node %d not present in the mesh!\n", id);
		return NULL;	
	}
	return M.nodes+i;
}

void SearchLink(mesh M, int id)
/* finds all links to a node */
{
	int i;
	for (i=0;i<M.Nn;i++)
	{
		if (IsInList(M.nodes[i].north, id))
			Print(DEBUG, "node %d linked by %d.north",id, M.nodes[i].id);
		if (IsInList(M.nodes[i].south, id))
			Print(DEBUG, "node %d linked by %d.south",id, M.nodes[i].id);	
		if (IsInList(M.nodes[i].west, id))
			Print(DEBUG, "node %d linked by %d.west",id, M.nodes[i].id);
		if (IsInList(M.nodes[i].east, id))
			Print(DEBUG, "node %d linked by %d.east",id, M.nodes[i].id);
	}
}

double R_Overlap(double a1, double b1, double a2, double b2)
{
	return (MIN(a2,b2)-MAX(a1,b1))/MIN(a2-a1,b2-b1);
}

void NodeConnector(node *N1, node *N2)
{
	/* first remove all connections */
	N1->north=RemoveFromList(N1->north, N2->id);
	N1->south=RemoveFromList(N1->south, N2->id);
	N1->west=RemoveFromList(N1->west, N2->id);
	N1->east=RemoveFromList(N1->east, N2->id);
	
	N2->north=RemoveFromList(N2->north, N1->id);
	N2->south=RemoveFromList(N2->south, N1->id);
	N2->west=RemoveFromList(N2->west, N1->id);
	N2->east=RemoveFromList(N2->east, N1->id);
	
	/* now connect as appropriate */
	if ((fabs(N1->y2-N2->y1)/MIN(N1->y2-N1->y1,N2->y2-N2->y1)<TINY) && (R_Overlap(N1->x1,N2->x1, N1->x2, N2->x2) > TINY))
	{
		N1->north=AddToList(N1->north, N2->id);
		N2->south=AddToList(N2->south, N1->id);
	}
	if ((fabs(N1->x2-N2->x1)/MIN(N1->x2-N1->x1,N2->x2-N2->x1)<TINY) && (R_Overlap(N1->y1,N2->y1, N1->y2, N2->y2) > TINY))
	{
		N1->east=AddToList(N1->east, N2->id);
		N2->west=AddToList(N2->west, N1->id);
	}
	if ((fabs(N1->y1-N2->y2)/MIN(N1->y2-N1->y1,N2->y2-N2->y1)<TINY) && (R_Overlap(N1->x1,N2->x1, N1->x2, N2->x2) > TINY))
	{
		N1->south=AddToList(N1->south, N2->id);
		N2->north=AddToList(N2->north, N1->id);
	}
	if ((fabs(N1->x1-N2->x2)/MIN(N1->x2-N1->x1,N2->x2-N2->x1)<TINY) && (R_Overlap(N1->y1,N2->y1, N1->y2, N2->y2) > TINY))
	{
		N1->west=AddToList(N1->west, N2->id);
		N2->east=AddToList(N2->east, N1->id);
	}
}

int * MeshOutline(mesh M)
/* create a list of nodes which lie along the edges of the mesh */
{
	int *list;
	int i;
	list=malloc(LISTBLOCK*sizeof(int));
	list[0]=0;	
	
	for (i=0;i<M.Nn;i++)
		if ((M.nodes[i].north[0]==0) || (M.nodes[i].south[0]==0) || (M.nodes[i].east[0]==0) || (M.nodes[i].west[0]==0))
			list=AddToList(list, M.nodes[i].id);
	return list;
}


void DuplicateNode(mesh M, node *d, int source_id)
/* duplicate a node */
{
	node *s;
	s=SearchNode(M, source_id);	
	d = memcpy(d, s, sizeof(node));	
	d->north=DuplicateList(s->north);	
	d->south=DuplicateList(s->south);
	d->west=DuplicateList(s->west);	
	d->east=DuplicateList(s->east);
}

void DuplicateProperties(mesh *M, local_prop *dest, local_prop *source)
/* duplicate local properties struct */
{
	int i, j;
	dest = memcpy(dest, source, sizeof(local_prop));	
	dest->name=malloc((strlen(source->name)+1)*sizeof(char));
	strncpy(dest->name, source->name, strlen(source->name)+1);
	dest->Rel=malloc((M->Nel+1)*sizeof(double));
	dest->Rvp=malloc((M->Nel+1)*sizeof(double));
	dest->Rvn=malloc((M->Nel+1)*sizeof(double));
	dest->conn=malloc((M->Nel)*sizeof(ElConn));
	for (j=0;j<M->Nel-1;j++)
	{
		dest->Rel[j]=source->Rel[j];		
		dest->Rvp[j]=source->Rvp[j];
		dest->Rvn[j]=source->Rvn[j];
		dest->conn[j].model=source->conn[j].model;
		dest->conn[j].ParSize=source->conn[j].ParSize;
		if (source->conn[j].ParStruct)
		{
			dest->conn[j].ParStruct=malloc(dest->conn[j].ParSize);
			memcpy(dest->conn[j].ParStruct, source->conn[j].ParStruct, source->conn[j].ParSize);
		}
		else
			dest->conn[j].ParStruct=NULL;
		dest->conn[j].N=source->conn[j].N;
		dest->conn[j].V=malloc((source->conn[j].N+1)*sizeof(double));
		dest->conn[j].J=malloc((source->conn[j].N+1)*sizeof(double));
		for (i=0;i<source->conn[j].N;i++)
		{
			dest->conn[j].V[i]=source->conn[j].V[i];
			dest->conn[j].J[i]=source->conn[j].J[i];
		}
	}
	dest->Rel[j]=source->Rel[j];		
	dest->Rvp[j]=source->Rvp[j];
	dest->Rvn[j]=source->Rvn[j];
}

void DuplicateResults(mesh M, results *dest, results *source)
/* duplicate results struct */
{
	int i, j, k;	
	dest = memcpy(dest, source, sizeof(results));	
	dest->Va=malloc((dest->Nva+1)*sizeof(double));
	dest->I=malloc((dest->Nva+1)*sizeof(double));
	dest->Vn=malloc((dest->Nva+1)*sizeof(double **));
	for (i=0;i<dest->Nva;i++)
	{
		dest->Vn[i]=malloc((M.Nel+1)*sizeof(double *));
		dest->Va[i]=source->Va[i];
		dest->I[i]=source->Va[i];
		for (j=0;j<M.Nel;j++)
		{
			dest->Vn[i][j]=malloc((M.Nn+1)*sizeof(double));
			for (k=0;k<M.Nn;k++)
				dest->Vn[i][j][k]=source->Vn[i][j][k];
		
		}
		
	}
}

mesh DuplicateMesh(mesh M)
/* duplicate a mesh */
{
	mesh res;
	int i;
	res.Nn=M.Nn;
	res.Nel=M.Nel;
	res.nodes=malloc(res.Nn*sizeof(node));
	for (i=0;i<res.Nn;i++)
		DuplicateNode(M,res.nodes+i, M.nodes[i].id);
	res.Na=M.Na;
	res.P=malloc((res.Na+1)*sizeof(local_prop));
	for (i=0;i<res.Na;i++)
		DuplicateProperties(&res, res.P+i, M.P+i);
	DuplicateResults(M, &(res.res), &(M.res));
	
	return res;

}

int FindProperties(mesh M, char *name)
{
	int i=0, l=0;
	
	l=strlen(name);
	while (i<M.Na)
	{
		if (strlen(M.P[i].name)==l)
			if (strncmp(M.P[i].name, name, l)==0)
				break;
		i++;
	}

	if (i<M.Na)
		return i;
	/* flag unavailable local properties */
	return -1;
}

void NewProperties(mesh *M, char *name)
{
	if (FindProperties(*M, name)<0)
	{
		M->P=realloc(M->P, (M->Na+1)*sizeof(local_prop));
		M->Na++;
		M->P[M->Na-1]=InitProperties (name, M->Nel);	
	}
	/* else properties exist, do nothing */
}



void ReInitResults(mesh *M)
/* re-initializes the results struct. To be used when the mesh changes */
{
	int i,j;
	for (i=0;i<M->res.Nva;i++)
	{
		for (j=0;j<M->Nel;j++)
			free(M->res.Vn[i][j]);
		free(M->res.Vn[i]);
	}
	free(M->res.Vn);
	free(M->res.Va);
	free(M->res.I);
	
	M->res.Nva=0;
	M->res.Va=malloc(sizeof(double));
	M->res.I=malloc(sizeof(double));
	M->res.Vn=malloc(sizeof(double **));
}

void AssignProperties(mesh *M, int *select, int P)
{
	int j;
	node *N;
	if ((P<0)||(P>=M->Na))
		Error("Area %d not defined\n", P);
	for (j=1;j<=select[0];j++)
	{
		N=SearchNode(*M, select[j]);
		N->P=P;
	} 
	ReInitResults(M);
}

void AssignPropertiesMesh(mesh *M, int P)
{
	int j;
	if ((P<0)||(P>=M->Na))
		Error("Area %d not defined\n", P);
	for (j=0;j<M->Nn;j++)
		M->nodes[j].P=P;
	ReInitResults(M);
}

void AddElectrode(mesh *M)
{
	int i;
	M->Nel++;
	for (i=0;i<M->Na;i++)
	{
		/* extend the properties arrays to add the new electrode */
		Print(NORMAL,"Please define the properties for electrode %d in area %s", M->Nel, M->P[i].name);	
		M->P[i].Rel=realloc(M->P[i].Rel, (M->Nel+1)*sizeof(double));
		M->P[i].Rvp=realloc(M->P[i].Rvp, (M->Nel+1)*sizeof(double));
		M->P[i].Rvn=realloc(M->P[i].Rvn, (M->Nel+1)*sizeof(double));	
		M->P[i].conn=realloc(M->P[i].conn, (M->Nel)*sizeof(ElConn));
	
		/* put in default values */
		M->P[i].Rel[M->Nel-1]=1;
		M->P[i].Rvp[M->Nel-1]=-1.0;
		M->P[i].Rvn[M->Nel-1]=-1.0;
		
		M->P[i].conn[M->Nel-2].model=JVD;
		M->P[i].conn[M->Nel-2].ParStruct=NULL;
		M->P[i].conn[M->Nel-2].ParSize=0;
			
		M->P[i].conn[M->Nel-2].V=malloc(2*sizeof(double));
		M->P[i].conn[M->Nel-2].J=malloc(2*sizeof(double));
		M->P[i].conn[M->Nel-2].V[0]=-1;
		M->P[i].conn[M->Nel-2].J[0]=-1;
		M->P[i].conn[M->Nel-2].V[1]=1;
		M->P[i].conn[M->Nel-2].J[1]=1;
		M->P[i].conn[M->Nel-2].N=2;			
	}		
	
	ReInitResults(M);
}

mesh JoinMeshes(mesh M1, mesh M2, double xoff, double yoff)
{
	mesh res;
	int *outl1, *outl2;
	int *prop_out;
	node *m1, *m2;
	int i,j;
	
	if (M1.Nel!=M2.Nel)
		Error("Cannot join meshes with unequal number of electrodes, M1:%d M2:%d\n",M1.Nel, M2.Nel);
	
	outl1=MeshOutline(M1);
	outl2=MeshOutline(M2);
	
	/* TODO check for overlap */
	
	res=DuplicateMesh(M1);
	res.nodes=realloc(res.nodes,(res.Nn+M2.Nn+1)*sizeof(node));
	res.Nn+=M2.Nn;
	
	/* merge local properties, if local properties have the same name we assume they are the same! */
	prop_out=malloc((M2.Na+1)*sizeof(int));
	for (i=0;i<M2.Na;i++)
	{
		if ((j=FindProperties(res, M2.P[i].name))<0)
		{
			res.Na++;
			res.P=realloc(res.P, (res.Na+1)*sizeof(local_prop));
			DuplicateProperties(&res, res.P+res.Na-1, M2.P+i);
			prop_out[i]=res.Na-1;		
		}
		else
			prop_out[i]=j;	
		
	}
	for (i=0;i<M2.Nn;i++)
	{
		DuplicateNode(M2,res.nodes+i+M1.Nn, M2.nodes[i].id);
		res.nodes[i+M1.Nn].id+=M1.Nn;
		res.nodes[i+M1.Nn].x1+=xoff;
		res.nodes[i+M1.Nn].x2+=xoff;
		res.nodes[i+M1.Nn].y1+=yoff;
		res.nodes[i+M1.Nn].y2+=yoff;
		res.nodes[i+M1.Nn].P=prop_out[res.nodes[i+M1.Nn].P];
		for (j=1;j<=res.nodes[i+M1.Nn].north[0];j++)
			res.nodes[i+M1.Nn].north[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].south[0];j++)
			res.nodes[i+M1.Nn].south[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].east[0];j++)
			res.nodes[i+M1.Nn].east[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].west[0];j++)
			res.nodes[i+M1.Nn].west[j]+=M1.Nn;
	}
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(res, outl1[i]);
		for (j=1;j<=outl2[0];j++)
		{
			m2=SearchNode(res, outl2[j]+M1.Nn);
			NodeConnector(m1, m2);
		}			
	}
	free(outl1);
	free(outl2);
	free(prop_out);
	ReInitResults(&res);

	return res;
}
mesh JoinMeshes_H(mesh M1, mesh M2, double yoff)
{
	mesh res;
	double xoff1, xoff2, xoff;
	int *outl1, *outl2;
	int *prop_out;
	node *m1, *m2;
	int i,j;
	
	if (M1.Nel!=M2.Nel)
		Error("Cannot join meshes with unequal number of electrodes, M1:%d M2:%d\n",M1.Nel, M2.Nel);
	outl1=MeshOutline(M1);
	outl2=MeshOutline(M2);
	xoff1=M1.nodes[0].x1;
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(M1, outl1[i]);
		if (xoff1<m1->x2)
			xoff1=m1->x2;
	}	
	xoff2=M2.nodes[0].x2;
	for (i=1;i<=outl2[0];i++)
	{
		m1=SearchNode(M2, outl2[i]);
		if (xoff2>m1->x1)
			xoff2=m1->x1;
	}	
	xoff=xoff1-xoff2;
	
	/* TODO check for overlap */
	
	res=DuplicateMesh(M1);
	res.nodes=realloc(res.nodes,(res.Nn+M2.Nn+1)*sizeof(node));
	res.Nn+=M2.Nn;
	
	/* merge local properties, if local properties have the same name we assume they are the same! */
	prop_out=malloc((M2.Na+1)*sizeof(int));
	for (i=0;i<M2.Na;i++)
	{
		if ((j=FindProperties(res, M2.P[i].name))<0)
		{
			res.Na++;
			res.P=realloc(res.P, (res.Na+1)*sizeof(local_prop));
			DuplicateProperties(&res, res.P+res.Na-1, M2.P+i);
			prop_out[i]=res.Na-1;		
		}
		else
			prop_out[i]=j;	
		
	}
	for (i=0;i<M2.Nn;i++)
	{
		DuplicateNode(M2,res.nodes+i+M1.Nn, M2.nodes[i].id);
		res.nodes[i+M1.Nn].id+=M1.Nn;
		res.nodes[i+M1.Nn].x1+=xoff;
		res.nodes[i+M1.Nn].x2+=xoff;
		res.nodes[i+M1.Nn].y1+=yoff;
		res.nodes[i+M1.Nn].y2+=yoff;
		res.nodes[i+M1.Nn].P=prop_out[res.nodes[i+M1.Nn].P];
		for (j=1;j<=res.nodes[i+M1.Nn].north[0];j++)
			res.nodes[i+M1.Nn].north[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].south[0];j++)
			res.nodes[i+M1.Nn].south[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].east[0];j++)
			res.nodes[i+M1.Nn].east[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].west[0];j++)
			res.nodes[i+M1.Nn].west[j]+=M1.Nn;
	}
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(res, outl1[i]);
		for (j=1;j<=outl2[0];j++)
		{
			m2=SearchNode(res, outl2[j]+M1.Nn);
			NodeConnector(m1, m2);
		}			
	}
	free(outl1);
	free(outl2);
	free(prop_out);
	ReInitResults(&res);

	return res;
}
mesh JoinMeshes_V(mesh M1, mesh M2, double xoff)
{
	mesh res;
	double yoff1, yoff2, yoff;
	int *outl1, *outl2;
	int *prop_out;
	node *m1, *m2;
	int i,j;
	
	if (M1.Nel!=M2.Nel)
		Error("Cannot join meshes with unequal number of electrodes, M1:%d M2:%d\n",M1.Nel, M2.Nel);
	outl1=MeshOutline(M1);
	outl2=MeshOutline(M2);
	yoff1=M1.nodes[0].y1;
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(M1, outl1[i]);
		if (yoff1<m1->y2)
			yoff1=m1->y2;
	}	
	yoff2=M2.nodes[0].y2;
	for (i=1;i<=outl2[0];i++)
	{
		m1=SearchNode(M2, outl2[i]);
		if (yoff2>m1->y1)
			yoff2=m1->y1;
	}	
	yoff=yoff1-yoff2;
	
	/* TODO check for overlap */
	
	res=DuplicateMesh(M1);
	res.nodes=realloc(res.nodes,(res.Nn+M2.Nn+1)*sizeof(node));
	res.Nn+=M2.Nn;
	
	/* merge local properties, if local properties have the same name we assume they are the same! */
	prop_out=malloc((M2.Na+1)*sizeof(int));
	for (i=0;i<M2.Na;i++)
	{
		if ((j=FindProperties(res, M2.P[i].name))<0)
		{
			res.Na++;
			res.P=realloc(res.P, (res.Na+1)*sizeof(local_prop));
			DuplicateProperties(&res, res.P+res.Na-1, M2.P+i);
			prop_out[i]=res.Na-1;		
		}
		else
			prop_out[i]=j;	
		
	}
	for (i=0;i<M2.Nn;i++)
	{
		DuplicateNode(M2,res.nodes+i+M1.Nn, M2.nodes[i].id);
		res.nodes[i+M1.Nn].id+=M1.Nn;
		res.nodes[i+M1.Nn].x1+=xoff;
		res.nodes[i+M1.Nn].x2+=xoff;
		res.nodes[i+M1.Nn].y1+=yoff;
		res.nodes[i+M1.Nn].y2+=yoff;
		res.nodes[i+M1.Nn].P=prop_out[res.nodes[i+M1.Nn].P];
		for (j=1;j<=res.nodes[i+M1.Nn].north[0];j++)
			res.nodes[i+M1.Nn].north[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].south[0];j++)
			res.nodes[i+M1.Nn].south[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].east[0];j++)
			res.nodes[i+M1.Nn].east[j]+=M1.Nn;
		for (j=1;j<=res.nodes[i+M1.Nn].west[0];j++)
			res.nodes[i+M1.Nn].west[j]+=M1.Nn;
	}
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(res, outl1[i]);
		for (j=1;j<=outl2[0];j++)
		{
			m2=SearchNode(res, outl2[j]+M1.Nn);
			NodeConnector(m1, m2);
		}			
	}
	free(outl1);
	free(outl2);
	free(prop_out);
	ReInitResults(&res);

	return res;
}

void SplitNodeX(int id, mesh *M)
{
	int newid, i, j, *list;
	node *old, *new, *a;
	/* add a node */
	newid=M->Nn;
	M->Nn++;
	M->nodes=realloc(M->nodes, (M->Nn+1)*sizeof(node));
	
	/* reallocate and copy node voltages */
	for (i=0;i<M->res.Nva;i++)
	{
		for (j=0;j<M->Nel;j++)
		{
			M->res.Vn[i][j]=realloc(M->res.Vn[i][j],(M->Nn+1)*sizeof(double));			
			M->res.Vn[i][j][M->Nn-1]=M->res.Vn[i][j][id];
		}
	}
	
	old=SearchNode(*M, id);
	new=M->nodes+newid;
	
	
	/* copy data */ 
	DuplicateNode(*M,new, id);
	/* set x coordinates of the two nodes */ 
	new->id=newid;
	new->x1=old->x1;
	new->x2=(old->x1+old->x2)/2;
	old->x1=(old->x1+old->x2)/2;
	
	old->west=realloc(old->west, LISTBLOCK*sizeof(int));
	old->west[1]=newid;
	old->west[0]=1;
	new->east=realloc(new->east, LISTBLOCK*sizeof(int));
	new->east[1]=id;
	new->east[0]=1;
	for (i=1;i<=new->west[0];i++)
	{
		a=SearchNode(*M, new->west[i]);
		a->east=AddToList(a->east, newid);
		a->east=RemoveFromList(a->east, id);
	}
	/* now search through the north and south nodes to see which are adjacent to the new nodes */
	list=DuplicateList(old->north);
	for (i=1;i<=list[0];i++)
	{
		a=SearchNode(*M, list[i]);
		if (R_Overlap(a->x1,new->x1, a->x2, new->x2) > 0)
			a->south=AddToList(a->south, newid);
		else
			new->north=RemoveFromList(new->north, a->id);
		if (R_Overlap(a->x1,old->x1, a->x2, old->x2) < 0)
		{
			old->north=RemoveFromList(old->north, a->id);
			a->south=RemoveFromList(a->south, id);
		}
	}
	free(list);
	list=DuplicateList(old->south);
	for (i=1;i<=list[0];i++)
	{
		a=SearchNode(*M, list[i]);
		if (R_Overlap(a->x1,new->x1, a->x2, new->x2) > 0)
			a->north=AddToList(a->north, newid);
		else
			new->south=RemoveFromList(new->south, a->id);
		if (R_Overlap(a->x1,old->x1, a->x2, old->x2) < 0)
		{
			old->south=RemoveFromList(old->south, a->id);
			a->north=RemoveFromList(a->north, id);
		}
	}
	free(list);
}
	
void SplitNodeY(int id, mesh *M)
{
	int newid, i, j, *list;
	node *old, *new, *a;
	/* add a node */
	newid=M->Nn;
	M->Nn++;
	M->nodes=realloc(M->nodes, (M->Nn+1)*sizeof(node));
	
	/* reallocate and copy node voltages */
	for (i=0;i<M->res.Nva;i++)
	{
		for (j=0;j<M->Nel;j++)
		{
			M->res.Vn[i][j]=realloc(M->res.Vn[i][j],(M->Nn+1)*sizeof(double));			
			M->res.Vn[i][j][M->Nn-1]=M->res.Vn[i][j][id];
		}
	}
	
	old=SearchNode(*M, id);
	new=M->nodes+newid;
	
	
	/* copy data */ 
	DuplicateNode(*M,new, id);
	/* set y coordinates of the two nodes */ 
	new->id=newid;
	new->y1=old->y1;
	new->y2=(old->y1+old->y2)/2;
	old->y1=(old->y1+old->y2)/2;
	
	old->south=realloc(old->south, LISTBLOCK*sizeof(int));
	old->south[1]=newid;
	old->south[0]=1;
	new->north=realloc(new->north, LISTBLOCK*sizeof(int));
	new->north[1]=id;
	new->north[0]=1;
	
	for (i=1;i<=new->south[0];i++)
	{
		a=SearchNode(*M, new->south[i]);
		a->north=AddToList(a->north, newid);
		a->north=RemoveFromList(a->north, id);
	}
	
	
	/* now search through the west and east nodes to see which are adjacent to two new nodes */
	list=DuplicateList(old->west);
	for (i=1;i<=list[0];i++)
	{
		a=SearchNode(*M, list[i]);
		if (R_Overlap(a->y1,new->y1, a->y2, new->y2) > 0)
			a->east=AddToList(a->east, newid);
		else
			new->west=RemoveFromList(new->west, a->id);
		if (R_Overlap(a->y1,old->y1, a->y2, old->y2) < 0)
		{
			old->west=RemoveFromList(old->west, a->id);
			a->east=RemoveFromList(a->east, id);
		}
	}
	free(list);
	list=DuplicateList(old->east);
	for (i=1;i<=list[0];i++)
	{
		a=SearchNode(*M, list[i]);
		if (R_Overlap(a->y1,new->y1, a->y2, new->y2) > 0)
			a->west=AddToList(a->west, newid);
		else
			new->east=RemoveFromList(new->east, a->id);
		if (R_Overlap(a->y1,old->y1, a->y2, old->y2) < 0)
		{
			old->east=RemoveFromList(old->east, a->id);
			a->west=RemoveFromList(a->west, id);
		}
	}
	free(list);
}	

void SplitNodeXY(int id, mesh *M)
{
	SplitNodeX(id, M);
	SplitNodeY(M->Nn-1, M);
	SplitNodeY(id, M);
}

void SplitMeshX(mesh *M)
{
	int i;
	for (i=M->Nn;i>0;i--)
		SplitNodeX(i-1, M);
}

void SplitMeshY(mesh *M)
{
	int i;
	for (i=M->Nn;i>0;i--)
		SplitNodeY(i-1, M);
}

void SplitMeshLong(mesh *M)
{
	int i;
	for (i=M->Nn;i>0;i--)
	{
		if ((M->nodes[i-1].x2-M->nodes[i-1].x1)>(M->nodes[i-1].y2-M->nodes[i-1].y1))
			SplitNodeX(M->nodes[i-1].id, M);
		else
			SplitNodeY(M->nodes[i-1].id, M);
	}
}

void SplitMeshXY(mesh *M)
{
	int i;
	for (i=M->Nn;i>0;i--)
		SplitNodeXY(i-1, M);
}

void SplitMeshWhileCoarse(mesh *M, double d)
{
	int i, *list, *newlist, c=0;	
	
	list=malloc(((M->Nn+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	for (i=0;i<M->Nn;i++)
		list[i+1]=i;
	list[0]=M->Nn;
	
	newlist=malloc(LISTBLOCK*sizeof(int));
	newlist[0]=0;
	
	do
	{
		Print(DEBUG, "%i elements to check %i",list[0], c);
		c++;
		for (i=list[0];i>0;i--) 
		{
			node *N;
			N=SearchNode(*M, list[i]);
			if (((N->x2-N->x1)>d)&&((N->y2-N->y1)>d))
			{
				SplitNodeXY(list[i], M);
				newlist=AddToList(newlist, list[i]);
				newlist=AddToList(newlist, M->Nn-1);
				newlist=AddToList(newlist, M->Nn-2);
				newlist=AddToList(newlist, M->Nn-3);
			}
			else if ((N->x2-N->x1)>d)
			{
				SplitNodeX(list[i], M);
				newlist=AddToList(newlist, list[i]);
				newlist=AddToList(newlist, M->Nn-1);
			}
			else if ((N->y2-N->y1)>d)
			{
				SplitNodeY(list[i], M);
				newlist=AddToList(newlist, list[i]);
				newlist=AddToList(newlist, M->Nn-1);
			}
		}
		free(list);
		list=newlist;		
		newlist=malloc(LISTBLOCK*sizeof(int));
		newlist[0]=0;
		
	}while (list[0]);
	free(list);
	free(newlist);
}

void SplitListX(mesh *M, int *list)
{
	int i;
	for (i=1;i<=list[0];i++)
		SplitNodeX(list[i], M);
}
void SplitListY(mesh *M, int *list)
{
	int i;
	for (i=1;i<=list[0];i++)
		SplitNodeY(list[i], M);
}


void SplitListXY(mesh *M, int *list)
{
	int i;
	for (i=1;i<=list[0];i++)
		SplitNodeXY(list[i], M);
}
void SplitListLong(mesh *M, int *list)
{
	int i;
	for (i=1;i<=list[0];i++)
	{
		node *N;
		N=SearchNode(*M, list[i]);
		if ((N->x2-N->x1)>(N->y2-N->y1))
			SplitNodeX(N->id, M);
		else
			SplitNodeY(N->id, M);
	}
}

void SplitListWhileCoarse(mesh *M, int *list, double d)
{
	int i;
	int *list_c, *newlist;
	
	list_c=DuplicateList(list);
	
	newlist=malloc(LISTBLOCK*sizeof(int));
	newlist[0]=0;
	do
	{
		Print(DEBUG, "%i elements to check",list_c[0]);
		for (i=list_c[0];i>0;i--) 
		/* Careful: This loop walks backward through the list. The reason is that any newly created node will end up at the 
		  back or the list (as the lists are sorted and newly created lists have indices at the end of the node list). By walking
		  backward we prevent the loop from interfering with itself as list_c[0] changes during the execuation of the loop.*/
		{
			node *N;
			N=SearchNode(*M, list_c[i]);
			if (((N->x2-N->x1)>d)&&((N->y2-N->y1)>d))
			{
				SplitNodeXY(list_c[i], M);
				newlist=AddToList(newlist, list_c[i]);
				newlist=AddToList(newlist, M->Nn-1);
				newlist=AddToList(newlist, M->Nn-2);
				newlist=AddToList(newlist, M->Nn-3);
			}
			else if ((N->x2-N->x1)>d)
			{
				SplitNodeX(list_c[i], M);
				newlist=AddToList(newlist, list_c[i]);
				newlist=AddToList(newlist, M->Nn-1);
			}
			else if ((N->y2-N->y1)>d)
			{
				SplitNodeY(list_c[i], M);
				newlist=AddToList(newlist, list_c[i]);
				newlist=AddToList(newlist, M->Nn-1);
			}
		}
		free(list_c);
		list_c=newlist;		
		newlist=malloc(LISTBLOCK*sizeof(int));
		newlist[0]=0;
	} while (list_c[0]);
	free(list_c);
	free(newlist);
}

/* apply basic geometrical transforms on meshes */
void ScaleMeshX(mesh *M, double f)
{
	int i;
	int *dummy;
	double xx;
	for (i=0;i<M->Nn;i++)
	{
		M->nodes[i].x1*=f;
		M->nodes[i].x2*=f;
		if (f<0) /* swap east and west */
		{
			/* adapt nodes to new orientation */ 
			/* swap east and west */
			dummy=M->nodes[i].east;
			M->nodes[i].east=M->nodes[i].west;
			M->nodes[i].west=dummy;
			xx=M->nodes[i].x1;
			M->nodes[i].x1=M->nodes[i].x2;
			M->nodes[i].x2=xx;
		}
	}
}
void ScaleMeshY(mesh *M, double f)
{
	int i;
	int *dummy;
	double yy;
	for (i=0;i<M->Nn;i++)
	{
		M->nodes[i].y1*=f;
		M->nodes[i].y2*=f;
		if (f<0) /* swap north and south */
		{
			/* adapt nodes to new orientation */ 
			/* swap north and south */
			dummy=M->nodes[i].north;
			M->nodes[i].north=M->nodes[i].south;
			M->nodes[i].south=dummy;
			yy=M->nodes[i].y1;
			M->nodes[i].y1=M->nodes[i].y2;
			M->nodes[i].y2=yy;
		
		}
	}
}
void ScaleMesh(mesh *M, double f)
{
	int i;
	int *dummy;
	double xx,yy;
	for (i=0;i<M->Nn;i++)
	{
		M->nodes[i].x1*=f;
		M->nodes[i].x2*=f;
		M->nodes[i].y1*=f;
		M->nodes[i].y2*=f;
		if (f<0)
		{
			/* adapt nodes to new orientation */ 
			/* swap east and west, and north and south */
			dummy=M->nodes[i].east;
			M->nodes[i].east=M->nodes[i].west;
			M->nodes[i].west=dummy;
		
			dummy=M->nodes[i].north;
			M->nodes[i].north=M->nodes[i].south;
			M->nodes[i].south=dummy;
			
			xx=M->nodes[i].x1;
			M->nodes[i].x1=M->nodes[i].x2;
			M->nodes[i].x2=xx;
			
			yy=M->nodes[i].y1;
			M->nodes[i].y1=M->nodes[i].y2;
			M->nodes[i].y2=yy;
		}
	}
}
void MoveMesh(mesh *M, double x, double y)
{
	int i;
	for (i=0;i<M->Nn;i++)
	{
		M->nodes[i].x1+=x;
		M->nodes[i].x2+=x;
		M->nodes[i].y1+=y;
		M->nodes[i].y2+=y;
	}
}

void RotateMesh(mesh *M, double x, double y, int d)
{
	int i;
	double c, s, xx, yy;
	int *dummy;
	
	/* rotate in multiples of 90 degrees */
	d/=90;
	while (d<0)
		d+=4;
	d%=4;
	switch (d)
	{
		case 1: /* rotate 90 degrees: cos(d)=0, sin(d)=1 */
			c=0.0;
			s=1.0;
			break;
		case 2: /* rotate 180 degrees: cos(d)=-1, sin(d)=0 */
			c=-1.0;
			s=0.0;
			break;
		case 3: /* rotate 270 degrees: cos(d)=0, sin(d)=-1 */
			c=0.0;
			s=-1.0;
			break;
		default: /* no rotation (d=0) */
			return;
	}
		
	for (i=0;i<M->Nn;i++)
	{
		/* move & rotate nodes */
		xx=(M->nodes[i].x1-x);
		yy=(M->nodes[i].y1-y);
		M->nodes[i].x1=x + xx*c - yy*s;
		M->nodes[i].y1=y + yy*c + xx*s;
		
		xx=(M->nodes[i].x2-x);
		yy=(M->nodes[i].y2-y);
		M->nodes[i].x2=x + xx*c - yy*s;
		M->nodes[i].y2=y + yy*c + xx*s;
		/* adapt nodes to new orientation */
		switch (d)
		{
			case 1: /* E->N, S->E, W->S, N->W */
				dummy=M->nodes[i].north;
				M->nodes[i].north=M->nodes[i].east;
				M->nodes[i].east=M->nodes[i].south;
				M->nodes[i].south=M->nodes[i].west;
				M->nodes[i].west=dummy;
				/* make sure the lower left corner and the upper right corner are exactly that */
				xx=M->nodes[i].x1;
				M->nodes[i].x1=M->nodes[i].x2;
				M->nodes[i].x2=xx;				
				break;
			case 2: /* E->W, S->N, W->E, N->S */
				dummy=M->nodes[i].north;
				M->nodes[i].north=M->nodes[i].south;
				M->nodes[i].south=dummy;				
				dummy=M->nodes[i].east;
				M->nodes[i].east=M->nodes[i].west;
				M->nodes[i].west=dummy;
				/* make sure the lower left corner and the upper right corner are exactly that */
				xx=M->nodes[i].x1;
				M->nodes[i].x1=M->nodes[i].x2;
				M->nodes[i].x2=xx;	
				yy=M->nodes[i].y1;
				M->nodes[i].y1=M->nodes[i].y2;
				M->nodes[i].y2=xx;	
				break;
			case 3:  /* E->S, S->W, W->N, N->E */
				dummy=M->nodes[i].south;
				M->nodes[i].south=M->nodes[i].east;
				M->nodes[i].east=M->nodes[i].north;
				M->nodes[i].north=M->nodes[i].west;
				M->nodes[i].west=dummy;	
				/* make sure the lower left corner and the upper right corner are exactly that */
				yy=M->nodes[i].y1;
				M->nodes[i].y1=M->nodes[i].y2;
				M->nodes[i].y2=xx;	
				break;
		}
		
	}
}


void GetMeshBB(mesh *M, double *x1, double *y1, double *x2, double *y2)
{
	/* compute bounding box of a mesh */
	int i;
	
	(*x1)=M->nodes[0].x1;
	(*y1)=M->nodes[0].y1;
	(*x2)=M->nodes[0].x2;
	(*y2)=M->nodes[0].y2;
	for (i=1;i<M->Nn;i++)
	{
		if (M->nodes[i].x1<(*x1))
			(*x1)=M->nodes[i].x1;
		if (M->nodes[i].y1<(*y1))
			(*y1)=M->nodes[i].y1;
		if (M->nodes[i].x2>(*x2))
			(*x2)=M->nodes[i].x2;
		if (M->nodes[i].y2>(*y2))
			(*y2)=M->nodes[i].y2;
	}
}

void SetMeshBB(mesh *M, double x1, double y1, double x2, double y2, int FixR)
{
	/* change bounding box of a mesh (using scaling and moving)  */
	double xx1, xx2, yy1, yy2;
	double fx,fy;
	GetMeshBB(M, &xx1, &yy1, &xx2, &yy2);
	
	fx=(x2-x1)/(xx2-xx1);
	fy=(y2-y1)/(yy2-yy1);
	
	if (FixR)
	{
		if (fx<fy)
			fy=fx;
		else
			fx=fy;
		
		ScaleMesh(M, fx);
	
	}
	else
	{	
		ScaleMeshX(M, fx);
		ScaleMeshY(M, fy);
	}
	
	xx1=(x1+x2)/2-fx*(xx1+xx2)/2;
	yy1=(y1+y2)/2-fy*(yy1+yy2)/2;
	
	MoveMesh(M, xx1, yy1);
}


void SortMesh(mesh *M)
{
	/* sort mesh after removing nodes */
	int i, j;
	node *N, *Nl;
			
	for (i=0;i<M->Nn;i++)
	{
		N=M->nodes+i;
		if (N->id!=i)
		{
			for (j=1;j<=N->north[0];j++)
			{
				Nl=SearchNode(*M,N->north[j]);
				Nl->south=RemoveFromList(Nl->south, N->id);
				Nl->south=AddToList(Nl->south, i);
			}
			for (j=1;j<=N->south[0];j++)	
			{
				Nl=SearchNode(*M,N->south[j]);
				Nl->north=RemoveFromList(Nl->north, N->id);
				Nl->north=AddToList(Nl->north, i);
			}
			for (j=1;j<=N->west[0];j++)
			{
				Nl=SearchNode(*M,N->west[j]);
				Nl->east=RemoveFromList(Nl->east, N->id);
				Nl->east=AddToList(Nl->east, i);
			}
			for (j=1;j<=N->east[0];j++)
			{
				Nl=SearchNode(*M,N->east[j]);
				Nl->west=RemoveFromList(Nl->west, N->id);
				Nl->west=AddToList(Nl->west, i);
			}
			N->id=i;			
		}			
	}	
}

/* valgrind trips over this(?) routine. I've been checking everything and I cannot find anything wrong with it.
   I tried mtrace which finds no leaks. Perhaps there is a bug, perhaps valgrind just is confused by this as it
   is kind of hard to keep track of where what is allocated and freed. 
   If you want to use valgrind you cannot simplify the mesh, it will segfault.
   */
void CleanUpMesh(mesh *M, int *merged)
{
	int i, j;
	for (i=1;i<=merged[0];i++)
	{
		node *N, *L;
		N=SearchNode(*M,merged[i]);
		for (j=1;j<=N->north[0];j++)
		{
			L=SearchNode(*M,N->north[j]);
			L->south=RemoveFromList(L->south, N->id);
		}
		for (j=1;j<=N->south[0];j++)
		{
			L=SearchNode(*M,N->south[j]);
			L->north=RemoveFromList(L->north, N->id);
		}
		for (j=1;j<=N->east[0];j++)
		{
			L=SearchNode(*M,N->east[j]);
			L->west=RemoveFromList(L->west, N->id);
		} 
		for (j=1;j<=N->west[0];j++)
		{
			L=SearchNode(*M,N->west[j]);
			L->east=RemoveFromList(L->east, N->id);
		}
		
		free(N->north);
		free(N->south);
		free(N->east);
		free(N->west);
	}
	
	j=0;
	i=0;
	while (i<M->Nn)
	{
		if(!IsInList(merged, M->nodes[i].id))
		{
			if (j!=i)
				M->nodes[j]=M->nodes[i];
			j++;
		}
		i++;	
	}
	M->Nn=j;
	M->nodes=realloc(M->nodes, (M->Nn+1)*sizeof(node));
	/* i=1;
 	while (i<=merged[0])
	{
 		SearchLink(*M, merged[i]);
		i++;
	} */
 	SortMesh(M);
	ReInitResults(M);

}

#define MaxR 2.0
/************************************************/

int *Chunkify_node(mesh *M, int id, int * merged)
{
	int j;
	int *list;
	int MRG;
	double R;	
	double newx1, newy1, newx2, newy2;
	node *N, *Nl;
	N=SearchNode(*M,id);	
	/* check east nodes to merge with */
	j=1;
	MRG=1;	
	while ((j<=N->east[0])&&(MRG==1))
	{
		if (IsInList(merged, N->east[j]))
		{
			MRG=0;
		}
		else
		{
			Nl=SearchNode(*M,N->east[j]);
			if (j==1)
				newx2=Nl->x2;
			newx2=MIN(newx2,Nl->x2);
				
			if (Nl->P!=N->P)
				MRG=0;											
			else if ((Nl->y2>N->y2+TINY)||(Nl->y1<N->y1-TINY))
				MRG=0;	
		}
		j++;
	}
	
	R=(newx2-N->x1)/(N->y2-N->y1);
	if (R>MaxR)
		MRG=0;
		
				
	if ((MRG)&&(N->east[0]>0))
	{
		int *a_nodes;
		a_nodes=malloc(LISTBLOCK*sizeof(int));
		a_nodes[0]=0;
		
		j=1;
		Nl=SearchNode(*M,N->east[j]);
		merged=AddToList(merged,Nl->id);	
		
		a_nodes=AddListToList(a_nodes, Nl->north);
		a_nodes=AddListToList(a_nodes, Nl->south);
		a_nodes=AddListToList(a_nodes, Nl->east);
		
		j++;
		while (j<=N->east[0])
		{
			Nl=SearchNode(*M,N->east[j]);
			merged=AddToList(merged,Nl->id);	
		
			a_nodes=AddListToList(a_nodes, Nl->north);
			a_nodes=AddListToList(a_nodes, Nl->south);
			a_nodes=AddListToList(a_nodes, Nl->east);
			
			j++;
		}
		N->x2=newx2;
		
		j=1;
		while (j<=N->east[0])
		{
			Nl=SearchNode(*M,N->east[j]);
			if (newx2<Nl->x2-TINY)
			{
				int k=1;
				node *L;
				merged=RemoveFromList(merged,Nl->id);
				Nl->x1=newx2;
				/* filter out the north and south nodes */
				list=DuplicateList(Nl->north);
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
				list=DuplicateList(Nl->south);
				k=1;
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
				
			}	
			j++;
		}
		j=1;
		while (j<=a_nodes[0])
		{
			Nl=SearchNode(*M,a_nodes[j]);
			NodeConnector(Nl, N);
			j++;
		}
		free(a_nodes);						
	}
	
	/*north*/	
	j=1;
	MRG=1;
	while ((j<=N->north[0])&&(MRG==1))
	{
		if (IsInList(merged, N->north[j]))
		{
			MRG=0;
		}
		else
		{
			Nl=SearchNode(*M,N->north[j]);
			if (j==1)
				newy2=Nl->y2;
			newy2=MIN(newy2,Nl->y2);
				
			if (Nl->P!=N->P)
				MRG=0;										
			else if ((Nl->x2>N->x2+TINY)||(Nl->x1<N->x1-TINY))
				MRG=0;	
		}
		j++;
	}
	R=(N->x2-N->x1)/(newy2-N->y1);
	if (R<1.0/MaxR)
		MRG=0;		
	if ((MRG)&&(N->north[0]>0))
	{
		int *a_nodes;
		a_nodes=malloc(LISTBLOCK*sizeof(int));
		a_nodes[0]=0;
		
		j=1;
		Nl=SearchNode(*M,N->north[j]);
		merged=AddToList(merged,Nl->id);	
		
		a_nodes=AddListToList(a_nodes, Nl->north);
		a_nodes=AddListToList(a_nodes, Nl->west);
		a_nodes=AddListToList(a_nodes, Nl->east);
		
		j++;
		while (j<=N->north[0])
		{
			Nl=SearchNode(*M,N->north[j]);
			merged=AddToList(merged,Nl->id);	
		
			a_nodes=AddListToList(a_nodes, Nl->north);
			a_nodes=AddListToList(a_nodes, Nl->west);
			a_nodes=AddListToList(a_nodes, Nl->east);
			
			j++;
		}
		N->y2=newy2;
		
		j=1;
		while (j<=N->north[0])
		{
			Nl=SearchNode(*M,N->north[j]);
			if (newy2<Nl->y2-TINY)
			{
				int k=1;
				node *L;
				merged=RemoveFromList(merged,Nl->id);
				Nl->y1=newy2;
				/* filter out the west and east nodes */
				list=DuplicateList(Nl->west);
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
				list=DuplicateList(Nl->east);
				k=1;
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
			}	
			j++;
		}
		j=1;
		while (j<=a_nodes[0])
		{
			Nl=SearchNode(*M,a_nodes[j]);
			NodeConnector(Nl, N);
			j++;
		}
		free(a_nodes);
	}
	
	/* west */	
	j=1;
	MRG=1;
	while ((j<=N->west[0])&&(MRG==1))
	{
		if (IsInList(merged, N->west[j]))
		{
			MRG=0;
		}
		else
		{
			Nl=SearchNode(*M,N->west[j]);
			if (j==1)
				newx1=Nl->x1;
			newx1=MAX(newx1,Nl->x1);
				
			if (Nl->P!=N->P)
				MRG=0;										
			else if ((Nl->y2>N->y2+TINY)||(Nl->y1<N->y1-TINY))
				MRG=0;
		}	
		j++;
	}
	
	R=(N->x2-newx1)/(N->y2-N->y1);
	if (R>MaxR)
		MRG=0;	
		
	if ((MRG)&&(N->west[0]>0))
	{
		int *a_nodes;
		a_nodes=malloc(LISTBLOCK*sizeof(int));
		a_nodes[0]=0;
		
		j=1;
		Nl=SearchNode(*M,N->west[j]);
		merged=AddToList(merged,Nl->id);	
		
		a_nodes=AddListToList(a_nodes, Nl->north);
		a_nodes=AddListToList(a_nodes, Nl->south);
		a_nodes=AddListToList(a_nodes, Nl->west);
		
		j++;
		while (j<=N->west[0])
		{
			Nl=SearchNode(*M,N->west[j]);
			merged=AddToList(merged,Nl->id);	
		
			a_nodes=AddListToList(a_nodes, Nl->north);
			a_nodes=AddListToList(a_nodes, Nl->south);
			a_nodes=AddListToList(a_nodes, Nl->west);
			
			j++;
		}
		N->x1=newx1;
		
		j=1;
		while (j<=N->west[0])
		{
			Nl=SearchNode(*M,N->west[j]);
			if (newx1>Nl->x1+TINY)
			{
				int k=1;
				node *L;
				merged=RemoveFromList(merged,Nl->id);
				Nl->x2=newx1;
				/*filter out the north and south nodes*/
				list=DuplicateList(Nl->north);
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
				list=DuplicateList(Nl->south);
				k=1;
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
				
			}	
			j++;
		}
		j=1;
		while (j<=a_nodes[0])
		{
			Nl=SearchNode(*M,a_nodes[j]);
			NodeConnector(Nl, N);
			j++;
		}
		free(a_nodes);
						
	}	
	/* south */
	j=1;
	MRG=1;
	while ((j<=N->south[0])&&(MRG==1))
	{
		if (IsInList(merged, N->south[j]))
		{
			MRG=0;
		}
		else
		{
			Nl=SearchNode(*M,N->south[j]);
			
			if (j==1)
				newy1=Nl->y1;
			newy1=MAX(newy1,Nl->y1);
			
			if (Nl->P!=N->P)
				MRG=0;						
			else if ((Nl->x2>N->x2+TINY)||(Nl->x1<N->x1-TINY))
				MRG=0;	
		}
		j++;
	}
	R=(N->x2-N->x1)/(N->y2-newy1);
	if (R<1.0/MaxR)
		MRG=0;		
	if ((MRG)&&(N->south[0]>0))
	{
		int *a_nodes;
		a_nodes=malloc(LISTBLOCK*sizeof(int));
		a_nodes[0]=0;
		
		j=1;
		Nl=SearchNode(*M,N->south[j]);
		merged=AddToList(merged,Nl->id);	
		
		a_nodes=AddListToList(a_nodes, Nl->south);
		a_nodes=AddListToList(a_nodes, Nl->west);
		a_nodes=AddListToList(a_nodes, Nl->east);
		
		j++;
		while (j<=N->south[0])
		{
			Nl=SearchNode(*M,N->south[j]);
			merged=AddToList(merged,Nl->id);	
		
			a_nodes=AddListToList(a_nodes, Nl->south);
			a_nodes=AddListToList(a_nodes, Nl->west);
			a_nodes=AddListToList(a_nodes, Nl->east);
			
			j++;
		}
		N->y1=newy1;
		
		j=1;
		while (j<=N->south[0])
		{
			Nl=SearchNode(*M,N->south[j]);
			if (newy1>Nl->y1+TINY)
			{
				int k=1;
				node *L;
				merged=RemoveFromList(merged,Nl->id);
				Nl->y2=newy1;
				/*filter out the west and east nodes*/
				list=DuplicateList(Nl->west);
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
				list=DuplicateList(Nl->east);
				k=1;
				while (k<=list[0])
				{
					L=SearchNode(*M,list[k]);
					NodeConnector(Nl, L);
					k++;
				}
				free(list);
			}	
			j++;
		}
		j=1;
		while (j<=a_nodes[0])
		{
			Nl=SearchNode(*M,a_nodes[j]);
			NodeConnector(Nl, N);
			j++;
		}
		free(a_nodes);
	}
	return merged;
}
/*
int Chunkify_nodes_(mesh *M, int skip, int offset)
{
	int i,Nold;
	int *merged;
	merged=malloc(LISTBLOCK*sizeof(int));
	merged[0]=0;
	Nold=M->Nn;
	i=0;
	
	while (i<Nold)
	{
		merged=Chunkify_node(M, i, merged);
		i+=skip;
		while (IsInList(merged, i)&&(i<M->Nn))
			i++;
	}
	CleanUpMesh(M, merged);
	free(merged);
	return Nold-M->Nn;
}
void Chunkify_(mesh *M, int J)
{
	int Nold, i=0, skip=1, offset=0;
	Nold=M->Nn+1; 
	
	while ((J>0)&&(skip<M->Nn))
	{
		skip*=2;
		J--;
	}
	J=skip;
	while ((M->Nn<Nold)&&(i<100))
	{
		Nold=M->Nn;
		skip=J;
		offset=J;
		Print(DEBUG,"Skip:%i Offset %i",skip,offset);
		Chunkify_nodes_(M, skip, 0);
		offset/=2;
		while (skip>1)
		{
			Print(DEBUG,"Skip:%i Offset %i",skip,offset);
			Chunkify_nodes_(M, skip, offset);
			offset/=2;
			skip/=2;
		}
		i++;
		Print(DEBUG,"Round %i, %i nodes left",i, M->Nn);
	}
}
*/
/* this time we add some true randomness to our mesh simplifier */
int Random(int rmin, int rmax)
{
	return rmin + (int) (1.0*(rmax-rmin+1) * rand()/(RAND_MAX+1.0) );
}


void InitRandom(void)
{
	time_t curtime;
	
	time( &curtime );
	srand( (unsigned int) curtime );
	Random(0, 1000);
}


void Shuffle (int *list)
/* shuffle integer array */
{
	int i;
	
	for( i = 2; i <list[0]; i++ ){
		int	where;
		int	temp;

		where = Random(0,32767) % (i-1) + 1;		
		temp = list[where];
		list[where] = list[i+1];
		list[i+1] = temp;
	}
}

int Chunkify_nodes(mesh *M, int *list, int skip, int offset)
{
	int i,Nold;
	int *merged;
	merged=malloc(LISTBLOCK*sizeof(int));
	merged[0]=0;
	Nold=M->Nn;
	i=offset;
	
	while (i<Nold)
	{
		merged=Chunkify_node(M, list[i+1], merged);
		i+=skip;
		while (IsInList(merged, list[i+1])&&(i<M->Nn))
			i++;
	}
	/* cleanup mesh, i.e. remove the merged nodes and sort the node id's */
	CleanUpMesh(M, merged);
	free(merged);
	return Nold-M->Nn;
}

void Chunkify(mesh *M)
{
	int Nold, i=0, j,J,skip=1, offset=0;
	int *list;
	Nold=M->Nn+1; 
	list=malloc((M->Nn+2)*sizeof(int));
	J=(int)(2*sqrt((double)M->Nn)/300);
	if (J>M->Nn/2)
		J=M->Nn/2;
	while (skip<J)
	{
		skip*=2;
	}
	J=skip;
	InitRandom();
	while ((M->Nn<Nold)&&(i<100))
	{
		/* we add a random component here to make the algorithm a bit more robust against special cases where a certain regular pattern may make everything turn bad */
		/* we therefore create a list of nodes and shuffle it to process the list in random order. This is only used to "seed" the "chunks" */
		Nold=M->Nn;
		list[0]=M->Nn;
		for (j=0;j<M->Nn;j++)
			list[j+1]=j;
		Shuffle (list);
		
		skip=J;
		offset=J;
		Print(DEBUG,"Skip:%i Offset %i %i",skip,offset,M->Nn);
		Chunkify_nodes(M,  list, skip, 0);
		
		/* here we proceed in normal order */
		list[0]=M->Nn;
		for (j=0;j<M->Nn;j++)
			list[j+1]=j;
		offset/=2;
		while (skip>1)
		{
			Print(DEBUG,"Skip:%i Offset %i %i",skip,offset,M->Nn);
			list[0]=M->Nn;
			Chunkify_nodes(M, list, skip, offset);
			offset/=2;
			skip/=2;
		}
		i++;
		Print(DEBUG,"Round %i, %i nodes left",i, M->Nn);
	}
}

/***************************************************************************
  The PVMOS binary file format is stupid. It just dumps the data structures 
  to a file with faily little stucture. Actually it is a fairly amorphous
  mess. In addition the sizes of data fiels may depend on compiler and 
  system, i.e. file formats may not be compatible for the same PVMOS 
  version compiled with different compilers. 
  
  I did put a few things in place to have some very basic error 
  checking while reading. The write and read functions are organized
  along the mesh data structures. As such I give each block and sub-block
  a 4 char long label. Upon reading I check whether we are still OK by 
  checking this label. In addition I use a signature of the mesh data 
  strutures to verify compatibility between the file format and the PVMOS 
  reading it (e.g. if your PVMOS uses a 2 byte int it will not attampt to 
  read a file written with a 4 byte int).
***************************************************************************/

enum {FF_NODE=0, FF_NODEARRAY, FF_PROP, FF_ELCONN, FF_PROPARRAY, FF_RESULTS};
#define LL 2 /* label length. As a rule of thumb every character adds 1.4% to the file size */
const char *FF_LABELS [] = {
	"no",
	"na",
	"pr",
	"co",
	"pa",
	"re"
};

/*dump a node in a file */
void WriteNode(FILE *f, node N)
{	
	fwrite(FF_LABELS[FF_NODE], sizeof(char), LL, f);
	fwrite(N.north, sizeof(int), N.north[0]+1, f);
	fwrite(N.south, sizeof(int), N.south[0]+1, f);
	fwrite(N.west, sizeof(int), N.west[0]+1, f);
	fwrite(N.east, sizeof(int), N.east[0]+1, f);
	fwrite(&N.id, sizeof(int), 1, f);
	fwrite(&N.x1, sizeof(double), 1, f);
	fwrite(&N.y1, sizeof(double), 1, f);
	fwrite(&N.x2, sizeof(double), 1, f);
	fwrite(&N.y2, sizeof(double), 1, f);
	fwrite(&N.P, sizeof(int), 1, f);
}

/* read a node, FIXME: error checking, in particular, end of file */
node ReadNode(FILE *f)
{
	int a;
	node N;
	char buf[LL];
	int i, c=0;
	
	if (!fread(buf, sizeof(char), LL, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_NODE], c++);
	for (i=0;i<LL;i++)
		if (buf[i]!=FF_LABELS[FF_NODE][i])
			Error("File format error, unexpected data field CODE: %s%02d\n", FF_LABELS[FF_NODE], c);	
	c++;
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	N.north=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.north[0]=a;
	if (fread(N.north+1, sizeof(int), N.north[0], f)<N.north[0])
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	N.south=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.south[0]=a;
	if (fread(N.south+1, sizeof(int), N.south[0], f)<N.south[0])
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	N.west=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.west[0]=a;
	if (fread(N.west+1, sizeof(int), N.west[0], f)<N.west[0])
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	N.east=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.east[0]=a;
	if (fread(N.east+1, sizeof(int), N.east[0], f)<N.east[0])
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	if (!fread(&N.id, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	if (!fread(&N.x1, sizeof(double), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	if (!fread(&N.y1, sizeof(double), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	if (!fread(&N.x2, sizeof(double), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	if (!fread(&N.y2, sizeof(double), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	if (!fread(&N.P, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE:  %s%02d\n", FF_LABELS[FF_NODE], c++);
	return N;
} 
void WriteNodeArray(FILE *f, mesh *M)
{
	int i;
	fwrite(FF_LABELS[FF_NODEARRAY], sizeof(char), LL, f);
	fwrite(&(M->Nn), sizeof(int), 1, f);
	for (i=0;i<M->Nn;i++)
		WriteNode(f, M->nodes[i]);
}

void ReadNodeArray(FILE *f, mesh *M)
{
	int i;
	char buf[LL];
	int c=0;
	
	if (!fread(buf, sizeof(char), LL, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_NODEARRAY], c);
	c++;
	for (i=0;i<LL;i++)
		if (buf[i]!=FF_LABELS[FF_NODEARRAY][i])
			Error("File format error, unexpected data field CODE: %s%02d\n", FF_LABELS[FF_NODEARRAY], c);
	c++;
	if (!fread(&i, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_NODEARRAY], c);
	M->Nn=i;
	M->nodes=malloc((M->Nn+1)*sizeof(node));
	for (i=0;i<M->Nn;i++)
		M->nodes[i]=ReadNode(f);
}
void WriteElConn(FILE *f, ElConn conn)
{
	fwrite(FF_LABELS[FF_ELCONN], sizeof(char), LL, f);
	fwrite(&conn.model, sizeof(diode_model), 1, f);
	fwrite(&conn.ParSize, sizeof(size_t), 1, f);
	if (conn.ParSize)
		fwrite(conn.ParStruct, conn.ParSize, 1, f);
	fwrite(&conn.N, sizeof(int), 1, f);
	fwrite(conn.V, sizeof(double), conn.N, f);
	fwrite(conn.J, sizeof(double), conn.N, f);
}

ElConn ReadElConn(FILE *f)
{
	ElConn conn;
	
	char buf[LL];
	int i, c=0;	
	if (!fread(buf, sizeof(char), LL, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	c++;
	for (i=0;i<LL;i++)
		if (buf[i]!=FF_LABELS[FF_ELCONN][i])
			Error("File format error, unexpected data field CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	c++;
	if (!fread(&conn.model, sizeof(diode_model), 1, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	c++;
	if (!fread(&conn.ParSize, sizeof(size_t), 1, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	
	c++;
	if (conn.ParSize)
	{		
		conn.ParStruct=malloc(conn.ParSize);	
		if (!fread(conn.ParStruct, conn.ParSize, 1, f))
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	}
	else
		conn.ParStruct=NULL;
	c++;		
	if (!fread(&conn.N, sizeof(int), 1, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	c++;
	conn.V=malloc((conn.N+1)*sizeof(double));	
	if (fread(conn.V, sizeof(double), conn.N, f)<conn.N)
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	c++;
	conn.J=malloc((conn.N+1)*sizeof(double));	
	if (fread(conn.J, sizeof(double), conn.N, f)<conn.N)
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_ELCONN], c);
	return conn;		
}

void WritePropertiesArray(FILE *f, mesh *M)
{
	int i,j;
	int len;
	fwrite(FF_LABELS[FF_PROPARRAY], sizeof(char), LL, f);
	fwrite(&(M->Na), sizeof(int), 1, f);
	fwrite(&(M->Nel), sizeof(int), 1, f);
	for (i=0;i<M->Na;i++)
	{
	
		fwrite(FF_LABELS[FF_PROP], sizeof(char), LL, f);
		len=strlen(M->P[i].name)+1;
		fwrite(&len, sizeof(int), 1, f);
		fwrite(M->P[i].name, sizeof(char), len, f);
		
		fwrite((M->P[i].Rel), sizeof(double), M->Nel, f);
		fwrite((M->P[i].Rvp), sizeof(double), M->Nel, f);
		fwrite((M->P[i].Rvn), sizeof(double), M->Nel, f);
		
		for (j=0;j<M->Nel-1;j++)
			WriteElConn(f, M->P[i].conn[j]);
		
		fwrite(&(M->P[i].T), sizeof(double), 1, f);
		fwrite(&(M->P[i].SplitX), sizeof(int), 1, f);
		fwrite(&(M->P[i].SplitY), sizeof(int), 1, f);	
	
	}
}
void ReadPropertiesArray(FILE *f, mesh *M)
{
	int i,j;
	char buf[LL];
	int c=0, c1;	
	if (!fread(buf, sizeof(char), LL, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROPARRAY], c++);
	for (i=0;i<LL;i++)
		if (buf[i]!=FF_LABELS[FF_PROPARRAY][i])
			Error("File format error, unexpected data field CODE: %s%02d\n", FF_LABELS[FF_PROPARRAY], c);
	c++;
	if (!fread(&i, sizeof(int), 1, f))
		Error("Premature end of file\n");
	M->Na=i;
	if (!fread(&i, sizeof(int), 1, f))
		Error("Premature end of file\n");
	M->Nel=i;
	M->P=malloc((M->Na+1)*sizeof(local_prop));
	c1=c;
	for (i=0;i<M->Na;i++)
	{
		c=c1;
		if (!fread(buf, sizeof(char), LL, f))
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		for (j=0;j<LL;j++)
			if (buf[j]!=FF_LABELS[FF_PROP][j])
				Error("File format error, unexpected data field CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
				
		if (!fread(&j, sizeof(int), 1, f))
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		M->P[i].name=malloc(MAXNAMELEN*sizeof(char));
		if (fread(M->P[i].name,sizeof(char), j, f)<j)
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
			
		M->P[i].Rel=malloc((M->Nel+1)*sizeof(double));
		M->P[i].Rvp=malloc((M->Nel+1)*sizeof(double));
		M->P[i].Rvn=malloc((M->Nel+1)*sizeof(double));
		if (fread(M->P[i].Rel, sizeof(double), M->Nel, f)<M->Nel)
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		if (fread(M->P[i].Rvp, sizeof(double), M->Nel, f)<M->Nel)
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		if (fread(M->P[i].Rvn, sizeof(double), M->Nel, f)<M->Nel)
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		
		M->P[i].conn=malloc((M->Nel)*sizeof(ElConn));	
		for (j=0;j<M->Nel-1;j++)
			M->P[i].conn[j]=ReadElConn(f);
		
		if (!fread(&(M->P[i].T), sizeof(double), 1, f))
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		if (!fread(&(M->P[i].SplitX), sizeof(int), 1, f))
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
		if (!fread(&(M->P[i].SplitY), sizeof(int), 1, f))
			Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_PROP], c);
		c++;
	}
}


void WriteResults(FILE *f, mesh *M)
{
	int i,j;
	
	fwrite(FF_LABELS[FF_RESULTS], sizeof(char), LL, f);
	fwrite(&(M->res.Nva), sizeof(int), 1, f);
	
	fwrite(M->res.Va, sizeof(double), M->res.Nva, f);
	fwrite(M->res.I, sizeof(double), M->res.Nva, f);
	for (i=0;i<M->res.Nva;i++)
		for (j=0;j<M->Nel;j++)
			fwrite(M->res.Vn[i][j], sizeof(double), M->Nn, f);
		
}

void ReadResults(FILE *f, mesh *M)
{
	int i,j;

	char buf[LL];
	int c=0;	
	if (!fread(buf, sizeof(char), LL, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_RESULTS], c++);
	for (i=0;i<LL;i++)
		if (buf[i]!=FF_LABELS[FF_RESULTS][i])
			Error("File format error, unexpected data field CODE: %s%02d\n", FF_LABELS[FF_RESULTS], c);
	c++;
	if (!fread(&(M->res.Nva), sizeof(int), 1, f))
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_RESULTS], c);
	c++;
	M->res.Va=malloc((M->res.Nva+1)*sizeof(double));
	M->res.I=malloc((M->res.Nva+1)*sizeof(double));
	M->res.Vn=malloc((M->res.Nva+1)*sizeof(double **));
	
	if (fread(M->res.Va, sizeof(double), M->res.Nva, f)<M->res.Nva)
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_RESULTS], c);
	c++;
	if (fread(M->res.I, sizeof(double), M->res.Nva, f)<M->res.Nva)
		Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_RESULTS], c);
	c++;
	for (i=0;i<M->res.Nva;i++)
	{
		M->res.Vn[i]=malloc((M->Nel+1)*sizeof(double *));
		for (j=0;j<M->Nel;j++)
		{
			M->res.Vn[i][j]=malloc((M->Nn+1)*sizeof(double));
			if (fread(M->res.Vn[i][j], sizeof(double), M->Nn, f)<M->Nn)
				Error("Premature end of mesh file CODE: %s%02d\n", FF_LABELS[FF_RESULTS], c);
		}
	}
		
}

/* write a complete mesh to a file */
void WriteMesh(char *fn, mesh *M)
{
	FILE *f;
	int i;
	
	if ((f=fopen(fn,"wb"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<NMESHHASH;i++)
		fwrite(MESHHASH+i, sizeof(char), 1, f); /* hash for compatibility verification */
	WriteNodeArray(f, M);
	WritePropertiesArray(f, M);
	WriteResults(f, M);
	fclose(f);
}

/* read a complete mesh from a file */
void ReadMesh(char *fn, mesh *M)
{
	int i;
	unsigned char c;
	FILE *f;	
	
	if ((f=fopen(fn,"rb"))==NULL)
		Error("Cannot open %s for reading\n", fn);
	for (i=0;i<NMESHHASH;i++) /* hash for compatibility verification */
	{
		if (!fread(&c, sizeof(char), 1, f))
			Error("File shorter than %d bytes!\n", NMESHHASH);
		if (c!=MESHHASH[i])
			Error("Error: Incompatible mesh format in file %s.\n",fn);
	}
	ReadNodeArray(f, M);
	ReadPropertiesArray(f, M);
	ReadResults(f, M);
	fclose(f);
}
