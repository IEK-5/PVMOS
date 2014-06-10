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
 *    Dr. Bart E. Pieters 2014                                   *
 *                                                               *             
 *****************************************************************/                                                                             

/*****************************************************************           
 * FUNCTION:                                                     *                
 * routines to create and modify meshes                          *      
 *                                                               *            
 *****************************************************************/     
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "mesh2d.h"
#include "main.h"
#include "list.h"
#include "utils.h"
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))
#define TINY 1e-12
#define TWOPI 6.28318530717959
#define MAXNAMELEN 255

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
	M.res.Vn=malloc(sizeof(double *));
	
	/* local properties struct */
	M.P=malloc((2)*sizeof(local_prop));
	M.Na=1;	
	M.P[0].name=malloc(MAXNAMELEN*sizeof(char));
	if (strlen(name)+1<MAXNAMELEN)
		strncpy(M.P[0].name, name,strlen(name)+1);
	else
	{
		strncpy(M.P[0].name, name,MAXNAMELEN-1);
		M.P[0].name[MAXNAMELEN]='\0';
	}
	M.P[0].Rp=1;
	M.P[0].Rn=1;
	M.P[0].Rpvp=-1.0;
	M.P[0].Rpvn=-1.0;
	M.P[0].Rnvp=-1.0;
	M.P[0].Rnvn=-1.0;
	
	M.P[0].model=JVD;
	M.P[0].J01=1e-12;
	M.P[0].J02=1e-8;
	M.P[0].Jph=0;
	M.P[0].nid1=1;
	M.P[0].nid2=2;
	M.P[0].Eg=1.12;
	M.P[0].T=300;
	M.P[0].Rs=1e-5;
	M.P[0].Rsh=1e4;

	M.P[0].V=malloc(2*sizeof(double));
	M.P[0].J=malloc(2*sizeof(double));
	M.P[0].V[0]=-1;
	M.P[0].J[0]=-1;
	M.P[0].V[1]=1;
	M.P[0].J[1]=1;
	M.P[0].N=2;
	M.P[0].SplitX=1;
	M.P[0].SplitY=1;
	
	/* nodes */
	Nn=Nx*Ny;
	M.Nn=Nn;
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

void FreeMesh(mesh *M)
{
	int i;
	for (i=0;i<M->Nn;i++)
	{
		free(M->nodes[i].north);
		free(M->nodes[i].south);
		free(M->nodes[i].east);
		free(M->nodes[i].west);
	}
	free(M->nodes);
	for (i=0;i<M->Na;i++)
	{
		free(M->P[i].name);
		free(M->P[i].V);
		free(M->P[i].J);
	}
	free(M->P);
	for (i=0;i<M->res.Nva;i++)
		free(M->res.Vn[i]);
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
	
	Print(DEBUG,"Chaos reighns the mesh!\nStarting a slow, exhaustive search for the requested node\n");	
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
			Print(DEBUG, "node %d linked by %d.north\n",id, M.nodes[i].id);
		if (IsInList(M.nodes[i].south, id))
			Print(DEBUG, "node %d linked by %d.south\n",id, M.nodes[i].id);	
		if (IsInList(M.nodes[i].west, id))
			Print(DEBUG, "node %d linked by %d.west\n",id, M.nodes[i].id);
		if (IsInList(M.nodes[i].east, id))
			Print(DEBUG, "node %d linked by %d.east\n",id, M.nodes[i].id);
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
/* create a list of nodes which lie along the outline of the mesh */
{
	int *list;
	node *N;
	int loop=0;
	
	list=malloc(LISTBLOCK*sizeof(int));
	list[0]=0;	
	/* first we walk north untill we cannot continue */
	N=SearchNode(M, 0);
			
	while(N->north[0]>0)
		N=SearchNode(M, N->north[1]);
	/* go north-east */
	list=AddToList(list, N->id);
	while(N->east[0]>0)
	{
		/* add to path */		
		N=SearchNode(M, N->east[1]);
		while(N->north[0]>0)
			N=SearchNode(M, N->north[1]);
		list=AddToList(list, N->id);
	}
	/* go south-east */
	while(N->south[0]>0)
	{
		/* add to path */
		N=SearchNode(M, N->south[1]);
		while(N->east[0]>0)
			N=SearchNode(M, N->east[1]);
		list=AddToList(list, N->id);
	}
	/* go south-west */
	while(N->west[0]>0)
	{
		/* add to path */
		N=SearchNode(M, N->west[1]);
		while(N->south[0]>0)
			N=SearchNode(M, N->south[1]);
		list=AddToList(list, N->id);
	}
	/* go north-west */
	while(N->north[0]>0)
	{
		/* add to path */
		N=SearchNode(M, N->north[1]);
		while(N->west[0]>0)
			N=SearchNode(M, N->west[1]);
		list=AddToList(list, N->id);
	}
	/* go north-east */
	while((N->east[0]>0)&&(loop==0))
	{
		/* add to path */
		N=SearchNode(M, N->east[1]);
		while(N->north[0]>0)
			N=SearchNode(M, N->north[1]);
		loop=IsInList(list, N->id);
		list=AddToList(list, N->id);
	}
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

void DuplicateProperties(local_prop *dest, local_prop *source)
/* duplicate local properties struct */
{
	int i;
	dest = memcpy(dest, source, sizeof(local_prop));	
	dest->name=malloc((strlen(source->name)+1)*sizeof(char));
	strncpy(dest->name, source->name, strlen(source->name)+1);
	dest->V=malloc((source->N+1)*sizeof(double));
	dest->J=malloc((source->N+1)*sizeof(double));
	for (i=0;i<source->N;i++)
	{
		dest->V[i]=source->V[i];
		dest->J[i]=source->J[i];
	}
}

void DuplicateResults(mesh M, results *dest, results *source)
/* duplicate results struct */
{
	int i, j;	
	dest = memcpy(dest, source, sizeof(results));	
	dest->Va=malloc((dest->Nva+1)*sizeof(double));
	dest->I=malloc((dest->Nva+1)*sizeof(double));
	dest->Vn=malloc((dest->Nva+1)*sizeof(double *));
	for (i=0;i<dest->Nva;i++)
	{
		dest->Vn[i]=malloc((2*M.Nn+1)*sizeof(double));
		dest->Va[i]=source->Va[i];
		dest->I[i]=source->Va[i];
		for (j=0;j<M.Nn*2;j++)
			dest->Vn[i][j]=source->Vn[i][j];
		
	}
}

mesh DuplicateMesh(mesh M)
/* duplicate a mesh */
{
	mesh res;
	int i;
	res.Nn=M.Nn;
	res.nodes=malloc(res.Nn*sizeof(node));
	for (i=0;i<res.Nn;i++)
		DuplicateNode(M,res.nodes+i, M.nodes[i].id);
	res.Na=M.Na;
	res.P=malloc((res.Na+1)*sizeof(local_prop));
	for (i=0;i<res.Na;i++)
		DuplicateProperties(res.P+i, M.P+i);
	DuplicateResults(M, &(res.res), &(M.res));
	
	return res;

}

void NewProperties(mesh *M, char *name)
{
	M->P=realloc(M->P, (M->Na+1)*sizeof(local_prop));
	M->Na++;
	M->P[M->Na-1].name=malloc(MAXNAMELEN*sizeof(char));
	strncpy(M->P[M->Na-1].name, name,strlen(name)+1);
	M->P[M->Na-1].Rp=1;
	M->P[M->Na-1].Rn=1;
	M->P[M->Na-1].Rpvp=-1.0;
	M->P[M->Na-1].Rpvn=-1.0;
	M->P[M->Na-1].Rnvp=-1.0;
	M->P[M->Na-1].Rnvn=-1.0;
	
	M->P[M->Na-1].model=JVD;
	M->P[M->Na-1].J01=1e-12;
	M->P[M->Na-1].J02=1e-8;
	M->P[M->Na-1].Jph=0;
	M->P[M->Na-1].nid1=1;
	M->P[M->Na-1].nid2=2;
	M->P[M->Na-1].Eg=1.12;
	M->P[M->Na-1].T=300;
	M->P[M->Na-1].Rs=1e-5;
	M->P[M->Na-1].Rsh=1e4;
	
	M->P[M->Na-1].V=malloc(2*sizeof(double));
	M->P[M->Na-1].J=malloc(2*sizeof(double));
	M->P[M->Na-1].V[0]=-1;
	M->P[M->Na-1].J[0]=-1;
	M->P[M->Na-1].V[1]=1;
	M->P[M->Na-1].J[1]=1;
	M->P[M->Na-1].N=2;
	M->P[M->Na-1].SplitX=1;
	M->P[M->Na-1].SplitY=1;
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
	for (j=0;j<M->res.Nva;j++)
		free(M->res.Vn[j]);
	free(M->res.Va);
	free(M->res.I);
	M->res.Nva=0;
	M->res.Va=malloc(sizeof(double));
	M->res.I=malloc(sizeof(double));
	M->res.Vn=malloc(sizeof(double *));
}

void AssignPropertiesMesh(mesh *M, int P)
{
	int j;
	if ((P<0)||(P>=M->Na))
		Error("Area %d not defined\n", P);
	for (j=0;j<M->Nn;j++)
		M->nodes[j].P=P;
	for (j=0;j<M->res.Nva;j++)
		free(M->res.Vn[j]);
	free(M->res.Va);
	free(M->res.I);
	M->res.Nva=0;
	M->res.Va=malloc(sizeof(double));
	M->res.I=malloc(sizeof(double));
	M->res.Vn=malloc(sizeof(double *));
}

mesh JoinMeshes(mesh M1, mesh M2, double xoff, double yoff)
{
	mesh res;
	int *outl1, *outl2;
	int *prop_out;
	node *m1, *m2;
	int i,j;
	
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
			DuplicateProperties(res.P+res.Na-1, M2.P+i);
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
	for (i=0;i<res.res.Nva;i++)
		free(res.res.Vn[i]);
	free(res.res.Va);
	free(res.res.I);
	res.res.Nva=0;
	res.res.Va=malloc(sizeof(double));
	res.res.I=malloc(sizeof(double));
	res.res.Vn=malloc(sizeof(double *));

	return res;
}
mesh JoinMeshes_H(mesh M1, mesh M2, double yoff)
{
	mesh res;
	double xoff;
	int *outl1, *outl2;
	int *prop_out;
	node *m1, *m2;
	int i,j;
	
	outl1=MeshOutline(M1);
	xoff=M1.nodes[0].x1;
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(M1, outl1[i]);
		if (xoff<m1->x2)
			xoff=m1->x2;
	}	
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
			DuplicateProperties(res.P+res.Na-1, M2.P+i);
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
	for (i=0;i<res.res.Nva;i++)
		free(res.res.Vn[i]);
	free(res.res.Va);
	free(res.res.I);
	res.res.Nva=0;
	res.res.Va=malloc(sizeof(double));
	res.res.I=malloc(sizeof(double));
	res.res.Vn=malloc(sizeof(double *));

	return res;
}
mesh JoinMeshes_V(mesh M1, mesh M2, double xoff)
{
	mesh res;
	double yoff;
	int *outl1, *outl2;
	int *prop_out;
	node *m1, *m2;
	int i,j;
	
	outl1=MeshOutline(M1);
	yoff=M1.nodes[0].y1;
	for (i=1;i<=outl1[0];i++)
	{
		m1=SearchNode(M1, outl1[i]);
		if (yoff<m1->y2)
			yoff=m1->y2;
	}	
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
			DuplicateProperties(res.P+res.Na-1, M2.P+i);
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
	for (i=0;i<res.res.Nva;i++)
		free(res.res.Vn[i]);
	free(res.res.Va);
	free(res.res.I);
	res.res.Nva=0;
	res.res.Va=malloc(sizeof(double));
	res.res.I=malloc(sizeof(double));
	res.res.Vn=malloc(sizeof(double *));

	return res;
}

void SplitNodeX(int id, mesh *M)
{
	int newid, i, *list;
	node *old, *new, *a;
	/* add a node */
	newid=M->Nn;
	M->Nn++;
	M->nodes=realloc(M->nodes, (M->Nn+1)*sizeof(node));
	
	/* reallocate and copy node voltages */
	for (i=0;i<M->res.Nva;i++)
	{
		int j;
		M->res.Vn[i]=realloc(M->res.Vn[i],(2*M->Nn+1)*sizeof(double));
		for (j=2*M->Nn-2;j>M->Nn-2;j--)
			M->res.Vn[i][j]=M->res.Vn[i][j-1];
		M->res.Vn[i][M->Nn-1]=M->res.Vn[i][id];
		M->res.Vn[i][2*M->Nn-1]=M->res.Vn[i][id+M->Nn];
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
		if (R_Overlap(a->x1,new->x1, a->x2, new->x2) > TINY)
			a->south=AddToList(a->south, newid);
		else
			new->north=RemoveFromList(new->north, a->id);
		if (R_Overlap(a->x1,old->x1, a->x2, old->x2) < TINY)
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
		if (R_Overlap(a->x1,new->x1, a->x2, new->x2) > TINY)
			a->north=AddToList(a->north, newid);
		else
			new->south=RemoveFromList(new->south, a->id);
		if (R_Overlap(a->x1,old->x1, a->x2, old->x2) < TINY)
		{
			old->south=RemoveFromList(old->south, a->id);
			a->north=RemoveFromList(a->north, id);
		}
	}
	free(list);
}
	
void SplitNodeY(int id, mesh *M)
{
	int newid, i, *list;
	node *old, *new, *a;
	/* add a node */
	newid=M->Nn;
	M->Nn++;
	M->nodes=realloc(M->nodes, (M->Nn+1)*sizeof(node));
	
	/* reallocate and copy node voltages */
	for (i=0;i<M->res.Nva;i++)
	{
		int j;
		M->res.Vn[i]=realloc(M->res.Vn[i],(2*M->Nn+1)*sizeof(double));
		for (j=2*M->Nn-2;j>M->Nn-2;j--)
			M->res.Vn[i][j]=M->res.Vn[i][j-1];
		M->res.Vn[i][2*M->Nn-1]=M->res.Vn[i][id+M->Nn];
		M->res.Vn[i][M->Nn-1]=M->res.Vn[i][id];
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
		if (R_Overlap(a->y1,new->y1, a->y2, new->y2) > TINY)
			a->east=AddToList(a->east, newid);
		else
			new->west=RemoveFromList(new->west, a->id);
		if (R_Overlap(a->y1,old->y1, a->y2, old->y2) < TINY)
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
		if (R_Overlap(a->y1,new->y1, a->y2, new->y2) > TINY)
			a->west=AddToList(a->west, newid);
		else
			new->east=RemoveFromList(new->east, a->id);
		if (R_Overlap(a->y1,old->y1, a->y2, old->y2) < TINY)
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
   I tried mtrace which finds no leaks. Either there is some hidden bug in here or valgrind gets confused by 
   this. If you want to use valgrind then do not simplify your mesh. */
void CleanUpMesh(mesh *M, int *merged)
{
	int i, j;
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
		else
		{
			node *L;
			int k;
			for (k=1;k<=M->nodes[i].north[0];k++)
			{
				L=SearchNode(*M,M->nodes[i].north[k]);
				L->south=RemoveFromList(L->south, M->nodes[i].id);
			}
			for (k=1;k<=M->nodes[i].south[0];k++)
			{
				L=SearchNode(*M,M->nodes[i].south[k]);
				L->north=RemoveFromList(L->north, M->nodes[i].id);
			}
			for (k=1;k<=M->nodes[i].east[0];k++)
			{
				L=SearchNode(*M,M->nodes[i].east[k]);
				L->west=RemoveFromList(L->west, M->nodes[i].id);
			} 
			for (k=1;k<=M->nodes[i].west[0];k++)
			{
				L=SearchNode(*M,M->nodes[i].west[k]);
				L->east=RemoveFromList(L->east, M->nodes[i].id);
			}
			
			free(M->nodes[i].north);
			free(M->nodes[i].south);
			free(M->nodes[i].east);
			free(M->nodes[i].west);
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
	for (i=0;i<M->res.Nva;i++)
		free(M->res.Vn[i]);
	free(M->res.Va);
	free(M->res.I);
	M->res.Nva=0;
	M->res.Va=malloc(sizeof(double));
	M->res.I=malloc(sizeof(double));
	M->res.Vn=malloc(sizeof(double *));

}

#define MaxR 2.0
void Chunkify_east(mesh *M)
{
	int i, j;
	int *merged;
	int *list;
	int MRG;
	node *N, *Nl;
	merged=malloc(LISTBLOCK*sizeof(int));
	merged[0]=0;
	Print(DEBUG,"Chunkify_east\n");
	for (i=0;i<M->Nn;i++)
	{
		if (!IsInList(merged, M->nodes[i].id))
		{
			/* check east nodes to merge with */
			double R;	
			double newx2;
			
			N=M->nodes+i;		
			j=1;
			MRG=1;	
			while ((j<=N->east[0])&&(MRG==1))
			{
				Nl=SearchNode(*M,N->east[j]);
				if (j==1)
					newx2=Nl->x2;
				newx2=MIN(newx2,Nl->x2);
					
				if (Nl->P!=N->P)
					MRG=0;											
				else if ((Nl->y2>N->y2+TINY)||(Nl->y1<N->y1-TINY))
					MRG=0;	
				j++;
			}
			
			R=(newx2-N->x1)/(N->y2-N->y1);
			if (R>MaxR)
				MRG=0;
				
						
			if ((MRG)&&(N->east[0]>0))
			{
				int *a_nodes;
				Print(DEBUG, "R: %e\n",R);
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
		}
		
	}
	/* cleanup mesh, i.e. remove the merged nodes and sort the node id's */
	CleanUpMesh(M, merged);
	free(merged);
}
void Chunkify_west(mesh *M)
{
	int i, j;
	int *merged;
	int *list;
	int MRG;
	node *N, *Nl;
	merged=malloc(LISTBLOCK*sizeof(int));
	merged[0]=0;
	Print(DEBUG,"Chunkify_west\n");
	
	for (i=0;i<M->Nn;i++)
	{
		if (!IsInList(merged, M->nodes[i].id))
		{
			/* check east nodes to merge with */	
			double R;
			double newx1;	
			N=M->nodes+i;		
			j=1;
			MRG=1;
			while ((j<=N->west[0])&&(MRG==1))
			{
				Nl=SearchNode(*M,N->west[j]);
				if (j==1)
					newx1=Nl->x1;
				newx1=MAX(newx1,Nl->x1);
					
				if (Nl->P!=N->P)
					MRG=0;										
				else if ((Nl->y2>N->y2+TINY)||(Nl->y1<N->y1-TINY))
					MRG=0;	
				j++;
			}
			
			R=(N->x2-newx1)/(N->y2-N->y1);
			if (R>MaxR)
				MRG=0;	
				
			if ((MRG)&&(N->west[0]>0))
			{
				int *a_nodes;
				Print(DEBUG, "R: %e\n",R);
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
		}
	}
	/* cleanup mesh, i.e. remove the merged nodes and sort the node id's */
	CleanUpMesh(M, merged);
	free(merged);
}
void Chunkify_north(mesh *M)
{
	int i, j;
	int *merged;
	int *list;
	int MRG;
	node *N, *Nl;
	merged=malloc(LISTBLOCK*sizeof(int));
	merged[0]=0;
	Print(DEBUG,"Chunkify_north\n");
	
	for (i=0;i<M->Nn;i++)
	{
		if (!IsInList(merged, M->nodes[i].id))
		{
			double R;
			double newy2;	
			N=M->nodes+i;		
			j=1;
			MRG=1;
			while ((j<=N->north[0])&&(MRG==1))
			{
				Nl=SearchNode(*M,N->north[j]);
				if (j==1)
					newy2=Nl->y2;
				newy2=MIN(newy2,Nl->y2);
					
				if (Nl->P!=N->P)
					MRG=0;										
				else if ((Nl->x2>N->x2+TINY)||(Nl->x1<N->x1-TINY))
					MRG=0;	
				j++;
			}
			R=(N->x2-N->x1)/(newy2-N->y1);
			if (R<1.0/MaxR)
				MRG=0;		
			if ((MRG)&&(N->north[0]>0))
			{
				int *a_nodes;
				Print(DEBUG, "R: %e\n",R);
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
		}
	}
	/* cleanup mesh, i.e. remove the merged nodes and sort the node id's */
	CleanUpMesh(M, merged);
	free(merged);
}
void Chunkify_south(mesh *M)
{
	int i, j;
	int *merged;
	int *list;
	int MRG;
	node *N, *Nl;
	merged=malloc(LISTBLOCK*sizeof(int));
	merged[0]=0;
	Print(DEBUG,"Chunkify_south\n");
	
	for (i=0;i<M->Nn;i++)
	{
		if (!IsInList(merged, M->nodes[i].id))
		{
			double R;
			double newy1;	
			N=M->nodes+i;		
			j=1;
			MRG=1;
			while ((j<=N->south[0])&&(MRG==1))
			{
				Nl=SearchNode(*M,N->south[j]);
				
				if (j==1)
					newy1=Nl->y1;
				newy1=MAX(newy1,Nl->y1);
				
				if (Nl->P!=N->P)
					MRG=0;						
				else if ((Nl->x2>N->x2+TINY)||(Nl->x1<N->x1-TINY))
					MRG=0;	
				j++;
			}
			R=(N->x2-N->x1)/(N->y2-newy1);
			if (R<1.0/MaxR)
				MRG=0;		
			if ((MRG)&&(N->south[0]>0))
			{
				int *a_nodes;
				Print(DEBUG, "R: %e\n",R);
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
		}
	}
	/* cleanup mesh, i.e. remove the merged nodes and sort the node id's */
	CleanUpMesh(M, merged);
	free(merged);
}
/************************************************/
int *Chunkify_node(mesh *M, int id, int * merged)
{
	int i, j;
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

int Chunkify_nodes(mesh *M, int skip, int offset)
{
	int i,j,Nold;
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
	/* cleanup mesh, i.e. remove the merged nodes and sort the node id's */
	CleanUpMesh(M, merged);
	return Nold-M->Nn;
}
void Chunkify(mesh *M, int J)
{
	int Nold, i=0, skip=1, offset=0;
	Nold=M->Nn+1; 
		
	/*
	while ((M->Nn<Nold)&&(i<100))
	{
		Nold=M->Nn;
		Chunkify_north(M);
		Chunkify_east(M);
		Chunkify_south(M);
		Chunkify_west(M);
		i++;	
	}*/
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
		Print(DEBUG,"Skip:%i Offset %i\n",skip,offset);
		Chunkify_nodes(M, skip, 0);
		offset/=2;
		while (skip>1)
		{
			Print(DEBUG,"Skip:%i Offset %i\n",skip,offset);
			Chunkify_nodes(M, skip, offset);
			offset/=2;
			skip/=2;
		}
		i++;
		Print(DEBUG,"Round %i, %i nodes left\n",i, M->Nn);
	}
}

/*dump a node in a file */
void WriteNode(FILE *f, node N)
{
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
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	N.north=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.north[0]=a;
	if (fread(N.north+1, sizeof(int), N.north[0], f)<N.north[0])
		Error("Premature end of mesh file\n");
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	N.south=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.south[0]=a;
	if (fread(N.south+1, sizeof(int), N.south[0], f)<N.south[0])
		Error("Premature end of mesh file\n");
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	N.west=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.west[0]=a;
	if (fread(N.west+1, sizeof(int), N.west[0], f)<N.west[0])
		Error("Premature end of mesh file\n");
	
	if (!fread(&a, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	N.east=malloc(((a+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	N.east[0]=a;
	if (fread(N.east+1, sizeof(int), N.east[0], f)<N.east[0])
		Error("Premature end of mesh file\n");
	if (!fread(&N.id, sizeof(int), 1, f))
		Error("Premature end of fmesh ile\n");
	if (!fread(&N.x1, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&N.y1, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&N.x2, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&N.y2, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&N.P, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	return N;
} 
void WriteNodeArray(FILE *f, mesh *M)
{
	int i;
	fwrite(&(M->Nn), sizeof(int), 1, f);
	for (i=0;i<M->Nn;i++)
		WriteNode(f, M->nodes[i]);
}

void ReadNodeArray(FILE *f, mesh *M)
{
	int i;
	if (!fread(&i, sizeof(int), 1, f))
		Error("Premature end of file\n");
	M->Nn=i;
	M->nodes=malloc((M->Nn+1)*sizeof(node));
	for (i=0;i<M->Nn;i++)
		M->nodes[i]=ReadNode(f);
}

/* write local properties array */
void WriteProperties(FILE *f, local_prop P)
{
	int len;
	len=strlen(P.name)+1;
	fwrite(&len, sizeof(int), 1, f);
	fwrite(P.name, sizeof(char), len, f);
	fwrite(&P.Rp, sizeof(double), 1, f);
	fwrite(&P.Rn, sizeof(double), 1, f);
	fwrite(&P.Rpvp, sizeof(double), 1, f);
	fwrite(&P.Rpvn, sizeof(double), 1, f);
	fwrite(&P.Rnvp, sizeof(double), 1, f);
	fwrite(&P.Rnvn, sizeof(double), 1, f);
	
	fwrite(&P.model, sizeof(diode_model), 1, f);
	fwrite(&P.J01, sizeof(double), 1, f);
	fwrite(&P.J02, sizeof(double), 1, f);
	fwrite(&P.Jph, sizeof(double), 1, f);
	fwrite(&P.nid1, sizeof(double), 1, f);
	fwrite(&P.nid2, sizeof(double), 1, f);
	fwrite(&P.Eg, sizeof(double), 1, f);
	fwrite(&P.T, sizeof(double), 1, f);
	fwrite(&P.Rs, sizeof(double), 1, f);
	fwrite(&P.Rsh, sizeof(double), 1, f);

	
	fwrite(&P.N, sizeof(int), 1, f);
	fwrite(P.V, sizeof(double), P.N, f);
	fwrite(P.J, sizeof(double), P.N, f);
	fwrite(&P.SplitX, sizeof(int), 1, f);
	fwrite(&P.SplitY, sizeof(int), 1, f);	
}
local_prop ReadProperties(FILE *f)
{
	local_prop P;
	int i;
	if (!fread(&i, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	P.name=malloc(MAXNAMELEN*sizeof(char));
	if (fread(P.name, sizeof(char), i, f)<i)
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rp, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rn, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rpvp, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rpvn, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rnvp, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rnvn, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	
	if (!fread(&P.model, sizeof(diode_model), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.J01, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.J02, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Jph, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.nid1, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.nid2, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Eg, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.T, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");	
	if (!fread(&P.Rs, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.Rsh, sizeof(double), 1, f))
		Error("Premature end of mesh file\n");	
	
	if (!fread(&P.N, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	P.V=malloc((P.N+1)*sizeof(double));	
	if (fread(P.V, sizeof(double), P.N, f)<P.N)
		Error("Premature end of mesh file\n");
	P.J=malloc((P.N+1)*sizeof(double));	
	if (fread(P.J, sizeof(double), P.N, f)<P.N)
		Error("Premature end of mesh file\n");
	if (!fread(&P.SplitX, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	if (!fread(&P.SplitY, sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	return P;		
}
void WritePropertiesArray(FILE *f, mesh *M)
{
	int i;
	fwrite(&(M->Na), sizeof(int), 1, f);
	for (i=0;i<M->Na;i++)
		WriteProperties(f, M->P[i]);
}
void ReadPropertiesArray(FILE *f, mesh *M)
{
	int i;
	if (!fread(&i, sizeof(int), 1, f))
		Error("Premature end of file\n");
	M->Na=i;
	M->P=malloc((M->Nn+1)*sizeof(local_prop));
	for (i=0;i<M->Na;i++)
		M->P[i]=ReadProperties(f);
}

void WriteResults(FILE *f, mesh *M)
{
	int i;
	
	fwrite(&(M->res.Nva), sizeof(int), 1, f);
	fwrite(M->res.Va, sizeof(double), M->res.Nva, f);
	fwrite(M->res.I, sizeof(double), M->res.Nva, f);
	for (i=0;i<M->res.Nva;i++)
		fwrite(M->res.Vn[i], sizeof(double), 2*M->Nn, f);
		
}

void ReadResults(FILE *f, mesh *M)
{
	int i;

	if (!fread(&(M->res.Nva), sizeof(int), 1, f))
		Error("Premature end of mesh file\n");
	
	M->res.Va=malloc((M->res.Nva+1)*sizeof(double));
	M->res.I=malloc((M->res.Nva+1)*sizeof(double));
	M->res.Vn=malloc((M->res.Nva+1)*sizeof(double *));
	
	if (fread(M->res.Va, sizeof(double), M->res.Nva, f)<M->res.Nva)
		Error("Premature end of mesh file\n");
	if (fread(M->res.I, sizeof(double), M->res.Nva, f)<M->res.Nva)
		Error("Premature end of mesh file\n");
	
	for (i=0;i<M->res.Nva;i++)
	{
		M->res.Vn[i]=malloc((2*M->Nn+1)*sizeof(double));
		if (fread(M->res.Vn[i], sizeof(double), 2*M->Nn, f)<2*M->Nn)
			Error("Premature end of mesh file\n");
	}
		
}

/* write a complete mesh to a file */
void WriteMesh(char *fn, mesh *M)
{
	FILE *f;
	if ((f=fopen(fn,"wb"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	WriteNodeArray(f, M);
	WritePropertiesArray(f, M);
	WriteResults(f, M);
	fclose(f);
}

/* read a complete mesh from a file */
void ReadMesh(char *fn, mesh *M)
{
	FILE *f;
	if ((f=fopen(fn,"rb"))==NULL)
		Error("Cannot open %s for reading\n", fn);
	ReadNodeArray(f, M);
	ReadPropertiesArray(f, M);
	ReadResults(f, M);
	fclose(f);
}