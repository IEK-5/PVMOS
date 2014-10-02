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
 * export data from a mesh                                       *
 *                                                               *            
 *****************************************************************/     
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mesh2d.h"
#include "main.h"
#include "list.h"
#include "utils.h"
#include "solve.h"
#include "select_nodes.h"
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))

void Jfield(mesh *M, int Vai, double **Jx, double **Jy, double **Ex, double **Ey)
/* compute current density and electric field in each node */
/* Input: Mesh M, Solution index Vai, two dimensional arrays for current densities and electric field in x an d y directions */
/* the two dimensions for these arrays are electrode index and node index */

{
	int i, j, k;
	node N1, N2;
	double J, *R;
	
	
	if ((Vai>=M->res.Nva)||(Vai<0))
		Error("No simulation with index %i available", Vai);
	
	R=malloc((M->Nel+1)*sizeof(double));
	
	/* remove any possible old data */
	for (i=0;i<M->Nn;i++)
		for (j=0;j<M->Nel;j++)
		{
			if (Jx)
				Jx[j][i]=0;
			if (Jy)
				Jy[j][i]=0;
			if (Ex)
				Ex[j][i]=0;
			if (Ey)
				Ey[j][i]=0;
		}
	
	for (i=0;i<M->Nn;i++)
	{
		N1=M->nodes[i];
		for (j=1;j<=N1.north[0];j++)
		{
			N2=(*SearchNode(*M, N1.north[j]));
			Resistance(*M, N1,N2, R);
			
			for (k=0;k<M->Nel;k++)
			{
				J=(M->res.Vn[Vai][k][N1.id]-M->res.Vn[Vai][k][N2.id])/R[k]/2;
				if (Jy)
				{
					Jy[k][N1.id]+=J/(N1.x2-N1.x1);
					Jy[k][N2.id]+=J/(N2.x2-N2.x1);
				
				}
				if (Ey)
				{
					Ey[k][N1.id]-=M->P[N1.P].Rel[k]*J/(N1.x2-N1.x1);
					Ey[k][N2.id]-=M->P[N2.P].Rel[k]*J/(N2.x2-N2.x1);
				}
			}
		}
		for (j=1;j<=N1.east[0];j++)
		{
			N2=(*SearchNode(*M, N1.east[j]));
			Resistance(*M, N1,N2, R);
			for (k=0;k<M->Nel;k++)
			{
				J=(M->res.Vn[Vai][k][N1.id]-M->res.Vn[Vai][k][N2.id])/R[k]/2;
				if (Jx)
				{
					Jx[k][N1.id]+=J/(N1.y2-N1.y1);
					Jx[k][N2.id]+=J/(N2.y2-N2.y1);
				
				}
				if (Ex)
				{
					Ex[k][N1.id]-=M->P[N1.P].Rel[k]*J/(N1.y2-N1.y1);
					Ex[k][N2.id]-=M->P[N2.P].Rel[k]*J/(N2.y2-N2.y1);
				}
			}
		}
	}
	free(R);
}

void LocalVoltage(mesh *M, int Vai, node N, double **Ex, double **Ey, double x, double y, double *V)
/* interpolate voltage within node using the electric field within the nodes */
{
	int i;
	for (i=0;i<M->Nel;i++)
		V[i]=M->res.Vn[Vai][i][N.id]+Ex[i][N.id]*(x-(N.x1+N.x2)/2)+Ey[i][N.id]*(y-(N.y1+N.y2)/2);
}

/*
void SurfVPlotNearest(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	x_step=(x2-x1)/((double)Nx);
	y_step=(y2-y1)/((double)Ny);
	x=x1;
	for (i=0;i<=Nx;i++)
	{
		y=y1;
		for (j=0;j<=Ny;j++)
		{	
			node N;
			double Vn, Vp;
			ln_y=FindPos(*M, ln_y, x, y);
			N=*SearchNode(*M,ln_y);
			fprintf(f,"%e %e %e %e\n", x, y, M->res.Vn[Vai][N.id], M->res.Vn[Vai][N.id+M->Nn]);
			if (j==0)
				ln_x=ln_y;
			y+=y_step;
			
		}
		fprintf(f,"\n");
		ln_y=ln_x;
		x+=x_step;
	}
}
*/

void SurfVPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, k, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double **Ex, **Ey, *V;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	Ex=malloc(M->Nel*sizeof(double *));
	Ey=malloc(M->Nel*sizeof(double *));
	V=malloc(M->Nel*sizeof(double));
	for (i=0;i<M->Nel;i++)
	{
		Ex[i]=malloc(M->Nn*sizeof(double));
		Ey[i]=malloc(M->Nn*sizeof(double));
	}
	/* compute the electric field in each node for each electrode */
	Jfield(M, Vai, NULL, NULL, Ex, Ey);
	
	/* scan a regular mesh */
	x_step=(x2-x1)/((double)Nx);
	y_step=(y2-y1)/((double)Ny);
	
	x=x1;
	for (i=0;i<=Nx;i++)
	{
		y=y1;
		for (j=0;j<=Ny;j++)
		{	
			node N;
			/* here I tried to optimize search performance. To this end I use the variables ln_y and ln_x */
			/* The FindPos routine searches a node by simply walking from the start node toward the desired coordinate */
			/* It is thus useful to try and choose a start node as close as possible to the desired coordinate */
			/* at the beginning of the routine we have no idea but after that we start always at a nearby node */
			ln_y=FindPos(*M, ln_y, x, y);
			N=*SearchNode(*M,ln_y);
			LocalVoltage(M, Vai, N, Ex, Ey, x, y, V);
			
			fprintf(f,"%e %e", x, y);
			for (k=0;k<M->Nel;k++)
				fprintf(f," %e",V[k]);
			fprintf(f,"\n");
				
			if (j==0)
				ln_x=ln_y;
			y+=y_step;
			
		}
		fprintf(f,"\n");
		ln_y=ln_x;
		x+=x_step;
	}
	for (i=0;i<M->Nel;i++)
	{
		free(Ex[i]);
		free(Ey[i]);
	}
	free(Ex);
	free(Ey);
	free(V);
	fclose(f);
}

void SurfPPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, k, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double **Ex, **Ey, **Jx, **Jy, *V;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
		
	Ex=malloc(M->Nn*sizeof(double *));
	Ey=malloc(M->Nn*sizeof(double *));
	Jx=malloc(M->Nn*sizeof(double *));
	Jy=malloc(M->Nn*sizeof(double *));
	V=malloc(M->Nel*sizeof(double));
	for (i=0;i<M->Nel;i++)
	{
		Ex[i]=malloc(M->Nn*sizeof(double));
		Ey[i]=malloc(M->Nn*sizeof(double));
		Jx[i]=malloc(M->Nn*sizeof(double));
		Jy[i]=malloc(M->Nn*sizeof(double));
	}
	Jfield(M, Vai, Jx, Jy, Ex, Ey);
	x_step=(x2-x1)/((double)Nx);
	y_step=(y2-y1)/((double)Ny);
	x=x1;
	for (i=0;i<=Nx;i++)
	{
		y=y1;
		for (j=0;j<=Ny;j++)
		{	
			node N;
			double P, Pj;
			double I; 
			ln_y=FindPos(*M, ln_y, x, y);
			N=*SearchNode(*M,ln_y);
			LocalVoltage(M, Vai, N, Ex, Ey, x, y, V);
			fprintf(f,"%e %e", x, y);
			for (k=0;k<M->Nel;k++)
			{
				P=sqrt((Jx[k][ln_y]*Jx[k][ln_y]+Jy[k][ln_y]*Jy[k][ln_y])*(Ex[k][ln_y]*Ex[k][ln_y]+Ey[k][ln_y]*Ey[k][ln_y]));
				if (k>0)
				{
					Diode(*M, N, k-1, V[k-1]-V[k], &I, NULL);
					Pj=(V[k-1]-V[k])*I/((N.x2-N.x1)*(N.y2-N.y1));
					fprintf(f," %e %e", Pj, P);
				}
				else
					fprintf(f," %e", P);
				
			}
			fprintf(f,"\n");
			
			if (j==0)
				ln_x=ln_y;
			y+=y_step;
			
		}
		fprintf(f,"\n");
		ln_y=ln_x;
		x+=x_step;
	}
	free(Ex);
	free(Ey);
	fclose(f);
}

void PrintMesh(char *fn, mesh *M)
{
	int i;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		fprintf(f,"%e %e %i\n", M->nodes[i].x1, M->nodes[i].y1, M->nodes[i].id);
		fprintf(f,"%e %e %i\n", M->nodes[i].x2, M->nodes[i].y1, M->nodes[i].id);		
		fprintf(f,"%e %e %i\n", M->nodes[i].x2, M->nodes[i].y2, M->nodes[i].id);
		fprintf(f,"%e %e %i\n", M->nodes[i].x1, M->nodes[i].y2, M->nodes[i].id);
		fprintf(f,"%e %e %i\n\n", M->nodes[i].x1, M->nodes[i].y1, M->nodes[i].id);
		fprintf(f,"\n");		
	}
	fclose(f);
}
void PrintSurfDef(char *fn, mesh *M)
{
	int i;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		fprintf(f,"%e %e %i %i\n", M->nodes[i].x1, M->nodes[i].y1, M->nodes[i].id, M->nodes[i].P);
		fprintf(f,"%e %e %i %i\n\n", M->nodes[i].x1, M->nodes[i].y2, M->nodes[i].id, M->nodes[i].P);
		fprintf(f,"%e %e %i %i\n", M->nodes[i].x2, M->nodes[i].y1, M->nodes[i].id, M->nodes[i].P);		
		fprintf(f,"%e %e %i %i\n\n", M->nodes[i].x2, M->nodes[i].y2, M->nodes[i].id, M->nodes[i].P);
		fprintf(f,"\n");
	}
	fclose(f);
}

void PrintSurfV(char *fn, mesh *M)
{
	int i, j, k;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		fprintf(f,"%e %e", M->nodes[i].x1, M->nodes[i].y1);
		for(j=0;j<M->res.Nva;j++)
			for(k=0;k<M->Nel;k++)
				fprintf(f," %e", M->res.Vn[j][k][i]);
		fprintf(f,"\n");	
					
		fprintf(f,"%e %e", M->nodes[i].x1, M->nodes[i].y2);
		for(j=0;j<M->res.Nva;j++)
			for(k=0;k<M->Nel;k++)
				fprintf(f," %e", M->res.Vn[j][k][i]);
		fprintf(f,"\n\n");	
			
			
		fprintf(f,"%e %e", M->nodes[i].x2, M->nodes[i].y1);
		for(j=0;j<M->res.Nva;j++)
			for(k=0;k<M->Nel;k++)
				fprintf(f," %e", M->res.Vn[j][k][i]);
		fprintf(f,"\n");	
			
					
		fprintf(f,"%e %e", M->nodes[i].x2, M->nodes[i].y2);
		for(j=0;j<M->res.Nva;j++)
			for(k=0;k<M->Nel;k++)
				fprintf(f," %e", M->res.Vn[j][k][i]);
		fprintf(f,"\n\n");	
		
		fprintf(f,"\n");
	}
	fclose(f);
}
void PrintConn(char *fn, mesh *M)
{
	int i, j;
	double xn, yn;
	double xnn, ynn;
	node * N;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		xn=(M->nodes[i].x1+M->nodes[i].x2)/2;
		yn=(M->nodes[i].y2+M->nodes[i].y1)/2;
		for (j=1;j<=M->nodes[i].north[0];j++)
		{
			N=SearchNode(*M, M->nodes[i].north[j]);
			xnn=(N->x1+N->x2)/2;
			ynn=(N->y2+N->y1)/2;
			fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].north[j],M->nodes[i].id);
		
		} 
		for (j=1;j<=M->nodes[i].south[0];j++)
		{
			N=SearchNode(*M, M->nodes[i].south[j]);
			xnn=(N->x1+N->x2)/2;
			ynn=(N->y2+N->y1)/2;
			fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].south[j],M->nodes[i].id);
		
		}
		for (j=1;j<=M->nodes[i].east[0];j++)
		{
			N=SearchNode(*M, M->nodes[i].east[j]);
			xnn=(N->x1+N->x2)/2;
			ynn=(N->y2+N->y1)/2;
			fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].east[j],M->nodes[i].id);
		
		}
		for (j=1;j<=M->nodes[i].west[0];j++)
		{
			N=SearchNode(*M, M->nodes[i].west[j]);
			xnn=(N->x1+N->x2)/2;
			ynn=(N->y2+N->y1)/2;
			fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].west[j],M->nodes[i].id );
		
		} 
	}
	fclose(f);
}
void PrintMeshSel(char *fn, mesh *M, double x1, double y1, double x2, double y2)
{
	int i;
	double xx1,yy1,xx2,yy2;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		xx1=MAX(M->nodes[i].x1,x1);
		xx2=MIN(M->nodes[i].x2,x2);
		yy1=MAX(M->nodes[i].y1,y1);
		yy2=MIN(M->nodes[i].y2,y2);		
		if ((M->nodes[i].x2>xx1) && (M->nodes[i].x1<xx2) && (M->nodes[i].y2>yy1) && (M->nodes[i].y1<yy2))
		{
			fprintf(f,"%e %e %i\n", xx1, yy1, M->nodes[i].id);
			fprintf(f,"%e %e %i\n", xx2, yy1, M->nodes[i].id);		
			fprintf(f,"%e %e %i\n", xx2, yy2, M->nodes[i].id);
			fprintf(f,"%e %e %i\n", xx1, yy2, M->nodes[i].id);
			fprintf(f,"%e %e %i\n\n", xx1, yy1, M->nodes[i].id);
			fprintf(f,"\n");
		}		
	}
	fclose(f);
}
void PrintSurfDefSel(char *fn, mesh *M, double x1, double y1, double x2, double y2)
{
	int i;
	double xx1,yy1,xx2,yy2;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		xx1=MAX(M->nodes[i].x1,x1);
		xx2=MIN(M->nodes[i].x2,x2);
		yy1=MAX(M->nodes[i].y1,y1);
		yy2=MIN(M->nodes[i].y2,y2);		
		if ((M->nodes[i].x2>xx1) && (M->nodes[i].x1<xx2) && (M->nodes[i].y2>yy1) && (M->nodes[i].y1<yy2))
		{
			fprintf(f,"%e %e %i %i\n", xx1, yy1, M->nodes[i].id, M->nodes[i].P);
			fprintf(f,"%e %e %i %i\n\n", xx1, yy2, M->nodes[i].id, M->nodes[i].P);
			fprintf(f,"%e %e %i %i\n", xx2, yy1, M->nodes[i].id, M->nodes[i].P);		
			fprintf(f,"%e %e %i %i\n\n", xx2, yy2, M->nodes[i].id, M->nodes[i].P);
			fprintf(f,"\n");
		}		
	}
	fclose(f);
}

void PrintSurfVSel(char *fn, mesh *M, double x1, double y1, double x2, double y2)
{
	int i, j, k;
	double xx1,yy1,xx2,yy2;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		xx1=MAX(M->nodes[i].x1,x1);
		xx2=MIN(M->nodes[i].x2,x2);
		yy1=MAX(M->nodes[i].y1,y1);
		yy2=MIN(M->nodes[i].y2,y2);		
		if ((M->nodes[i].x2>xx1) && (M->nodes[i].x1<xx2) && (M->nodes[i].y2>yy1) && (M->nodes[i].y1<yy2))
		{
			fprintf(f,"%e %e", xx1, yy1);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][i]);
			fprintf(f,"\n");	
						
			fprintf(f,"%e %e", xx1, yy2);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][i]);
			fprintf(f,"\n\n");	
				
				
			fprintf(f,"%e %e", xx2, yy1);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][i]);
			fprintf(f,"\n");	
				
						
			fprintf(f,"%e %e", xx2, yy2);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][i]);
			fprintf(f,"\n\n");	
			
			fprintf(f,"\n");
		}
	}
	fclose(f);
}
void PrintConnSel(char *fn, mesh *M, double x1, double y1, double x2, double y2)
{
	int i, j;
	double xx1,yy1,xx2,yy2;
	double xn, yn;
	double xnn, ynn;
	node * N;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		xx1=MAX(M->nodes[i].x1,x1);
		xx2=MIN(M->nodes[i].x2,x2);
		yy1=MAX(M->nodes[i].y1,y1);
		yy2=MIN(M->nodes[i].y2,y2);		
		if ((M->nodes[i].x2>xx1) && (M->nodes[i].x1<xx2) && (M->nodes[i].y2>yy1) && (M->nodes[i].y1<yy2))
		{
			xn=(xx1+xx2)/2;
			yn=(yy2+yy1)/2;
			for (j=1;j<=M->nodes[i].north[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].north[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].north[j],M->nodes[i].id);
			
			} 
			for (j=1;j<=M->nodes[i].south[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].south[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].south[j],M->nodes[i].id);
			
			}
			for (j=1;j<=M->nodes[i].east[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].east[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].east[j],M->nodes[i].id);
			
			}
			for (j=1;j<=M->nodes[i].west[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].west[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].west[j],M->nodes[i].id );
			
			} 
		}
	}	
	fclose(f);
}
void PrintPars(char *fn, mesh *M)
{
	int i, j, k;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Na;i++)
	{
		fprintf(f,"*****************Parameters for area %s\n", M->P[i].name);
		fprintf(f,"id:     %i\n", i);
		for (k=0;k<M->Nel;k++)
		{
			fprintf(f,"Rel %i: %e\n", k, M->P[i].Rel[k]);
			fprintf(f,"Rvp %i: %e\tRvn %i: %e\n", k, M->P[i].Rvp[k], k, M->P[i].Rvn[k]);
			if (k>0)
			{
				fprintf(f,"Electrode Connection %i to %i:\n", k-1, k);
				switch (M->P[i].conn[k-1].model)
				{
					case JVD:
						fprintf(f,"V [V]\t\tJ [A/cm2]\n");
						for (j=0;j<M->P[i].conn[k-1].N;j++)
							fprintf(f,"%e\t%e\n", M->P[i].conn[k-1].V[j], M->P[i].conn[k-1].J[j]);
						break;
					case ONED:
						fprintf(f,"J0:     %e\tnid:    %e\n", M->P[i].conn[k-1].J01, M->P[i].conn[k-1].nid1);
						fprintf(f,"Jph:    %e\n", M->P[i].conn[k-1].Jph);
						fprintf(f,"Rs:     %e\tRsh:    %e\n", M->P[i].conn[k-1].Rs, M->P[i].conn[k-1].Rsh);
						fprintf(f,"Eg:     %e\n", M->P[i].conn[k-1].Eg);
						break;
					case TWOD:
						fprintf(f,"J01:    %e\tnid1:   %e\n", M->P[i].conn[k-1].J01, 1.0);
						fprintf(f,"J02:    %e\tnid2:   %e\n", M->P[i].conn[k-1].J02, 2.0);
						fprintf(f,"Jph:    %e\n", M->P[i].conn[k-1].Jph);
						fprintf(f,"Rs:     %e\tRsh:    %e\n", M->P[i].conn[k-1].Rs, M->P[i].conn[k-1].Rsh);
						fprintf(f,"Eg:     %e\n", M->P[i].conn[k-1].Eg);
						break;
				}
			}
		}
		fprintf(f,"T:      %e\n", M->P[i].T);
		fprintf(f,"SplitX: %d\t\tSplitY: %d\n", M->P[i].SplitX, M->P[i].SplitY);
		fprintf(f,"\n");
	}
	fclose(f);
}
void PrintIV(char *fn, mesh *M)
{
	int i;
	FILE *f;
	double *V, *I;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	
	V=malloc((M->res.Nva+1)*sizeof(double));
	I=malloc((M->res.Nva+1)*sizeof(double));
	for (i=0;i<M->res.Nva;i++)
	{
		V[i]=M->res.Va[i];
		I[i]=M->res.I[i];
	}
	BubbleSortJV(M->res.Nva, V, I);
	for (i=0;i<M->res.Nva;i++)
		fprintf(f,"%e %e\n", V[i], I[i]);
	free(V);
	free(I);
	fclose(f);
}

