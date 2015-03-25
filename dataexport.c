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
	PrintFileHeader(f);
	fprintf(f, "# Simulated Potentials mapped to a regular mesh\n");
	fprintf(f, "# V(i):      Potential in the i-th electrode\n");
	fprintf(f, "# x [cm]\ty [cm]\tU(i) [V]\tU(i+1) [V]...\n");
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

void SurfVjPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, k, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double **Ex, **Ey, *V;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	fprintf(f, "# Simulated Junction voltages mapped to a regular mesh\n");
	fprintf(f, "# Vj(i+0.5):      Junction voltage between the i-th and i+1 the electrode electrode\n");
	fprintf(f, "# x [cm]\ty [cm]\tVj(i+0.5) [V]\tVj(i+1.5) [V]...\n");
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
			double Vj;
			ln_y=FindPos(*M, ln_y, x, y);
			N=*SearchNode(*M,ln_y);
			LocalVoltage(M, Vai, N, Ex, Ey, x, y, V);
			fprintf(f,"%e %e", x, y);
			for (k=0;k<M->Nel;k++)
			{
				if (k>0)
				{
					Diode(*M, N, k-1, V[k-1]-V[k], NULL, NULL, &Vj);
					fprintf(f," %e", Vj);
				}
				
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
		
	PrintFileHeader(f);
	fprintf(f, "# Simulated Power Density mapped to a regular mesh\n");
	fprintf(f, "# P(i):      power dissipation in the i-th electrode\n");
	fprintf(f, "# P(i+0.5):  power dissipation in the connection beteen the i-th and i+1-th electrodes\n");	
	fprintf(f, "# x [cm]\ty [cm]\tP(i) [W/cm^2]\tP(i+0.5) [W/cm^2]\tP(i+1) [W/cm^2]...\n");
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
					Diode(*M, N, k-1, V[k-1]-V[k], &I, NULL, NULL);
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
	free(Jx);
	free(Jy);
	fclose(f);
}
void SurfJPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, k, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double **Ex, **Ey, **Jx, **Jy, *V;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
		
	PrintFileHeader(f);
	fprintf(f, "# Simulated Current Densities mapped to a regular mesh\n");
	fprintf(f, "# Jx(i) Jy(i): Current densities in x an y direction in the i-th electrode\n");
	fprintf(f, "# Jz(i+0.5):   Current densities in the connection beteen the i-th and i+1-th electrodes\n");
	fprintf(f, "# Take note of the units!\n");	
	fprintf(f, "# x [cm]\ty [cm]\tJx(i) [A/cm]\tJy(i) [A/cm]\tJz(i+0.5) [A/cm^2]...\n");
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
			double Jz;
			ln_y=FindPos(*M, ln_y, x, y);
			N=*SearchNode(*M,ln_y);
			LocalVoltage(M, Vai, N, Ex, Ey, x, y, V);
			fprintf(f,"%e %e", x, y);
			for (k=0;k<M->Nel;k++)
			{
				if (k>0)
				{
					Diode(*M, N, k-1, V[k-1]-V[k], &Jz, NULL, NULL);
					Jz/=((N.x2-N.x1)*(N.y2-N.y1));
					fprintf(f," %e %e %e", Jz, Jx[k][ln_y], Jy[k][ln_y]);
				}
				else
					fprintf(f," %e %e", Jx[k][ln_y], Jy[k][ln_y]);
				
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
	free(Jx);
	free(Jy);
	fclose(f);
}
void SurfEPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, k, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double **Ex, **Ey;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
		
	PrintFileHeader(f);
	fprintf(f, "# Simulated Electric Fields mapped to a regular mesh\n");
	fprintf(f, "# Ex(i) Ey(i):     Electric field in x an y direction in the i-th electrode\n");
	fprintf(f, "# x [cm]\ty [cm]\tEx(i) [V/cm]\tEy(i) [V/cm]...\n");
	Ex=malloc(M->Nn*sizeof(double *));
	Ey=malloc(M->Nn*sizeof(double *));
	for (i=0;i<M->Nel;i++)
	{
		Ex[i]=malloc(M->Nn*sizeof(double));
		Ey[i]=malloc(M->Nn*sizeof(double));
	}
	Jfield(M, Vai, NULL, NULL, Ex, Ey);
	x_step=(x2-x1)/((double)Nx);
	y_step=(y2-y1)/((double)Ny);
	x=x1;
	for (i=0;i<=Nx;i++)
	{
		y=y1;
		for (j=0;j<=Ny;j++)
		{	
			ln_y=FindPos(*M, ln_y, x, y);
			fprintf(f,"%e %e", x, y);
			for (k=0;k<M->Nel;k++)
				fprintf(f," %e %e", Ex[k][ln_y], Ey[k][ln_y]);
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

void PrintMesh(char *fn, mesh *M, int *selected)
{
	int i;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	if (selected[0]==0)
		fprintf(f, "# Mesh element borders\n");
	else
		fprintf(f, "# Selected mesh element borders\n");
	fprintf(f, "# gnuplot line plot\n");
	fprintf(f, "# x [cm]\ty [cm]\tElement-id\n");
	if (selected[0]==0)
		for (i=0;i<M->Nn;i++)
		{
			fprintf(f,"%e %e %i\n", M->nodes[i].x1, M->nodes[i].y1, M->nodes[i].id);
			fprintf(f,"%e %e %i\n", M->nodes[i].x2, M->nodes[i].y1, M->nodes[i].id);		
			fprintf(f,"%e %e %i\n", M->nodes[i].x2, M->nodes[i].y2, M->nodes[i].id);
			fprintf(f,"%e %e %i\n", M->nodes[i].x1, M->nodes[i].y2, M->nodes[i].id);
			fprintf(f,"%e %e %i\n", M->nodes[i].x1, M->nodes[i].y1, M->nodes[i].id);
			fprintf(f,"\n");		
		}
	else
	{
		node *N;
		for (i=1;i<=selected[0];i++)
		{
			
			N=SearchNode(*M, selected[i]);
			fprintf(f,"%e %e %i\n", N->x1, N->y1, N->id);
			fprintf(f,"%e %e %i\n", N->x2, N->y1, N->id);		
			fprintf(f,"%e %e %i\n", N->x2, N->y2, N->id);
			fprintf(f,"%e %e %i\n", N->x1, N->y2, N->id);
			fprintf(f,"%e %e %i\n", N->x1, N->y1, N->id);
			fprintf(f,"\n");		
		}
	}
		
	fclose(f);
}
void PrintSurfDef(char *fn, mesh *M, int *selected)
{
	int i;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	if (selected[0]==0)
		fprintf(f, "# Mesh element area definition\n");
	else
		fprintf(f, "# Selected mesh element area definition\n");
	fprintf(f, "# gnuplot surface plot\n");
	fprintf(f, "# x [cm]\ty [cm]\tElement-id\tArea-id\n");
	if (selected[0]==0)
		for (i=0;i<M->Nn;i++)
		{
			fprintf(f,"%e %e %i %i\n", M->nodes[i].x1, M->nodes[i].y1, M->nodes[i].id, M->nodes[i].P);
			fprintf(f,"%e %e %i %i\n\n", M->nodes[i].x1, M->nodes[i].y2, M->nodes[i].id, M->nodes[i].P);
			fprintf(f,"%e %e %i %i\n", M->nodes[i].x2, M->nodes[i].y1, M->nodes[i].id, M->nodes[i].P);		
			fprintf(f,"%e %e %i %i\n\n", M->nodes[i].x2, M->nodes[i].y2, M->nodes[i].id, M->nodes[i].P);
			fprintf(f,"\n");
		}
	else
	{
		node *N;
		for (i=1;i<=selected[0];i++)
		{			
			N=SearchNode(*M, selected[i]);			
			fprintf(f,"%e %e %i %i\n", N->x1, N->y1, N->id, N->P);
			fprintf(f,"%e %e %i %i\n\n", N->x1, N->y2, N->id, N->P);
			fprintf(f,"%e %e %i %i\n", N->x2, N->y1, N->id, N->P);		
			fprintf(f,"%e %e %i %i\n\n", N->x2, N->y2, N->id, N->P);
			fprintf(f,"\n");
		}
	}
	fclose(f);
}

void PrintSurfV(char *fn, mesh *M, int *selected)
{
	int i, j, k;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	if (selected[0]==0)
		fprintf(f, "# Simulated Potentials per element\n");
	else
		fprintf(f, "# Simulated Potentials per selected element\n");
	fprintf(f, "# Simulated Potentials per element\n");
	fprintf(f, "# gnuplot surface plot\n");
	fprintf(f, "# V(i):      Potential in the i-th electrode\n");
	fprintf(f, "# x [cm]\ty [cm]\tU(i) [V]\tU(i+1) [V]...\n");
	if (selected[0]==0)
		for (i=0;i<M->Nn;i++)
		{
			fprintf(f,"%e %e", M->nodes[i].x1, M->nodes[i].y1);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][M->nodes[i].id]);
			fprintf(f,"\n");	
						
			fprintf(f,"%e %e", M->nodes[i].x1, M->nodes[i].y2);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][M->nodes[i].id]);
			fprintf(f,"\n\n");	
				
				
			fprintf(f,"%e %e", M->nodes[i].x2, M->nodes[i].y1);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][M->nodes[i].id]);
			fprintf(f,"\n");	
				
						
			fprintf(f,"%e %e", M->nodes[i].x2, M->nodes[i].y2);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][M->nodes[i].id]);
			fprintf(f,"\n\n");	
			
			fprintf(f,"\n");
		}
	else
	{
		node *N;
		for (i=1;i<=selected[0];i++)
		{			
			N=SearchNode(*M, selected[i]);	
	
			fprintf(f,"%e %e", N->x1, N->y1);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][N->id]);
			fprintf(f,"\n");	
						
			fprintf(f,"%e %e", N->x1, N->y2);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][N->id]);
			fprintf(f,"\n\n");	
				
				
			fprintf(f,"%e %e", N->x2, N->y1);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][N->id]);
			fprintf(f,"\n");	
				
						
			fprintf(f,"%e %e", N->x2, N->y2);
			for(j=0;j<M->res.Nva;j++)
				for(k=0;k<M->Nel;k++)
					fprintf(f," %e", M->res.Vn[j][k][N->id]);
			fprintf(f,"\n\n");	
			
			fprintf(f,"\n");
		}
	}
	fclose(f);
}
void PrintConn(char *fn, mesh *M, int *selected)
{
	int i, j;
	double xn, yn;
	double xnn, ynn;
	node * N;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	if (selected[0]==0)
		fprintf(f, "# Connections between elements\n");
	else
		fprintf(f, "# Connections between selected elements\n");
	fprintf(f, "# gnuplot vector plot\n");
	fprintf(f, "# x [cm]\ty [cm]\tdx [cm]\tdy [cm]\telement_id1\t element_id2\n");
	if (selected[0]==0)
		for (i=0;i<M->Nn;i++)
		{
			xn=(M->nodes[i].x1+M->nodes[i].x2)/2;
			yn=(M->nodes[i].y2+M->nodes[i].y1)/2;
			fprintf(f, "# north of %d\n", M->nodes[i].id);
			for (j=1;j<=M->nodes[i].north[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].north[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].north[j],M->nodes[i].id);
			
			} 
			fprintf(f, "# south of %d\n", M->nodes[i].id);
			for (j=1;j<=M->nodes[i].south[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].south[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].south[j],M->nodes[i].id);
			
			}
			fprintf(f, "# east of %d\n", M->nodes[i].id);
			for (j=1;j<=M->nodes[i].east[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].east[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].east[j],M->nodes[i].id);
			
			}
			fprintf(f, "# west of %d\n", M->nodes[i].id);
			for (j=1;j<=M->nodes[i].west[0];j++)
			{
				N=SearchNode(*M, M->nodes[i].west[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, M->nodes[i].west[j],M->nodes[i].id );
			
			} 
		}
	else
	{
		for (i=1;i<=selected[0];i++)
		{			
			N=SearchNode(*M, selected[i]);
			xn=(N->x1+N->x2)/2;
			yn=(N->y2+N->y1)/2;
			for (j=1;j<=N->north[0];j++)
			{
				N=SearchNode(*M, N->north[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, N->north[j],N->id);
			
			} 
			for (j=1;j<=N->south[0];j++)
			{
				N=SearchNode(*M, N->south[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, N->south[j],N->id);
			
			}
			for (j=1;j<=N->east[0];j++)
			{
				N=SearchNode(*M, N->east[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, N->east[j],N->id);
			
			}
			for (j=1;j<=N->west[0];j++)
			{
				N=SearchNode(*M, N->west[j]);
				xnn=(N->x1+N->x2)/2;
				ynn=(N->y2+N->y1)/2;
				fprintf(f,"%e %e %e %e %i %i\n", xn,yn,xnn-xn, ynn-yn, N->west[j],N->id );
			
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
	PrintFileHeader(f);
	fprintf(f, "# Simulation parameters per area definition\n");
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
	PrintFileHeader(f);
	fprintf(f, "# Simulated Current-Voltage pairs\n");
	fprintf(f, "# U [V]\tI [A]\n");
	
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

void PrintProbe(char *fn, mesh *M, double x, double y)
{
	int i,k;
	double **Ex, **Ey, *V;
	node N;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	fprintf(f, "# Simulated Potentials at coordinate (%e %e)\n", x, y);
	fprintf(f, "# Va [V]\tVprobe [V]\n");
	Ex=malloc(M->Nel*sizeof(double *));
	Ey=malloc(M->Nel*sizeof(double *));
	V=malloc(M->Nel*sizeof(double));
	for (i=0;i<M->Nel;i++)
	{
		Ex[i]=malloc(M->Nn*sizeof(double));
		Ey[i]=malloc(M->Nn*sizeof(double));
	}
	                              
	N=*SearchNode(*M,FindPos(*M, 0, x, y));
	
	for (i=0;i<M->res.Nva;i++)
	{
		/* compute the electric field in each node for each electrode */
		Jfield(M, i, NULL, NULL, Ex, Ey);
		LocalVoltage(M, i, N, Ex, Ey, x, y, V);
		fprintf(f,"%e",M->res.Va[i]);
		for (k=0;k<M->Nel;k++)
			fprintf(f," %e",V[k]);
		fprintf(f,"\n");
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

void BubbleSortIIIV(int n, double *V, double *I1, double *I2, double *I3)
{
	int a, b, s;
	double c;
	for(a=0;a<(n-1);a++)
	{
		s=1;
		for(b=0;b<n-a-1;b++)
		{
			if (V[b]>V[b+1])
			{
				c=V[b];
				V[b]=V[b+1];
				V[b+1]=c;
				
				c=I1[b];
				I1[b]=I1[b+1];
				I1[b+1]=c;
				
				c=I2[b];
				I2[b]=I2[b+1];
				I2[b+1]=c;
				
				c=I3[b];
				I3[b]=I3[b+1];
				I3[b+1]=c;
				s=0;
			}
		}
		if (s)
			break;
	}
}

#define Area (N->x2-N->x1)*(N->y2-N->y1)
void PrintInIp(char *fn, mesh *M, int *selected)
{
	int i, j, k;
	FILE *f;
	node *N;
	double *V, *I, *In, *Ip;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	PrintFileHeader(f);
	fprintf(f, "# Total contact currents for selected elements\n");
	fprintf(f, "# U [V]\tIp [A]\tIn [A]\tI [A]\n");
	V=malloc((M->res.Nva+1)*sizeof(double));
	I=malloc((M->res.Nva+1)*sizeof(double));
	Ip=malloc((M->res.Nva+1)*sizeof(double));
	In=malloc((M->res.Nva+1)*sizeof(double));
	for (i=0;i<M->res.Nva;i++)
	{
		V[i]=M->res.Va[i];
		I[i]=M->res.I[i];
		Ip[i]=0;
		In[i]=0;
	}
	if (selected[0]==0)
		for (j=0;j<M->Nn;j++)
		{			
			N=M->nodes+j;
			
			for(i=0;i<M->res.Nva;i++)
			{
				for(k=0;k<M->Nel;k++)
				{
					if (M->P[N->P].Rvp[k]>0)
						Ip[i]+=Area*(V[i]-M->res.Vn[i][k][N->id])/M->P[N->P].Rvp[k];
					if (M->P[N->P].Rvn[k]>0)
						In[i]+=Area*M->res.Vn[i][k][N->id]/M->P[N->P].Rvn[k];
				
				}
			}
		}
	else
		for (j=1;j<=selected[0];j++)
		{			
			N=SearchNode(*M, selected[j]);
			
			for(i=0;i<M->res.Nva;i++)
			{
				for(k=0;k<M->Nel;k++)
				{
					if (M->P[N->P].Rvp[k]>0)
						Ip[i]+=Area*(V[i]-M->res.Vn[i][k][N->id])/M->P[N->P].Rvp[k];
					if (M->P[N->P].Rvn[k]>0)
						In[i]+=Area*M->res.Vn[i][k][N->id]/M->P[N->P].Rvn[k];
				
				}
			}
		}
	BubbleSortIIIV(M->res.Nva, V, In , Ip, I);
	
	for (i=0;i<M->res.Nva;i++)
		fprintf(f,"%e %e %e %e\n", V[i], Ip[i], In[i], I[i]);
	fclose(f);
	free(V);
	free(I);
	free(Ip);
	free(In);
}
