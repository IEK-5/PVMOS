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

void Jfield(mesh *M, int Vai, double *Jx, double *Jy, double *Ex, double *Ey)
{
	int i, j;
	node N1, N2;
	double J, Rp, Rn;
	double *Vn, *Vp, *Jxn, *Jxp, *Jyn, *Jyp, *Exn, *Exp, *Eyn, *Eyp;
	if (Vai>=M->res.Nva)
		Error("No simulation with index %i available", Vai);
	Vp=M->res.Vn[Vai];
	Vn=M->res.Vn[Vai]+M->Nn;
	for (i=0;i<M->Nn;i++)
	{
		if (Jx)
		{
			Jx[i]=0;
			Jx[i+M->Nn]=0;
		}
		if (Jy)
		{
			Jy[i]=0;
			Jy[i+M->Nn]=0;
		}
		if (Ex)
		{
			Ex[i]=0;
			Ex[i+M->Nn]=0;
		}
		if (Ey)
		{
			Ey[i]=0;
			Ey[i+M->Nn]=0;
		}
	}
	if (Jx)
	{
		Jxp=Jx;
		Jxn=Jx+M->Nn;
	}
	if (Jy)
	{
		Jyp=Jy;
		Jyn=Jy+M->Nn;
	}
	if (Ex)
	{
		Exp=Ex;
		Exn=Ex+M->Nn;
	}	
	if (Ey)
	{
		Eyp=Ey;
		Eyn=Ey+M->Nn;
	}
	
	for (i=0;i<M->Nn;i++)
	{
		N1=M->nodes[i];
		for (j=1;j<=N1.north[0];j++)
		{
			N2=(*SearchNode(*M, N1.north[j]));
			Resistance(*M, N1,N2, &Rp, &Rn);
			J=(Vn[N1.id]-Vn[N2.id])/Rn/2;
			if (Jy)
			{
				Jyn[N1.id]+=J/(N1.x2-N1.x1);
				Jyn[N2.id]+=J/(N2.x2-N2.x1);
			}
			if (Ey)
			{
				Eyn[N1.id]-=M->P[N1.P].Rn*J/(N1.x2-N1.x1);
				Eyn[N2.id]-=M->P[N2.P].Rn*J/(N2.x2-N2.x1);
			}
			J=(Vp[N1.id]-Vp[N2.id])/Rp/2;
			if (Jy)
			{
				Jyp[N1.id]+=J/(N1.x2-N1.x1);
				Jyp[N2.id]+=J/(N2.x2-N2.x1);
			}
			if (Ey)
			{
				Eyp[N1.id]-=M->P[N1.P].Rp*J/(N1.x2-N1.x1);
				Eyp[N2.id]-=M->P[N2.P].Rp*J/(N2.x2-N2.x1);
			}			
		}
		for (j=1;j<=N1.east[0];j++)
		{
			N2=(*SearchNode(*M, N1.east[j]));
			Resistance(*M, N1,N2, &Rp, &Rn);
			J=(Vn[N1.id]-Vn[N2.id])/Rn/2;
			if (Jx)
			{
				Jxn[N1.id]+=J/(N1.y2-N1.y1);
				Jxn[N2.id]+=J/(N2.y2-N2.y1);
			}
			if (Ex)
			{
				Exn[N1.id]-=M->P[N1.P].Rn*J/(N1.y2-N1.y1);
				Exn[N2.id]-=M->P[N2.P].Rn*J/(N2.y2-N2.y1);
			}
			J=(Vp[N1.id]-Vp[N2.id])/Rp/2;
			if (Jx)
			{
				Jxp[N1.id]+=J/(N1.y2-N1.y1);
				Jxp[N2.id]+=J/(N2.y2-N2.y1);
			}
			if (Ex)
			{
				Exp[N1.id]-=M->P[N1.P].Rp*J/(N1.y2-N1.y1);
				Exp[N2.id]-=M->P[N2.P].Rp*J/(N2.y2-N2.y1);
			}			
		}
	}
}

void LocalVoltage(mesh *M, int Vai, node N, double *Ex, double *Ey, double x, double y, double *Vn, double *Vp)
{
	(*Vp)=M->res.Vn[Vai][N.id]+Ex[N.id]*(x-(N.x1+N.x2)/2)+Ey[N.id]*(y-(N.y1+N.y2)/2);
	(*Vn)=M->res.Vn[Vai][N.id+M->Nn]+Ex[N.id+M->Nn]*(x-(N.x1+N.x2)/2)+Ey[N.id+M->Nn]*(y-(N.y1+N.y2)/2);
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
	int i,j, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double *Ex, *Ey;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	Ex=malloc(2*M->Nn*sizeof(double));
	Ey=malloc(2*M->Nn*sizeof(double));
	Jfield(M, Vai, NULL, NULL, Ex, Ey);
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
			LocalVoltage(M, Vai, N, Ex, Ey, x, y, &Vn, &Vp);
			fprintf(f,"%e %e %e %e\n", x, y, Vp, Vn);
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
}

void SurfPPlot(char *fn, mesh *M, int Vai, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	int i,j, ln_y=0, ln_x=0;
	double x,y, x_step, y_step;
	double *Ex, *Ey, *Jx, *Jy;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	Ex=malloc(2*M->Nn*sizeof(double));
	Ey=malloc(2*M->Nn*sizeof(double));
	Jx=malloc(2*M->Nn*sizeof(double));
	Jy=malloc(2*M->Nn*sizeof(double));
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
			double Pp, Pn, Pj;
			double Vn, Vp, I; 
			ln_y=FindPos(*M, ln_y, x, y);
			N=*SearchNode(*M,ln_y);
			Pp=sqrt((Jx[ln_y]*Jx[ln_y]+Jy[ln_y]*Jy[ln_y])*(Ex[ln_y]*Ex[ln_y]+Ey[ln_y]*Ey[ln_y]));
			Pn=sqrt((Jx[ln_y+M->Nn]*Jx[ln_y+M->Nn]+Jy[ln_y+M->Nn]*Jy[ln_y+M->Nn])*(Ex[ln_y+M->Nn]*Ex[ln_y+M->Nn]+Ey[ln_y+M->Nn]*Ey[ln_y+M->Nn]));
			LocalVoltage(M, Vai, N, Ex, Ey, x, y, &Vn, &Vp);
			Diode(*M, N, Vp-Vn, &I, NULL);
			Pj=(Vp-Vn)*I/((N.x2-N.x1)*(N.y2-N.y1));
			fprintf(f,"%e %e %e %e %e\n", x, y,Pp, Pn, Pj);
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
	int i, j;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Nn;i++)
	{
		fprintf(f,"%e %e", M->nodes[i].x1, M->nodes[i].y1);
		for(j=0;j<M->res.Nva;j++)
			fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
		fprintf(f,"\n");	
					
		fprintf(f,"%e %e", M->nodes[i].x1, M->nodes[i].y2);
		for(j=0;j<M->res.Nva;j++)
			fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
		fprintf(f,"\n\n");	
			
			
		fprintf(f,"%e %e", M->nodes[i].x2, M->nodes[i].y1);
		for(j=0;j<M->res.Nva;j++)
			fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
		fprintf(f,"\n");	
			
					
		fprintf(f,"%e %e", M->nodes[i].x2, M->nodes[i].y2);
		for(j=0;j<M->res.Nva;j++)
			fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
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
	int i, j;
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
				fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
			fprintf(f,"\n");	
						
			fprintf(f,"%e %e", xx1, yy2);
			for(j=0;j<M->res.Nva;j++)
				fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
			fprintf(f,"\n\n");	
				
				
			fprintf(f,"%e %e", xx2, yy1);
			for(j=0;j<M->res.Nva;j++)
				fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
			fprintf(f,"\n");	
				
						
			fprintf(f,"%e %e", xx2, yy2);
			for(j=0;j<M->res.Nva;j++)
				fprintf(f," %e %e", M->res.Vn[j][i], M->res.Vn[j][M->Nn+i]);
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
	int i, j;
	FILE *f;
	if ((f=fopen(fn,"w"))==NULL)
		Error("Cannot open %s for writing\n", fn);
	for (i=0;i<M->Na;i++)
	{
		fprintf(f,"*****************Parameters for area %s\n", M->P[i].name);
		fprintf(f,"id:     %i\n", i);
		fprintf(f,"Rp:     %e\tRn:     %e\n", M->P[i].Rp, M->P[i].Rn);
		fprintf(f,"Rpvp:   %e\tRpvn:   %e\n", M->P[i].Rpvp, M->P[i].Rpvn);
		fprintf(f,"Rnvp:   %e\tRnvn:   %e\n", M->P[i].Rnvp, M->P[i].Rnvn);
		switch (M->P[i].model)
		{
			case JVD:
				fprintf(f,"V [V]\t\tJ [A/cm2]\n");
				for (j=0;j<M->P[i].N;j++)
					fprintf(f,"%e\t%e\n", M->P[i].V[j], M->P[i].J[j]);
				break;
			case ONED:
				fprintf(f,"J0:     %e\tnid:    %e\n", M->P[i].J01, M->P[i].nid1);
				fprintf(f,"Rs:     %e\tRsh:    %e\n", M->P[i].Rs, M->P[i].Rsh);
				fprintf(f,"T:      %e\tEg:     %e\n", M->P[i].T, M->P[i].Eg);
				break;
			case TWOD:
				fprintf(f,"J01:    %e\tnid1:   %e\n", M->P[i].J01, 1.0);
				fprintf(f,"J02:    %e\tnid2:   %e\n", M->P[i].J02, 2.0);
				fprintf(f,"Rs:     %e\tRsh:    %e\n", M->P[i].Rs, M->P[i].Rsh);
				fprintf(f,"T:      %e\tEg:     %e\n", M->P[i].T, M->P[i].Eg);
				break;
				
		}
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

