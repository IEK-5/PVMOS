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
 * routines to select nodes in a mesh                            *
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
#include "select_nodes.h"
#define TWOPI 6.28318530717959
#define TINY 1e-12

int FindPos(mesh M, int id, double x, double y)
/* find the node within which the coordinate x,y falls */
{
	int n_id=0, i, *list;
	double mind=-1;
	node *N;
	
	list=malloc(LISTBLOCK*sizeof(int));
	list[0]=0;
	while(!IsInList(list, id))
	{
newid:		list=AddToList(list,id);
		N=SearchNode(M, id);
		Print(DEBUG,"%e %e %i\n",(N->x1+N->x2)/2, (N->y1+N->y2)/2, N->id);
		if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY)&&(y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
		{
			free(list);
			return N->id;
		}
			
		if ((mind<0)||(mind>((N->x1+N->x2)/2-x)*((N->x1+N->x2)/2-x)+((N->y1+N->y2)/2-y)*((N->y1+N->y2)/2-y)))
		{
			mind=((N->x1+N->x2)/2-x)*((N->x1+N->x2)/2-x)+((N->y1+N->y2)/2-y)*((N->y1+N->y2)/2-y);
			n_id=id;
		}
		
		if (N->x2<x)
			for (i=1;i<=N->east[0];i++)
				if (!IsInList(list, N->east[i]))
				{
					id=N->east[i];
					goto newid;
				}		
		if (N->x1>x)
			for (i=1;i<=N->west[0];i++)
				if (!IsInList(list, N->west[i]))
				{
					id=N->west[i];
					goto newid;
				}
		if (N->y2<y)
			for (i=1;i<=N->north[0];i++)
				if (!IsInList(list, N->north[i]))
				{
					id=N->north[i];
					goto newid;
				}		
		if (N->y1>y)
			for (i=1;i<=N->south[0];i++)
				if (!IsInList(list, N->south[i]))
				{
					id=N->south[i];
					goto newid;
				}
	}
	free(list);
	/* fallback option, check all nodes 
	for (id=1;id<M.Nn;id++)
	{
		N=SearchNode(M, id);
		d=((N->x1+N->x2)/2-x)*((N->x1+N->x2)/2-x)+((N->y1+N->y2)/2-y)*((N->y1+N->y2)/2-y);
		if (d<mind)
		{
			mind=d;
			n_id=id;
		}
		if ((x>=N->x1)&&(x<=N->x2)&&(y>=N->y1)&&(y<=N->y2))
			return N->id;
	} */
	Warning("Warning: No node at coordinate (%e,%e)\n", x, y);
	return n_id;
}

double Angle (double y, double x)
{
	double a;
	a=atan2(y,x);
	if (a<0)
		a+=TWOPI;
	return a;
}
	
double dAngle(double x1, double x2, double y1, double y2, double x, double y)
{
	/* change in angle over a line segmentn*/
	double da; 
	
	double a1=Angle(y1-y, x1-x);
	double a2=Angle(y2-y,x2-x);
	
	da=a2-a1;
	/* if the angle is larger than pi rad, choose shortest angular change */
	if (da<-TWOPI/2)
		da+=TWOPI;
	
	if (da>TWOPI/2)
		da-=TWOPI;
	
	return da;
}


int IsInPolygon(polygon P, double x, double y)
{
	int i;
	
	/* if the point is enclosed by the Contour the total change in angle should be 2*pi */
	double sa=0;
	
	for (i=0;i<P.N-1;i++)
		sa+=dAngle(P.x[i], P.x[i+1],P.y[i], P.y[i+1],x, y);
	sa+=dAngle(P.x[P.N-1], P.x[0],P.y[P.N-1], P.y[0],x, y);	
	if (fabs(sa)<TWOPI-1e-6)
		return 0;
	else
		return 1;
}

int * PolySelectNodes(polygon P, mesh M, int *sel_nodes) 
{
	int i;
	
	if (sel_nodes[0]>0)
	{
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));		
	}
	for (i=0;i<M.Nn;i++)
		if (IsInPolygon(P, (M.nodes[i].x1+M.nodes[i].x2)/2, (M.nodes[i].y1+M.nodes[i].y2)/2))
			sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
	if (sel_nodes[0]==0)
		Warning("No nodes selected in PolySelectNodes\n");
	return sel_nodes;		
}

int * CircSelectNodes(double x, double y, double r, mesh M, int *sel_nodes)
{
	int i;
	double xn, yn;
	if (sel_nodes[0]>0)
	{
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));		
	}
	r=r*r;
	for (i=0;i<M.Nn;i++)
	{
		xn=x-(M.nodes[i].x1+M.nodes[i].x2)/2;
		yn=y-(M.nodes[i].y1+M.nodes[i].y2)/2;
		if (xn*xn+yn*yn<=r)
			sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in CircSelectNodes\n");	
	return sel_nodes;		
}
int * RectSelectNodes(double x1, double y1, double x2, double y2, mesh M, int *sel_nodes)
{
	int i;
	double xn, yn;
	if (sel_nodes[0]>0)
	{
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));		
	}
	for (i=0;i<M.Nn;i++)
	{
		xn=(M.nodes[i].x1+M.nodes[i].x2)/2;
		yn=(M.nodes[i].y1+M.nodes[i].y2)/2;
		if ((x1<=xn)&&(x2>xn)&&(y1<=yn)&&(y2>yn))
			sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
	}	
	if (sel_nodes[0]==0)
		Warning("No nodes selected in RectSelectNodes\n");
	return sel_nodes;		
}

polygon ReadPoly(char *fn)
{
	polygon res;
	char c, *line;
	int k, Na=50;
	FILE *f;
	if ((f=fopen(fn,"r"))==NULL)
		Error("Cannot open %s for reading\n", fn);
		
	line=malloc(MAXSTRLEN*sizeof(char));
    	fgets(line, MAXSTRLEN-1, f);
	res.N=0;
	res.x=malloc(Na*sizeof(double));
	res.y=malloc(Na*sizeof(double));
	while(feof(f)==0)
	{
		
    		k=sscanf(line, " %c", &c);
		if((k!=-1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le", res.x+res.N,res.y+res.N);
			if(k!=-1)
			{
				res.N++;
				if (Na-1==res.N)
				{
					Na+=50;
					res.x=realloc(res.x, Na*sizeof(double));	
					res.y=realloc(res.y, Na*sizeof(double));					
				}
			}
			
		}
    		fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	res.x=realloc(res.x, (res.N+1)*sizeof(double));	
	res.y=realloc(res.y, (res.N+1)*sizeof(double));
	fclose(f);
	return res;
}
