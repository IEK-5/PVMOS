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
#include "polygon.h"
#define TWOPI 6.28318530717959
#define FOURPI 12.5663706143592
#define TINY 1e-12

node *NextNodeNE(mesh M, node *N, double a, double b)
{
	/* we are in node N and follow the line a*x+b, in NE direction */
	/* this function returns the next node */
	double x,y;
	node *N2;
	int i;
	if (a<0)
		Error("Cannot follow line y=%e x + %e in north-east direction\n", a,b);
	x=(N->y2-b)/a;
	if ((x>=N->x1)&&(x<=N->x2))
		for (i=1;i<=N->north[0];i++)
		{
			N2=SearchNode(M, N->north[i]);
			if ((x>=N2->x1)&&(x<=N2->x2))
				return N2;			
		}
	y=a*N->x2+b;
	if ((y>=N->y1)&&(y<=N->y2))
		for (i=1;i<=N->east[0];i++)
		{
			N2=SearchNode(M, N->east[i]);
			if ((y>=N2->y1)&&(y<=N2->y2))
				return N2;			
		}
	
	return NULL;
}

node *NextNodeSW(mesh M, node *N, double a, double b)
{
	/* we are in node N and follow the line a*x+b, in SW direction */
	/* this function returns the next node */
	double x,y;
	node *N2;
	int i;
	if (a<0)
		Error("Cannot follow line y=%e x + %e in north-east direction\n", a,b);
	x=(N->y1-b)/a;
	if ((x>=N->x1)&&(x<=N->x2))
		for (i=1;i<=N->south[0];i++)
		{
			N2=SearchNode(M, N->south[i]);
			if ((x>=N2->x1)&&(x<=N2->x2))
				return N2;			
		}
	y=a*N->x1+b;
	if ((y>=N->y1)&&(y<=N->y2))
		for (i=1;i<=N->west[0];i++)
		{
			N2=SearchNode(M, N->west[i]);
			if ((y>=N2->y1)&&(y<=N2->y2))
				return N2;			
		}
	
	return NULL;
}

node *NextNodeNW(mesh M, node *N, double a, double b)
{
	/* we are in node N and follow the line a*x+b, in NW direction */
	/* this function returns the next node */
	double x,y;
	node *N2;
	int i;
	if (a>0)
		Error("Cannot follow line y=%e x + %e in north-west direction\n", a,b);
	x=(N->y2-b)/a;
	if ((x>=N->x1)&&(x<=N->x2))
		for (i=1;i<=N->north[0];i++)
		{
			N2=SearchNode(M, N->north[i]);
			if ((x>=N2->x1)&&(x<=N2->x2))
				return N2;			
		}
	y=a*N->x1+b;
	if ((y>=N->y1)&&(y<=N->y2))
		for (i=1;i<=N->west[0];i++)
		{
			N2=SearchNode(M, N->west[i]);
			if ((y>=N2->y1)&&(y<=N2->y2))
				return N2;			
		}
	
	return NULL;
}

node *NextNodeSE(mesh M, node *N, double a, double b)
{
	/* we are in node N and follow the line a*x+b, in SE direction */
	/* this function returns the next node */
	double x,y;
	node *N2;
	int i;
	if (a>0)
		Error("Cannot follow line y=%e x + %e in north-west direction\n", a,b);
	x=(N->y1-b)/a;
	if ((x>=N->x1)&&(x<=N->x2))
		for (i=1;i<=N->south[0];i++)
		{
			N2=SearchNode(M, N->south[i]);
			if ((x>=N2->x1)&&(x<=N2->x2))
				return N2;			
		}
	y=a*N->x2+b;
	if ((y>=N->y1)&&(y<=N->y2))
		for (i=1;i<=N->east[0];i++)
		{
			N2=SearchNode(M, N->east[i]);
			if ((y>=N2->y1)&&(y<=N2->y2))
				return N2;			
		}
	
	return NULL;
}


int FindPos(mesh M, int id, double x, double y, int *NOTINMESH)
/* find the node within which the coordinate x,y falls */
{
	double a,b, xx, yy;
	int n_id;
	node *N;
	*NOTINMESH=0;
	N=SearchNode(M, id);
	n_id=N->id;
	/* parameterize the line from current node center to desired point */
	xx=(N->x1+N->x2)/2;
	yy=(N->y1+N->y2)/2;
	b=(x-xx);
	if (fabs(b)<TINY)
	{
		if (b<0)
			b=-TINY;
		else
			b=TINY;
	}
	a=(y-yy)/b;
	b=y-a*x;
	
	/* In the next loop we try to walk in a straight line from the current node to the desired coordinate */
	
	while(N!=NULL)
	{
		if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY)&&(y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
			return N->id;
		if ((x-xx>=0)&&(y-yy>=0))
			N=NextNodeNE(M, N, a, b);			
		else if ((x-xx<0)&&(y-yy>=0))
			N=NextNodeNW(M, N, a, b);			
		else if ((x-xx>=0)&&(y-yy<=0))
			N=NextNodeSE(M, N, a, b);			
		else if ((x-xx<0)&&(y-yy<=0))
			N=NextNodeSW(M, N, a, b);
		else
			N=NULL;	
	}
	/* That did not work */
	/* It could be that the coordinate is not within the mesh or the shape of the mesh is such that there is no straight line */
	/* between start and end which connects the two. Here we do a brute force search */
	/* fallback option, check all nodes */
	for (id=0;id<M.Nn;id++)
	{
		N=SearchNode(M, id);
		if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY)&&(y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
			return N->id;
	}
	/* position is not in the mesh */
	*NOTINMESH=1;
	return n_id;
}



#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))
inline int Overlap(double a1, double a2, double b1, double b2)
{
	if (a2<a1)
	{
		double d;
		d=a1;
		a1=a2;
		a2=d;
	} 
	if (b2<b1)
	{
		double d;
		d=b1;
		b1=b2;
		b2=d;
	}
	return ((MIN(a2,b2)-MAX(a1,b1))>-TINY);
}


int PolygonCrossNode(polygon P, node *N, int loop)
{
	int i;
	int j;
	int F, L;
	double a, b, x0, y0, t1, t2, x, y;
	L=-1;
	for (j=1;j<=P.BR[0];j++)
	{
		F=L+1;
		L=P.BR[j]-1;
		for (i=F;i<L;i++)
		{
			if (Overlap(P.x[i], P.x[i+1], N->x1, N->x2)&&Overlap(P.y[i], P.y[i+1],N->y1, N->y2))
			{
				/* check endpoints */
				if ((P.x[i+1]-N->x1>=-TINY)&&(N->x2-P.x[i+1]>=-TINY)&&(P.y[i+1]-N->y1>=-TINY)&&(N->y2-P.y[i+1]>=-TINY))
					return 1;
				if ((P.x[i]-N->x1>=-TINY)&&(N->x2-P.x[i]>=-TINY)&&(P.y[i]-N->y1>=-TINY)&&(N->y2-P.y[i]>=-TINY))
					return 1;
				
				
				/* parameterize line */		
				a=(P.x[i+1]-P.x[i]);
				b=(P.y[i+1]-P.y[i]);
				x0=P.x[i];
				y0=P.y[i];
				
				t1=(N->y1-y0)/b;
				if ((t1>=0)&&(t1<=1.0))
				{
					x=a*t1+x0;
					if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY))
						return 1;
				}
				t2=(N->y2-y0)/b;
				if ((t2>=0)&&(t2<=1.0))
				{
					x=a*t2+x0;
					if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY))
						return 1;
				}	
				
				t1=(N->x1-x0)/a;
				if ((t1>=0)&&(t1<=1.0))
				{
					y=b*t1+y0;
					if ((y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
						return 1;
				}
				t2=(N->x2-x0)/a;
			
				if ((t2>=0)&&(t2<=1.0))
				{
					y=b*t2+y0;
					if ((y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
						return 1;
				}
			}
			
		}
		if (loop)
		{
			if (Overlap(P.x[F], P.x[L],N->x1, N->x2)&&Overlap(P.y[F], P.y[L], N->y1, N->y2))
			{
				/* check endpoints */
				if ((P.x[F]-N->x1>=-TINY)&&(N->x2-P.x[F]>=-TINY)&&(P.y[F]-N->y1>=-TINY)&&(N->y2-P.y[F]>=-TINY))
					return 1;
				if ((P.x[L]-N->x1>=-TINY)&&(N->x2-P.x[L]>=-TINY)&&(P.y[L]-N->y1>=-TINY)&&(N->y2-P.y[L]>=-TINY))
					return 1;
				/* parameterize line */		
				a=(P.x[F]-P.x[L]);
				b=(P.y[F]-P.y[L]);
				x0=P.x[L];
				y0=P.y[L];
				
				t1=(N->y1-y0)/b;
				if ((t1>=0)&&(t1<=1.0))
				{
					x=a*t1+x0;
					if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY))
						return 1;
				}
				t2=(N->y2-y0)/b;
				if ((t2>=0)&&(t2<=1.0))
				{
					x=a*t2+x0;
					if ((x-N->x1>=-TINY)&&(N->x2-x>=-TINY))
						return 1;
				}
	
				t1=(N->x1-x0)/a;
				if ((t1>=0)&&(t1<=1.0))
				{
					y=b*t1+y0;
					if ((y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
						return 1;
				}
				t2=(N->x2-x0)/a;
				if ((t2>=0)&&(t2<=1.0))
				{
					y=b*t2+y0;
					if ((y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
						return 1;
				}
			}
		}
	}
	return 0;
}

void ResolvSplit(polygon P, mesh *M, int id, int loop, double Pxmax, double Pymax, double Pxmin, double Pymin, double D)
{
	node *N;
	/* we try to avoid running the PolygonCrossNode routine as it is expensive */
	N=SearchNode(*M, id);
	/* resolution criterion fulfilled */
	if ((N->x2-N->x1<D)&&(N->y2-N->y1<D))
		return;
	/* possibility of crossing, overlap with entire polygon?  */
	if (Overlap(Pxmin, Pxmax, N->x1, N->x2)&&Overlap(Pymin, Pymax, N->y1, N->y2))
		if (PolygonCrossNode(P, N, loop))
		{
			/* recurse  */
			if ((N->x2-N->x1)>(N->y2-N->y1))
			{
				SplitNodeX(id, M);
				ResolvSplit(P, M, M->Nn-1, loop, Pxmax, Pymax, Pxmin, Pymin, D);
				ResolvSplit(P, M, id, loop, Pxmax, Pymax, Pxmin, Pymin, D);
			} 
			else
			{
				SplitNodeY(id, M);
				ResolvSplit(P, M, M->Nn-1, loop, Pxmax, Pymax, Pxmin, Pymin, D);
				ResolvSplit(P, M, id, loop, Pxmax, Pymax, Pxmin, Pymin, D);
			}
		}
}

void ResolvContour(polygon P, mesh *M, int loop, double D)
{
	int k, NN;
	double Pxmax, Pymax, Pxmin, Pymin;
	if (!P.N)
		return;
		
	NN=M->Nn;
	
	Pxmax=P.x[0];
	Pxmin=P.x[0];
	Pymax=P.y[0];
	Pymin=P.y[0];
	for (k=1;k<P.N;k++)
	{
		Pxmax=MAX(Pxmax,P.x[k]);
		Pymax=MAX(Pymax,P.y[k]);
		Pxmin=MIN(Pxmin,P.x[k]);
		Pymin=MIN(Pymin,P.y[k]);
	}
	for (k=0;k<NN;k++)
		ResolvSplit(P, M, M->nodes[k].id, loop, Pxmax, Pymax, Pxmin, Pymin, D);
	
}


int * PolySelectNodes(polygon P, mesh M, int *sel_nodes) 
/* selects nodes within a polygon */
{
	int i, *old_sel;
	
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
		{
			node *N;
			N=SearchNode(M, old_sel[i]);			
			if (IsInPolygon(P, (N->x1+N->x2)/2, (N->y1+N->y2)/2))
				sel_nodes=AddToList(sel_nodes, N->id);
		}		
		free(old_sel);			
	}
	else
	{
		/* make new selection */
		for (i=0;i<M.Nn;i++)
		{
			Print(DEBUG, "node %i of %i\n",i, M.Nn);
			if (IsInPolygon(P, (M.nodes[i].x1+M.nodes[i].x2)/2, (M.nodes[i].y1+M.nodes[i].y2)/2))
				sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
		}
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in PolySelectNodes\n");
	return sel_nodes;		
}

int * PolyContourSelectNodes(double d, polygon P, int loop, mesh M, int *sel_nodes) 
/* selects nodes around the contour described by a polygon */
{
	int i, *old_sel;
	
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
		{
			node *N;
			N=SearchNode(M, old_sel[i]);			
			if (IsNearPolygon(P, (N->x1+N->x2)/2, (N->y1+N->y2)/2, d,loop))
				sel_nodes=AddToList(sel_nodes, N->id);
		}		
		free(old_sel);			
	}
	else
	{
		/* make new selection */
		for (i=0;i<M.Nn;i++)
			if (IsNearPolygon(P, (M.nodes[i].x1+M.nodes[i].x2)/2, (M.nodes[i].y1+M.nodes[i].y2)/2, d, loop))
				sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in PolyContourSelectNodes\n");
	return sel_nodes;		
}

int * CircSelectNodes(double x, double y, double r, mesh M, int *sel_nodes)
{
	int i, *old_sel;
	double xn, yn;
	r=r*r;
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
		{
			node *N;
			N=SearchNode(M, old_sel[i]);
			xn=x-(N->x1+N->x2)/2;
			yn=y-(N->y1+N->y2)/2;
			if (xn*xn+yn*yn<=r)
				sel_nodes=AddToList(sel_nodes, N->id);
		}
		
		free(old_sel);	
	}
	else
	{
		/* make new selection */
		for (i=0;i<M.Nn;i++)
		{
			xn=x-(M.nodes[i].x1+M.nodes[i].x2)/2;
			yn=y-(M.nodes[i].y1+M.nodes[i].y2)/2;
			if (xn*xn+yn*yn<=r)
				sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
		}
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in CircSelectNodes\n");	
	return sel_nodes;		
}

int CircThroughNode(node *N, double x, double y, double r2)
{
	double t, a, b, c, d;
	
	/* check line (x1,y1) to (x2,y1) */
	a=(N->x2-N->x1);
	b=N->x1-x;
	c=N->y1-y;
	c=c*c;
	if (r2-c>0)
	{
		d=sqrt(r2-c);
		t=-(d+b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;
		t=(d-b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;		
	}
	/* check line (x1,y2) to (x2,y2) */
	c=N->y2-y;
	c=c*c;
	if (r2-c>0)
	{
		d=sqrt(r2-c);
		t=-(d+b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;
		t=(d-b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;		
	}
	/* check line (x1,y1) to (x1,y2) */
	a=(N->y2-N->y1);
	b=N->y1-y;
	c=N->x1-x;
	c=c*c;
	if (r2-c>0)
	{
		d=sqrt(r2-c);
		t=-(d+b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;
		t=(d-b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;		
	}
	/* check line (x2,y1) to (x2,y2) */
	c=N->x2-x;
	c=c*c;
	if (r2-c>0)
	{
		d=sqrt(r2-c);
		t=-(d+b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;
		t=(d-b)/a;
		if ((t>=0)&&(t<=1.0))
			return 1;		
	}
	/* check if the circle is contained entirely by the node */
	r2=x+sqrt(r2);
	if ((r2-N->x1>=-TINY)&&(N->x2-r2>=-TINY)&&(y-N->y1>=-TINY)&&(N->y2-y>=-TINY))
		return 1;
	return 0;
	
}


void ResolvCircSplit(double x, double y, double r, mesh *M, int id, double D)
{
	node *N;
	/* we try to avoid running the PolygonCrossNode routine as it is expensive */
	N=SearchNode(*M, id);
	/* resolution criterion fulfilled */
	if ((N->x2-N->x1<D)&&(N->y2-N->y1<D))
		return;
	/* possibility of crossing, overlap with entire polygon?  */
	if (Overlap(x-r, x+r, N->x1, N->x2)&&Overlap(y-r, y+r, N->y1, N->y2))
		if (CircThroughNode(N, x, y, r*r))
		{
			/* recurse  */
			if ((N->x2-N->x1)>(N->y2-N->y1))
			{
				SplitNodeX(id, M);
				ResolvCircSplit(x, y, r, M, M->Nn-1, D);
				ResolvCircSplit(x, y, r, M, id, D);
			} 
			else
			{
				SplitNodeY(id, M);
				ResolvCircSplit(x, y, r, M, M->Nn-1, D);
				ResolvCircSplit(x, y, r, M, id, D);
			}
		}
}

void ResolvCircle(double x, double y, double r, mesh *M, double D)
{
	int k, NN;		
	NN=M->Nn;
	for (k=0;k<NN;k++)
		ResolvCircSplit(x,y,r, M, M->nodes[k].id, D);
	
}


int * CircContourSelectNodes(double x, double y, double r, double d, mesh M, int *sel_nodes)
{
	int i, *old_sel;
	double xn, yn, r1, r2;
	r1=(r-d)*(r-d);
	r2=(r+d)*(r+d);
	
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
		{
			node *N;
			N=SearchNode(M, old_sel[i]);
			xn=x-(N->x1+N->x2)/2;
			yn=y-(N->y1+N->y2)/2;
			if ((xn*xn+yn*yn<=r2)&&(xn*xn+yn*yn>=r1))
				sel_nodes=AddToList(sel_nodes, N->id);
		}
		
		free(old_sel);	
	}
	else
	{
		/* make new selection */
		for (i=0;i<M.Nn;i++)
		{
			xn=x-(M.nodes[i].x1+M.nodes[i].x2)/2;
			yn=y-(M.nodes[i].y1+M.nodes[i].y2)/2;
			if ((xn*xn+yn*yn<=r2)&&(xn*xn+yn*yn>=r1))
				sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
		}
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in CircContourSelectNodes\n");	
	return sel_nodes;		
}

int * RectSelectNodes(double x1, double y1, double x2, double y2, mesh M, int *sel_nodes)
{
	int i, *old_sel;
	double xn, yn;
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
		{
			node *N;
			N=SearchNode(M, old_sel[i]);
			xn=(N->x1+N->x2)/2;
			yn=(N->y1+N->y2)/2;
			if ((x1<=xn)&&(x2>xn)&&(y1<=yn)&&(y2>yn))
				sel_nodes=AddToList(sel_nodes, N->id);
		}
		free(old_sel);						
	}
	else
	{
		/* make new selection */
		for (i=0;i<M.Nn;i++)
		{
			xn=(M.nodes[i].x1+M.nodes[i].x2)/2;
			yn=(M.nodes[i].y1+M.nodes[i].y2)/2;
			if ((x1<=xn)&&(x2>xn)&&(y1<=yn)&&(y2>yn))
				sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
		}
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in RectSelectNodes\n");
	return sel_nodes;		
}

int * RectContourSelectNodes(double x1, double y1, double x2, double y2, double d, mesh M, int *sel_nodes)
{
	int i, *old_sel;
	double xn, yn;
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
		{
			node *N;
			N=SearchNode(M, old_sel[i]);
			xn=(N->x1+N->x2)/2;
			yn=(N->y1+N->y2)/2;
			if ((yn>=y1-d)&&(yn<=y2+d))
			{
				/* select around vertical edges */
				if ((xn>=x1-d)&&(xn<=x1+d))
					sel_nodes=AddToList(sel_nodes, N->id);					
				if ((xn>=x2-d)&&(xn<=x2+d))
					sel_nodes=AddToList(sel_nodes, N->id);
			}
			if ((xn>=x1-d)&&(xn<=x2+d))
			{
				/* select around horizontal edges */
				if ((yn>=y1-d)&&(yn<=y1+d))
					sel_nodes=AddToList(sel_nodes, N->id);					
				if ((yn>=y2-d)&&(yn<=y2+d))
					sel_nodes=AddToList(sel_nodes, N->id);
			}
		}
		free(old_sel);						
	}
	else
	{
		/* make new selection */
		for (i=0;i<M.Nn;i++)
		{
			xn=(M.nodes[i].x1+M.nodes[i].x2)/2;
			yn=(M.nodes[i].y1+M.nodes[i].y2)/2;
			if ((yn>=y1-d)&&(yn<=y2+d))
			{
				/* select around vertical edges */
				if ((xn>=x1-d)&&(xn<=x1+d))
					sel_nodes=AddToList(sel_nodes, M.nodes[i].id);					
				if ((xn>=x2-d)&&(xn<=x2+d))
					sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
			}
			if ((xn>=x1-d)&&(xn<=x2+d))
			{
				/* select around horizontal edges */
				if ((yn>=y1-d)&&(yn<=y1+d))
					sel_nodes=AddToList(sel_nodes, M.nodes[i].id);					
				if ((yn>=y2-d)&&(yn<=y2+d))
					sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
			}
		}
	}
	if (sel_nodes[0]==0)
		Warning("No nodes selected in RectContourSelectNodes\n");
	return sel_nodes;		
}

void ResolvRect(double x1, double y1, double x2, double y2, mesh *M, double D)
{
	polygon P;
	P.N=4;
	P.x=malloc(5*sizeof(double));
	P.y=malloc(5*sizeof(double));
	P.BR=malloc(3*sizeof(int));
	P.BR[0]=1;
	P.BR[1]=4;
	P.x[0]=x1;
	P.x[1]=x2;
	P.x[2]=x2;
	P.x[3]=x1;
	P.y[0]=y1;
	P.y[1]=y1;
	P.y[2]=y2;
	P.y[3]=y2;
	
	ResolvContour(P, M, 1, D);

	free(P.x);
	free(P.y);
	free(P.BR);

}


int * SelectArea(mesh M, int *sel_nodes, char *name)
{
	int i, *old_sel, P;
	
	P=FindProperties(M, name);
	if (P<0)
		Error("Area %s does not exist\n", name);
	
	if (sel_nodes[0]>0)
	{
		/* make selection within selection */
		old_sel=DuplicateList(sel_nodes);
		sel_nodes[0]=0;
		sel_nodes=realloc(sel_nodes,LISTBLOCK*sizeof(int));
		for (i=1;i<=old_sel[0];i++)
			if (SearchNode(M, old_sel[i])->P==P)
				sel_nodes=AddToList(sel_nodes, old_sel[i]);
		free(old_sel);				
	}
	else
	{
		for (i=0;i<M.Nn;i++)
			if (M.nodes[i].P==P)
				sel_nodes=AddToList(sel_nodes, M.nodes[i].id);
	}	
	if (sel_nodes[0]==0)
		Warning("No nodes selected in RectSelectNodes\n");
	return sel_nodes;		
}
