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
#include "main.h"
#include "list.h"
#include "utils.h"
#include "polygon.h"
#define TWOPI 6.28318530717959
#define FOURPI 12.5663706143592
#define TINY 1e-12
#define RADPERDEG 0.017453292519943

void Poly_Rotate(polygon P, double xc, double yc, double deg)
{
	int i;
	double xx, yy, d;
	d=RADPERDEG*deg;
	for (i=0;i<P.N;i++)
	{
		xx=(P.x[i]-xc);
		yy=(P.y[i]-yc);
		P.x[i]=xc+xx*cos(d)-yy*sin(d);
		P.y[i]=yc+yy*cos(d)+xx*sin(d);
	}
}


void Poly_ScaleMove(polygon P, double Fx, double Fy, double Dx, double Dy)
/* fisrt scale then move all points in a polygon */
{
	int i;
	for (i=0;i<P.N;i++)
	{
		P.x[i]=Fx*P.x[i]+Dx;
		P.y[i]=Fy*P.y[i]+Dy;
	}
}

void Poly_GetBoundingBox(polygon P,double *x1, double *y1, double *x2, double *y2)
{
	int i;
	if (P.N==0)
		return;
	*x1=P.x[1];
	*y1=P.y[1];
	*x2=P.x[1];
	*y2=P.y[1];
	for (i=1;i<P.N;i++)
	{
		if (*x1>P.x[i])
			*x1=P.x[i];
		if (*y1>P.y[i])
			*y1=P.y[i];
		if (*x2<P.x[i])
			*x2=P.x[i];
		if (*y2<P.y[i])
			*y2=P.y[i];
			
	}
}

void Poly_SetBoundingBox(polygon P, double x1, double y1, double x2, double y2, int fixR, double *Fx, double *Fy, double *Dx, double *Dy)
{
	double xx1, yy1, xx2, yy2;
	
	Poly_GetBoundingBox(P, &xx1, &yy1, &xx2, &yy2);
	*Fx=(x2-x1)/(xx2-xx1);
	*Fy=(y2-y1)/(yy2-yy1);
	
	if (fixR)
	{
		if (*Fx<*Fy)
			*Fy=*Fx;
		else
			*Fx=*Fy;
	}
	
	*Dx=(x1+x2)/2-(*Fx)*(xx1+xx2)/2;
	*Dy=(y1+y2)/2-(*Fy)*(yy1+yy2)/2;
	Poly_ScaleMove(P, *Fx, *Fy, *Dx, *Dy);
}

void Poly_FlipY(polygon P)
{
	double xx1, yy1, xx2, yy2;
	Poly_GetBoundingBox(P, &xx1, &yy1, &xx2, &yy2);	
	Poly_ScaleMove(P, 1.0, -1.0, 0, (yy1+yy2));
}

void Poly_FlipX(polygon P)
{
	double xx1, yy1, xx2, yy2;
	Poly_GetBoundingBox(P, &xx1, &yy1, &xx2, &yy2);	
	Poly_ScaleMove(P, -1.0, 1.0, (xx1+xx2), 0);
}


/* use Jordan curve theorem with a horizontal line through the point of interest*/
int IsInPolygon(polygon P, double x, double y)
{
	int i=0, c=0;
	int j;
	int F, L;
	L=-1;
	for (j=1;j<=P.BR[0];j++)
	{
		F=L+1;
		L=P.BR[j]-1;
		for (i=F;i<L;i++)
    			if ( ((P.y[i+1]>y) != (P.y[i]>y)) && (x < (P.x[i]-P.x[i+1]) * (y-P.y[i+1]) / (P.y[i]-P.y[i+1]) + P.x[i+1]) )
       				c = !c;
    		if ( ((P.y[F]>y) != (P.y[L]>y)) && (x < (P.x[L]-P.x[F]) * (y-P.y[F]) / (P.y[L]-P.y[F]) + P.x[F]) )
       			c = !c;
	}
	return c;
}

int IsNearPolygon(polygon P, double x, double y, double D, int loop)
{
	int i;
	int j;
	int F, L;
	double dx, dy, xx, yy, d;
	D=D*D;
	
	L=-1;
	for (j=1;j<=P.BR[0];j++)
	{
		F=L+1;
		L=P.BR[j]-1;
		for (i=F;i<L;i++)
		{
			/* check endpoints */
			d=(x-P.x[i+1])*(x-P.x[i+1])+(y-P.y[i+1])*(y-P.y[i+1]);
			if (d<D)
				return 1;
			d=(x-P.x[i])*(x-P.x[i])+(y-P.y[i])*(y-P.y[i]);
			if (d<D)
				return 1;
				
			/* compute distance point to line */
			/* line dx and dy values */
			dx=P.x[i+1]-P.x[i];
			dy=P.y[i+1]-P.y[i];
			/* if the line is actually a point we can ignore it as we already checked the distance to he endpoints */
			if ((dy*dy+dx*dx)>TINY)
			{
				/* Intersectipon point, perhaps not so readable but trust me, it is correct :) */
				xx=(dx*dy*y-dx*dy*P.y[i]+dx*dx*x+dy*dy*P.x[i])/(dy*dy+dx*dx);
				yy=(dy*dy*y+dx*dx*P.y[i]+dx*dy*x-dx*dy*P.x[i])/(dy*dy+dx*dx);
				
				/* perhaps not so readable but trust me, it is correct :) */
				if (((P.x[i]-xx)*(P.x[i+1]-xx)<=TINY)&&((P.y[i]-yy)*(P.y[i+1]-yy)<=TINY))
				{
					d=(x-xx)*(x-xx)+(y-yy)*(y-yy);
					if (d<D)
						return 1;
				}
			}
		}
		if (loop)
		{
			/* check endpoints */
			d=(x-P.x[F])*(x-P.x[F])+(y-P.y[F])*(y-P.y[F]);
			if (d<D)
				return 1;
			d=(x-P.x[L])*(x-P.x[L])+(y-P.y[L])*(y-P.y[L]);
			if (d<D)
				return 1;
				
			/* compute distance point to line */
			/* line dx and dy values */
			dx=P.x[F]-P.x[L];
			dy=P.y[F]-P.y[L];
			/* if the line is actually a point we can ignore it as we already checked the distance to he endpoints */
			if ((dy*dy+dx*dx)>TINY)
			{
				/* Intersectipon point, perhaps not so readable but trust me, it is correct :) */
				xx=(dx*dy*y-dx*dy*P.y[L]+dx*dx*x+dy*dy*P.x[L])/(dy*dy+dx*dx);
				yy=(dy*dy*y+dx*dx*P.y[L]+dx*dy*x-dx*dy*P.x[L])/(dy*dy+dx*dx);
				
				/* perhaps not so readable but trust me, it is correct :) */
				if (((P.x[L]-xx)*(P.x[F]-xx)<=TINY)&&((P.y[L]-yy)*(P.y[F]-yy)<=TINY))
				{
					d=(x-xx)*(x-xx)+(y-yy)*(y-yy);
					if (d<D)
						return 1;
				}
			}
		}
		
	}
	return 0;

}

/* extendable to polygons with holes if we add the possibility to track breaks, i.e. an empty line to indicate no connection between two susequent nodes
   This would allow clean holes cut out of a polygon shape. See comment at IsInPoly */
polygon ReadPoly(char *fn)
{
	polygon res;
	char c, *line;
	int k, Na=50;
	FILE *f;
	if ((f=fopen(fn,"rb"))==NULL)
		Error("Cannot open %s for reading\n", fn);
		
	line=malloc(MAXSTRLEN*sizeof(char));
    	fgets(line, MAXSTRLEN-1, f);
	res.N=0;
	res.x=malloc(Na*sizeof(double));
	res.y=malloc(Na*sizeof(double));
	res.BR=malloc(LISTBLOCK*sizeof(int));
	res.BR[0]=0;
	while(feof(f)==0)
	{
		
    		k=sscanf(line, " %c", &c);
		if((k==1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le", res.x+res.N,res.y+res.N);
			if(k==2)
			{
				res.N++;
				if (Na-1==res.N)
				{
					Na+=50;
					res.x=realloc(res.x, Na*sizeof(double));	
					res.y=realloc(res.y, Na*sizeof(double));					
				}
			}
			else if (res.N)
				res.BR=AddToList(res.BR, res.N);
			
		}
		else if (res.N)
			res.BR=AddToList(res.BR, res.N);
    		fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	res.x=realloc(res.x, (res.N+1)*sizeof(double));	
	res.y=realloc(res.y, (res.N+1)*sizeof(double));
	res.BR=AddToList(res.BR, res.N);
	fclose(f);
	return res;
}
