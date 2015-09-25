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
 * solve the non-linear set of kirchoff current law equations    *                 
 * for a given mesh. Also provides an adaptive meshing routine   * 
 *                                                               *            
 *****************************************************************/     
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <suitesparse/cholmod.h>

#include "mesh2d.h"
#include "list.h"
#include "main.h"
#include "solve.h"
#include "utils.h"
#include "diode.h"
#include "phototransistor.h"
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))

NumSet Numeric_Settings={
	1e-5, 1e-5, 1e-5, 1e-5,
	25, 10,
	6,0,
	1e-2,
	10,
	0.0
};                                                                          

void Resistance(mesh M, node N1, node N2, double *R)
{
	/* computes the resistance in the electrodes between two elements */
	int n,s,w,e, i;
	double L1,L2, W;
	n=IsInList(N1.north, N2.id);
	s=IsInList(N1.south, N2.id);
	w=IsInList(N1.west, N2.id);
	e=IsInList(N1.east, N2.id);
	/* compute the resistances between nodes */
	if (!(n+s+w+e))
		Error("Node %d is not connected to node %d\n", N1.id, N2.id);
	if (n||s)
	{
		/* length is diff in y coordinate, w is overlap in x coordinate */
		L1=(N1.y2-N1.y1)/2;
		L2=(N2.y2-N2.y1)/2;
		W=MIN(N1.x2, N2.x2)-MAX(N1.x1, N2.x1);
	}
	else
	{
		/* length is diff in y coordinate, w is overlap in x coordinate */
		L1=(N1.x2-N1.x1)/2;
		L2=(N2.x2-N2.x1)/2;
		W=MIN(N1.y2, N2.y2)-MAX(N1.y1, N2.y1);	
	}
	if (W<=0)
		W=(L1+L2)*1e-100;
	for (i=0;i<M.Nel;i++)
	{
		R[i]=(L1*M.P[N1.P].Rel[i]+L2*M.P[N2.P].Rel[i])/W;
	}
}

void BubbleSortJV(int n, double *V, double *J)
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
				c=J[b];
				J[b]=J[b+1];
				J[b+1]=c;
				s=0;
			}
		}
		if (s)
			break;
	}
}


#define vv M.P[P].conn[inter_index].V
#define jj M.P[P].conn[inter_index].J
void Diode_JVD(mesh M, int P, int inter_index, double V, double *J, double *dJdV)
/* returns diode current density (J) or the derivative of current versus voltage (dJdV)
   for a given area (indexed with P) and voltage (V). */
{
	int min=0, max, i;
	max=M.P[P].conn[inter_index].N-1;
	if (V<vv[min])
	{
		if (J)
			(*J)=(jj[min]+(V-vv[min])*(jj[min+1]-jj[min])/(vv[min+1]-vv[min]));
		if (dJdV)
			(*dJdV)=(jj[min+1]-jj[min])/(vv[min+1]-vv[min]);
		return;
	}
	if (V>vv[max])
	{
		if (J)
			(*J)=(jj[max]+(V-vv[max])*(jj[max]-jj[max-1])/(vv[max]-vv[max-1]));
		if (dJdV)
			(*dJdV)=(jj[max]-jj[max-1])/(vv[max]-vv[max-1]);
		return;
	}
	while(max-min>1)
	{
		i=(max+min)/2;
		if (V<vv[i])
			max=i;
		else
			min=i;
	}
	
	if (J)
		(*J)=((V-vv[min])*jj[max]+(vv[max]-V)*jj[min])/(vv[max]-vv[min]);
	if (dJdV)
		(*dJdV)=(jj[max]-jj[min])/(vv[max]-vv[min]);
}
#undef vv
#undef jj
void Diode_J(mesh M, int P, int inter_index, double V, double *J, double *dJdV, double *Vj)
{
	switch (M.P[P].conn[inter_index].model)
	{
		case JVD:
			Diode_JVD(M, P, inter_index, V, J, dJdV);
			if (Vj)
				(*Vj)=0; /* unknown junction voltage */

			break;
		case ONED:
			OneDiode(V, *((OneTwoDiode *) M.P[P].conn[inter_index].ParStruct), M.P[P].T, J, dJdV, Vj);
			break;
		case TWOD:
			TwoDiode(V, *((OneTwoDiode *) M.P[P].conn[inter_index].ParStruct), M.P[P].T, J, dJdV, Vj);
			break;
		case PHOTOT:
		{
			double v2;
			Phototransistor(V, *((PhotoTransistor *) M.P[P].conn[inter_index].ParStruct), M.P[P].T, Vj, &v2, J, dJdV);
			break;
		}
		default:
			Error("Unknown diode model in function Diode_J\n");
	}
}

#define Area (N.x2-N.x1)*(N.y2-N.y1)
void Diode(mesh M, node N, int inter_index, double V, double *I, double *dIdV, double *Vj)
{
	Diode_J(M, N.P, inter_index, V, I, dIdV, Vj);
	if (I)
		(*I)*=Area;
	if (dIdV)
		(*dIdV)*=Area;
}

int *CollectColumn(int *list, node N)
{
	int j;
	list[0]=0;
	list=AddToList(list, N.id);
	for (j=1;j<=N.north[0];j++)
		list=AddToList(list, N.north[j]);
	for (j=1;j<=N.south[0];j++)
		list=AddToList(list, N.south[j]);
	for (j=1;j<=N.east[0];j++)
		list=AddToList(list, N.east[j]);
	for (j=1;j<=N.west[0];j++)
		list=AddToList(list, N.west[j]);
	return list;
}

cholmod_sparse * SystemMatrix(mesh M, cholmod_common *c)
/* creates the sparce system matrix with the Kirchhoff Current Law equations for each node */
{
	cholmod_sparse *S;
	int n, i, j, k;
	int *ii=NULL, *pp=NULL;
	double *xx=NULL;
	int *list;
	int nnz=0, nnz_a; /* counter for the number of non zeros and the number of allocated non zeros in the matrix S */
	double *R, *diag;
	
	R=malloc((M.Nel+1)*sizeof(double));
	diag=malloc((M.Nel+1)*sizeof(double));
	
	list=malloc(LISTBLOCK*sizeof(int));
	
	n=M.Nel*M.Nn;
	nnz_a=6*n; /* the ratio between the number of nodes and nnz depends on the mesh geometry, for a regular mesh it should be a little less than 5, for irregular meshes this can be larger */
	S=cholmod_allocate_sparse(n, n, nnz_a, 1, 1, 1, CHOLMOD_REAL, c);
	
	
	pp=(int*)S->p;
	/* we scan the mesh to determine the structure of the system matrix */
	for (i=0;i<M.Nn;i++)
	{
		pp[i]=nnz;
		list=CollectColumn(list, M.nodes[i]);
			
		ii=((int*)S->i)+nnz;
		nnz+=list[0];
		if (nnz>=nnz_a-1)
		{
			nnz_a=nnz*1.1+1;
			cholmod_reallocate_sparse(nnz_a, S, c);
		}
		for (j=0;j<list[0];j++)
			ii[j]=list[j+1];
	}
	for (k=0;k<M.Nel;k++)
		pp[(k+1)*M.Nn]=(k+1)*nnz;
	free(list);
	/* at this point we know the final size as the parts for the other electrodes have the exact same structure as for the first electrode */
	/* we reallocate the space for the matrix */
	nnz_a=M.Nel*nnz;
	
	cholmod_reallocate_sparse(nnz_a+1, S, c);
	/* we now enter a loop which fills in the matrix values and structure for all electrodes. */
	
	xx=(double *)S->x;
	ii=(int*)S->i;
	for (i=0;i<M.Nel*nnz;i++)
		xx[i]=0;
	
	for (i=0;i<M.Nn;i++)
	{
		int D=0;
		for (k=0;k<M.Nel;k++)
		{
			if (k>0)
				pp[i+k*M.Nn+1]=pp[i+k*M.Nn]+pp[i+1]-pp[i];
			diag[k]=0;
		}
		for (j=pp[i];j<pp[i+1];j++)
		{
			node N;
			N=*SearchNode(M, i);
			
			for (k=1;k<M.Nel;k++)
				ii[j+k*nnz]=ii[j]+k*M.Nn;
			
			if (ii[j]!=i)
			{
				/* non-diagonal element */
				Resistance(M, N, *SearchNode(M, ii[j]), R);
				for (k=0;k<M.Nel;k++)
				{
					diag[k]+=1/R[k];
					xx[j+k*nnz]-=1/R[k];
				}
			}
			else
			{
			/* diagonal element */
				for (k=0;k<M.Nel;k++)
				{
					if (M.P[N.P].Rvp[k]>0)
						diag[k]+=Area/M.P[N.P].Rvp[k];
					if (M.P[N.P].Rvn[k]>0)
						diag[k]+=Area/M.P[N.P].Rvn[k];
					D=j;
				}
			}
		}
		for (k=0;k<M.Nel;k++)
			xx[D+k*nnz]=diag[k];
		
	}
	/*for (i=0;i<M.Nel*M.Nn;i++)
		printf("%i %i\n",i,pp[i]);
		printf("%i %i\n",i,pp[i]);*/
	/*
	for (i=0;i<M.Nel*M.Nn;i++)
		for (j=pp[i];j<pp[i+1];j++)
			printf("%i %i %e\n",i,ii[j],xx[j]);
	cholmod_print_sparse (S, "S", c);*/
	free(R);
	free(diag);
	return S;
} 



cholmod_sparse * JacobiMatrix(mesh M, double *V, double Gmin, cholmod_sparse *S, cholmod_common *c)
/* Takes the system mathix and adds the linearized (non-linear) inter-electrode connection */
{
	cholmod_sparse *J, *jj;
	int n, i;
	int *ii=NULL, *pp=NULL;
	double *xx=NULL;
	double a[2]={1,0}, b[2]={1,0};
	
	n=M.Nel*M.Nn;
	jj=cholmod_allocate_sparse(n, n, n+2*(n-M.Nn), 1, 1, 1, CHOLMOD_REAL, c);
	
	
	pp=(int*)jj->p;
	ii=(int*)jj->i;
	xx=(double *)jj->x;
	
	
			
	for (i=0;i<n+2*(n-M.Nn);i++)
		xx[i]=0;
	
	pp[0]=0;
	for (i=0;i<n;i++)
	{
		node *N;
		double dIdV;
		if (i<M.Nn)
		{
			/* first element on the diagonal, 2 elements per column/row */
			ii[pp[i]]=i;
			/* second element is Nn off diagonal */
			ii[pp[i]+1]=i+M.Nn;
			pp[i+1]=pp[i]+2;
			
		}
		else if (i<(M.Nel-1)*M.Nn)
		{
			/* first element off diagonal, 3 elements per column/row */
			ii[pp[i]]=i-M.Nn;
			ii[pp[i]+1]=i;	
			ii[pp[i]+2]=i+M.Nn;
			pp[i+1]=pp[i]+3;
			N=SearchNode(M, i%M.Nn);
			Diode(M, *N, i/M.Nn -1, V[i-M.Nn]-V[i], NULL, &dIdV, NULL);
			dIdV+=(N->x2-N->x1)*(N->y2-N->y1)*Gmin;
			/* diagonal */
			xx[pp[i]+1]+=dIdV;
			/* off diagonal */
			xx[pp[i]]-=dIdV;
			
			/* previous electrode block */
			if (i-M.Nn<M.Nn)
			{
				/* diagonal */
				xx[pp[i-M.Nn]]+=dIdV;
				/* off diagonal */
				xx[pp[i-M.Nn]+1]-=dIdV;
			}
			else
			{
				/* diagonal */
				xx[pp[i-M.Nn]+1]+=dIdV;	
				/* off diagonal */
				xx[pp[i-M.Nn]+2]-=dIdV;	
			}			
		}
		else
		{
			/* first element off diagonal, 2 elements per column/row */
			ii[pp[i]]=i-M.Nn;
			ii[pp[i]+1]=i;	
			pp[i+1]=pp[i]+2;
			
			N=SearchNode(M, i%M.Nn);
			Diode(M, *N, i/M.Nn -1, V[i-M.Nn]-V[i], NULL, &dIdV, NULL);
			dIdV+=(N->x2-N->x1)*(N->y2-N->y1)*Gmin;
			/* diagonal */
			xx[pp[i]+1]+=dIdV;
			/* off diagonal */
			xx[pp[i]]-=dIdV;
			
			/* previous electrode block */
			if (i-M.Nn<M.Nn)
			{
				/* diagonal */
				xx[pp[i-M.Nn]]+=dIdV;
				/* off diagonal */
				xx[pp[i-M.Nn]+1]-=dIdV;
			}
			else
			{
				/* diagonal */
				xx[pp[i-M.Nn]+1]+=dIdV;	
				/* off diagonal */
				xx[pp[i-M.Nn]+2]-=dIdV;	
			}			
		
		}
	}
	
	J=cholmod_add(S,jj,a,b,1,1,c);
	/*
	cholmod_print_sparse (jj, "jj", c);
	cholmod_print_sparse (S, "S", c);
	cholmod_print_sparse (J, "J", c);
	
	pp=(int*)J->p;
	ii=(int*)J->i;
	xx=(double *)J->x;
	for (i=0;i<M.Nel*M.Nn;i++)
	{
		for (j=pp[i];j<pp[i+1];j++)
			printf("%i %i %e\n", i, ii[j], xx[j]);
	}
	fflush(stdout);*/
	cholmod_free_sparse(&jj, c);
	return J;
}

void Residual(mesh M, double *V, double Va, double *I, double *E, double *Erel, cholmod_sparse *S, cholmod_dense *res, cholmod_common *c)
/* Compute Residuals:
   input:
   	M	The mesh
	V	voltages for each node
	Va	Applied voltage
	I	if not NULL the value of (*I) will be the total current afzter completion of this routine
	E	if not NULL the value of (*E) will be theabsolute RMS error
	S	system matrix
	res	allocated dense vector residuals, we sould have that S V=I(V), where I(V) are the (non-linear) current sources
		thus we compute the residula as res=S V-I(V)
	c	cholmod common
*/	
{
	int i, k;
	double Id;
	double *bb, *vv;
	double alpha[2]={1.0,0.0}, beta[2]={0.0,0.0};
	cholmod_dense *v;
	clock_t start, end;
	
     	start = clock();
	
	v=cholmod_allocate_dense(M.Nel*M.Nn,1,M.Nel*M.Nn,CHOLMOD_REAL, c);
	
	vv=(double *)v->x;
	vv=memcpy(vv, V, M.Nel*M.Nn*sizeof(double));
		
	/* b = S V */
	cholmod_sdmult (S, 0, alpha, beta, v, res, c);
	
	bb=(double *)res->x;
	if (I)
		(*I)=0;
	if (Erel)
		(*Erel)=1e-10;
	for (i=0;i<M.Nn;i++)
	{
		node N;
		N=*SearchNode(M, i);
		
		for (k=0;k<M.Nel;k++)
		{
			if (k>0)
			{
				/* connection between electrodes */
				Diode(M, N, k-1, V[i+(k-1)*M.Nn]-V[i+k*M.Nn], &Id, NULL, NULL);
				bb[i+(k-1)*M.Nn]+=Id;
				bb[i+k*M.Nn]-=Id;
				if (Erel)
					(*Erel)+=fabs(Id);
			}
			/* connection of electrodes to applied bias node */
			if (M.P[N.P].Rvp[k]>0)
			{
				bb[i+k*M.Nn]-=Va*Area/M.P[N.P].Rvp[k];
				if (I)
					(*I)+=(Va-V[i+k*M.Nn])*Area/M.P[N.P].Rvp[k];
			}
		}
	}
	if (E)
	{
		(*E)=0;
		for (i=0;i<M.Nel*M.Nn;i++)
			(*E)+=bb[i]*bb[i];
		(*E)=sqrt((*E)/((double)(M.Nel*M.Nn)));
		
	}
	if (Erel)
	{
		(*Erel)=(*Erel)/((double)((M.Nel-1)*M.Nn));
		(*Erel)=(*E)/(*Erel);
	}
	cholmod_free_dense(&v, c);
     	end = clock();
	cpu_time_rhs+=((double) (end - start)) / CLOCKS_PER_SEC;
}

#undef Area

int NewtonStep(mesh M, double *Vin, double *Vout, double Va, double *I, double *Ekcl, double *Ekcl_rel, double *Ev, int N_lin_search, double Gmin, cholmod_sparse *S, cholmod_common *c)
/* Do one newton iteration  */
{
	cholmod_sparse *J;
	cholmod_factor *F ;
	cholmod_dense *res;
	cholmod_dense *dv;
	double *dvv, a=1.0, E0;
	double *vv;
	int i=0, j, conv=0;
	clock_t start, end;
	
	res=cholmod_allocate_dense(M.Nel*M.Nn,1,M.Nel*M.Nn,CHOLMOD_REAL, c);
	
	/* compute right hand side and error in KCL (E0) */
	Residual(M, Vin, Va, I, &E0, NULL, S, res, c);
	
     	start = clock();	
	/* Determine Jacobi matrix */
	J=JacobiMatrix(M, Vin, Gmin, S, c);
     	end = clock();
	cpu_time_jacobi+=((double) (end - start)) / CLOCKS_PER_SEC;
	
     	start = clock();
	/* Solve system */
	F = cholmod_analyze( J, c) ;
	cholmod_factorize (J, F, c) ;
	if (c->nmethods >1)
	{
		Print(DEBUG, "Selected ordering: %i", c->selected);
		c->nmethods = 1 ;
		c->method [0].ordering =  c->method [c->selected].ordering;
	}
	dv = cholmod_solve ( CHOLMOD_A, F, res, c) ;
     	end = clock();
	cpu_time_cholmod+=((double) (end - start)) / CLOCKS_PER_SEC;
	
	
	if (c->status==CHOLMOD_NOT_POSDEF)
	{
		Print(NORMAL, "Please check for floating nodes");
		conv=1;
	}
	/* Process voltage step */
	/* determine magnitude of voltage step */
	dvv=(double *)dv->x;
	(*Ev)=0;
	for (j=0;j<M.Nel*M.Nn;j++)
		(*Ev)+=dvv[j]*dvv[j];
	(*Ev)=sqrt((*Ev)/((double)(M.Nel*M.Nn)));
	if (!isfinite(*Ev))
	{
		conv=1;
		Print(VERBOSE, "Voltage error not finite.");
	}
	/* compute new KCL Error */
	(*Ekcl)=E0+1;
	
	if (!Vout)
		vv=malloc(M.Nel*M.Nn*sizeof(double));
	else
		vv=Vout;
	/* if new KCL Error is larger than the old one, do a linear seach for a voltage step with a smaller KCL error */
	/* We simply divide the adaption by 2 repeatedly untill either the new KCL error is smaller than the old one or a (hard-coded) maximum number of 7 steps (a=0.00781250)*/
	for (j=0;j<M.Nel*M.Nn;j++)
		vv[j]=Vin[j];
		
	while (((*Ekcl)>E0)&&(i<N_lin_search)&&(!conv))
	{
		for (j=0;j<M.Nel*M.Nn;j++)
		{
			vv[j]=Vin[j]-a*dvv[j];
			if (isfinite(vv[j])==0)
			{
				/*do note step voltage, trigger gmin stepping */
				for (j=0;j<M.Nel*M.Nn;j++)
					vv[j]=Vin[j];
				conv=1;
				Print(VERBOSE, "Voltage vector not finite.");
					
			}
		}
		Residual(M, vv, Va, I, Ekcl, Ekcl_rel, S, res, c);
		a/=2;	
		i++;	
	}
	if (!conv&&(i>=N_lin_search)&&(N_lin_search>1))
	{
		conv=1;	
		for (j=0;j<M.Nel*M.Nn;j++)
			vv[j]=Vin[j];
		Residual(M, vv, Va, I, Ekcl, Ekcl_rel, S, res, c);
		Print(VERBOSE, "Linear search in Newton direction failed.");
	}
	if ((i>1)&&(!conv))
		Print(VERBOSE,"Step size reduced by a factor of %e",a*2);
	if (!Vout)
		free(vv);
	cholmod_free_dense(&res, c);
	cholmod_free_dense(&dv, c);
	cholmod_free_sparse(&J, c);
	cholmod_free_factor(&F, c);
	return conv;
	
}

int FindVa(double Va, double *list, int Nva)
/* returns a pointer to the closest simulated element voltages */
{
	double min;
	int i, imin;
	if (Nva==0)
		return -1;
	
	imin=0;
	min=fabs(Va-list[0]);
	for (i=Nva-1;i>0;i--)
		if (min>fabs(Va-list[i]))
		{
			min=fabs(Va-list[i]);
			imin=i;		
		}
	return imin;
}



void SolveVa(mesh *M, double Vstart, double Vend, int Nstep)
/* Do an IV sweep */
{
	cholmod_common c ;
	cholmod_sparse *S;
	double Ekcl, Ekcl_rel, Ev, Va, Vaf=1, Gmin;
	double *V, *Vout, *vv;
	int i, j, k, GminSteps=0;
	clock_t start, end;
	
	Gmin=Numeric_Settings.Gmin;
	Print(NORMAL, "________________________________________________________________");
	Print(NORMAL, "Mesh with %d layers, each with %d elements", M->Nel, M->Nn);
	cholmod_start (&c);
	c.nmethods = 2 ; /* seems that metis leads to less memory (less than 10% improvement) but it comes at a considerable speed penalty (up to a factor of 2!) */ 
	                 /* this avoids metis (if you have it enabled in your cholmod version */
	/* c.method [0].ordering =  CHOLMOD_AMD; */
	c.postorder = 1 ; 
     	start = clock();
	S=SystemMatrix(*M, &c);	
     	end = clock();
	cpu_time_system+=((double) (end - start)) / CLOCKS_PER_SEC;
	
	Print(NORMAL, "Va          iter    Ev          Ev_rel      Ekcl        Ekcl_rel");
	Print(NORMAL, "----------------------------------------------------------------");
	
	V=calloc((M->Nel*M->Nn+1),sizeof(double));
	Vout=calloc((M->Nel*M->Nn+1),sizeof(double));
	M->res.Vn=realloc(M->res.Vn, (M->res.Nva+Nstep+1)*sizeof(double **));	
	for (k=0;k<Nstep;k++)
	{
		int conv;
		if (Nstep>1)
			Va=Vstart+(double)k*(Vend-Vstart)/((double)Nstep-1.0);
		else
			Va=Vstart;
		i=FindVa(Va, M->res.Va, M->res.Nva);
		
		M->res.Va=realloc(M->res.Va, (M->res.Nva+1)*sizeof(double));
		M->res.I=realloc(M->res.I, (M->res.Nva+1)*sizeof(double));
		M->res.Va[M->res.Nva]=Va;
		M->res.Vn[M->res.Nva]=malloc((M->Nel+1)*sizeof(double *));	
		for (j=0;j<M->Nel;j++)	
			M->res.Vn[M->res.Nva][j]=calloc(M->Nn+1,sizeof(double));
		
		/* do a linear extrapolation to the new applied voltage 
		   (i.e.  scale all node voltages with a factor). We limit this
		   to a sensible range to avoid problems. */
		if (i>=0)
		{
			if (fabs(M->res.Va[i])>1e-4)
				Vaf=Va/M->res.Va[i];
			if (fabs(Vaf)>100)
				Vaf=1;
			for (j=0;j<M->Nel*M->Nn;j++)				
				V[j]=M->res.Vn[i][j/M->Nn][j%M->Nn]*Vaf;
		}
			
		M->res.Nva++;
				
		i=0;
		do
		{
			if (i<Numeric_Settings.GminStart)
				conv=NewtonStep(*M, V, Vout,Va, &(M->res.I[M->res.Nva-1]), &Ekcl, &Ekcl_rel, &Ev,Numeric_Settings.N_lin_search, Numeric_Settings.GminMax, S, &c);
			else
				conv=NewtonStep(*M, V, Vout,Va, &(M->res.I[M->res.Nva-1]), &Ekcl, &Ekcl_rel, &Ev,Numeric_Settings.N_lin_search, Gmin, S, &c);
			if ((GminSteps>1)||(i<Numeric_Settings.GminStart))
				Print(VERBOSE, "%-12.2e*%-7d%-12.2e%-12.2e%-12.2e%-8.2e",Va, i+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
			else
				Print(VERBOSE, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-8.2e",Va, i+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);

			if (conv)
			{
				/* do gmin stepping */
				if (!GminSteps)
				{
					GminSteps=Numeric_Settings.GminStep+1;
					Gmin=Numeric_Settings.GminMax*Numeric_Settings.GminFac;
				}
				else
					Gmin*=(1.5*Numeric_Settings.GminFac);
			}
			else
			{
				/* swap input and output arrays */
				vv=Vout;
				Vout=V;
				V=vv;
			}
			i++;
			if (GminSteps)
			{
				j=i;
				Gmin/=Numeric_Settings.GminFac;
				GminSteps--;
				if (GminSteps==1)
					Gmin=Numeric_Settings.Gmin;
			}
		} while ((GminSteps||(i<=Numeric_Settings.GminStart))||((i<Numeric_Settings.max_iter)&&(((Ekcl>Numeric_Settings.tol_kcl_abs)&&(Ekcl_rel>Numeric_Settings.tol_kcl_rel))||((Ev>Numeric_Settings.tol_v_abs)&&(Ev/(fabs(Va)+1e-10)>Numeric_Settings.tol_v_rel)))));
		if (verbose<VERBOSE)
			Print(NORMAL, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-8.2e",Va, i,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
	
		for (j=0;j<M->Nel*M->Nn;j++)				
			M->res.Vn[M->res.Nva-1][j/M->Nn][j%M->Nn]=V[j];	
	}
	Print(NORMAL, "----------------------------------------------------------------");
	free(V);
	free(Vout);
	cholmod_free_sparse(&S, &c);	
	cholmod_finish(&c);
}

/* returns indexes to datapoints in the IV characteristics around points of interest, namely Short circuit, Open circuit and mpp */
/* This serves either to refine those points or to determine the solar cell parameters */
/* As we can set the applied voltage to 0 we dpo not need to refine Isc, thus I return just the index of the minimal absolute voltage */
/* You may want to make sure it is 0 volts */
/* To refine Voc I need one point above and one point below, thus I return two indexes */
/* To refine Pmpp I need Pmax, one below this and one abovce it, i.e. I look for three indexes */
int SolPar(mesh *M, int *isc, int *imp_m, int *imp, int *imp_p, int *ioc_m, int *ioc_p)
{
	double Vsc;
	double Ioc_m, Ioc_p;
	
	double Vmp_m;
	double Vmp_p;
	double Vmp;
	
	double Pmax;
	
	int i;
	(*isc)=-1;
	(*imp_m)=-1;
	(*imp)=-1;
	(*imp_p)=-1;
	(*ioc_m)=-1;
	(*ioc_p)=-1;
	
	if (M->res.Nva>0)
	{
		Vsc=M->res.Va[0];
		(*isc)=0;
		if (M->res.I[0]<0)
		{
			(*ioc_m)=0;
			Ioc_m=M->res.I[0];
		}
		else
		{
			(*ioc_p)=0;
			Ioc_p=M->res.I[0];
		}
		Pmax=-(M->res.I[i]*M->res.Va[i]);
		Vmp=Vsc;
		(*imp)=0;
						
	}
	
	for (i=1;i<M->res.Nva;i++)
	{
		if (fabs(Vsc)>fabs(M->res.Va[i]))
		{
			Vsc=M->res.Va[i];
			(*isc)=i;		
		}
		
		if ((M->res.I[i]<0)&&((Ioc_m<M->res.I[i])||((*ioc_m)<0)))
		{
			(*ioc_m)=i;
			Ioc_m=M->res.I[i];			
		}
			
		if ((M->res.I[i]>=0)&&((Ioc_p>M->res.I[i])||((*ioc_p)<0)))
		{
			(*ioc_p)=i;
			Ioc_p=M->res.I[i];			
		}
		
		if (Pmax<-(M->res.I[i]*M->res.Va[i]))		
		{
			if (Vmp<M->res.Va[i])
			{
				(*imp_m)=(*imp);
				Vmp_m=Vmp;				
			}
			else
			{
				(*imp_p)=(*imp);
				Vmp_p=Vmp;				
			}
			Pmax=-(M->res.I[i]*M->res.Va[i]);
			Vmp=M->res.Va[i];
			(*imp)=i;
		} 
		else if (((*imp_p)<0)&&(Vmp>M->res.Va[i]))
		{
			(*imp_m)=i;
			Vmp_m=M->res.Va[i];		
		} 
		else if (((*imp_p)<0)&&(Vmp<M->res.Va[i]))
		{
			(*imp_p)=i;
			Vmp_p=M->res.Va[i];		
		} 
		else if ((Vmp<M->res.Va[i])&&(Vmp_p>M->res.Va[i]))
		{
			(*imp_p)=i;
			Vmp_p=M->res.Va[i];		
		} 
		else if ((Vmp>M->res.Va[i])&&(Vmp_m<M->res.Va[i]))
		{
			(*imp_m)=i;
			Vmp_m=M->res.Va[i];		
		}
			
	}
	if (M->res.Nva>0)
	{
		if (fabs(M->res.Va[(*isc)])>1e-3)
			Warning("Warning: Short circuit conditions were not simulated\n\t-->Using current at %e V instead\n", M->res.Va[(*isc)]);
	}
	if (((*isc)<0)||((*imp_m)<0)||((*imp)<0)||((*imp_p)<0)||((*ioc_m)<0)||((*ioc_p)<0))
		return 1;
	return 0;
}

#define rat_bisect 1e-2
#define sign(a) (((a)>=0)?(1.0):(-1.0))
void RefineOC(mesh *M, double tol_i, double tol_v, int Niter)
{
	double x1,x2,x3,x4;
	double f1,f2,f3,f4;
	int isc, imp_m, imp, imp_p, ioc_m, ioc_p, iter=0;
	SolPar(M, &isc, &imp_m, &imp, &imp_p, &ioc_m, &ioc_p);
	if ((ioc_m<0)||(ioc_p<0))
	{
		Warning("Cannot refine Voc, I need at least one simulated bias above and one below Voc\n");
		return;
	}
	
	/* Ridder's method */	
	
	x1=M->res.Va[ioc_m];
	f1=M->res.I[ioc_m];
	
	x2=M->res.Va[ioc_p];
	f2=M->res.I[ioc_p];
	
	while ((fabs(x2-x1)>tol_v)&&(iter<Niter))
	{
		x3=(x1+x2)/2;
		SolveVa(M, x3, x3, 1);
		f3=M->res.I[M->res.Nva-1];
		if (fabs(f3)<tol_i)
			return;
		iter++;
		if (iter==Niter)
			return;
		x4=x3+(x3-x1)*sign(f1-f2)*f3/(sqrt(f3*f3-f2*f1));
		SolveVa(M, x4, x4, 1);
		f4=M->res.I[M->res.Nva-1];
		if (fabs(f4)<tol_i)
			return;
		if (sign(f4)!=sign(f3))
		{
			x1=x3;
			f1=f3;
			x2=x4;
			f2=f4;
		}
		else if (sign(f1)!=sign(f4))
		{
			x2=x4;
			f2=f4;
		} 
		else if  (sign(f2)!=sign(f4))
		{
			x1=x4;
			f1=f4;
		} 
		else
		{
			Warning("Cannot find Voc, aborting\n");
			return;
		}
		iter++;	
	}
}
void RefineMPP(mesh *M, double tol_i, double tol_v, int Niter)
{
	double Pm,Pmin,Pmax;
	double Vm,Vmin,Vmax;
	double Im,Imin,Imax;
	int isc, imp_m, imp, imp_p, ioc_m, ioc_p, iter=0;
	SolPar(M, &isc, &imp_m, &imp, &imp_p, &ioc_m, &ioc_p);
	if ((imp_m<0)||(imp_p<0)||(imp<0))
	{
		Warning("Cannot refine maximum powerpoint, I need at least three bias points\nspanning a range around mpp\n");
		return;
	}
	
	Pmin=M->res.Va[imp_m]*M->res.I[imp_m];
	Pmax=M->res.Va[imp_p]*M->res.I[imp_p];
	Pm=M->res.Va[imp]*M->res.I[imp];
	
	Vmin=M->res.Va[imp_m];
	Vmax=M->res.Va[imp_p];
	Vm=M->res.Va[imp];
	
	Imin=M->res.I[imp_m];
	Imax=M->res.I[imp_p];
	Im=M->res.I[imp];
	
	while ((Vmax-Vmin>tol_v)&&((Im-Imin>tol_i)||(Imax-Im>tol_i))&&(iter<Niter))
	{
		double v1, v2, vv, p1, p2, pp;
		v1=(Vmin+Vm)/2;
		v2=(Vmax+Vm)/2;
		
		/* the derivative of the power is a somewhat exponential function, taking the log makes it a tad more linear and thus faster to optimize */
		p1=log(1+(Pmin-Pm)/(Vm-Vmin));
		p2=log(1+(Pmax-Pm)/(Vmax-Vm));
		
		/* attempt a sort of false position scheme to maximizing the power */
		/* i.e. use false position type of derivative calculation between Vmin Vm and between Vm and Vmax.*/
		vv=(p2*v1+p1*v2)/(p1+p2);			
		if (fabs(Vm-vv)/(Vmax-Vmin)<rat_bisect) /* we modify the false positions using some sort of half-baked Illinois algorithm */
		{
			if (Vmax-Vm>Vm-Vmin)
				vv=(0.5*p2*v1+p1*v2)/(p1+0.5*p2);
			else
				vv=(p2*v1+0.5*p1*v2)/(0.5*p1+p2);
		}
				
		
		SolveVa(M, vv, vv, 1);

		pp=vv*M->res.I[M->res.Nva-1];
		if (pp<Pm)
		{			
			if (vv<Vm)
			{
				Vmax=Vm;
				Imax=Im;
				Pmax=Pm;
			}
			else
			{
				Vmin=Vm;
				Imin=Im;
				Pmin=Pm;
			}
			
			Pm=pp;
			Vm=vv;	
			Im=M->res.I[M->res.Nva-1];	
		}
		else
		{
			if (vv<Vm)
			{
				Vmin=vv;
				Imin=Im;
				Pmin=pp;
			}
			else
		 	{
				Vmax=vv;
				Imax=Im;
				Pmax=pp;	
			}
		}
		iter++;	
	}

	
}



#define MAXRATIO 4
void AdaptMesh(mesh *M, int Vai, double rel_threshold)
{
	int i,j, k, Nno;
	double *dVx,*dVy, dvmax=0, ddv;
	node *N1;
	dVx=calloc(M->Nn, sizeof(double));
	dVy=calloc(M->Nn, sizeof(double));
	if ((M->res.Nva<=Vai)||(Vai<0))
		Error("Requested solution to adapt the mesh to is not available\n");
		
	for (i=0;i<M->Nn;i++)
	{
		N1=SearchNode(*M, i);
		for (j=1;j<=N1->north[0];j++)
		{
			if ((M->P[N1->P].SplitY) && (M->P[(SearchNode(*M, N1->north[j]))->P].SplitY))
			{
				ddv=0;
				for (k=0;k<M->Nel;k++)
					ddv+=fabs(M->res.Vn[Vai][k][i]-M->res.Vn[Vai][k][N1->north[j]]);
				dVy[i]+=ddv;
				dVy[N1->north[j]]+=ddv;
				if (dVy[i]>dvmax)
					dvmax=dVy[i];
				if (dVy[N1->north[j]]>dvmax)
					dvmax=dVy[N1->north[j]];
			}
		}
		for (j=1;j<=N1->west[0];j++)
		{
			if ((M->P[N1->P].SplitX) && (M->P[(SearchNode(*M, N1->west[j]))->P].SplitX))
			{
				ddv=0;
				for (k=0;k<M->Nel;k++)
					ddv+=fabs(M->res.Vn[Vai][k][i]-M->res.Vn[Vai][k][N1->west[j]]);
				dVx[i]+=ddv;
				dVx[N1->west[j]]+=ddv;
				if (dVx[i]>dvmax)
					dvmax=dVx[i];
				if (dVx[N1->west[j]]>dvmax)
					dvmax=dVx[N1->west[j]];
			}
		}
	}
	dvmax*=rel_threshold;
	Print(NORMAL,"Split threshold set to %e V", dvmax);
	Nno=M->Nn;
	for (i=0;i<Nno;i++)
	{
		if ((dVx[i]>dvmax)&&(dVy[i]>dvmax))
			SplitNodeXY(i, M);
		else if (dVx[i]>dvmax)
		{
			/* in general it is not a good idea to get nodes which are very long and thin,
			   as it *can*, depending on the situation, lead to considerable errors in your 
			   solution. In turn these errors often lead to a situation where nodes are being 
			   split in one direction over and over again, where the solution gets no more
			   accurate. For this reason I put a maximum to the ratio. Note that you are still 
			   free to screw up your mesh by disallowing nodes to be split in one or another 
			   direction. */
			N1=SearchNode(*M, i);
			if (((N1->y2-N1->y1)/(N1->x2-N1->x1)>MAXRATIO)&&(M->P[N1->P].SplitY))
				SplitNodeY(i, M);
			else
				SplitNodeX(i, M);
		}
		else if (dVy[i]>dvmax)
		{
			N1=SearchNode(*M, i);
			if (((N1->x2-N1->x1)/(N1->y2-N1->y1)>MAXRATIO)&&(M->P[N1->P].SplitX))
				SplitNodeX(i, M);
			else
				SplitNodeY(i, M);
		}
	}
	
	free(dVx);
	free(dVy);
}


void AdaptiveSolveVa(mesh *M, double Va, double rel_threshold, int N)
{
	int i, j;
	double Ev=0, Ei=0;
	SolveVa(M, Va, Va, 1);
	
	/* clean up old data */
	
	for (i=0;i<M->res.Nva-1;i++)
	{
		for (j=0;j<M->Nel;j++)
			free(M->res.Vn[i][j]);
		free(M->res.Vn[i]);
	}
	M->res.Vn[0]=M->res.Vn[M->res.Nva-1];
	M->res.I[0]=M->res.I[M->res.Nva-1];
	M->res.Va[0]=Va;
	
	M->res.Nva=1;
	
	for (i=0;i<N;i++)
	{
		Print(NORMAL, "Adapting mesh iteration %d",i+1);
		fflush(stdout);
		AdaptMesh(M, M->res.Nva-1, rel_threshold);
		Print(NORMAL, "Solving System");
		fflush(stdout);
		SolveVa(M, Va, Va, 1);
		Ev=0;
		for (j=0;j<M->Nel*M->Nn;j++)
			Ev+=(M->res.Vn[M->res.Nva-1][j/M->Nn][j%M->Nn]-M->res.Vn[M->res.Nva-2][j/M->Nn][j%M->Nn])*(M->res.Vn[M->res.Nva-1][j/M->Nn][j%M->Nn]-M->res.Vn[M->res.Nva-2][j/M->Nn][j%M->Nn]);
		Ev/=(double)(M->Nel*M->Nn);
		Ev=sqrt(Ev);
		Ei=(M->res.I[M->res.Nva-1]-M->res.I[M->res.Nva-2]);
		
		Print(NORMAL, "Adapt Mesh Errors: %e V and %e A",Ev, Ei);
		
		for (j=0;j<M->Nel;j++)
		{
			free(M->res.Vn[M->res.Nva-2][j]);
			M->res.Vn[M->res.Nva-2][j]=M->res.Vn[M->res.Nva-1][j];
		}
		M->res.I[M->res.Nva-2]=M->res.I[M->res.Nva-1];
		M->res.Nva--;
	}
}

#define Njv 101
/* the locally collected current is purely for lateral transport, i.e. it is assumed that the current flows out of the diode, hence it does not consider the internal series resistance */ 
/* differential mode: compute small signal current collection efficiency
   non-differential mode: large signal response, take local short circuit current and monitor the non-linear response to removing that current. This produces the locally collected currentdensity */
double *LocalyCollectedCurrent(mesh *M, double Va, int diode_index, int *nodes, int diff, int Ri)
{
	cholmod_common c ;
	cholmod_sparse *S;
	int i, j, Na_old;
	double Ev=0;
	double Ekcl, Ekcl_rel;
	double *res;
	double *V, *Vnew, *Vref, *vv;
	double Ilocal, Vstart=0, Vend=0;
	int pc, pcl=0;
	clock_t start, end;
	
	if (!nodes[0])
		return NULL;
		
	res=malloc((nodes[0]+1)*sizeof(double));
	Print(NORMAL, "Simulating the Locally Collected Current");
	Print(NORMAL, "Doing reference calculation");
	/* solva system */
	SolveVa(M, Va, Va, 1);
	V=calloc((M->Nel*M->Nn+1),sizeof(double));
	Vref=calloc((M->Nel*M->Nn+1),sizeof(double));
	Vnew=calloc((M->Nel*M->Nn+1),sizeof(double));
	
	for (j=0;j<M->Nel*M->Nn;j++)	
	{		
		Vref[j]=M->res.Vn[M->res.Nva-1][j/M->Nn][j%M->Nn];
		V[j]=Vref[j];
	}
	for (j=0;j<M->Nn;j++)	
	{	
		double v;	
		if (IsInList(nodes, M->nodes[j].id))
		{
			v=V[M->nodes[j].id+diode_index*M->Nn]-V[M->nodes[j].id+(diode_index+1)*M->Nn];		
			if (Vstart>v)
				Vstart=v;
			if (Vend<v)
				Vend=v;
		}
	}
	Vstart-=0.5;
	Vend+=0.5;
	
	/* for each area definition, create a new area with diode removed */
	Na_old=M->Na;
	for (i=0;i<Na_old;i++)
	{
		M->Na++;
		M->P=realloc(M->P, (M->Na+1)*sizeof(local_prop));
		DuplicateProperties(M, M->P+M->Na-1, M->P+i);
				
		if (!Ri)
		{
			if (M->P[Na_old+i].conn[diode_index].V)
				free(M->P[Na_old+i].conn[diode_index].V);
			if (M->P[Na_old+i].conn[diode_index].J)
				free(M->P[Na_old+i].conn[diode_index].J);
			if (diff)
			{	
				/* only need to force one operating point to be shifted, i.e. put in one current value at thge operating point */
				/* this means we need to know the operating point of each individual element in this area before we can fill in the table */				
				M->P[Na_old+i].conn[diode_index].V=calloc((3),sizeof(double));
				M->P[Na_old+i].conn[diode_index].J=calloc((3),sizeof(double));
				M->P[Na_old+i].conn[diode_index].N=2;
				M->P[Na_old+i].conn[diode_index].model=JVD;
			
			}
			else
			{
				double J;
				/* need a complete JV for the non-linear response. As I only want the lateral collection I force superposition to hold 
				   i.e. I shift the entire JV characteristics over the short circuit current */
				Diode_J(*M, i, diode_index, 0, &Ilocal,NULL,NULL);
				
				M->P[Na_old+i].conn[diode_index].V=calloc((Njv+1),sizeof(double));
				M->P[Na_old+i].conn[diode_index].J=calloc((Njv+1),sizeof(double));
				M->P[Na_old+i].conn[diode_index].N=Njv;
				M->P[Na_old+i].conn[diode_index].model=JVD;
				/* voltage range 1 V more than the maximum voltage difference found for the selected nodes */
				for (j=0;j<Njv;j++)
				{
					M->P[Na_old+i].conn[diode_index].V[j]=Vstart+j*(Vend-Vstart)/(Njv-1);
					Diode_J(*M, i, diode_index, M->P[Na_old+i].conn[diode_index].V[j], &J,NULL,NULL);
					M->P[Na_old+i].conn[diode_index].J[j]=J-Ilocal;
				}
				
			}
		}
		else
		{
			if (!diff)
			{
				switch(M->P[Na_old+i].conn[diode_index].model)
				{
					case JVD:
						{
							Diode_J(*M, i, diode_index, 0, &Ilocal,NULL,NULL);								
							for (j=0;j<M->P[Na_old+i].conn[diode_index].N;j++)
								M->P[Na_old+i].conn[diode_index].J[j]=M->P[i].conn[diode_index].J[j]-Ilocal;					
						}
						break;
					case ONED:
					case TWOD:
						((OneTwoDiode *) M->P[Na_old+i].conn[diode_index].ParStruct)->Jph=0;
						break;
					default:
						Warning("Unknown or unsupported diode model in function LocallyCollectedCurrent\n");
				}
			}
			else
			{
				if (M->P[Na_old+i].conn[diode_index].V)
					free(M->P[Na_old+i].conn[diode_index].V);
				if (M->P[Na_old+i].conn[diode_index].J)
					free(M->P[Na_old+i].conn[diode_index].J);			
				M->P[Na_old+i].conn[diode_index].V=calloc((3),sizeof(double));
				M->P[Na_old+i].conn[diode_index].J=calloc((3),sizeof(double));
				M->P[Na_old+i].conn[diode_index].N=2;
				M->P[Na_old+i].conn[diode_index].model=JVD;
			}
			
		}
		
	}
	
	
	Print(NORMAL, "Running through selected elements");
	cholmod_start (&c);
	c.nmethods = 2 ;
	c.postorder = 1 ; 
     	start = clock();
	S=SystemMatrix(*M, &c);	
     	end = clock();
	cpu_time_system+=((double) (end - start)) / CLOCKS_PER_SEC;
	
	
	if (diff)
		diff=1;
	else
		diff=0;
	if (Ri)
		Ri=1;
	else 
		Ri=0;
	
	switch (diff+(Ri<<1))
	{
		case 0:
			Print(NORMAL, "Simulating Lateral Current Density Collection");
			break;
		case 1:
			Print(NORMAL, "Simulating Differential Lateral Current Collection Efficiency");
			break;	
		case 2:
			Print(NORMAL, "Simulating Photo-Current Density Collection (including internal series resistance)");
			break;
		case 3:
			Print(NORMAL, "Simulating Differential Photo-Current Collection Efficiency (including internal series resistance)");
			break;
		default:
			Warning("Bit's are a flippin %d'\n",diff+(Ri<<1));
	}
	
	/* for each element in the mesh, move it into the new area copy, simulate the current and move it back */
	printf("\n");
	for (i=1;i<=nodes[0];i++)
	{
		double Inew, A;
		node *N;
		N=SearchNode(*M, nodes[i]);
		/* store local current */
		A=((N->x2-N->x1)*(N->y2-N->y1));
		
		if (diff)
		{
			double dIdV, v;
			v=Vref[N->id+diode_index*M->Nn]-Vref[N->id+(diode_index+1)*M->Nn];
			/* linearize and shift the operating point */
			Diode(*M, *N, diode_index, v, &Ilocal, &dIdV, NULL); /* current operating point */
			
			if (!Ri)
			{
				M->P[N->P+Na_old].conn[diode_index].V[0]=v-1;
				M->P[N->P+Na_old].conn[diode_index].J[0]=(Ilocal-fabs(Ilocal)-1000*Numeric_Settings.tol_kcl_abs-dIdV)/A;
				M->P[N->P+Na_old].conn[diode_index].V[1]=v+1;
				M->P[N->P+Na_old].conn[diode_index].J[1]=(Ilocal-fabs(Ilocal)-1000*Numeric_Settings.tol_kcl_abs+dIdV)/A;				
				Ilocal=fabs(Ilocal)+1000*Numeric_Settings.tol_kcl_abs;
			}
			else
			{
				switch(M->P[N->P].conn[diode_index].model)
				{
					case JVD:
						M->P[N->P+Na_old].conn[diode_index].V[0]=v-1;
						M->P[N->P+Na_old].conn[diode_index].J[0]=(Ilocal-fabs(Ilocal)-1000*Numeric_Settings.tol_kcl_abs-dIdV)/A;
						M->P[N->P+Na_old].conn[diode_index].V[1]=v+1;
						M->P[N->P+Na_old].conn[diode_index].J[1]=(Ilocal-fabs(Ilocal)-1000*Numeric_Settings.tol_kcl_abs+dIdV)/A;				
						Ilocal=fabs(Ilocal)+1000*Numeric_Settings.tol_kcl_abs;
					case ONED:
					case TWOD:
					{
						double rj, is, Rs;
					/* here we have to be careful, it must be a differential value. The jaccobian makes that if we take tabular data it does not make any difference whether we 
					take a small or large change in current as we have the linearized system. However, in this case we want to compute the impact of the series resistance
					this means we must truely work with small values otherwise the diode model goes wrong! as it is not included in the Jacobian so to say */
						/*
						Ilocal=MAX(M->P[N->P].conn[diode_index].Jph/20,0.001);
						Ilocal*=A;
						M->P[Na_old+N->P].conn[diode_index].Jph=M->P[N->P].conn[diode_index].Jph+Ilocal/A;
						*/
						/* for the one and two diode models we know 1/dIdV-Rs=differential resistance of parallel diodes and shunt */	
									
						Ilocal=fabs(Ilocal)+1000*Numeric_Settings.tol_kcl_abs;
						Rs=((OneTwoDiode *) M->P[Na_old+i].conn[diode_index].ParStruct)->Rs;
						rj=1/dIdV-Rs/A;
						is=-rj*dIdV*Ilocal;
						M->P[N->P+Na_old].conn[diode_index].V[0]=v-1;
						M->P[N->P+Na_old].conn[diode_index].J[0]=(is-dIdV)/A;
						M->P[N->P+Na_old].conn[diode_index].V[1]=v+1;
						M->P[N->P+Na_old].conn[diode_index].J[1]=(is+dIdV)/A;	
						
						/*Print(NORMAL,"Jph %e %e %e %e",M->P[Na_old+N->P].conn[diode_index].Jph, M->P[N->P].conn[diode_index].Jph, Ilocal, A);*/
						break;
					}
					default:
						Warning("Unknown or unsupported diode model in function LocallyCollectedCurrent\n");
				}
			}
		
		}
			
			
		for (j=0;j<M->Nel*M->Nn;j++)	
			V[j]=Vref[j];
		
		/* move node into modified area */
		N->P+=Na_old;
		
		/* simulate */
		if (diff)
		{
			/* small signal, only need one iteration with the current jaccobnian */
			NewtonStep(*M, V, Vnew,Va, &Inew, &Ekcl, &Ekcl_rel, &Ev,1, 0, S, &c);
			/* swap input and output arrays */
			vv=Vnew;
			Vnew=V;
			V=vv;
			Print(VERBOSE, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-8.2e",Va, j+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
		}
		else
		{ 
			j=0;
			do
			{
				NewtonStep(*M, V, Vnew,Va, &Inew, &Ekcl, &Ekcl_rel, &Ev,Numeric_Settings.N_lin_search, 0, S, &c);
				/* swap input and output arrays */
				vv=Vnew;
				Vnew=V;
				V=vv;
				Print(VERBOSE, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-8.2e",Va, j+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
			} while ((j<Numeric_Settings.max_iter)&&(((Ekcl>Numeric_Settings.tol_kcl_abs)&&(Ekcl_rel>Numeric_Settings.tol_kcl_rel))||((Ev>Numeric_Settings.tol_v_abs)&&(Ev/(fabs(Va)+1e-10)>Numeric_Settings.tol_v_rel))));
		}
		
		
		if (diff)
			res[i]=(M->res.I[M->res.Nva-1]-Inew)/Ilocal;
		else
			res[i]=(Inew-M->res.I[M->res.Nva-1])/A;
		/* move node back */
		N->P-=Na_old;
		pc=floor(1000*(double)i/((double)nodes[0])+0.5);
		if (pcl<pc)
		{
			printf("\r%.1f %% completed   ", 100*(double)i/((double)nodes[0]));
			pcl=pc;
		}
		fflush(stdout);
		
	}
	printf("\n");

	Print(NORMAL, "----------------------------------------------------------------");
	/* clean up modifications to the mesh */
	for (i=Na_old;i<M->Na;i++)
		FreeProperties(M->P+i, M->Nel);
	M->Na=Na_old;
	
	free(V);
	free(Vnew);
	free(Vref);
	cholmod_free_sparse(&S, &c);	
	cholmod_finish(&c);
	return res;
}
