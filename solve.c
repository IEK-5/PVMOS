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
#include "utils.h"
#include "diode.h"
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))

void Resistance(mesh M, node N1, node N2, double *Rp, double *Rn)
{
	/* computes the resistance in the back and front electrode between two nodes */
	int n,s,w,e;
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
	(*Rp)=(L1*M.P[N1.P].Rp+L2*M.P[N2.P].Rp)/W;
	(*Rn)=(L1*M.P[N1.P].Rn+L2*M.P[N2.P].Rn)/W;	
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


#define vv M.P[N.P].V
#define jj M.P[N.P].J
#define Area (N.x2-N.x1)*(N.y2-N.y1)
void Diode_JVD(mesh M, node N, double V, double *I, double *dIdV)
/* returns diode current (I) or the derivative of current versus voltage (dIdV)
   for a given node (n) and voltage (V).
   It takes the JV characteristics stored in the mesh for thie given node.
   As the JV characteritics are ordered with increasing voltage
   we try to search the requested voltage efficiently */
{
	int min=0, max, i;
	max=M.P[N.P].N-1;
	if (V<vv[min])
	{
		if (I)
			(*I)=Area*(jj[min]+(V-vv[min])*(jj[min+1]-jj[min])/(vv[min+1]-vv[min]));
		if (dIdV)
			(*dIdV)=Area*(jj[min+1]-jj[min])/(vv[min+1]-vv[min]);
		return;
	}
	if (V>vv[max])
	{
		if (I)
			(*I)=Area*(jj[max]+(V-vv[max])*(jj[max]-jj[max-1])/(vv[max]-vv[min-1]));
		if (dIdV)
			(*dIdV)=Area*(jj[max]-jj[max-1])/(vv[max]-vv[max-1]);
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
	
	if (I)
		(*I)=Area*((V-vv[min])*jj[max]+(vv[max]-V)*jj[min])/(vv[max]-vv[min]);
	if (dIdV)
		(*dIdV)=Area*(jj[max]-jj[min])/(vv[max]-vv[min]);
}
#undef vv
#undef jj

void Diode(mesh M, node N, double V, double *I, double *dIdV)
{
	switch (M.P[N.P].model)
	{
		case JVD:
			Diode_JVD(M, N, V, I, dIdV);
			break;
		case ONED:
			OneDiode(V, M.P[N.P].J01, M.P[N.P].nid1, M.P[N.P].Eg, M.P[N.P].T, M.P[N.P].Jph, M.P[N.P].Rs, M.P[N.P].Rsh, I, dIdV);
			if (I)
				(*I)*=Area;
			if (dIdV)
				(*dIdV)*=Area;
			break;
		case TWOD:
			TwoDiode(V, M.P[N.P].J01, M.P[N.P].J02, M.P[N.P].Eg, M.P[N.P].T, M.P[N.P].Jph, M.P[N.P].Rs, M.P[N.P].Rsh, I, dIdV);
			if (I)
				(*I)*=Area;
			if (dIdV)
				(*dIdV)*=Area;
			break;
		default:
			Error("Unknown diode model in function Diode\n");
	}
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
	int n, i, j;
	int *ii=NULL, *pp=NULL;
	double *xx=NULL;
	int *list;
	int nnz=0, nnz_a; /* counter for the number of non zeros and the number of allocated non zeros in the matrix S */
	
	
	list=malloc(LISTBLOCK*sizeof(int));
	
	n=2*M.Nn;
	nnz_a=6*n; /* the ratio between the number of nodes and nnz depends on the mesh geometry, for a regular mesh it should be a little less than 5, for irregular meshes this can be larger */
	S=cholmod_allocate_sparse(n, n, nnz_a, 1, 1, 1, CHOLMOD_REAL, c);
	
	
	pp=(int*)S->p;
	/* we first write the part for the positive electrode */
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
	pp[M.Nn]=nnz;
	free(list);
	/* at this point we know the final size as the part for the bottom electrode has the exact same structure as for the top electrode */
	/* we reallocate the space for the matrix */
	nnz_a=2*nnz;
	
	cholmod_reallocate_sparse(nnz_a+1, S, c);
	/* we now enter a loop which fills in the matrix values and fills in the structure of the matrix for the negative electrode part. */
	xx=(double *)S->x;
	ii=(int*)S->i;
	for (i=0;i<2*nnz;i++)
		xx[i]=0;
	
	for (i=0;i<M.Nn;i++)
	{
		int j, D;
		double Rp, Rn, diagp, diagn;
		pp[i+M.Nn+1]=pp[i+M.Nn]+pp[i+1]-pp[i];
		diagp=0;
		diagn=0;
		for (j=pp[i];j<pp[i+1];j++)
		{
			node N;
			N=*SearchNode(M, i);
			ii[j+nnz]=ii[j]+M.Nn;
			if (ii[j]!=i)
			{
				/* non-diagonal element */
				Resistance(M, N, *SearchNode(M, ii[j]), &Rp, &Rn);
				diagp+=1/Rp;
				diagn+=1/Rn;
				xx[j]-=1/Rp;
				xx[j+nnz]-=1/Rn;
			}
			else
			{
				/* diagonal element */
				if (M.P[N.P].Rpvp>0)
					diagp+=Area/M.P[N.P].Rpvp;
				if (M.P[N.P].Rpvn>0)
					diagp+=Area/M.P[N.P].Rpvn;
				if (M.P[N.P].Rnvp>0)
					diagn+=Area/M.P[N.P].Rnvp;
				if (M.P[N.P].Rnvn>0)
					diagn+=Area/M.P[N.P].Rnvn;
				D=j;
			}
		}
		xx[D]=diagp;
		xx[D+nnz]=diagn;
	}/*
	for (i=0;i<2*M.Nn;i++)
		for (j=pp[i];j<pp[i+1];j++)
			printf("%i %i %e\n",i,ii[j],xx[j]);
	cholmod_print_sparse (S, "S", c);*/
		
	return S;
} 



cholmod_sparse * JacobiMatrix(mesh M, double *V, cholmod_sparse *S, cholmod_common *c)
{
	cholmod_sparse *J, *jj;
	int n, i;
	int *ii=NULL, *pp=NULL;
	double *xx=NULL;
	double a[2]={1,0}, b[2]={1,0};
	
	n=2*M.Nn;
	jj=cholmod_allocate_sparse(n, n, 4*M.Nn, 1, 1, 1, CHOLMOD_REAL, c);
	
	
	pp=(int*)jj->p;
	ii=(int*)jj->i;
	xx=(double *)jj->x;
	
	pp[0]=0;
	pp[M.Nn]=2*M.Nn;	
	for (i=0;i<M.Nn;i++)
	{
		double dIdV;	
		ii[2*i]=i;
		ii[2*i+1]=i+M.Nn;
		
		ii[2*(i+M.Nn)]=i;
		ii[2*(i+M.Nn)+1]=i+M.Nn;
		
		pp[i+1]=pp[i]+2;
		pp[i+M.Nn+1]=pp[i+M.Nn]+2;
		
		Diode(M, *SearchNode(M, i), V[i]-V[i+M.Nn], NULL, &dIdV);
		/* diagonal elements */
		xx[2*i]=dIdV;	
		xx[2*(i+M.Nn)+1]=dIdV;
		
		/* off-diagobals */
		xx[2*i+1]=-dIdV;
		xx[2*(i+M.Nn)]=-dIdV;		
	}/*
	pp=(int*)jj->p;
	ii=(int*)jj->i;
	xx=(double *)jj->x;
	for (i=0;i<2*M.Nn;i++)
	{
		for (j=pp[i];j<pp[i+1];j++)
			printf("%i %i %e\n", i, ii[j], xx[j]);
	}*/
	
	J=cholmod_add(S,jj,a,b,1,1,c);
	/*cholmod_print_sparse (J, "J", c);*/
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
	int i;
	double Id;
	double *bb, *vv;
	double alpha[2]={1.0,0.0}, beta[2]={0.0,0.0};
	cholmod_dense *v;
	
	v=cholmod_allocate_dense(2*M.Nn,1,2*M.Nn,CHOLMOD_REAL, c);
	
	vv=(double *)v->x;
	vv=memcpy(vv, V, 2*M.Nn*sizeof(double));
		
	/* b = S V */
	cholmod_sdmult (S, 0, alpha, beta, v, res, c);
	
	bb=(double *)res->x;
	if (*I)
		(*I)=0;
	if (Erel)
		(*Erel)=1e-10;
	for (i=0;i<M.Nn;i++)
	{
		node N;
		N=*SearchNode(M, i);
		Diode(M, N, V[i]-V[i+M.Nn], &Id, NULL);
		bb[i]+=Id;
		bb[i+M.Nn]-=Id;
		
		if (Erel)
			(*Erel)+=fabs(Id);
			
		if (M.P[N.P].Rpvp>0)
		{
			bb[i]-=Va*Area/M.P[N.P].Rpvp;
			if (I)
				(*I)+=(Va-V[i])*Area/M.P[N.P].Rpvp;
		}
		if (M.P[N.P].Rnvp>0)
		{
			bb[i+M.Nn]-=Va*Area/M.P[N.P].Rnvp;
			if (I)
				(*I)+=(Va-V[i+M.Nn])*Area/M.P[N.P].Rnvp;
		}
	}
	if (E)
	{
		(*E)=0;
		for (i=0;i<2*M.Nn;i++)
			(*E)+=bb[i]*bb[i];
		(*E)=sqrt((*E)/((double)M.Nn*2));
		
	}
	if (Erel)
	{
		(*Erel)=(*Erel)/((double)M.Nn);
		(*Erel)=(*E)/(*Erel);
	}
	cholmod_free_dense(&v, c);
}

#undef Area

void NewtonStep(mesh M, double *Vin, double *Vout, double Va, double *I, double *Ekcl, double *Ekcl_rel, double *Ev, cholmod_sparse *S, cholmod_common *c)
{
	cholmod_sparse *J;
	cholmod_factor *F ;
	cholmod_dense *res;
	cholmod_dense *dv;
	double *dvv, a=1.0, E0;
	double *vv;
	int i=0, j;
	clock_t start, end;
	
	res=cholmod_allocate_dense(2*M.Nn,1,2*M.Nn,CHOLMOD_REAL, c);
	
     	start = clock();
	Residual(M, Vin, Va, I, &E0, NULL, S, res, c);
     	end = clock();
	cpu_time_rhs+=((double) (end - start)) / CLOCKS_PER_SEC;
	
     	start = clock();
	J=JacobiMatrix(M, Vin, S, c);
     	end = clock();
	cpu_time_jacobi+=((double) (end - start)) / CLOCKS_PER_SEC;
	
     	start = clock();
	F = cholmod_analyze( J, c) ;
	cholmod_factorize (J, F, c) ;
	if (c->nmethods >1)
	{
		Print(DEBUG, "Selected ordering: %i\n", c->selected);
		c->nmethods = 1 ;
		c->method [0].ordering =  c->method [c->selected].ordering;
	}
	dv = cholmod_solve ( CHOLMOD_A, F, res, c) ;
     	end = clock();
	cpu_time_cholmod+=((double) (end - start)) / CLOCKS_PER_SEC;
	dvv=(double *)dv->x;
	(*Ev)=0;
	for (j=0;j<2*M.Nn;j++)
		(*Ev)+=dvv[j]*dvv[j];
	(*Ev)=sqrt((*Ev)/((double)M.Nn*2.0));
	
	(*Ekcl)=E0+1;
	
	if (!Vout)
		vv=malloc(2*M.Nn*sizeof(double));
	else
		vv=Vout;
	while (((*Ekcl)>E0)&&(i<8))
	{
		for (j=0;j<2*M.Nn;j++)
			vv[j]=Vin[j]-a*dvv[j];
		Residual(M, vv, Va, I, Ekcl, Ekcl_rel, S, res, c);
		a/=2;	
		i++;	
	}
	if (i>1)
		Print(NORMAL,"Step size reduced by a factor of %e\n",a*2);
	if (!Vout)
		free(vv);
	cholmod_free_dense(&res, c);
	cholmod_free_dense(&dv, c);
	cholmod_free_sparse(&J, c);
	cholmod_free_factor(&F, c);
	
}

int FindVa(double Va, double *list, int Nva)
/* returns a pointer to the closest simulated node voltages */
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



void SolveVa(mesh *M, double Vstart, double Vend, int Nstep, double tol_kcl_abs, double tol_kcl_rel, double tol_v_abs, double tol_v_rel, int max_iter)
{
	cholmod_common c ;
	cholmod_sparse *S;
	double Ekcl, Ekcl_rel, Ev, Va, Vaf=1;
	int i, j, k;
	clock_t start, end;
	Print(NORMAL, "________________________________________________________________\n");
	Print(NORMAL, "Solving using a Mesh with %d elements\n", M->Nn);
	cholmod_start (&c);
	c.nmethods = 2 ; /* seems that metis leads to less memory (less than 10% improvement) but it comes at a considerable speed penalty (up to a factor of 2!) */ 
	                 /* this avoids metis */
	/* c.method [0].ordering =  CHOLMOD_AMD; */
	c.postorder = 1 ; 
     	start = clock();
	S=SystemMatrix(*M, &c);	
     	end = clock();
	cpu_time_system+=((double) (end - start)) / CLOCKS_PER_SEC;
	
	Print(NORMAL, "Va          iter    Ev          Ev_rel      Ekcl        Ekcl_rel\n");
	Print(NORMAL, "----------------------------------------------------------------\n");
	
	for (k=0;k<Nstep;k++)
	{
		if (Nstep>1)
			Va=Vstart+(double)k*(Vend-Vstart)/((double)Nstep-1.0);
		else
			Va=Vstart;
		i=FindVa(Va, M->res.Va, M->res.Nva);
		
		M->res.Va=realloc(M->res.Va, (M->res.Nva+1)*sizeof(double));
		M->res.I=realloc(M->res.I, (M->res.Nva+1)*sizeof(double));
		M->res.Va[M->res.Nva]=Va;
		M->res.Vn=realloc(M->res.Vn, (M->res.Nva+1)*sizeof(double *));		
		M->res.Vn[M->res.Nva]=calloc(2*M->Nn,sizeof(double));
		
		/* do a linear extrapolation to the new applied voltage 
		   (i.e.  scale all node voltages with a factor). We limit this
		   to a sensible range to avoid problems. */
		if (fabs(M->res.Va[i])>1e-4)
			Vaf=Va/M->res.Va[i];
		if (fabs(Vaf)>100)
			Vaf=1;
		if (i>=0)
			for (j=0;j<2*M->Nn;j++)				
				M->res.Vn[M->res.Nva][j]=M->res.Vn[i][j]*Vaf;
		M->res.Nva++;
				
		i=0;
		do
		{
			NewtonStep(*M, M->res.Vn[M->res.Nva-1], M->res.Vn[M->res.Nva-1],Va, &(M->res.I[M->res.Nva-1]), &Ekcl, &Ekcl_rel, &Ev,S, &c);
			Print(VERBOSE, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-12.2e\n",Va, i+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
			i++;
		} while ((i<max_iter)&&(((Ekcl>tol_kcl_abs)&&(Ekcl_rel>tol_kcl_rel))||((Ev>tol_v_abs)&&(Ev/(fabs(Va)+1e-10)>tol_v_rel))));
		if (verbose<VERBOSE)
			Print(NORMAL, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-12.2e\n",Va, i,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
	}
	Print(NORMAL, "----------------------------------------------------------------\n");
	cholmod_free_sparse(&S, &c);	
	cholmod_finish(&c);
}
#define MAXRATIO 4
void AdaptMesh(mesh *M, int Vai, double rel_threshold)
{
	int i,j, Nno;
	double *dVx,*dVy, dvmax=0, ddv, *V;
	node *N1;
	dVx=calloc(M->Nn, sizeof(double));
	dVy=calloc(M->Nn, sizeof(double));
	if ((M->res.Nva<=Vai)||(Vai<0))
		Error("Requested solution to adapt the mesh to is not available\n");
	V=M->res.Vn[Vai];
	for (i=0;i<M->Nn;i++)
	{
		N1=SearchNode(*M, i);
		for (j=1;j<=N1->north[0];j++)
		{
			if ((M->P[N1->P].SplitY) && (M->P[(SearchNode(*M, N1->north[j]))->P].SplitY))
			{
				ddv=fabs(V[i]-V[N1->north[j]])+fabs(V[i+M->Nn]-V[N1->north[j]+M->Nn]);
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
				ddv=fabs(V[i]-V[N1->west[j]])+fabs(V[i+M->Nn]-V[N1->west[j]+M->Nn]);
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
	Print(NORMAL,"Split threshold set to %e V\n", dvmax);
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
			   split in one direction over and over again, where the solution only gets less
			   accurate. My solution is to put a maximum to the ratio. Note that you are still 
			   free to screw up your mesh by disallowing nodes to be split in one or another 
			   direction. Any powerful tool should alow you to screw up, right? */
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


double AdaptiveSolveVa(mesh *M, double Va, double rel_threshold, int N, double tol_kcl_abs, double tol_kcl_rel, double tol_v_abs, double tol_v_rel, int max_iter)
{
	int i, j;
	double E;
	SolveVa(M, Va, Va, 1, tol_kcl_abs, tol_kcl_rel, tol_v_abs, tol_v_rel, max_iter);
	
	/* clean up old data */
	for (i=0;i<M->res.Nva-1;i++)
	{
		free(M->res.Vn[i]);
	}
	M->res.Vn[0]=M->res.Vn[M->res.Nva-1];
	M->res.Nva=1;
	
	for (i=0;i<N;i++)
	{
		Print(NORMAL, "Adapting mesh iteration %d\n",i+1);
		fflush(stdout);
		AdaptMesh(M, M->res.Nva-1, rel_threshold);
		Print(NORMAL, "Solving System\n");
		fflush(stdout);
		SolveVa(M, Va, Va, 1, tol_kcl_abs, tol_kcl_rel, tol_v_abs, tol_v_rel, max_iter);
		E=0;
		for (j=0;j<2*M->Nn;j++)
			E+=(M->res.Vn[M->res.Nva-1][j]-M->res.Vn[M->res.Nva-2][j])*(M->res.Vn[M->res.Nva-1][j]-M->res.Vn[M->res.Nva-2][j]);
		E/=(double)2*M->Nn;
		E=sqrt(E);
		Print(NORMAL, "Adapt Mesh Error: %e\n\n",E);
		free(M->res.Vn[M->res.Nva-2]);
		M->res.Vn[M->res.Nva-2]=M->res.Vn[M->res.Nva-1];
		M->res.Nva--;
	}
	return E;
}

/*
double *CollectionEfficiency(mesh *M, int *list, double Va, double dJph, double tol_kcl_abs, double tol_kcl_rel, double tol_v_abs, double tol_v_rel, int max_iter)
{
	cholmod_common c ;
	cholmod_sparse *S;
	double Ekcl, Ekcl_rel, Ev, Va, I, *V;
	int i, j, Na;
	

	Print(NORMAL, "________________________________________________________________\n");
	Print(NORMAL, "Solving operating point using a Mesh with %d elements\n", M->Nn);
	cholmod_start (&c);
	S=SystemMatrix(*M, &c);	
	
	Print(NORMAL, "Va          iter    Ev          Ev_rel      Ekcl        Ekcl_rel\n");
	Print(NORMAL, "----------------------------------------------------------------\n");
	
	i=FindVa(Va, M->res.Va, M->res.Nva);
		
	M->res.Va=realloc(M->res.Va, (M->res.Nva+1)*sizeof(double));
	M->res.I=realloc(M->res.I, (M->res.Nva+1)*sizeof(double));
	M->res.Va[M->res.Nva]=Va;
	M->res.Vn=realloc(M->res.Vn, (M->res.Nva+1)*sizeof(double *));		
	M->res.Vn[M->res.Nva]=calloc(2*M->Nn,sizeof(double));
		
	if (i>=0)
		for (j=0;j<2*M->Nn;j++)
			M->res.Vn[M->res.Nva][j]=M->res.Vn[i][j];
	M->res.Nva++;
			
	i=0;
	do
	{
		NewtonStep(*M, M->res.Vn[M->res.Nva-1], M->res.Vn[M->res.Nva-1],Va, &(M->res.I[M->res.Nva-1]), &Ekcl, &Ekcl_rel, &Ev,S, &c);
		Print(VERBOSE, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-12.2e\n",Va, i+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
		i++;
	} while ((i<max_iter)&&(((Ekcl>tol_kcl_abs)&&(Ekcl_rel>tol_kcl_rel))||((Ev>tol_v_abs)&&(Ev/(fabs(Va)+1e-10)>tol_v_rel))));
	if (verbose<VERBOSE)
		Print(NORMAL, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-12.2e\n",Va, i,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
	
	Print(NORMAL, "----------------------------------------------------------------\n");



	Print(NORMAL, "Computing collection efficiency\n");
	
	M->P=realloc(M->P, (2*M->Na+1)*sizeof(local_prop));
	for (i=0;i<M->Na;i++)
		DuplicateProperties(M->P+i+M->Na, M->P+i);
	Na=M->Na;
	M->Na*=2;
	
	// change photocurrent here
	
	V=malloc((2*M->Nn+1)*sizeof(double));
	if (lost[0]==0)
		for (i=0;i<M->Nn;i++)
		{
			M->nodes[i].P+=Na;
			
			do
			{
				NewtonStep(*M, M->res.Vn[M->res.Nva-1], V, Va, &I, &Ekcl, &Ekcl_rel, &Ev,S, &c);
				Print(VERBOSE, "%-12.2e%-8d%-12.2e%-12.2e%-12.2e%-12.2e\n",Va, i+1,Ev, Ev/(fabs(Va)+1e-10), Ekcl, Ekcl_rel);
				i++;
			} while ((i<max_iter)&&(((Ekcl>tol_kcl_abs)&&(Ekcl_rel>tol_kcl_rel))||((Ev>tol_v_abs)&&(Ev/(fabs(Va)+1e-10)>tol_v_rel))));
			
			M->nodes[i].P-=Na;
		}
	else
		for (i=1;i<=list[0];i++)
		{
			node *N;
			N=SearchNode(M, list[i]);
			N->P+=Na;
			
			N->P-=Na;
		}
		
	cholmod_free_sparse(&S, &c);	
	cholmod_finish(&c);
	
}
*/

