#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#define RADPERDEG 0.017453292519943

/* Quick 'n Dirty code to do basic transformations on the polygons I use for PVMOS. */
/* My usecase:
   I have an extension for inkscape that lets me export the coordinates of a path.
   This is cool as I can then load those in PVMOS. Thus with inkscape I have a 
   comfortable editor for my geometries. Problem is only that inkscape has more 
   or less arbitrary units. This is where this tool helps me to get the coordinates 
   in the right units for PVMOS */

typedef struct POS{
	double x;
	double y;
} POS;

typedef struct ARGS{
	double f;
	POS P0;
	POS P1;
} ARGS;
typedef enum TRANS {MOVE, ROT, FLIPY, FLIPX, SCALE, SCALEX, SCALEY, BB, BBFIXR} TRANS;

void Move(POS *P, int N, POS P0)
{
	int i;
	for (i=0;i<N;i++)
	{
		P[i].x+=P0.x;
		P[i].y+=P0.y;	
	}
}

void Rotate(POS *P, int N, POS P0, double deg)
{
	int i;
	double x, y, d;
	d=RADPERDEG*deg;
	for (i=0;i<N;i++)
	{
		x=(P[i].x-P0.x);
		y=(P[i].y-P0.y);
		P[i].x=P0.x+x*cos(d)-y*sin(d);
		P[i].y=P0.y+y*cos(d)+x*sin(d);
	}
}

void FlipY(POS *P, int N)
{
	int i;
	for (i=0;i<N;i++)
		P[i].y=-P[i].y;	
}

void FlipX(POS *P, int N)
{
	int i;
	for (i=0;i<N;i++)
		P[i].x=-P[i].x;	
}

void ScaleX(POS *P, int N, double f)
{
	int i;
	for (i=0;i<N;i++)
		P[i].x=f*P[i].x;
}
void ScaleY(POS *P, int N, double f)
{
	int i;
	for (i=0;i<N;i++)
		P[i].y=f*P[i].y;
}

void Scale(POS *P, int N, double f)
{
	int i;
	for (i=0;i<N;i++)
	{
		P[i].x=f*P[i].x;
		P[i].y=f*P[i].y;
	}
}

void BoundingBox(POS *P, int N, POS *LL, POS *UR)
{
	int i;
	LL->x=0;
	LL->y=0;
	UR->x=0;
	UR->y=0;
	if (N==0)
		return;
	LL->x=P[1].x;
	LL->y=P[1].y;
	UR->x=P[1].x;
	UR->y=P[1].y;
	for (i=1;i<N;i++)
	{
		if (LL->x>P[i].x)
			LL->x=P[i].x;
		if (LL->y>P[i].y)
			LL->y=P[i].y;
		if (UR->x<P[i].x)
			UR->x=P[i].x;
		if (UR->y<P[i].y)
			UR->y=P[i].y;
			
	}
}

void MoveBB(POS *P, int N, POS LL0, POS UR0)
{
	int i;
	POS UR;
	POS LL;
	POS S;
	double fx,fy;
	
	BoundingBox(P, N, &LL, &UR);
	fx=(UR0.x-LL0.x)/(UR.x-LL.x);
	fy=(UR0.y-LL0.y)/(UR.y-LL.y);
	S.x=(LL0.x+UR0.x)/2-fx*(LL.x+UR.x)/2;
	S.y=(LL0.y+UR0.y)/2-fy*(LL.y+UR.y)/2;
	fprintf(stderr,"scalex %.12e scaley %.12e\nmove %.12e %.12e\n", fx, fy,S.x, S.y);
	ScaleX(P, N, fx);
	ScaleY(P, N, fy);
	Move(P, N, S);
}

void MoveBBFixR(POS *P, int N, POS LL0, POS UR0)
{
	int i;
	POS UR;
	POS LL;
	POS S;
	double fx,fy;
	
	BoundingBox(P, N, &LL, &UR);
	fx=(UR0.x-LL0.x)/(UR.x-LL.x);
	fy=(UR0.y-LL0.y)/(UR.y-LL.y);
		
	if (fx<fy)
	{
		S.x=(LL0.x+UR0.x)/2-fx*(LL.x+UR.x)/2;
		S.y=(LL0.y+UR0.y)/2-fx*(LL.y+UR.y)/2;
		fprintf(stderr,"scale %.12e move %.12e %.12e\n", fx, S.x, S.y);
		Scale(P, N, fx);
	}
	else
	{
		S.x=(LL0.x+UR0.x)/2-fy*(LL.x+UR.x)/2;
		S.y=(LL0.y+UR0.y)/2-fy*(LL.y+UR.y)/2;
		fprintf(stderr,"scale %.12e move %.12e %.12e\n", fy, S.x, S.y);
		Scale(P, N, fy);
	}
	
	Move(P, N, S);
}

void TransForm(POS *P, int Np, TRANS *T, ARGS *A, int Nc)
{
	int i;
	for (i=0;i<Nc;i++)
	{
		switch (T[i])
		{
			case MOVE:
				Move(P, Np, A[i].P0);
				break;
			case ROT:
				Rotate(P, Np, A[i].P0, A[i].f);
				break;
			case FLIPY:
				FlipY(P, Np);
				break;
			case FLIPX:
				FlipX(P, Np);
				break;
			case SCALE:
				Scale(P, Np, A[i].f);
				break;	
			case SCALEX:
				ScaleX(P, Np, A[i].f);
				break;	
			case SCALEY:
				ScaleY(P, Np, A[i].f);
				break;	
			case BB:
				MoveBB(P, Np, A[i].P0, A[i].P1);
				break;	
			case BBFIXR:
				MoveBBFixR(P, Np, A[i].P0, A[i].P1);
				break;			
		}
	
	}
}

#define MAXSTRLEN 2000
POS * ReadPoly(char * fn, int **BR, int *Nbr, int *N)
{
	FILE *f;
	POS *P;
	int Na=16, Nabr=16;
	char *line, c;
	int k;
	if ((f=fopen(fn,"r"))==NULL)
	{
		fprintf(stderr, "Cannot open %s for reading\n", fn);
		exit(1);
	}
	(*Nbr)=0;
	(*N)=0;
	
	P=malloc(Na*sizeof(POS));
	(*BR)=malloc(Nabr*sizeof(int));	
	line=malloc(MAXSTRLEN*sizeof(char));
    	fgets(line, MAXSTRLEN-1, f);
	while(feof(f)==0)
	{
		
    		k=sscanf(line, " %c", &c);
		if((k!=-1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le", &P[(*N)].x,&P[(*N)].y);
			if(k!=-1)
				(*N)++;
			else
			{
				(*BR)[(*Nbr)]=(*N);
				if ((*Nbr)==0)
					(*Nbr)++;
				else if ((*BR)[(*Nbr)]!=(*BR)[(*Nbr)-1])
					(*Nbr)++;
			}
			if ((*Nbr)>Nabr-1)
			{
				Nabr*=2;
				(*BR)=realloc((*BR), Nabr*sizeof(int));
			}
			if ((*N)>Na-1)
			{
				Na*=2;
				P=realloc(P, Na*sizeof(POS));
			}
			
		}
		else
		{
			(*BR)[(*Nbr)]=(*N);
			if ((*Nbr)==0)
				(*Nbr)++;
			else if ((*BR)[(*Nbr)]!=(*BR)[(*Nbr)-1])
				(*Nbr)++;
			if ((*Nbr)>Nabr-1)
			{
				Nabr*=2;
				(*BR)=realloc((*BR), Nabr*sizeof(int));
			}
		}
		
    		fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	(*BR)=realloc((*BR), ((*Nbr)+1)*sizeof(int));
	P=realloc(P, ((*N)+1)*sizeof(POS));	
	return P;	
}

double LengthPoly(POS *P, int Np, int *BR, int Nbr)
{
	int i, j=0;
	double L=0;
	double x,y;
	for (i=1;i<Np;i++)
	{
	
		if (i==BR[j])
		{
			if (j==0)
			{
				x=P[0].x;
				y=P[0].y;
			}
			else
			{
				x=P[BR[j-1]].x;
				y=P[BR[j-1]].y;
			
			}
			L+=sqrt((P[i-1].x - x)*(P[i-1].x - x)+ (P[i-1].y - y)*(P[i-1].y - y));
			j++;	
		}
		else
			L+=sqrt((P[i].x - P[i-1].x)*(P[i].x - P[i-1].x)+ (P[i].y - P[i-1].y)*(P[i].y - P[i-1].y));
	
	}
	if (j==0)
	{
		x=P[0].x;
		y=P[0].y;
	}
	else
	{
		x=P[BR[j-1]].x;
		y=P[BR[j-1]].y;	
	}
	L+=sqrt((P[i-1].x - x)*(P[i-1].x - x)+ (P[i-1].y - y)*(P[i-1].y - y));
	return L;
}

void PrintPoly( POS *P, int Np, int *BR, int Nbr)
{
	int i, j=0;
	for (i=0;i<Np;i++)
	{
		printf("%.12e\t%.12e\n",P[i].x, P[i].y);
		if (j<Nbr)
			if (i+1==BR[j])
			{
				printf("\n");
				j++;				
			}
	
	}
	printf("\n");
}

int main(int argc, char **argv)
{
	TRANS *T;
	ARGS *A;
	int Nc=0;
	int Np=0;
	int NewT=0, i;
	POS *P;
	int *BR, Nbr;
	if (argc<2) 
	{
        	fprintf(stderr," usage: polytransform <file> [trans 1] [trans 2] ... \n");
        	return 1;
    	}
	T=malloc(argc*sizeof(TRANS));
	A=malloc(argc*sizeof(ARGS));
	i=2;
	while (i<argc)
	{
		if (NewT==0)
		{
			if (strncmp(argv[i],"move",4)==0)
			{
				T[Nc]=MOVE;
				Nc++;
				NewT++;
			}
			if (strncmp(argv[i],"rotate",6)==0)
			{
				T[Nc]=ROT;
				Nc++;
				NewT++;
			}
			else if (strncmp(argv[i],"flipx",5)==0)
			{
				T[Nc]=FLIPX;
				Nc++;
			}
			else if (strncmp(argv[i],"flipy",5)==0)
			{
				T[Nc]=FLIPY;
				Nc++;
			}
			else if (strncmp(argv[i],"scalex",6)==0)
			{
				T[Nc]=SCALEX;
				Nc++;
				NewT++;
			}
			else if (strncmp(argv[i],"scaley",6)==0)
			{
				T[Nc]=SCALEY;
				Nc++;
				NewT++;
			}
			else if (strncmp(argv[i],"scale",5)==0)
			{
				T[Nc]=SCALE;
				Nc++;
				NewT++;
			}
			else if (strncmp(argv[i],"bb_fixR",7)==0)
			{
				T[Nc]=BBFIXR;
				Nc++;
				NewT++;
			} else	if (strncmp(argv[i],"bb",2)==0)
			{
				T[Nc]=BB;
				Nc++;
				NewT++;
			}
				
		}
		else
		{
			if (T[Nc-1]==MOVE)
			{
				if (NewT==1)
				{
					A[Nc-1].P0.x=atof(argv[i]);
					NewT++;
				}
				else
				{
					A[Nc-1].P0.y=atof(argv[i]);
					NewT=0;
				}
			}
			if (T[Nc-1]==ROT)
			{
				if (NewT==1)
				{
					A[Nc-1].P0.x=atof(argv[i]);
					NewT++;
				}
				else if (NewT==2)
				{
					A[Nc-1].P0.y=atof(argv[i]);
					NewT++;
				}
				else if (NewT==3)
				{
					A[Nc-1].f=atof(argv[i]);
					NewT=0;
				}
			}
			if ((T[Nc-1]==SCALE)||(T[Nc-1]==SCALEX)||(T[Nc-1]==SCALEY))
			{
				A[Nc-1].f=atof(argv[i]);
				NewT=0;
			}
			if ((T[Nc-1]==BB)||(T[Nc-1]==BBFIXR))
			{
				if (NewT==1)
				{
					A[Nc-1].P0.x=atof(argv[i]);
					NewT++;
				}
				else if (NewT==2)
				{
					A[Nc-1].P0.y=atof(argv[i]);
					NewT++;
				}
				else if (NewT==3)
				{
					A[Nc-1].P1.x=atof(argv[i]);
					NewT++;
				}
				else
				{
					A[Nc-1].P1.y=atof(argv[i]);
					NewT=0;
				}
			}
		}
		i++;			
	}
	P=ReadPoly(argv[1], &BR, &Nbr, &Np);
	TransForm(P, Np, T, A, Nc);
	PrintPoly(P, Np, BR, Nbr);
	fprintf(stderr,"Poly Length: %e\n", LengthPoly(P, Np, BR, Nbr));
 	free(T);
 	free(A);
 	free(P);
 	free(BR);
}
