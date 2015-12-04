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
 * Main routine:                                                 * 
 * Sort out command line option and call appropriate functions   *              
 *                                                               *            
 *****************************************************************/     

#ifdef __MINGW32__ 
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mesh2d.h"
#include "main.h"
#include "parse.h"
#include "utils.h"
#include "expr.h"
#include <unistd.h>
#ifndef __MINGW32__ 
	#include <sys/resource.h>
#endif
/* #define DEB */
#ifdef DEB
#include <mcheck.h>
#endif 
VERB_LEVEL verbose=NORMAL;

typedef enum{VAR, VARFILE} STATE; 

int fixverb=0;
double cpu_time=0;
double cpu_time_cholmod=0;
double cpu_time_system=0;
double cpu_time_jacobi=0;
double cpu_time_rhs=0;
double cpu_time_rest=0;
int peak_mem;
void printenv(char *name)
{
	printf ("%s = \"%s\"\n", name, getenv (name) ? getenv (name) : "<null>");
}

/* if we are linked to OpenBLAS we better disable multithreading */
#ifdef OPENBLAS
extern void openblas_set_num_threads(int num_threads);
#endif


int main(int argc, char **argv)
{
	int i;
	clock_t start, end;
	STATE arg=VAR;
#ifndef __MINGW32__ 
	struct rusage rusage;
#endif 	

#ifdef DEB
	mtrace();
#endif 
#ifdef __MINGW32__ 
	putenv("PRINTF_EXPONENT_DIGITS=2");
#endif

/* if we are linked to OpenBLAS we better disable multithreading */
#ifdef OPENBLAS
	openblas_set_num_threads(1);
#endif
	/* printenv("OMP_NUM_THREADS");*/
	
	InitExprEval();
	if (argc<2)
	{
		PrintHeader();
		fprintf(stderr,"USAGE:\n%s [options/variable] <inputfile>\n", argv[0]);
		fprintf(stderr,"optional argument can be:\n");
		fprintf(stderr,"       -q   			-   quiet\n");
		fprintf(stderr,"       -n   			-   normal\n");
		fprintf(stderr,"       -v   			-   verbose\n");
		fprintf(stderr,"       -db  			-   debug\n");
		fprintf(stderr,"       -vf  			-   next argument is a variable file\n");
		fprintf(stderr,"       <variable>=<value>   	-   define a variable\n");
		return 1;
	}
	
	for(i=1;i<argc-1;i++)
	{
    		switch (arg)
		{
			case VAR: /* either read in a variable or decite the argument is not a variable */
				if (strncmp(argv[i], "-vf", 3)==0) /* this argument is no variable, next argument is a variable file */
					arg=VARFILE;
				else if (strncmp(argv[i], "-q", 2)==0) /* this argument is no variable */
					verbose=QUIET;
				else if (strncmp(argv[i], "-v", 2)==0) /* this argument is no variable */
					verbose=VERBOSE;
				else if (strncmp(argv[i], "-db", 3)==0) /* this argument is no variable */
					verbose=DEBUG;
				else if (strncmp(argv[i], "-n", 3)==0) /* this argument is no variable */
					verbose=NORMAL;
				else /* this argument is a variable */
				{
					char *name, *value;
					GetNameValue(argv[i], &(name), &(value));
					if (name && value)
						DefineVar(name, atof(value));
					else
					{
						PrintHeader();
						fprintf(stderr,"USAGE:\n%s [options/variable] <inputfile>\n", argv[0]);
						fprintf(stderr,"optional argument can be:\n");
						fprintf(stderr,"       -q   			-   quiet\n");
						fprintf(stderr,"       -n   			-   normal\n");
						fprintf(stderr,"       -v   			-   verbose\n");
						fprintf(stderr,"       -db  			-   debug\n");
						fprintf(stderr,"       -vf  			-   next argument is a variable file\n");
						fprintf(stderr,"       <variable>=<value>   	-   define a variable\n");
						return 1;
					}
				}
				break;
			case VARFILE:
			{
				char *name, *value, *line; 
				FILE *varfile;
			
				printf("Opening Variable File  : %s\n", argv[i]);
                       		if ((varfile=fopen(argv[i],"r"))==NULL)
				{	
					fprintf(stderr,"Couldn't open variable file %s\n",argv[i]);
					exit(1);
				}
				while(feof(varfile)==0)
				{
				
					line=GetLine(varfile);
					GetNameValue(line, &(name), &(value));
					if (name && value)
						DefineVar(name, atof(value));
					free(line);	
				}
				fclose(varfile);
				arg=VAR;
			}
			break;
			default:
				PrintHeader();
				fprintf(stderr,"USAGE:\n%s [options/variable] <inputfile>\n", argv[0]);
				fprintf(stderr,"optional argument can be:\n");
				fprintf(stderr,"       -q   			-   quiet\n");
				fprintf(stderr,"       -n   			-   normal\n");
				fprintf(stderr,"       -v   			-   verbose\n");
				fprintf(stderr,"       -db  			-   debug\n");
				fprintf(stderr,"       -vf  			-   next argument is a variable file\n");
				fprintf(stderr,"       <variable>=<value>   	-   define a variable\n");
				return 1;
		}
	}
	/* last argument is the input file */
	if( access(argv[argc-1], F_OK ) == -1 ) 
	{
		PrintHeader();
     		Print(NORMAL, "* Input file \"%s\" does not exist",argv[argc-1]);
		fprintf(stderr,"USAGE:\n%s [options/variable] <inputfile>\n", argv[0]);
		fprintf(stderr,"optional argument can be:\n");
		fprintf(stderr,"       -q   			-   quiet\n");
		fprintf(stderr,"       -n   			-   normal\n");
		fprintf(stderr,"       -v   			-   verbose\n");
		fprintf(stderr,"       -db  			-   debug\n");
		fprintf(stderr,"       -vf  			-   next argument is a variable file\n");
		fprintf(stderr,"       <variable>=<value>   	-   define a variable\n");
		DestroyExprEval();
		return 1;	
	}	
	PrintHeader();
     	Print(NORMAL, "* Parsing file \"%s\"",argv[argc-1]);
	if (var_count)
	{
     		Print(NORMAL, "* The following variables are defined:");
		for(i=0;i<var_count;i++)
     			Print(NORMAL, "* %s\t=\t%e", var_names[i], var_values[i]);
	}
     	start = clock();
	Parse (argv[argc-1]);
     	end = clock();
	cpu_time=((double) (end - start)) / CLOCKS_PER_SEC;
	cpu_time_rest=cpu_time-cpu_time_rhs-cpu_time_jacobi-cpu_time_system-cpu_time_cholmod;
	#ifndef __MINGW32__ 
	getrusage( RUSAGE_SELF, &rusage );
	#endif 
     	Print(NORMAL, "________________________________________________________________");
     	Print(NORMAL, "Resource usage:\n");
     	Print(NORMAL, "----------------------------------------------------------------");
	#ifndef __MINGW32__ 
     	Print(NORMAL, "Maximum resident set size:                         %10.0f MB",(double)rusage.ru_maxrss/1000.0);
	#endif 
	Print(NORMAL, "Total CPU time                                     %10.2f  s",cpu_time);	
     	Print(NORMAL, "CPU time summary                            time (s)    time (%%)");
     	Print(NORMAL, "----------------------------------------------------------------");
     	Print(NORMAL, "cholmod_solve                               %-10.2f       %3.0f",cpu_time_cholmod, 100.0*cpu_time_cholmod/cpu_time);
     	Print(NORMAL, "build system                                %-10.2f       %3.0f",cpu_time_system, 100.0*cpu_time_system/cpu_time);
     	Print(NORMAL, "build Jacobi                                %-10.2f       %3.0f",cpu_time_jacobi, 100.0*cpu_time_jacobi/cpu_time);
     	Print(NORMAL, "build RHS                                   %-10.2f       %3.0f",cpu_time_rhs, 100.0*cpu_time_rhs/cpu_time);	
     	Print(NORMAL, "rest (mesh operations, parsing, etc.)       %-10.2f       %3.0f",cpu_time_rest, 100.0*cpu_time_rest/cpu_time);	
     	Print(NORMAL, "----------------------------------------------------------------");
#ifdef DEB
	muntrace();
#endif 
	DestroyExprEval();
	return 0;	
}
