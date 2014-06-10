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
 * Main routine:                                                 * 
 * Sort out command line option and call appropriate functions   *              
 *                                                               *            
 *****************************************************************/     

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mesh2d.h"
#include "main.h"
#include "parse.h"
#include "utils.h"
#include <unistd.h>
#ifndef __MINGW32__ 
	#include <sys/resource.h>
#endif
/* #define DEB */
#ifdef DEB
#include <mcheck.h>
#endif 
VERB_LEVEL verbose=NORMAL;
int fixverb=0;
double cpu_time=0;
double cpu_time_cholmod=0;
double cpu_time_system=0;
double cpu_time_jacobi=0;
double cpu_time_rhs=0;
double cpu_time_rest=0;
int peak_mem;

int main(int argc, char **argv)
{
	int f=1;
	clock_t start, end;
#ifdef DEB
	mtrace();
#endif 
	if(argc==2||argc==3)
	{
		#ifndef __MINGW32__ 
		struct rusage rusage;
		#endif 
		if (argc==3)
		{
			fixverb=1;
			if (strncmp(argv[1],"-q",2)==0)
				verbose=QUIET;
			else if (strncmp(argv[1],"-v",2)==0)
				verbose=VERBOSE;
			else if (strncmp(argv[1],"-db",3)==0)
				verbose=DEBUG;
			else if (strncmp(argv[1],"-n",2)==0)
				verbose=NORMAL;
			else
			{
				PrintHeader();
				fprintf(stderr,"USAGE:\n%s [verbose-level] <inputfile>\n", argv[0]);
				fprintf(stderr,"optional verbose-level argument can be:\n");
				fprintf(stderr,"       -q    -   quiet\n");
				fprintf(stderr,"       -n    -   normal\n");
				fprintf(stderr,"       -v    -   verbose\n");
				fprintf(stderr,"       -db   -   debug\n");
				return 1;
			}
			f=2;				
		}
		PrintHeader();
     		Print(NORMAL, "* Parsing file \"%s\"\n",argv[f]);
     		start = clock();
		Parse (argv[f]);
     		end = clock();
		cpu_time=((double) (end - start)) / CLOCKS_PER_SEC;
		cpu_time_rest=cpu_time-cpu_time_rhs-cpu_time_jacobi-cpu_time_system-cpu_time_cholmod;
		#ifndef __MINGW32__ 
		getrusage( RUSAGE_SELF, &rusage );
		#endif 
     		Print(NORMAL, "________________________________________________________________\n");
     		Print(NORMAL, "Resource usage:\n");
     		Print(NORMAL, "----------------------------------------------------------------\n");
		#ifndef __MINGW32__ 
     		Print(NORMAL, "Maximum resident set size:                         %10.0f MB\n",(double)rusage.ru_maxrss/1000.0);
		#endif 
		Print(NORMAL, "Total CPU time                                     %10.2f  s\n",cpu_time);	
     		Print(NORMAL, "CPU time summary                            time (s)    time (%%)\n");
     		Print(NORMAL, "----------------------------------------------------------------\n");
     		Print(NORMAL, "cholmod_solve                               %-10.2f       %3.0f\n",cpu_time_cholmod, 100.0*cpu_time_cholmod/cpu_time);
     		Print(NORMAL, "build system                                %-10.2f       %3.0f\n",cpu_time_system, 100.0*cpu_time_system/cpu_time);
     		Print(NORMAL, "build Jacobi                                %-10.2f       %3.0f\n",cpu_time_jacobi, 100.0*cpu_time_jacobi/cpu_time);
     		Print(NORMAL, "build RHS                                   %-10.2f       %3.0f\n",cpu_time_rhs, 100.0*cpu_time_rhs/cpu_time);	
     		Print(NORMAL, "rest (mesh operations, parsing, etc.)       %-10.2f       %3.0f\n",cpu_time_rest, 100.0*cpu_time_rest/cpu_time);	
     		Print(NORMAL, "----------------------------------------------------------------\n");
	}
	else
	{
		PrintHeader();
		fprintf(stderr,"USAGE:\n%s [verbose-level] <inputfile>\n", argv[0]);
		fprintf(stderr,"optional verbose-level argument can be:\n");
		fprintf(stderr,"       -q        -   quiet\n");
		fprintf(stderr,"       -n        -   normal\n");
		fprintf(stderr,"       -v        -   verbose\n");
		fprintf(stderr,"       -db       -   debug\n");
		Disclamer();
		return 1;
	}
#ifdef DEB
	muntrace();
#endif 
	return 0;	
}
