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
 * some small utilitiy routines, mostly regarding printing       *
 * output with respect to the set verbosity level                *         
 *                                                               *            
 *****************************************************************/     
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "main.h"
#include "mesh2d.h"
#include "utils.h"

void Print(VERB_LEVEL v_level, const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	if (v_level<=verbose)
		vprintf(format_str, ap);
}

void Error( const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	vfprintf(stderr,format_str, ap); 
	exit(1);
}


void Warning( const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	vfprintf(stderr,format_str, ap); 
}

int ProgressBar(int pcn, int pco, int len, int tics)
/* pcn: new percentage complete */
/* pco: old percentage complete */
/* len: length of the progress bar */
/* tics: number of tics */
{
	int pc_n=(pcn*len)/100;
	int pc=pco, tic;
	int i;
	pco=(pco*len)/100;
	if (pco==len)
		return pc;
	
	tic=len/tics;
	
	for (i=pco+1;i<=pc_n;i++)
	{
		if(i%tic == 0)
		{
			if (i==0)
				printf("|");
			else
			{
				printf("\b\b\b\b");
				printf("|");
			}
		}
		else
		{
			printf("\b\b\b\b");
			printf("=");
		}
		pc=pcn;
		printf("%3i%%",pc);
	}
	if (pcn>pc)
	{
		printf("\b\b\b\b");
		printf("%3i%%",pcn);
		pc=pcn;	
	}
		
	if (pc_n==len)
		printf("\n");	
	fflush(stdout);
	
	return pc;
}

void PrintHeader()
{
	if (NORMAL>verbose)
		return;
	printf("  _____   __   __  __  ___  ___                                  \n");
	printf(" | _ \\ \\ / /__|  \\/  |/ _ \\/ __|                                 \n");
	printf(" |  _/\\ V /___| |\\/| | (_) \\__ \\                                 \n");
	printf(" |_|   \\_/    |_|  |_|\\___/|___/    PV-Module Simulator          \n");
	printf("________________________________________________________________\n");
	printf("PV-MOS Version %s   %s   B.E. Pieters \n", VERSION, __DATE__); 
	printf("IEK-5 Photovoltaik Forschungszentrum Juelich, Germany\n"); 
	printf("\n");
}

void Disclamer()
{                                            
	printf("\n");
	printf("DISCLAMER:\n"); 
	printf("PVMOS  Copyright (C) 2014  B. E. Pieters\n");
	printf("This program comes with ABSOLUTELY NO WARRANTY\n");
	printf("This is free software, and you are welcome to redistribute\n");
	printf("under certain conditions. You should have received a copy \n");
	printf("of the GNU General Public License along with this program. \n");
	printf("If not, see <http://www.gnu.org/licenses/>.\n"); 
	printf("\n");
}
