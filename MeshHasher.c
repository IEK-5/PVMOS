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
 * dummy program to create the MESHHASH used to check            *
 * compatibility between the mesh creator and user.              *                     
 * If a mesh is saved to file the MESHHASH is stored in the file *                         
 * If a mesh is read the MESHHASH ins the file and of the        *                      
 * current PVMOS version are compared.                           *     
 *                                                               *            
 *****************************************************************/     
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <openssl/md5.h>
#include "main.h"
#include "mesh2d.h"
#include "utils.h"

VERB_LEVEL verbose=QUIET;
/* 
	Rationale: 
		I want PVMOS to recognize whether a saved meshes is incompatible with the current version.
		I could check for version number but that changes more often than the actual mesh structure.
		This means I need to check whether the mesh changed, i.e. I need a mesh version number. The
		cloolest would be if the mesh version is computed from the mesh data structures themselves 
		so that the mesh version is automatically changed when I change the mesh data structures.
	Solution:
		mesh2d gets a meshhash variable, stored in meshhash.h. This function generates the meshhash.
		The procedure is as follows:
		This function generates a standard and simple mesh using the InitMesh function. The mesh is 
		saved to a file. Note that mesh2d.o stores the meshhash variable in the saved mesh. We 
		determine the md5sum of the saved mesh but we exclude the stored hash. This way we get a 
		hash for the simple mesh data structure. I expect it to be highly improbable that the hash 
		is not changed where the data structures are. Note that some changes to data initialization 
		will also change the hash wher not all such changes are in fact incompatible.
*/

void WriteHexSignature(FILE *f, unsigned const char c[])
{
	int i;
	fprintf(f, " ***********************************************************************\n");
	fprintf(f, " * Current PVMOS mesh data structure signature:                        *\n");
	fprintf(f, " *---------------------------------------------------------------------*\n");
	fprintf(f, " *                                                                     *\n");
	fprintf(f, " *");
	for(i = 0; i < 34-MD5_DIGEST_LENGTH; i++)
		fprintf(f, " ");
	for(i = 0; i < MD5_DIGEST_LENGTH; i++)
		fprintf(f, "%02x", c[i]);
	for(i = 0; i < 35-MD5_DIGEST_LENGTH; i++)
		fprintf(f, " ");
	fprintf(f, "*\n");
	fprintf(f, " *                                                                     *\n");
	fprintf(f, " * *********************************************************************");

}

void WriteHeaderFile(FILE *f, unsigned const char c[])
{
	int i;
	
	fprintf(f, "/* *********************************************************************\n");
	fprintf(f, " * PVMOS Header file generated by MeshHasher, DO NOT EDIT              *\n");
	fprintf(f, " *                                                                     *\n");
	fprintf(f, " * What is this?                                                       *\n");
	fprintf(f, " * PVMOS binary mesh format is not guaranteed to be compatable with    *\n");
	fprintf(f, " * future or past PVMOS versions. For simplicity the format is a       *\n");
	fprintf(f, " * straightforward binary dump of the mesh data structure. Obviosously *\n");
	fprintf(f, " * if anything changes in the data structures this format must change. *\n");
	fprintf(f, " * To verify compatibility we thus need some sort of data structure    *\n");
	fprintf(f, " * versioning. The coolest way is obviously to generate a data         *\n");
	fprintf(f, " * structure signature so PVMOS can automatically detect incompatible  *\n");
	fprintf(f, " * changes.                                                            *\n");
	fprintf(f, " * MeshHasher generates an md5 sum of a simple mesh and stores the     *\n");
	fprintf(f, " * hash in a header file for inclusion in PVMOS. The hash is included  *\n");
	fprintf(f, " * in the binary mesh files to verify compatibility.                   *\n");
	fprintf(f, " *                                                                     *\n");
	WriteHexSignature(f, c);
	fprintf(f, "/\n");	
	fprintf(f, "#define _HAS_MESHHASH\n");	
	fprintf(f, "int NMESHHASH=%d;\n",MD5_DIGEST_LENGTH);
	fprintf(f, "unsigned char MESHHASH[] = {\n\t0x%02x",c[0]);
	for(i = 1; i < MD5_DIGEST_LENGTH; i++)
		fprintf(f, ",\n\t0x%02x", c[i]);
	fprintf(f, "\n};\n");
}


int main (int argc, char **argv)
{
	unsigned char c[MD5_DIGEST_LENGTH];
	char *fn;
	double x1=0, x2=1, y1=0, y2=1;
	char *name="meshhash";
	int Nx=3, Ny=3;
	int i;
	int bytes;
   	MD5_CTX mdContext;
	unsigned char data[1024];
	FILE *f;	
	mesh M;
	
	/* temporary filename */
	fn=malloc(sizeof(char));
	fn=tmpnam (fn);
	
	/* create and save mesh */
	M=InitMesh(name, x1, x2, y1, y2, Nx, Ny);
	WriteMesh(fn,&M);
	FreeMesh(&M);
	
	
	/* create hash */	
	f = fopen (fn, "rb");
	if (f == NULL) 
	{
		printf ("meshfile could not be opened.\n");
		exit(1);
	}

	MD5_Init (&mdContext);
	
	for (i=0;i<NMESHHASH;i++) /* discard the hash itself to prevent a circular dependency */
	{
		unsigned char d;
		fread (&d, sizeof(char), 1, f);
	}
	
	while ((bytes = fread (data, 1, 1024, f)) != 0)
		MD5_Update (&mdContext, data, bytes);
	fclose (f);
	remove(fn);
	
	MD5_Final (c,&mdContext);
	
	
	/* output results */
	WriteHexSignature(stdout, c);
	printf("\n");
	
	/* create header file */
	if (argc==2)
	{
		f = fopen (argv[1], "w");
		if (f == NULL) 
		{
			printf ("File %s could not be opened.\n", argv[1]);
			exit(1);
		}
		WriteHeaderFile(f, c);
	}
	return 0;
	
}
