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
 * Parse commands to create edit and modify meshes               *                     
 * and input/output data                                         *         
 *                                                               *            
 *****************************************************************/     

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "main.h"
#include "mesh2d.h"
#include "list.h"
#include "select_nodes.h"
#include "dataexport.h"
#include "utils.h"
#include "solve.h"
#include "parsedef.h"
#include "parse.h"

#define MAXSTRLEN 256
#define TWOPI 6.28318530717959

/* parsing utility functions */

/* lookup a word (name) in the keyword table (AKeyTable) */
static PRSDEF LookupKey (char *name,  const KeyWord * AKeyTable)
{
      unsigned int len;
      len=strlen(name);
      for (; AKeyTable->name; AKeyTable++)
      {
      	    if(strlen(AKeyTable->name)==len)
	    	if (strncmp (name, AKeyTable->name, len) == 0)
			break;
      }
      return AKeyTable->PAR;
}

/* Find the beginning of the next word in a string (line) */
/* i.e. skip whitespace until we hit a non whitespace character or reach the end of the string */
/* returns either a pointer to the first character of the next word or NULL if we've reached the end of the string  */
static char * Begin (char *line)
{
	/* returns a pointer to the beginning of a word, NULL if end is reached */
	while(isspace((*line)) && *line)
		line++; 
	if (*line)
		return line;
	else
		return NULL;
}
/* Find the end of the current word in a string (line) */
/* i.e. skip non-whitespace character untill we reach whitespace of the end of the string */
/* returns a pointer to the first white space character or null-character */
static char * End (char *line)
{
	/* returns a pointer to the end of a word */
	while(!isspace((*line)) && *line)
		line++;
	return line;
}

/* copies one word from begin to word */
/* begin should point to the next word in a string */
/* the function returns a poiinter to the next word in the string  */
static char *GetWord (char *begin, char *word)
{
	char *end;
	if(!begin)
	{
		*word='\0';
		return NULL;
	}
	end=End(begin);	
	word=strncpy(word,begin,end-begin);
	word[end-begin]='\0';
	begin=Begin(end);
	return begin;
}

/* looks up a mesh by its name in the meshes array */
/* a mesh vbariable consists of
	1: a name
	2: a mesh
	3: list of selected nodes in the mesh
	see parse.h for its definition. 
*/
/* input:
         name: name of the mesh
	 meshes: an array of defined mesh variables
	 Nm: the number of meshes in meshes */
/* If a mesh variable with the same name exists the function return a pointer to the mesh varaible
   otherwise it returns NULL (i.e. there is no mesh with the requested name) */
meshvar *LookupMesh (char *name,  meshvar * meshes, int Nm)
{
      unsigned int len;
      int i;
      
      len=strlen(name);
      for (i=0; i<Nm; i++)
      {
      	    if(strlen(meshes[i].name)==len)
	    	if (strncmp (name, meshes[i].name, len) == 0)
			break;
      }
      if (i==Nm)
      	return NULL;
      return meshes+i;
}

/* looks up a mesh by its name in the meshes array */
/* does more or less the same as the function above except that it returns a pointer to the mesh in the mesh variable
   rather than to the mesh variable (i.e. mesh+name+selected nodes)*/
mesh *FetchMesh(char *name, meshvar * meshes, int Nm)
{
	meshvar *MV;
	
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("Mesh %s does not exist\n",name);
	return &(MV->M);
}	


/* Adds a mesh (M) to the meshes mesh variable array */	
/* Input:
	M: the mesh to add
	name: pointer to an allocated string with the name for the mesh
	meshes: mesh variable array	
	Nm: number of mesh variables in meshes */					
void AddMeshVar (mesh M, char **name,  meshvar ** meshes, int *Nm)
{
	if (strchr(*name, '.'))
		Error("Mesh name %s invalid, no dots (.) allowed.\n", *name);

	if (LookupMesh (*name,  *meshes, *Nm))
		Error("Mesh already exists\n");
	(*meshes)[*Nm].M=M;
	(*meshes)[*Nm].name=*name;
	(*meshes)[*Nm].nodes=malloc(LISTBLOCK*sizeof(int));
	(*meshes)[*Nm].nodes[0]=0;
	(*Nm)++;
	(*meshes)=realloc((*meshes),((*Nm)+1)*sizeof(meshvar));	
}

/* removes the mesh with name from the mesh variable array */	
/* Input:
	name: name of the mesh to remove
	meshes: mesh variable array	
	Nm: number of mesh variables in meshes */					
void RemoveMeshVar (char **name,  meshvar ** meshes, int *Nm)
{
      	unsigned int len;
      	int i;
      
	if (strchr(*name, '.'))
		Error("Mesh name %s invalid, no dots (.) allowed.\n", *name);
      
      	len=strlen((*name));
      	for (i=0; i<(*Nm); i++)
      	{
      	    	if(strlen((*meshes)[i].name)==len)
			if (strncmp ((*name), (*meshes)[i].name, len) == 0)
				break;
      	}
      	if (i==(*Nm))
		Error("Mesh does not exist exists\n");
	
	FreeMesh(&((*meshes)[i].M));
	free((*meshes)[i].name);
	free((*meshes)[i].nodes);
      	for (; i<(*Nm)-1; i++)
	{
		(*meshes)[i].M=(*meshes)[i+1].M;
		(*meshes)[i].name=(*meshes)[i+1].name;
		(*meshes)[i].nodes=(*meshes)[i+1].nodes;
	}
	
	(*Nm)--;
	(*meshes)=realloc((*meshes),((*Nm)+1)*sizeof(meshvar));	
}


/* select rectangular area in a mesh variable */
/* Input: 
	name: name of the mesh
	meshes: mesh variable array
	Nm: number of meshes in the mesh variable array
	x1: x-coordinate of the lower left corner to the to be selected rectangle
	y1: y-coordinate of the lower left corner to the to be selected rectangle 
	x2: x-coordinate of the upper right corner to the to be selected rectangle
	y2: y-coordinate of the upper right corner to the to be selected rectangle 

  The routing sets the selected nodes list for the mesh variable with the name "name"	
*/	
void SelectRectNodes (char *name,  meshvar * meshes, int Nm, double x1, double y1, double x2, double y2)
{
	meshvar *MV;
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("Mesh %s is not defined\n", name);
	MV->nodes=RectSelectNodes(x1, y1, x2, y2, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected\n",MV->nodes[0]);
}
/* select circular area in a mesh variable */
/* Input: 
	name: name of the mesh
	meshes: mesh variable array
	Nm: number of meshes in the mesh variable array
	x: x-coordinate center of the to be selected circle
	y: y-coordinate center of the to be selected circle
	r: radius of the to be selected circle

  The routing sets the selected nodes list for the mesh variable with the name "name"	
*/	
void SelectCircNodes (char *name,  meshvar * meshes,  int Nm, double x, double y, double r)
{
	meshvar *MV;
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("Mesh %s is not defined\n", name);
	MV->nodes=CircSelectNodes(x, y, r, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected\n",MV->nodes[0]);
}
/* select area enclosed by a polygon in a mesh variable */
/* Input: 
	name: name of the mesh
	meshes: mesh variable array
	Nm: number of meshes in the mesh variable array
	P: polygon (see select_nodes.h for the polygon struct definition)
	
  The routing sets the selected nodes list for the mesh variable with the name "name"	
*/	
void SelectPolyNodes (char *name,  meshvar * meshes, int Nm, polygon P)
{
	meshvar *MV;
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("Mesh %s is not defined\n", name);
	MV->nodes=PolySelectNodes(P, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected\n",MV->nodes[0]);
}

/* looks up an area within a mesh by its name in the meshes array */
/* input:
         name: name of the area in a mesh. general form: <mesh-name>.<area-name>
	 meshes: an array of defined mesh variables
	 Nm: the number of meshes in meshes 
	 P: pointer to the index of the area in the mesh variable
 The function sets the P pointer and retuns the mesh variable */
meshvar *LookupMeshArea (char *name,  meshvar * meshes, int Nm, int *P)
{
	meshvar *MV;

	int len, i;
	i=0;
	len=strlen(name);
	while ((i<len)&&(name[i]!='.'))
		i++;
	if (i==len)
		Error("Expected a string of the form <mesh_name>.<area_name>, got %s\n", name);
		
	name[i]='\0';
	
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("Mesh %s is not defined\n", name);
	(*P)=FindProperties(MV->M, name+i+1);
	name[i]='.';
	return MV;
}

/* the parse routine, parses the commands in a file */
void Parse (char *file)
{
	FILE *f;
	char *line, *word, c;
	char *begin;
	meshvar * Meshes;
	int Nm=0, k;
	int line_nr=1;
	int MaxIter=25;
	double TolV=1e-5, RelTolV=1e-5, TolKcl=1e-5, RelTolKcl=1e-5;
	PRSDEF key;
	polygon P;
	P.N=0;
	
	if ((f=fopen(file,"r"))==NULL)
		Error("Cannot open file %s\n", file);
	

	line=malloc(MAXSTRLEN*sizeof(char));
	word=malloc(MAXSTRLEN*sizeof(char));
	
	Meshes=malloc((Nm+1)*sizeof(meshvar));
	
    	fgets(line, MAXSTRLEN-1, f);
	
	while(feof(f)==0)
	{
	
    		k=sscanf(line, " %c", &c);
		if((k!=-1)&&(c!='#'))
		{
			/* read a word, separate at whitespace, check what the word means then read its value or execute */
			begin=Begin(line);
			while(begin)
			{
				begin=GetWord (begin, word);
				key=LookupKey (word,  KeyTable);
				switch(key)
				{
					/* parameters */
					/********************************* Secion Meshing */
					case NEWMESH:
					{
						double x1, x2, y1, y2;
						char *name;
						int Nx, Ny, len;
						/* read next word and process, if no next word trow an error */
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						x1=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						y1=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						x2=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						y2=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						Nx=atoi(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						Ny=atoi(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						len=strlen(word);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,word,len+1);	
						
						Print(NORMAL,"* line %3d: Creating new mesh %s\n", line_nr, word);
						AddMeshVar (InitMesh(word, x1, x2, y1, y2, Nx, Ny), &name,  &Meshes, &Nm);
						break;
					}
					case JOINMESH:
					{
						double xoff, yoff;
						char *name;
						int len;
						mesh *M1, *M2;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						xoff=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						yoff=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M1=FetchMesh(word, Meshes, Nm);
						
						Print(NORMAL,"* line %3d: Joining mesh %s ", line_nr,word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						M2=FetchMesh(word, Meshes, Nm);	
						
						Print(NORMAL,"and mesh %s ", word);		
								
						begin=GetWord (begin, word);				
						if(word[0]=='\0')
							goto premature_end;
							
						len=strlen(word);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,word,len+1);	
										
						Print(NORMAL,"to %s\n", word);
						AddMeshVar (JoinMeshes(*M1, *M2, xoff, yoff), &name,  &Meshes, &Nm);
						break;
					}
					case JOINMESH_H:
					{
						double yoff;
						char *name;
						int len, i;
						mesh *M1, *M2;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						yoff=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M1=FetchMesh(word, Meshes, Nm);
						
						Print(NORMAL,"* line %3d: Joining mesh %s ",line_nr, word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						M2=FetchMesh(word, Meshes, Nm);	
						
						Print(NORMAL,"and mesh %s ", word);		
								
						begin=GetWord (begin, word);				
						if(word[0]=='\0')
							goto premature_end;
							
						len=strlen(word);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,word,len+1);	
										
						Print(NORMAL,"to %s\n", word);
						AddMeshVar (JoinMeshes_H(*M1, *M2, yoff), &name,  &Meshes, &Nm);
						break;
					
					}
					case JOINMESH_V:
					{
						double xoff;
						char *name;
						int len;
						mesh *M1, *M2;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						xoff=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M1=FetchMesh(word, Meshes, Nm);
						
						Print(NORMAL,"* line %3d: Joining mesh %s ",line_nr, word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						M2=FetchMesh(word, Meshes, Nm);	
						
						Print(NORMAL,"and mesh %s ", word);		
								
						begin=GetWord (begin, word);				
						if(word[0]=='\0')
							goto premature_end;
							
						len=strlen(word);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,word,len+1);	
										
						Print(NORMAL,"to %s\n", word);
						AddMeshVar (JoinMeshes_V(*M1, *M2, xoff), &name,  &Meshes, &Nm);
						break;
					
					}
					case SPLITX:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("Mesh %s does not exist\n",word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in x-direction\n",line_nr, word);
							fflush(stdout);
							SplitMeshX(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in x-direction\n",line_nr, word);
							fflush(stdout);
							SplitListX(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}					
						break;
					}
					case SPLITY:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("Mesh %s does not exist\n",word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in y-direction\n",line_nr, word);
							fflush(stdout);
							SplitMeshY(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in y-direction\n",line_nr, word);
							fflush(stdout);
							SplitListY(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}					
						break;
					}
					case SPLITXY:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("Mesh %s does not exist\n",word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in x- and y-direction\n",line_nr, word);
							fflush(stdout);
							SplitMeshXY(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in x- and y-direction\n",line_nr, word);
							fflush(stdout);
							SplitListXY(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}					
						break;
					}
					case SPLITLONG:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("Mesh %s does not exist\n",word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in the longest direction\n",line_nr, word);
							fflush(stdout);
							SplitMeshLong(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in the longest direction\n",line_nr, word);
							fflush(stdout);
							SplitListLong(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}					
						break;
					}
					case SIMPLIFY_MESH:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						Print(NORMAL,"* line %3d: Simplifying mesh of %d nodes\n",line_nr, MV->M.Nn);
						fflush(stdout);		
						Chunkify(&(MV->M), 4);
						Print(NORMAL,"            -->  Resulting mesh is %d nodes large\n", MV->M.Nn);
						MV->nodes[0]=0;
						break;
					}
					case LOADMESH:
					{
						int len;
						char *name;
						mesh Mnew;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL, "* line %3d: Reading mesh from file %s\n",line_nr, word);
						ReadMesh(word, &Mnew);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						len=strlen(word);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,word,len+1);							
						
						AddMeshVar (Mnew, &name,  &Meshes, &Nm);
						Print(NORMAL, "* line %3d: Asigning %s as a name for the new mesh\n",line_nr,word);
						break;
					}
					case SAVEMESH:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Saving mesh %s ",line_nr,word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL, "to file %s\n",word);
						WriteMesh(word,M);
						
						break;
					}
					/********************************* Secion Output */
					/********************************* SubSecion mesh data */
					case PRINTMESH:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print nodes of mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						
						Print(NORMAL, "to file %s\n",word);
						PrintMesh(word,M);
						break;
					}
					case PRINTCONN:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print connections in mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL, "to file %s\n",word);
						
						PrintConn(word,M);
						break;
					}
					case PRINTSURF:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print area definition per node in mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL, "to file %s\n",word);	
						
						PrintSurfDef(word,M);
						break;
					}
					case PRINTPOT:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print node potentials in mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						
						Print(NORMAL, "to file %s\n",word);
						PrintSurfV(word,M);
						break;
					}
					case PRINTMESHSEL:
					{
						double x1,y1,x2,y2;
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print selected nodes of mesh %s ",line_nr,word);
									
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x1=atof(word);		
							
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y1=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x2=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y2=atof(word);								
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						
						Print(NORMAL, "to file %s\n",word);
						PrintMeshSel(word,M, x1, y1, x2, y2);
						break;
					}
					case PRINTCONNSEL:
					{
						double x1,y1,x2,y2;
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
									
						Print(NORMAL, "* line %3d: Print selected connections in mesh %s ",line_nr,word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x1=atof(word);		
							
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y1=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x2=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y2=atof(word);								

												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						
						Print(NORMAL, "to file %s\n",word);
						PrintConnSel(word,M, x1, y1, x2, y2);
						break;
					}
					case PRINTSURFSEL:
					{
						double x1,y1,x2,y2;
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print area definition per selected node in mesh %s ",line_nr,word);
									
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x1=atof(word);		
							
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y1=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x2=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y2=atof(word);								

												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						
						Print(NORMAL, "to file %s\n",word);
						PrintSurfDefSel(word,M, x1, y1, x2, y2);
						break;
					}
					case PRINTPOTSEL:
					{
						double x1,y1,x2,y2;
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print selected node potentials in mesh %s ",line_nr,word);
									
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x1=atof(word);		
							
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y1=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x2=atof(word);		
			
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y2=atof(word);								

												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						
						Print(NORMAL, "to file %s\n",word);
						PrintSurfVSel(word,M, x1, y1, x2, y2);
						break;
					}
					case PRINTIV:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print simulated I-V pairs of mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						
						Print(NORMAL, "to file %s\n",word);
						PrintIV(word,M);
						break;
					}
					case PRINTPARS:
					{
						mesh *M;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Print parameters per area in mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						
						Print(NORMAL, "to file %s\n",word);
						PrintPars(word, M);
						break;
					}
					case SURFVPLOT:
					{
						mesh *M;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Export potentials from mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						x1=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						y1=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						x2=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						y2=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Nx=atoi(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Ny=atoi(word);	
									
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Va=atof(word);						
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Vai=FindVa(Va, M->res.Va, M->res.Nva);
						Print(NORMAL, "to file %s\n",word);
						Print(NORMAL,"            -->  Using simulation at %e V\n", M->res.Va[Vai]);
						SurfVPlot(word, M, Vai, x1, y1, x2, y2, Nx, Ny);
						break;
					}
					case SURFPPLOT:
					{
						mesh *M;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL, "* line %3d: Export power density from mesh %s ",line_nr,word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						x1=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						y1=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						x2=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						y2=atof(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Nx=atoi(word);
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Ny=atoi(word);	
									
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Va=atof(word);						
												
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						
						Vai=FindVa(Va, M->res.Va, M->res.Nva);
						Print(NORMAL, "to file %s\n",word);
						Print(NORMAL,"            -->  Using simulation at %e V\n", M->res.Va[Vai]);
						SurfPPlot(word, M, Vai, x1, y1, x2, y2, Nx, Ny);
						break;
					}
					/********************************* Secion Node Selection */
					case LOAD_POLY:
					{		
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL, "* line %3d: Load polygon from file %s ",line_nr,word);
						P=ReadPoly(word);
						break;
					}	
					case SELECT_RECT:
					{
						double x1,x2,y1,y2;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x1=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y1=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x2=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y2=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Make rectangular selection of nodes\n", line_nr);
						fflush(stdout);
						SelectRectNodes (word,  Meshes, Nm, x1,y1,x2,y2);
						break;
					}	
					case SELECT_CIRC:
					{
						double x,y,r;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						x=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						y=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						r=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Make circular selection of nodes\n", line_nr);
						SelectCircNodes (word,  Meshes, Nm, x,y,r);
						break;
					}
					case SELECT_POLY:
					{
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Make poly-selection of nodes\n", line_nr);
						fflush(stdout);
						SelectPolyNodes (word,  Meshes, Nm, P);
						break;
					}
					
					/***************** Section Local Properties */
					/*
						double Rp, Rn; 
						double Rpvp, Rpvn, Rnvp, Rnvn; 
						double *V, *J; 
						int N;
						int SplitX, SplitY;
					*/
					case ASSIGN_PROP:
					{
						meshvar *MV;
						int P;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						
						if (P<0)
							Error("Area %s is not defined\n", word);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Assigning all nodes in the mesh to area %s\n",line_nr, word);
							fflush(stdout);
							AssignPropertiesMesh(&(MV->M), P);
						}
						else
						{
							Print(NORMAL,"* line %3d: Assigning selected nodes to area %s\n",line_nr, word);
							fflush(stdout);
							AssignProperties(&(MV->M), MV->nodes, P);
						}					
						break;
					}
					case SET_RP:
					{
						meshvar *MV;
						int P;
						double R;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						Print(NORMAL,"* line %3d: Setting Rp of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"to %s\n",word);								
						R=atof(word);
						MV->M.P[P].Rp=R;
						break;
					}
					case SET_RN:
					{
						meshvar *MV;
						int P;
						double R;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting Rn of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"to %s\n",word);								
						R=atof(word);
						MV->M.P[P].Rn=R;
						break;
					}
					case SET_RPVP:
					{
						meshvar *MV;
						int P;
						double R;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting Rpvp of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"to %s\n",word);										
						R=atof(word);
						MV->M.P[P].Rpvp=R;
						break;
					}
					case SET_RNVP:
					{
						meshvar *MV;
						int P;
						double R;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);
						}	
						
						Print(NORMAL,"* line %3d: Setting Rnvp of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"to %s\n",word);								
						R=atof(word);
						MV->M.P[P].Rnvp=R;
						break;
					}
					case SET_RPVN:
					{
						meshvar *MV;
						int P;
						double R;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);
						}	
						
						Print(NORMAL,"* line %3d: Setting Rpvn of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"to %s\n",word);									
						R=atof(word);
						MV->M.P[P].Rpvn=R;
						break;
					}
					case SET_RNVN:
					{
						meshvar *MV;
						int P;
						double R;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting Rnvn of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"to %s\n",word);								
						R=atof(word);
						MV->M.P[P].Rnvn=R;
						break;
					}
					case SET_JV:
					{
						meshvar *MV;
						int P;
						polygon JV;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting tabular JV of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"to data from file %s\n",word);	
						JV=ReadPoly(word);
						if (JV.N<2)
							Error("JV characteristic from file %s is too short\nAt least two voltage current pairs are requires\n", word);
						free(MV->M.P[P].V);
						free(MV->M.P[P].J);
						/* ensure monotoneous increasing voltages */
						BubbleSortJV(JV.N, JV.x, JV.y);
						MV->M.P[P].V=JV.x;
						MV->M.P[P].J=JV.y;
						MV->M.P[P].N=JV.N;
						MV->M.P[P].model=JVD;						
						break;
					}
					case SET_2DJV:
					{
						meshvar *MV;
						int P;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, area);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting 2 diode model parameters of %s\n",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].J01=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].J02=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Jph=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Rs=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Rsh=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].T=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Eg=atof(word);
						MV->M.P[P].nid1=1.0;
						MV->M.P[P].nid2=2.0;							
						MV->M.P[P].model=TWOD;	
						break;
					}
					case SET_1DJV:
					{
						meshvar *MV;
						int P;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting 1 diode model parameters of %s\n",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].J01=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].nid1=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Jph=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Rs=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Rsh=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].T=atof(word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						MV->M.P[P].Eg=atof(word);					
						MV->M.P[P].model=ONED;							
						break;
					}
					case SET_R:
					{
						meshvar *MV;
						int P;
						double R;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Setting R of %s ",line_nr, word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"to %s\n",word);	
						R=atof(word);
						
						MV->M.P[P].V=malloc(3*sizeof(double));
						MV->M.P[P].J=malloc(3*sizeof(double));
						MV->M.P[P].N=2;	
						MV->M.P[P].V[0]=-1.0;
						MV->M.P[P].J[0]=-1.0/R;
						MV->M.P[P].V[1]=1.0;
						MV->M.P[P].J[1]=1.0/R;
						break;
					}
					case SET_SPLITY:
					{
						meshvar *MV;
						int P;
						double S;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						Print(NORMAL,"* line %3d: Toggle split-y parameter of %s\n",line_nr, word);
						if (MV->M.P[P].SplitY)
							MV->M.P[P].SplitY=0;
						else
							MV->M.P[P].SplitY=1;
						break;
					}					
					case SET_SPLITX:
					{
						meshvar *MV;
						int P;
						double S;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=word;
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s\n",line_nr, word);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);
						}	
						
						Print(NORMAL,"* line %3d: Toggle split-x parameter of %s\n",line_nr, word);
						if (MV->M.P[P].SplitX)
							MV->M.P[P].SplitX=0;
						else
							MV->M.P[P].SplitX=1;
						break;
					}
					/********************************* Solving*/		
					case MAXITER:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting maximum number of iterations to %s\n",line_nr, word);							
						MaxIter=atoi(word);
						break;		
					case TOLV:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting absolute voltage totelance to %s\n",line_nr, word);							
						TolV=atof(word);
						break;		
					case RELTOLV:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting relative voltage totelance to %s\n",line_nr, word);							
						RelTolV=atof(word);
						break;		
					case TOLKCL:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting absolute KCL totelance to %s\n",line_nr, word);							
						TolKcl=atof(word);
						break;		
					case RELTOLKCL:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting relative KCL totelance to %s\n",line_nr, word);								
						RelTolKcl=atof(word);
						break;		
					case SOLVE:
					{
						mesh *M;
						double Vstart, Vend;
						int Nstep;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL,"* line %3d: Solving potentials in mesh %s\n",line_nr, word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Vstart=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Vend=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Nstep=atoi(word);
						
						SolveVa(M, Vstart, Vend, Nstep, TolKcl, RelTolKcl, TolV, RelTolV, MaxIter);
						break;
					}
					case ADAPTIVE_SOLVE:
					{
						mesh *M;
						double Va, rel_th;
						int Na;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;							
						M=FetchMesh (word,  Meshes, Nm);
						Print(NORMAL,"* line %3d: Adaptive solving of potentials in mesh %s\n",line_nr, word);
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Va=atof(word);
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						rel_th=atof(word);	
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Na=atoi(word);	
						AdaptiveSolveVa(M, Va, rel_th, Na, TolKcl, RelTolKcl, TolV, RelTolV, MaxIter);
						break;
					
					}
					/********************************* Secion verbosity settings*/
					case _QUIET:
						if(!fixverb)
							verbose=QUIET;
						break;
					case _NORMAL:
						if(!fixverb)
							verbose=NORMAL;
						break;
					case _VERBOSE:
						if(!fixverb)
							verbose=VERBOSE;
						break;
					case _DEBUG:
						if(!fixverb)
							verbose=DEBUG;
						break;
					case NONE:
						Warning("* line %3d: Warning: Word \"%s\" is not recognized\n", line_nr, word);			
					default:
						break;
premature_end:
						Error("Premature end of input\n");
						exit(1);
				
				}
			}
			
		}
    		fgets(line, MAXSTRLEN-1, f);
		line_nr++;
	}
	fclose(f);
	free(line);
	free(word);
	if (P.N)
	{
		free(P.x);
		free(P.y);
	}
	for (k=0;k<Nm;k++)
	{
		FreeMesh(&(Meshes[k].M));
		free(Meshes[k].name);
		free(Meshes[k].nodes);
	}
	free(Meshes);
	return;
}
