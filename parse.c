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


static char **GetArgs (char **begin, int Nw)
{
	char ** res;
	int i;
	res=malloc((Nw+1)*sizeof(char *));
	for (i=0;i<Nw;i++)
	{
		res[i]=malloc((MAXSTRLEN)*sizeof(char));
		(*begin)=GetWord ((*begin), res[i]);
		if(res[i][0]=='\0')
			return NULL;	
	}
	return res;
}

void FreeArgs (char **res, int Nw)
{
	int i;
	for (i=0;i<Nw;i++)
		free(res[i]);
	free(res);
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
		return NULL;
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
		Error("In AddMeshVar: Mesh name %s invalid, no dots (.) allowed.\n", *name);

	if (LookupMesh (*name,  *meshes, *Nm))
		Error("In AddMeshVar: Mesh already exists\n");
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
		Error("In RemoveMeshVar: Mesh name %s invalid, no dots (.) allowed.\n", *name);
      
      	len=strlen((*name));
      	for (i=0; i<(*Nm); i++)
      	{
      	    	if(strlen((*meshes)[i].name)==len)
			if (strncmp ((*name), (*meshes)[i].name, len) == 0)
				break;
      	}
      	if (i==(*Nm))
		Error("In RemoveMeshVar: Mesh does not exist exists\n");
	
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
		Error("In SelectRectNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection");
	MV->nodes=RectSelectNodes(x1, y1, x2, y2, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("In SelectRectNodes: No nodes selected, try selecting a larger area or refine the mesh first");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
}
void SelectRectContourNodes (char *name,  meshvar * meshes, int Nm, double x1, double y1, double x2, double y2, double d)
{
	meshvar *MV;
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("In SelectRectContourNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection");
	MV->nodes=RectContourSelectNodes(x1, y1, x2, y2, d, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("In SelectRectContourNodes: No nodes selected, try selecting a larger area or refine the mesh first");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
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
		Error("In SelectCircNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection");
	MV->nodes=CircSelectNodes(x, y, r, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("In SelectCircNodes: No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
}
void SelectCircContourNodes (char *name,  meshvar * meshes,  int Nm, double x, double y, double r,double d)
{
	meshvar *MV;
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("In SelectCircContourNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection");
	MV->nodes=CircContourSelectNodes(x, y, r, d, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("In SelectCircContourNodes: No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
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
		Error("In SelectPolyNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection");
	MV->nodes=PolySelectNodes(P, MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("In SelectPolyNodes: No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
}
void SelectPolyContourNodes (char *name,  meshvar * meshes, int Nm, polygon P, double d, int loop)
{
	meshvar *MV;
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("In SelectPolyContourNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection");
	MV->nodes=PolyContourSelectNodes(d, P, loop,MV->M, MV->nodes);
	if (MV->nodes[0]==0)
		Error("In SelectPolyContourNodes: No nodes selected, try selecting a larger area or refine the mesh first\n");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
}

void SelectAreaNodes (char *name,  meshvar * meshes, int Nm)
{
	meshvar *MV;
	int len, i;
	i=0;
	len=strlen(name);
	while ((i<len)&&(name[i]!='.'))
		i++;
	if (i==len)
		Error("In LookupMeshArea: Expected a string of the form <mesh_name>.<area_name>, got %s\n", name);
		
	name[i]='\0';
	
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("In SelectAreaNodes: Mesh %s is not defined\n", name);
	if (MV->nodes[0]>0)	
		Print(NORMAL,"            -->  Making sub-selection %s", name+i+1);
	MV->nodes=SelectArea(MV->M, MV->nodes, name+i+1);
	if (MV->nodes[0]==0)
		Error("In SelectAreaNodes: No nodes selected.\n");
	Print(NORMAL,"            -->  %d nodes selected",MV->nodes[0]);
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
		Error("In LookupMeshArea: Expected a string of the form <mesh_name>.<area_name>, got %s\n", name);
		
	name[i]='\0';
	
	MV=LookupMesh (name,  meshes, Nm);
	if (!MV)
		Error("In LookupMeshArea: Mesh %s is not defined\n", name);
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
		Error("In Parse: Cannot open file %s\n", file);
	

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
						char **args;
						char *name;
						int Nx, Ny, len;
						args=GetArgs (&begin, 7);
						if (args==NULL)
							goto premature_end;	
						Print(NORMAL,"* line %3d: Creating new mesh %s", line_nr, args[6]);
						
								
						x1=atof(args[0]);		
						y1=atof(args[1]);	
						x2=atof(args[2]);		
						y2=atof(args[3]);								
						Nx=atoi(args[4]);								
						Ny=atoi(args[5]);	
						len=strlen(args[6]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[6],len+1);					
						AddMeshVar (InitMesh(args[6], x1, x2, y1, y2, Nx, Ny), &name,  &Meshes, &Nm);
						FreeArgs (args, 7);
						break;
					}
					case RMMESH:
					{
						char *name;
						int len;
						/* read next word and process, if no next word trow an error */
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Removing mesh %s", line_nr, word);
						
								
						len=strlen(word);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,word,len+1);	
								
						RemoveMeshVar (&name,  &Meshes, &Nm);
						free(name);
						break;
					}
					case JOINMESH:
					{
						double xoff, yoff;
						char *name;
						char **args;
						int len;
						mesh *M1, *M2;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 5);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Joining mesh %s and mesh %s to %s", line_nr,args[2],args[3],args[4]);	
						
														
						xoff=atof(args[0]);								
						yoff=atof(args[1]);						
						M1=FetchMesh(args[2], Meshes, Nm);
						if (!M1)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[2]);
						M2=FetchMesh(args[3], Meshes, Nm);
						if (!M2)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[3]);	
						
						len=strlen(args[4]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[4],len+1);
						AddMeshVar (JoinMeshes(*M1, *M2, xoff, yoff), &name,  &Meshes, &Nm);
						FreeArgs (args, 5);
						break;
					}
					case JOINMESH_H:
					{
						double yoff;
						char *name;
						char **args;
						int len;
						mesh *M1, *M2;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 4);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Joining mesh %s and mesh %s to %s", line_nr,args[1],args[2],args[3]);	
						
											
						yoff=atof(args[0]);						
						M1=FetchMesh(args[1], Meshes, Nm);
						if (!M1)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						M2=FetchMesh(args[2], Meshes, Nm);
						if (!M2)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[2]);	
						
						len=strlen(args[3]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[3],len+1);
						AddMeshVar (JoinMeshes_H(*M1, *M2, yoff), &name,  &Meshes, &Nm);
						FreeArgs (args, 4);
						break;
					
					}
					case JOINMESH_V:
					{
						double xoff;
						char *name;
						char **args;
						int len;
						mesh *M1, *M2;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 4);
						if (args==NULL)
							goto premature_end;	
						Print(NORMAL,"* line %3d: Joining mesh %s and mesh %s to %s", line_nr,args[1],args[2],args[3]);	
						
										
						xoff=atof(args[0]);						
						M1=FetchMesh(args[1], Meshes, Nm);
						if (!M1)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						M2=FetchMesh(args[2], Meshes, Nm);
						if (!M2)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[2]);	
						
						len=strlen(args[3]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[3],len+1);
						AddMeshVar (JoinMeshes_V(*M1, *M2, xoff), &name,  &Meshes, &Nm);
						FreeArgs (args, 4);
						break;
					
					}
					case DUPMESH:
					{
						int len;
						char *name;
						mesh Mnew;
						char **args;
						mesh *M;
						
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Duplicating mesh %s and storing duplicate in %s",line_nr, args[0], args[1]);
						
						M=FetchMesh (args[0],  Meshes, Nm);							
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);											
						Mnew=DuplicateMesh((*M));						
						len=strlen(args[1]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[1],len+1);						
						AddMeshVar (Mnew, &name,  &Meshes, &Nm);
						FreeArgs (args, 2);
						break;					
					}
					case ADDEL:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Adding an electrode to mesh %s",line_nr, word);
						Print(NORMAL,"            Please define the appropriate properties");	
													
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr, word);		
						
						AddElectrode(&(MV->M));			
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
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in x-direction",line_nr, word);
							fflush(stdout);
							SplitMeshX(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in x-direction",line_nr, word);
							fflush(stdout);
							SplitListX(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d nodes",MV->name, MV->M.Nn);				
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
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in y-direction",line_nr, word);
							fflush(stdout);
							SplitMeshY(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in y-direction",line_nr, word);
							fflush(stdout);
							SplitListY(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d nodes",MV->name, MV->M.Nn);					
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
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in x- and y-direction",line_nr, word);
							fflush(stdout);
							SplitMeshXY(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in x- and y-direction",line_nr, word);
							fflush(stdout);
							SplitListXY(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d nodes",MV->name, MV->M.Nn);					
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
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s in the longest direction",line_nr, word);
							fflush(stdout);
							SplitMeshLong(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s in the longest direction",line_nr, word);
							fflush(stdout);
							SplitListLong(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}
						Print(NORMAL,"            ---> Mesh %s consists of %d nodes",MV->name, MV->M.Nn);
						break;
					}
					case SPLITCOARSE:
					{
						meshvar *MV;
						double d;
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
							
						d=atof(args[1]);
											
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all nodes in %s until all edges shorter than %e",line_nr, MV->name, d);
							fflush(stdout);
							SplitMeshWhileCoarse(&(MV->M), d);
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected nodes in %s until all edges shorter than %e",line_nr, MV->name, d);
							fflush(stdout);
							SplitListWhileCoarse(&(MV->M), MV->nodes, d);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d nodes",MV->name, MV->M.Nn);
						FreeArgs (args, 2);											
						break;
					}
					case SIMPLIFY_MESH:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);							
						Print(NORMAL,"* line %3d: Simplifying mesh of %d nodes",line_nr, MV->M.Nn);
						fflush(stdout);		
						Chunkify(&(MV->M));
						MV->nodes[0]=0;
						Print(NORMAL,"            ---> Mesh %s consists of %d nodes",MV->name, MV->M.Nn);	
						break;
					}
					case LOADMESH:
					{
						int len;
						char *name;
						mesh Mnew;
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
							
						
						Print(NORMAL, "* line %3d: Reading mesh from file %s into new mesh %s",line_nr, args[0],args[1]);
						ReadMesh(args[0], &Mnew);
						
						len=strlen(args[1]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[1],len+1);							
						
						AddMeshVar (Mnew, &name,  &Meshes, &Nm);
						FreeArgs (args, 2);								
						break;
					}
					case SAVEMESH:
					{
						mesh *M;						
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL, "* line %3d: Saving mesh %s to file %s",line_nr,args[0], args[1]);
														
						M=FetchMesh (args[0],  Meshes, Nm);
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						WriteMesh(args[1],M);
						FreeArgs (args, 2);						
						break;
					}
					/********************************* Secion Output */
					/********************************* SubSecion mesh data */
					case PRINTMESH:
					{
						meshvar *MV;						
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL, "* line %3d: Print elements of mesh %s to file %s",line_nr,args[0],args[1]);
						
											
						MV=LookupMesh (args[0],  Meshes, Nm);	
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
								
						PrintMesh(args[1],&(MV->M), MV->nodes);
						FreeArgs (args, 2);	
						break;
					}
					case PRINTCONN:
					{
						meshvar *MV;					
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print connections in mesh %s to file %s",line_nr,args[0], args[1]);
													
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
															
						PrintConn(args[1],&(MV->M), MV->nodes);
						FreeArgs (args, 2);	
						break;
					}
					case PRINTSURF:
					{
						meshvar *MV;					
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print area definition per node in mesh %s to file %s",line_nr,args[0], args[1]);
							
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
							
						PrintSurfDef(args[1],&(MV->M), MV->nodes);
						FreeArgs (args, 2);	
						break;
					}
					case PRINTPOT:
					{
						meshvar *MV;					
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print node potentials in mesh %s to file %s",line_nr,args[0], args[1]);
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
							
						PrintSurfV(args[1],&(MV->M), MV->nodes);
						FreeArgs (args, 2);	
						break;
					}
					case PRINTIV:
					{
						mesh *M;				
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print simulated I-V pairs of mesh %s to file %s",line_nr,args[0], args[1]);
											
						M=FetchMesh (args[0],  Meshes, Nm);
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						PrintIV(args[1],M);
						FreeArgs (args, 2);	
						break;
					}
					case PRINTPARS:
					{
						mesh *M;				
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print parameters per area in mesh %s to file %s",line_nr,args[0], args[1]);
											
						M=FetchMesh (args[0],  Meshes, Nm);
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
								
						PrintPars(args[1], M);
						FreeArgs (args, 2);	
						break;
					}
					case SURFVPLOT:
					{
						mesh *M;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export potentials from mesh %s to file %s",line_nr,args[0], args[8]);
							
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);	
						Va=atof(args[7]);
																	
						M=FetchMesh (args[0],  Meshes, Nm);
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						Vai=FindVa(Va, M->res.Va, M->res.Nva);
						if (Vai>=0)
						{
							Print(NORMAL,"            -->  Using simulation at %e V", M->res.Va[Vai]);
							SurfVPlot(args[8], M, Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present.\n", line_nr);	
							
						FreeArgs (args, 9);	
						break;
					}
					case SURFPPLOT:
					{
						mesh *M;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export power density from mesh %s  to file %s",line_nr,args[0], args[8]);
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);	
						Va=atof(args[7]);
													
						M=FetchMesh (args[0],  Meshes, Nm);
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						
						Vai=FindVa(Va, M->res.Va, M->res.Nva);
						if (Vai>=0)
						{
							Print(NORMAL,"            -->  Using simulation at %e V", M->res.Va[Vai]);
							SurfPPlot(args[8], M, Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present\n", line_nr);	
						FreeArgs (args, 9);	
						break;
					}
					/********************************* Secion Node Selection */
					case LOAD_POLY:
					{		
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL, "* line %3d: Load polygon from file %s",line_nr,word);
						P=ReadPoly(word);
						break;
					}	
					case SELECT_RECT:
					{
						double x1,x2,y1,y2;			
						char **args;
						args=GetArgs (&begin, 5);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Make rectangular selection of elements in mesh %s", line_nr, args[4]);
						
						x1=atof(args[0]);
						y1=atof(args[1]);
						x2=atof(args[2]);
						y2=atof(args[3]);
						SelectRectNodes (args[4],  Meshes, Nm, x1,y1,x2,y2);
						FreeArgs (args, 5);	
						break;
					}		
					case SELECT_RECT_CONTOUR:
					{
						double x1,x2,y1,y2, d;			
						char **args;
						args=GetArgs (&begin, 6);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Make selection of elements near a rectangular contour in mesh %s", line_nr, args[5]);
						
						x1=atof(args[0]);
						y1=atof(args[1]);
						x2=atof(args[2]);
						y2=atof(args[3]);
						d=atof(args[4]);
						SelectRectContourNodes (args[5],  Meshes, Nm, x1,y1,x2,y2,d);
						FreeArgs (args, 6);	
						break;
					}	
					case SELECT_CIRC:
					{
						double x,y,r;			
						char **args;
						args=GetArgs (&begin, 4);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Make circular selection of elements in mesh %s", line_nr, args[3]);
						x=atof(args[0]);
						y=atof(args[1]);
						r=atof(args[2]);

						SelectCircNodes (args[3],  Meshes, Nm, x,y,r);
						FreeArgs (args, 4);	
						break;
					}	
					case SELECT_CIRC_CONTOUR:
					{
						double x,y,r, d;			
						char **args;
						args=GetArgs (&begin, 5);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Make selection of elements near a circular contour in mesh %s", line_nr, args[4]);
						x=atof(args[0]);
						y=atof(args[1]);
						r=atof(args[2]);
						d=atof(args[3]);

						SelectCircContourNodes (args[4],  Meshes, Nm, x,y,r,d);
						FreeArgs (args, 5);	
						break;
					}
					case SELECT_POLY:
					{
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Make poly-selection of elements in mesh %s", line_nr, word);
						SelectPolyNodes (word,  Meshes, Nm, P);
						break;
					}
					case SELECT_POLY_CONTOUR:
					{
						double d;	
						int loop=0;		
						char **args;
						args=GetArgs (&begin, 3);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL,"* line %3d: Make poly-contour selection of elements in mesh %s", line_nr, args[1]);
						d=atof(args[0]);
						loop=atoi(args[1]);
						SelectPolyContourNodes (args[2],  Meshes, Nm, P, d, loop);
						FreeArgs (args, 3);	
						break;
					}
					case SELECT_AREA:
					{						

						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Make area selection of elements in mesh %s", line_nr, word);
						SelectAreaNodes (word,  Meshes, Nm);
						break;
					}	
					case DESELECT:
					{
						meshvar *MV;
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);	
						Print(NORMAL,"* line %3d: Deselecting selection in mesh %s", line_nr, MV->name);	
						MV->nodes[0]=0;
						break;
					}
					
					/***************** Section Local Properties */
					case ASSIGN_PROP:
					{
						meshvar *MV;
						int P;
						
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						
						if (P<0)
							Error("* line %3d: Area %s is not defined\n", line_nr, word);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Assigning all nodes in the mesh to area %s",line_nr, word);
							fflush(stdout);
							AssignPropertiesMesh(&(MV->M), P);
						}
						else
						{
							Print(NORMAL,"* line %3d: Assigning selected nodes to area %s",line_nr, word);
							fflush(stdout);
							AssignProperties(&(MV->M), MV->nodes, P);
						}					
						break;
					}
					case SET_REL:
					{
						meshvar *MV;
						int P, el;
						double R;		
						char **args;
						args=GetArgs (&begin, 3);
						if (args==NULL)
							goto premature_end;
										
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
							
						el=atoi(args[1]);	
						
						
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);
						Print(NORMAL,"* line %3d: Setting R for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[P].name, args[2]);
												
						R=atof(args[2]);
						MV->M.P[P].Rel[el]=R;
						FreeArgs (args, 3);	
						break;
					}
					case SET_RVP:
					{
						meshvar *MV;
						int P, el;
						double R;		
						char **args;
						args=GetArgs (&begin, 3);
						if (args==NULL)
							goto premature_end;
													
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);
						Print(NORMAL,"* line %3d: Setting Rvp for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[P].name, args[2]);
															
						R=atof(args[2]);
						MV->M.P[P].Rvp[el]=R;
						FreeArgs (args, 3);	
						break;
					}
					case SET_RVN:
					{
						meshvar *MV;
						int P, el;
						double R;	
						char **args;
						args=GetArgs (&begin, 3);
						if (args==NULL)
							goto premature_end;
													
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);
						}	
						
						el=atoi(args[1]);
						
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);
							
						Print(NORMAL,"* line %3d: Setting Rvn for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[P].name, args[2]);
														
						R=atof(args[2]);
						MV->M.P[P].Rvn[el]=R;
						FreeArgs (args, 3);	
						break;
					}
					case SET_JV:
					{
						meshvar *MV;
						int P, el;
						polygon JV;							
						char **args;
						args=GetArgs (&begin, 3);
						if (args==NULL)
							goto premature_end;
												
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}	
						
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter-electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
							
						Print(NORMAL,"* line %3d: Setting tabular JV between electrodes %d and %d in %s.%s to data from file %s",line_nr, el, el+1, MV->name, MV->M.P[P].name, args[2]);
						
						JV=ReadPoly(args[2]);
						if (JV.N<2)
							Error("JV characteristic from file %s is too short\nAt least two voltage current pairs are requires\n", word);
						free(MV->M.P[P].conn[el].V);
						free(MV->M.P[P].conn[el].J);
						/* ensure monotoneous increasing voltages */
						BubbleSortJV(JV.N, JV.x, JV.y);
						MV->M.P[P].conn[el].V=JV.x;
						MV->M.P[P].conn[el].J=JV.y;
						MV->M.P[P].conn[el].N=JV.N;
						MV->M.P[P].conn[el].model=JVD;	
						FreeArgs (args, 3);						
						break;
					}
					case SET_2DJV:
					{
						meshvar *MV;
						int P, el;						
						char **args;
						args=GetArgs (&begin, 8);
						if (args==NULL)
							goto premature_end;
							
													
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, area);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter-electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
							
						Print(NORMAL,"* line %3d: Setting 2 diode model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[P].name);
						
						MV->M.P[P].conn[el].J01=atof(args[2]);
						MV->M.P[P].conn[el].J02=atof(args[3]);
						MV->M.P[P].conn[el].Jph=atof(args[4]);
						MV->M.P[P].conn[el].Rs=atof(args[5]);
						MV->M.P[P].conn[el].Rsh=atof(args[6]);
						MV->M.P[P].conn[el].Eg=atof(args[7]);
						MV->M.P[P].conn[el].nid1=1.0;
						MV->M.P[P].conn[el].nid2=2.0;							
						MV->M.P[P].conn[el].model=TWOD;	
						FreeArgs (args, 8);		
						break;
					}
					case SET_1DJV:
					{
						meshvar *MV;
						int P, el;						
						char **args;
						args=GetArgs (&begin, 8);
						if (args==NULL)
							goto premature_end;
								
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, area);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter-electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
							
						Print(NORMAL,"* line %3d: Setting 1 diode model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[P].name);
						
						MV->M.P[P].conn[el].J01=atof(args[2]);
						MV->M.P[P].conn[el].nid1=atof(args[3]);
						MV->M.P[P].conn[el].Jph=atof(args[4]);
						MV->M.P[P].conn[el].Rs=atof(args[5]);
						MV->M.P[P].conn[el].Rsh=atof(args[6]);
						MV->M.P[P].conn[el].Eg=atof(args[7]);					
						MV->M.P[P].conn[el].model=ONED;	
						FreeArgs (args, 8);								
						break;
					}
					case SET_R:
					{
						meshvar *MV;
						int P, el;
						double R;						
						char **args;
						args=GetArgs (&begin, 3);
						if (args==NULL)
							goto premature_end;
								
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, area);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter-electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
							
						Print(NORMAL,"* line %3d: Setting R between electrodes %d and %d in %s.%s to %s",line_nr, el, el+1, MV->name, MV->M.P[P].name, args[2]);
						R=atof(args[2]);
						if (MV->M.P[P].conn[el].V)
							free(MV->M.P[P].conn[el].V);
						if (MV->M.P[P].conn[el].J)
							free(MV->M.P[P].conn[el].J);
						MV->M.P[P].conn[el].V=malloc(3*sizeof(double));
						MV->M.P[P].conn[el].J=malloc(3*sizeof(double));
						MV->M.P[P].conn[el].N=2;	
						MV->M.P[P].conn[el].V[0]=-1.0;
						MV->M.P[P].conn[el].J[0]=-1.0/R;
						MV->M.P[P].conn[el].V[1]=1.0;
						MV->M.P[P].conn[el].J[1]=1.0/R;
						MV->M.P[P].conn[el].model=JVD;	
						FreeArgs (args, 3);		
						break;
					}
					case SET_T:
					{
						meshvar *MV;
						int P;
						double T;						
						char **args;
						args=GetArgs (&begin, 2);
						if (args==NULL)
							goto premature_end;
												
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
														
						Print(NORMAL,"* line %3d: Setting initial temperature in %s to %s",line_nr, args[0], args[1]);
						T=atof(args[1]);
						
						MV->M.P[P].T=T;
						FreeArgs (args, 2);		
						break;
					}
					case SET_SPLITY:
					{
						meshvar *MV;
						int P;						
						char **args;
						args=GetArgs (&begin, 1);
						if (args==NULL)
							goto premature_end;
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
						
						Print(NORMAL,"* line %3d: Toggle split-y parameter of %s",line_nr, args[0]);
						if (MV->M.P[P].SplitY)
							MV->M.P[P].SplitY=0;
						else
							MV->M.P[P].SplitY=1;
						break;
						FreeArgs (args, 1);
					}					
					case SET_SPLITX:
					{
						meshvar *MV;
						int P;					
						char **args;
						args=GetArgs (&begin, 1);
						if (args==NULL)
							goto premature_end;
						MV=LookupMeshArea (args[0],  Meshes, Nm, &P);
						if (P<0)
						{
							char *area;
							area=args[0];
							while ((*area)!='.')
								area++;
							area++;
							Print(NORMAL,"* line %3d: Creating new area definition %s",line_nr, args[0]);
							P=MV->M.Na;
							NewProperties(&(MV->M), area);	
						}
						
						Print(NORMAL,"* line %3d: Toggle split-x parameter of %s",line_nr, args[0]);
						if (MV->M.P[P].SplitX)
							MV->M.P[P].SplitX=0;
						else
							MV->M.P[P].SplitX=1;
						FreeArgs (args, 1);
						break;
					}
					/********************************* Solving*/		
					case MAXITER:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting maximum number of iterations to %s",line_nr, word);							
						MaxIter=atoi(word);
						break;		
					case TOLV:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting absolute voltage totelance to %s",line_nr, word);							
						TolV=atof(word);
						break;		
					case RELTOLV:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting relative voltage totelance to %s",line_nr, word);							
						RelTolV=atof(word);
						break;		
					case TOLKCL:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting absolute KCL totelance to %s",line_nr, word);							
						TolKcl=atof(word);
						break;		
					case RELTOLKCL:
						begin=GetWord (begin, word);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting relative KCL totelance to %s",line_nr, word);								
						RelTolKcl=atof(word);
						break;		
					case SOLVE:
					{
						mesh *M;
						double Vstart, Vend;
						int Nstep;				
						char **args;
						args=GetArgs (&begin, 4);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL,"* line %3d: Solving potentials in mesh %s",line_nr, args[0]);						
						M=FetchMesh (args[0],  Meshes, Nm);
						if (!M)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						
						Vstart=atof(args[1]);
						Vend=atof(args[2]);
						Nstep=atoi(args[3]);
						
						SolveVa(M, Vstart, Vend, Nstep, TolKcl, RelTolKcl, TolV, RelTolV, MaxIter);
						FreeArgs (args, 4);
						break;
					}
					case ADAPTIVE_SOLVE:
					{
						meshvar *MV;
						double Va, rel_th;
						int Na;				
						char **args;
						args=GetArgs (&begin, 4);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL,"* line %3d: Adaptive solving of potentials in mesh %s",line_nr, args[0]);					
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						Va=atof(args[1]);
						rel_th=atof(args[2]);
						Na=atoi(args[3]);
						MV->nodes[0]=0;
						AdaptiveSolveVa(&(MV->M), Va, rel_th, Na, TolKcl, RelTolKcl, TolV, RelTolV, MaxIter);
						FreeArgs (args, 4);
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
						Error("* line %3d: Premature end of input\n", line_nr);
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
