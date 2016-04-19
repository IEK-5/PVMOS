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
 * Parse commands to create edit and modify meshes               *                     
 * and input/output data                                         *         
 *                                                               *            
 *****************************************************************/     

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include "main.h"
#include "mesh2d.h"
#include "list.h"
#include "polygon.h"
#include "select_nodes.h"
#include "dataexport.h"
#include "utils.h"
#include "diode.h"
#include "phototransistor.h"
#include "solve.h"
#include "parsedef.h"
#include "parse.h"
#include "expr.h"

#define MAXSTRLEN 256
#define TWOPI 6.28318530717959

/* parsing utility functions */

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
	(*meshes)[*Nm].setsel=0;
	GetMeshBB(&M, &((*meshes)[*Nm].x1), &((*meshes)[*Nm].y1), &((*meshes)[*Nm].x2), &((*meshes)[*Nm].y2));
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
		(*meshes)[i].setsel=(*meshes)[i+1].setsel;
		(*meshes)[i].x1=(*meshes)[i+1].x1;
		(*meshes)[i].y1=(*meshes)[i+1].y1;
		(*meshes)[i].x2=(*meshes)[i+1].x2;
		(*meshes)[i].y2=(*meshes)[i+1].y2;
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
		Warning("In SelectRectNodes: No elements selected, try selecting a larger area or refine the mesh first");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
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
		Warning("In SelectRectContourNodes: No elements selected, try selecting a larger area or refine the mesh first");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
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
		Warning("In SelectCircNodes: No elements selected, try selecting a larger area or refine the mesh first\n");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
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
		Warning("In SelectCircContourNodes: No elements selected, try selecting a larger area or refine the mesh first\n");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
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
		Warning("In SelectPolyNodes: No elements selected, try selecting a larger area or refine the mesh first\n");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
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
		Warning("In SelectPolyContourNodes: No elements selected, try selecting a larger area or refine the mesh first\n");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
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
		Warning("In SelectAreaNodes: No elements selected.\n");
	MV->setsel=0;
	Print(NORMAL,"            -->  %d elements selected",MV->nodes[0]);
	
	
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

static PRSMESHVAR LookupMeshVar (char *name,  const MV_KeyWord * AKeyTable)
{
      unsigned int len;
      for (; AKeyTable->name; AKeyTable++)
      {
      	    len=strlen(AKeyTable->name);
      	    if(len<MAXSTRLEN)
	    	if (strncmp (name, AKeyTable->name, len) == 0)
			break;
      }
      return AKeyTable->PAR;
}
static char **GetMeshVarArgs (char *expr, int Nw)
{
	char ** res;
	char *p, *q, c;
	int i;
	res=malloc((Nw+1)*sizeof(char *));
	p=expr;
	while ((*p)&&(*p!='('))
		p++;
	if (*p!='(')
		Error("In EvalMeshVar: Expected a string of the form <mesh_name>.<var_name>(..,..), got %s\n", expr);
	for (i=0;i<Nw;i++)
	{
		res[i]=malloc((MAXSTRLEN)*sizeof(char));
		p++;	
		q=p;
		while ((*p)&&(*p!=',')&&(*p!=')'))
			p++;
		
		if (*p)
		{
			c=*p;
			*p='\0';
			strncpy(res[i],q,(p-q+1));
			*p=c;
		}
		else
			Error("In EvalMeshVar: Expected a string of the form <mesh_name>.<var_name>(..,..), got %s\n", expr);
		
	}
	return res;
}


void FreeArgs (char **res, int Nw);
static int EvalMeshVar(char *expr, char *res, meshvar * meshes, int Nm)
{
	char * p;
	mesh *M;
	meshvar *MV;
	PRSMESHVAR var;
	p=expr;
	while((*p)&&(*p!='.'))
		p++;	
	if (*p!='.')
		Error("In EvalMeshVar: Expected a string of the form <mesh_name>.<var_name>, got %s\n", expr);
	*p='\0';
	MV=LookupMesh (expr,  meshes, Nm);
	if (!MV)
		Error("In EvalMeshVar: mesh %s does not exist\n",expr);
	M=&(MV->M);
	*p='.';
	p++;	
	var=LookupMeshVar (p,  MV_KeyTable);
	switch(var)
	{
		case MV_NELEC:
			/* number of electrodes */
			snprintf(res, MAXSTRLEN-1, "%d", M->Nel);
			return 0;
		case MV_NEL:
			/* number of elements */
			snprintf(res, MAXSTRLEN-1, "%d", M->Nn);
			return 0;
		case MV_NAREA:		
			/* number of areas */
			snprintf(res, MAXSTRLEN-1, "%d", M->Na);
			return 0;
		case MV_NVA:
			/* number of applied voltages simulated */
			snprintf(res, MAXSTRLEN-1, "%d", M->res.Nva);
			return 0;
		case MV_NSEL:
			/* number of selected elements */
			snprintf(res, MAXSTRLEN-1, "%d", MV->nodes[0]);
			return 0;	
		case MV_BB_X1:	
			/* bounding box lower x-coordinate */
			snprintf(res, MAXSTRLEN-1, "%e", MV->x1);
			return 0;	
		case MV_BB_X2:	
			/* bounding box upper x-coordinate */
			snprintf(res, MAXSTRLEN-1, "%e", MV->x2);
			return 0;	
		case MV_BB_Y1:	
			/* bounding box lower y-coordinate */
			snprintf(res, MAXSTRLEN-1, "%e", MV->y1);
			return 0;	
		case MV_BB_Y2:	
			/* bounding box upper y-coordinate */
			snprintf(res, MAXSTRLEN-1, "%e", MV->y2);
			return 0;	
		case MV_ISC:
		{
			/* short circuit current */
			int isc, imp_m, imp, imp_p, ioc_m, ioc_p;
						
			SolPar(M, &isc, &imp_m, &imp, &imp_p, &ioc_m, &ioc_p);
			if (isc>=0)					
				snprintf(res, MAXSTRLEN-1, "%e", M->res.I[isc]);
			else
			{
				Warning("In EvalMeshVar: Cannot determine Isc, setting it to 0\n");
				snprintf(res, MAXSTRLEN-1, "0.0");
			}
			return 0;
		}
		case MV_VOC:
		{
			/* open circuit voltage */
			int isc, imp_m, imp, imp_p, ioc_m, ioc_p;
						
			SolPar(M, &isc, &imp_m, &imp, &imp_p, &ioc_m, &ioc_p);			
			if ((ioc_m>=0)||(ioc_p>=0))
			{
				
				if ((ioc_m>=0)&&(ioc_p>=0))
				{
					double Voc;
					Voc=(MV->M.res.Va[ioc_m]*fabs(MV->M.res.I[ioc_p])+MV->M.res.Va[ioc_p]*fabs(MV->M.res.I[ioc_m]))/(fabs(MV->M.res.I[ioc_p])+fabs(MV->M.res.I[ioc_m]));
					snprintf(res, MAXSTRLEN-1, "%e", Voc);
				}
				else
				{
					if (ioc_m>=0)
						snprintf(res, MAXSTRLEN-1, "%e", MV->M.res.Va[ioc_m]);
					else 
						snprintf(res, MAXSTRLEN-1, "%e", MV->M.res.Va[ioc_p]);
				}
			}
			else
			{
				Warning("In EvalMeshVar: Cannot determine Voc, setting it to 0\n");
				snprintf(res, MAXSTRLEN-1, "0.0");
			}
			return 0;
		}
		case MV_IMPP:
		{
			/* Maximum Power Point Current */
			int isc, imp_m, imp, imp_p, ioc_m, ioc_p;
						
			SolPar(M, &isc, &imp_m, &imp, &imp_p, &ioc_m, &ioc_p);		
			if (imp>=0)
				snprintf(res, MAXSTRLEN-1, "%e", MV->M.res.I[imp]);
			else
			{
				Warning("In EvalMeshVar: Cannot determine Impp, setting it to 0\n");
				snprintf(res, MAXSTRLEN-1, "0.0");
			}
			return 0;
		}
		case MV_VMPP:
		{
			/* Maximum Power Point Voltage */
			int isc, imp_m, imp, imp_p, ioc_m, ioc_p;
						
			SolPar(M, &isc, &imp_m, &imp, &imp_p, &ioc_m, &ioc_p);		
			if (imp>=0)
				snprintf(res, MAXSTRLEN-1, "%e", MV->M.res.Va[imp]);
			else
			{
				Warning("In EvalMeshVar: Cannot determine Vmpp, setting it to 0\n");
				snprintf(res, MAXSTRLEN-1, "0.0");
			}	
			return 0;
		}
		case MV_EL_XY:
		{
			/* return the element id at (x,y) */
			char **args;
			double x, y;
			int notinmesh, id;
			
			
			args=GetMeshVarArgs (expr, 2);
			
			x=atof(args[0]);
			y=atof(args[1]);
			id=FindPos(*M, 0, x, y, &notinmesh);
			
			if (notinmesh)
			{
				Warning("In EvalMeshVar: Position (%e,%e) is not in the mesh, returning -1\n", x, y);
				snprintf(res, MAXSTRLEN-1, "-1");
			}
			else
				snprintf(res, MAXSTRLEN-1, "%d", id);
			FreeArgs(args,2);
			return 0;
		}
		case MV_X_EL:		
		{
			/* return the x coordinate of element id*/
			char **args;
			int i;
			
			
			args=GetMeshVarArgs (expr, 1);
			
			i=atoi(args[0]);
			
			if ((i<0)||(i>=M->Nn))
			{
				Warning("In EvalMeshVar: Element index %d out of range (0..%d), returning NaN\n", i, M->Nn);
				snprintf(res, MAXSTRLEN-1, "NaN");
				
			}
			else		
				snprintf(res, MAXSTRLEN-1, "%e", (M->nodes[i].x1+M->nodes[i].x2)/2);
			FreeArgs(args,1);
			return 0;
		}
		case MV_Y_EL:		
		{
			/* return the element id at (x,y) */
			char **args;
			int i;
			
			
			args=GetMeshVarArgs (expr, 1);
			
			i=atoi(args[0]);
			
			if ((i<0)||(i>=M->Nn))
			{
				Warning("In EvalMeshVar: Element index %d out of range (0..%d), returning NaN\n", i, M->Nn);
				snprintf(res, MAXSTRLEN-1, "NaN");
			}
			else		
				snprintf(res, MAXSTRLEN-1, "%e", (M->nodes[i].y1+M->nodes[i].y2)/2);
			FreeArgs(args,1);
			return 0;
		}				
		case MV_AREA_I:		
		{
			/* return the area name for index i */
			char **args;
			int id;
			
			
			args=GetMeshVarArgs (expr, 1);
			
			id=atoi(args[0]);
			if ((id>=0)&&(id<M->Na))
				snprintf(res, MAXSTRLEN-1, "%s", M->P[id].name);	
			else
			{
				Warning("In EvalMeshVar: area index %d out of range (0..%d), returning _AREA_NOT_DEFINED_\n", id, M->Na-1);
				snprintf(res, MAXSTRLEN-1, "_AREA_NOT_DEFINED_");
			}
			FreeArgs(args,1);
			return 0;
		}
		case MV_I_V:
		{
			/* return the i-th applied voltage */
			char **args;
			int Vai;
			double Va;
			
			
			args=GetMeshVarArgs (expr, 1);
			
			Va=atof(args[0]);
			Vai=FindVa(Va, M->res.Va, M->res.Nva);
			if (M->res.Nva)
				snprintf(res, MAXSTRLEN-1, "%d", Vai);	
			else
			{
				Warning("In EvalMeshVar: could not determine a simulation index, no simulation results available, returning -1\n");
				snprintf(res, MAXSTRLEN-1, "-1");
			}
			FreeArgs(args,1);
			return 0;
		}	
				
		case MV_V_I:			
		{
			/* return the i-th applied voltage */
			char **args;
			int id;
			
			
			args=GetMeshVarArgs (expr, 1);
			
			id=atoi(args[0]);
			if ((id>=0)&&(id<M->res.Nva))
				snprintf(res, MAXSTRLEN-1, "%e", M->res.Va[id]);	
			else
			{
				Warning("In EvalMeshVar: simulation index %d out of range (0..%d), returning NaN\n", id, M->res.Nva-1);
				snprintf(res, MAXSTRLEN-1, "NaN");
			}
			FreeArgs(args,1);
			return 0;
		}	
		case MV_I_I:		
		{
			/* return the i-th applied voltage */
			char **args;
			int id;
			
			
			args=GetMeshVarArgs (expr, 1);
			
			id=atoi(args[0]);
			if ((id>=0)&&(id<M->res.Nva))
				snprintf(res, MAXSTRLEN-1, "%e", M->res.I[id]);	
			else
			{
				Warning("In EvalMeshVar: simulation index %d out of range (0..%d), returning NaN\n", id, M->res.Nva-1);
				snprintf(res, MAXSTRLEN-1, "NaN");
			}
			FreeArgs(args,1);
			return 0;
		}
		case MV_VI_I:
		{
			/* return the i-th applied voltage */
			char **args;
			double x;
			double *V, *I;
			int i;
			
			args=GetMeshVarArgs (expr, 1);
			x=atof(args[0]);
			if (M->res.Nva==0)
			{
				Warning("In EvalMeshVar: no simulated IV data available, cannot estimate voltage @ %e A, returning NaN\n\n", x);
				snprintf(res, MAXSTRLEN-1, "NaN");
				FreeArgs(args,1);
				return 0;
			
			}
			
			
			if (M->res.Nva==1)
			{
			
				snprintf(res, MAXSTRLEN-1, "%e", M->res.Va[0]);				
				FreeArgs(args,1);
				return 0;
			}
			
			V=malloc((M->res.Nva+1)*sizeof(double));
			I=malloc((M->res.Nva+1)*sizeof(double));
			
			for (i=0;i<M->res.Nva;i++)
			{
				V[i]=M->res.Va[i];
				I[i]=M->res.I[i];
			}
			
			BubbleSortJV(M->res.Nva, I, V);
			if (I[0]>x)
				snprintf(res, MAXSTRLEN-1, "%e", V[0]+(x-I[0])*(V[1]-V[0])/(I[1]-I[0]));	
			else
			{
				i=0;
				while(i<M->res.Nva-1)
				{
					if (x==I[i])
						break;
					if ((x-I[i])*(x-I[i+1])<0)
						break;
					i++;				
				}
				if (x==I[i])
					snprintf(res, MAXSTRLEN-1, "%e", V[i]);
				else
					snprintf(res, MAXSTRLEN-1, "%e", V[i]+(x-I[i])*(V[i+1]-V[i])/(I[i+1]-I[i]));
				
			}
			free(V);
			free(I);
			
			FreeArgs(args,1);
			return 0;
		
		}
		case MV_IV_V:	
		{
			/* return the i-th applied voltage */
			char **args;
			double x;
			double *V, *I;
			int i;
			
			args=GetMeshVarArgs (expr, 1);
			x=atof(args[0]);
			if (M->res.Nva==0)
			{
				Warning("In EvalMeshVar: no simulated IV data available, cannot estimate current @ %e V, returning NaN\n", x);
				snprintf(res, MAXSTRLEN-1, "NaN");
				FreeArgs(args,1);
				return 0;
			
			}
			
			
			if (MV->M.res.Nva==1)
			{
			
				snprintf(res, MAXSTRLEN-1, "%e", MV->M.res.Va[0]);				
				FreeArgs(args,1);
				return 0;
			}
			
			V=malloc((M->res.Nva+1)*sizeof(double));
			I=malloc((M->res.Nva+1)*sizeof(double));
			
			for (i=0;i<M->res.Nva;i++)
			{
				V[i]=M->res.Va[i];
				I[i]=M->res.I[i];
			}
			
			BubbleSortJV(M->res.Nva, V, I);
			if (V[0]>x)
				snprintf(res, MAXSTRLEN-1, "%e", I[0]+(x-V[0])*(I[1]-I[0])/(V[1]-V[0]));	
			else
			{
				i=0;
				while(i<M->res.Nva-1)
				{
					if (x==V[i])
						break;
					if ((x-V[i])*(x-V[i+1])<0)
						break;
					i++;				
				}
				if (x==V[i])
					snprintf(res, MAXSTRLEN-1, "%e", I[i]);
				else
					snprintf(res, MAXSTRLEN-1, "%e", I[i]+(x-V[i])*(I[i+1]-I[i])/(V[i+1]-V[i]));
				
			}
			free(V);
			free(I);
			
			FreeArgs(args,1);
			return 0;
		
		}
		case MV_VEL_I:
		{
			/* return the voltage for element i, electrode j and simulation index k (i,j,k) */
			char **args;
			int i,j,k;
			
			
			args=GetMeshVarArgs (expr, 3);
			
			i=atoi(args[0]);
			j=atoi(args[1]);
			k=atoi(args[2]);
			
			if ((i>=0)&&(i<M->Nn)&&(j>=0)&&(j<M->Nel)&&(k>=0)&&(k<M->res.Nva))
				snprintf(res, MAXSTRLEN-1, "%e",  M->res.Vn[k][j][i]);
			else
			{
				if ((i<0)||(i>=M->Nn))
					Warning("In EvalMeshVar: Element index %d out of range (0..%d), returning NaN\n", i, M->Nn);				
				if ((j<0)||(j>=M->Nel))
					Warning("In EvalMeshVar: Electrode index %d out of range (0..%d), returning NaN\n", j, M->Nel);				
				if ((k<0)||(k>=M->res.Nva))
					Warning("In EvalMeshVar: simulation index %d out of range (0..%d), returning NaN\n", k, M->res.Nva-1);
				snprintf(res, MAXSTRLEN-1, "NaN");	
			}
				
			FreeArgs(args,3);
			return 0;
		}
		case MV_VVEL_VEL:
		{
		
			/* return the i-th applied voltage */
			char **args;
			double x;
			double *V, *Vel;
			int i,j,k;
			
			args=GetMeshVarArgs (expr, 3);
			i=atoi(args[0]);
			j=atoi(args[1]);
			x=atof(args[2]);
			if (M->res.Nva==0)
			{
				Warning("In EvalMeshVar: no simulated IV data available, cannot estimate applied voltage @ element %d voltage %e V in electrode %d\n", i, x, j);
				*res='\0';			
				FreeArgs(args,1);
				return 0;
			
			}
			
			if (!((i>=0)&&(i<M->Nn)&&(j>=0)&&(j<M->Nel)))
			{
				if ((i<0)||(i>=M->Nn))
					Warning("In EvalMeshVar: Element index %d out of range (0..%d), returning NaN\n", i, M->Nn);				
				if ((j<0)||(j>=M->Nel))
					Warning("In EvalMeshVar: Electrode index %d out of range (0..%d), returning NaN\n", j, M->Nel);
				snprintf(res, MAXSTRLEN-1, "NaN");
			}
			
			if (M->res.Nva==1)
			{
			
				snprintf(res, MAXSTRLEN-1, "%e", M->res.Va[0]);				
				FreeArgs(args,1);
				return 0;
			}
			
			V=malloc((M->res.Nva+1)*sizeof(double));
			Vel=malloc((M->res.Nva+1)*sizeof(double));
			
			for (k=0;k<M->res.Nva;k++)
			{
				V[k]=M->res.Va[k];
				Vel[k]=M->res.Vn[k][j][i];
			}
			
			BubbleSortJV(M->res.Nva, Vel, V);
			if (Vel[0]>x)
				snprintf(res, MAXSTRLEN-1, "%e", V[0]+(x-Vel[0])*(V[1]-V[0])/(Vel[1]-Vel[0]));	
			else
			{
				i=0;
				while(i<M->res.Nva-1)
				{
					if (x==Vel[i])
						break;
					if ((x-Vel[i])*(x-Vel[i+1])<0)
						break;
					i++;				
				}
				if (x==Vel[i])
					snprintf(res, MAXSTRLEN-1, "%e", V[i]);
				else
					snprintf(res, MAXSTRLEN-1, "%e", V[i]+(x-Vel[i])*(V[i+1]-V[i])/(Vel[i+1]-Vel[i]));
				
			}
			free(V);
			free(Vel);
			
			FreeArgs(args,3);
			return 0;
		}
		case MV_V4POINT:
		{
			char **args;
			double x;
			double *V, *Vel;
			int i,j,k,l,n;
			
			args=GetMeshVarArgs (expr, 5);
			i=atoi(args[0]);	/*element i*/
			k=atoi(args[1]);	/*electrode k*/
			j=atoi(args[2]);	/*element j*/
			l=atoi(args[3]);	/*electrode l*/
			x=atof(args[4]);	/*voltage difference*/
			if (M->res.Nva==0)
			{
				Warning("In EvalMeshVar: no simulated IV data available, cannot estimate applied voltage @ element %d voltage %e V in electrode %d\n", i, x, j);
				*res='\0';			
				FreeArgs(args,1);
				return 0;
			
			}
			
			if (!((i>=0)&&(i<M->Nn)&&(j>=0)&&(j<M->Nn)&&(k>=0)&&(k<M->Nel)&&(l>=0)&&(l<M->Nel)))
			{
				if ((i<0)||(i>=M->Nn))
					Warning("In EvalMeshVar: Element index %d out of range (0..%d), returning NaN\n", i, M->Nn);	
				if ((j<0)||(j>=M->Nn))
					Warning("In EvalMeshVar: Element index %d out of range (0..%d), returning NaN\n", j, M->Nn);				
				if ((k<0)||(k>=M->Nel))
					Warning("In EvalMeshVar: Electrode index %d out of range (0..%d), returning NaN\n", k, M->Nel);			
				if ((l<0)||(l>=M->Nel))
					Warning("In EvalMeshVar: Electrode index %d out of range (0..%d), returning NaN\n", l, M->Nel);
				snprintf(res, MAXSTRLEN-1, "NaN");
			}
			
			if (M->res.Nva==1)
			{
			
				snprintf(res, MAXSTRLEN-1, "%e", M->res.Va[0]);				
				FreeArgs(args,1);
				return 0;
			}
			
			V=malloc((M->res.Nva+1)*sizeof(double));
			Vel=malloc((M->res.Nva+1)*sizeof(double));
			
			for (n=0;n<M->res.Nva;n++)
			{
				V[n]=M->res.Va[n];
				Vel[n]=M->res.Vn[n][k][i]-M->res.Vn[n][l][j];
			}
			
			BubbleSortJV(M->res.Nva, Vel, V);
			if (Vel[0]>x)
				snprintf(res, MAXSTRLEN-1, "%e", V[0]+(x-Vel[0])*(V[1]-V[0])/(Vel[1]-Vel[0]));	
			else
			{
				i=0;
				while(i<M->res.Nva-1)
				{
					if (x==Vel[i])
						break;
					if ((x-Vel[i])*(x-Vel[i+1])<0)
						break;
					i++;				
				}
				if (x==Vel[i])
					snprintf(res, MAXSTRLEN-1, "%e", V[i]);
				else
					snprintf(res, MAXSTRLEN-1, "%e", V[i]+(x-Vel[i])*(V[i+1]-V[i])/(Vel[i+1]-Vel[i]));
				
			}
			free(V);
			free(Vel);
			
			FreeArgs(args,5);
			return 0;
		}
		case MV_SURFAREA:
		{
			char **args;
			double S=0;
			int i, P;
			
			args=GetMeshVarArgs (expr, 1);
			
			P=FindProperties(MV->M, args[0]);
			if (P<0)
				snprintf(res, MAXSTRLEN-1, "0.0");
			else
			{
				for (i=0;i<M->Nn;i++)
					if (M->nodes[i].P==P)
						S+=(M->nodes[i].x2-M->nodes[i].x1)*(M->nodes[i].y2-M->nodes[i].y1);
				snprintf(res, MAXSTRLEN-1, "%e",S);				
			}			
			FreeArgs(args,1);
			return 0;
			
		}
		case MV_NONE:
		default:
			Error("In EvalMeshVar: Variable name \"%s\"  not recognized\n",expr);
			return 1;
		
	}

}
/* copies one word from begin to word */
/* begin should point to the next word in a string */
/* the function returns a poiinter to the next word in the string  */
static char *GetWord (char *begin, char *word, meshvar * meshes, int Nm)
{
	char *end, *be;
	int parseExp=1;
	if(!begin)
	{
		*word='\0';
		return NULL;
	}
	end=End(begin);	
	word=strncpy(word,begin,end-begin);
	word[end-begin]='\0';
	begin=Begin(end);
	
	/* check for expressions as part of the word */
	/* we support nested expressions */
	while (parseExp)
	{
		/* we search backwards from the end for the last '[' */
		be=word+strlen(word);
		
		while ((be>=word)&&(*be!='[')&&(*be!='{'))
			be--;
		if (be>=word)
		{
			/* we have a '[', find the matching ']' and substitute the expression in the word with its result */
			char *ee, *expr, *res, *p, *q;
			int fail;
			ee=be;
			/* in the ASCII table the matching bracket is always 2 places further, hence (*be)+2 */
			while ((*ee)&&(*ee!=(*be)+2))
				ee++;
			if (!(*ee))
				Error("No matching \"%c\" found\n", (*be)+2);
			
			expr=malloc(MAXSTRLEN*sizeof(char));
			res=malloc(MAXSTRLEN*sizeof(char));
			expr=strncpy(expr,be+1,ee-be-1);
			expr[ee-be-1]='\0';
			if (*be=='[')
				fail=ExprEval(expr, res);
			else
				fail=EvalMeshVar(expr, res, meshes, Nm);
				
			if (fail)
				Error("Failed to evaluate expression: \"%s\"\n", expr);
			else
			{
				p=res+strlen(res);
				q=ee+1;
				while ((*q)&&(p-res+be-word<MAXSTRLEN))
				{
					*p=*q;
					q++;
					p++;
				}
				*p='\0';
				p=res;
				while (*p)
				{
					*be=*p;
					be++;
					p++;
				}
				*be='\0';
			}
			free(expr);
			free(res);
		}
		else
			/* we have no '[', no more expression in the word */
			parseExp=0;
	}
	return begin;
}

static char *GetWord_noexp (char *begin, char *word)
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


static char **GetArgs (char **begin, int Nw, meshvar * meshes, int Nm)
{
	char ** res;
	int i;
	res=malloc((Nw+1)*sizeof(char *));
	for (i=0;i<Nw;i++)
	{
		res[i]=malloc((MAXSTRLEN)*sizeof(char));
		(*begin)=GetWord ((*begin), res[i], meshes, Nm);
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




#define LOOPALLOCBLOCK 5
/* the parse routine, parses the commands in a file */
void Parse (char *file)
{
	FILE *f;
	char *line, *word, c;
	char *begin;
	meshvar * Meshes;
	int Nm=0, k, ExitPVMOS=1;
	int line_nr=1;
	int *loop_line_nr;
	long int *loop_stack, last_pos; 
	int loop=0, loop_allocated=LOOPALLOCBLOCK;
	clock_t tic, toc;
	PRSDEF key;
	polygon P;
	
	P.N=0;
	
	if ((f=fopen(file,"rb"))==NULL)
		Error("In Parse: Cannot open file %s\n", file);
	

	line=malloc(MAXSTRLEN*sizeof(char));
	word=malloc(MAXSTRLEN*sizeof(char));
	
	Meshes=malloc((Nm+1)*sizeof(meshvar));
	
	loop_stack=malloc(loop_allocated*sizeof(long int));
	loop_line_nr=malloc(loop_allocated*sizeof(int));
	
	last_pos=ftell(f);
    	fgets(line, MAXSTRLEN-1, f);
	
	tic = clock();
	while((feof(f)==0)&&ExitPVMOS)
	{
    		k=sscanf(line, " %c", &c);
		if((k!=-1)&&(c!='#'))
		{
			/* read a word, separate at whitespace, check what the word means then read its value or execute */
			begin=Begin(line);
			while(begin)
			{
				begin=GetWord (begin, word, Meshes, Nm);
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
						args=GetArgs (&begin, 7, Meshes, Nm);
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
					case ADDCOL_E:
					{
						double dx;
						char **args;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Adding column to the right side of mesh %s", line_nr,args[1]);	
						
														
						dx=atof(args[0]);
						if (dx<0)
							Error("* line %3d: Negative length is not allowed\n",line_nr);						
						MV=LookupMesh(args[1], Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						AddColEast(&(MV->M), dx);
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 2);
						break;
					}
					case ADDROW_N:
					{
						double dy;
						char **args;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Adding row to the top side of mesh %s", line_nr,args[1]);	
						
														
						dy=atof(args[0]);
						if (dy<0)
							Error("* line %3d: Negative length is not allowed\n",line_nr);							
						MV=LookupMesh(args[1], Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						AddRowNorth(&(MV->M), dy);
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 2);
						break;
					}
					case ADDCOL_W:
					{
						double dx;
						char **args;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Adding column to the left side of mesh %s", line_nr,args[1]);	
						
														
						dx=atof(args[0]);
						if (dx<0)
							Error("* line %3d: Negative length is not allowed\n",line_nr);							
						MV=LookupMesh(args[1], Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						AddColWest(&(MV->M), dx);
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 2);
						break;
					}
					case ADDROW_S:
					{
						double dy;
						char **args;
						meshvar *MV;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Adding row to the bottom side of mesh %s", line_nr,args[1]);	
						
														
						dy=atof(args[0]);
						if (dy<0)
							Error("* line %3d: Negative length is not allowed\n",line_nr);							
						MV=LookupMesh(args[1], Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						AddRowSouth(&(MV->M), dy);
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 2);
						break;
					}
					case RMMESH:
					{
						char *name;
						int len;
						/* read next word and process, if no next word trow an error */
						
						begin=GetWord (begin, word, Meshes, Nm);
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
						meshvar *MV1, *MV2;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 5, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Joining mesh %s and mesh %s to %s", line_nr,args[2],args[3],args[4]);	
						
														
						xoff=atof(args[0]);								
						yoff=atof(args[1]);						
						MV1=LookupMesh(args[2], Meshes, Nm);
						if (!MV1)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[2]);
						MV2=LookupMesh(args[3], Meshes, Nm);
						if (!MV2)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[3]);	
						
						len=strlen(args[4]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[4],len+1);
						AddMeshVar (JoinMeshes(MV1->M, MV2->M, xoff, yoff), &name,  &Meshes, &Nm);	
						MV1=LookupMesh(name, Meshes, Nm);						
						GetMeshBB(&(MV1->M), &(MV1->x1), &(MV1->y1), &(MV1->x2), &(MV1->y2));
						FreeArgs (args, 5);
						break;
					}
					case JOINMESH_H:
					{
						double yoff;
						char *name;
						char **args;
						int len;
						meshvar *MV1, *MV2;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: Joining mesh %s and mesh %s to %s", line_nr,args[1],args[2],args[3]);	
						
											
						yoff=atof(args[0]);					
						MV1=LookupMesh(args[1], Meshes, Nm);
						if (!MV1)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						MV2=LookupMesh(args[2], Meshes, Nm);
						if (!MV2)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[2]);	
						
						
						len=strlen(args[3]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[3],len+1);
						AddMeshVar (JoinMeshes_H(MV1->M, MV2->M, yoff), &name,  &Meshes, &Nm);	
						MV1=LookupMesh(name, Meshes, Nm);						
						GetMeshBB(&(MV1->M), &(MV1->x1), &(MV1->y1), &(MV1->x2), &(MV1->y2));
						FreeArgs (args, 4);
						break;
					
					}
					case JOINMESH_V:
					{
						double xoff;
						char *name;
						char **args;
						int len;
						meshvar *MV1, *MV2;
						/* read next word and process, if no next word trow an error */
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;	
						Print(NORMAL,"* line %3d: Joining mesh %s and mesh %s to %s", line_nr,args[1],args[2],args[3]);	
						
										
						xoff=atof(args[0]);				
						MV1=LookupMesh(args[1], Meshes, Nm);
						if (!MV1)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[1]);
						MV2=LookupMesh(args[2], Meshes, Nm);
						if (!MV2)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[2]);	
						
						len=strlen(args[3]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[3],len+1);
						AddMeshVar (JoinMeshes_V(MV1->M, MV2->M, xoff), &name,  &Meshes, &Nm);	
						MV1=LookupMesh(name, Meshes, Nm);						
						GetMeshBB(&(MV1->M), &(MV1->x1), &(MV1->y1), &(MV1->x2), &(MV1->y2));
						FreeArgs (args, 4);
						break;
					
					}
					case DUPMESH:
					{
						int len;
						char *name;
						mesh Mnew;
						char **args;
						meshvar *MV;
						
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Duplicating mesh %s and storing duplicate in %s",line_nr, args[0], args[1]);
						
						MV=LookupMesh (args[0],  Meshes, Nm);							
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);											
						Mnew=DuplicateMesh(MV->M);						
						len=strlen(args[1]);
						name=malloc((len+2)*sizeof(char));
						name=strncpy(name,args[1],len+1);						
						AddMeshVar (Mnew, &name,  &Meshes, &Nm);
						MV=LookupMesh(name, Meshes, Nm);						
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 2);
						break;					
					}
					case CLEANMESH:
					{
						int nc;
						char **args;
						meshvar *MV;
						
						args=GetArgs (&begin, 1, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Cleaning up unused areas in mesh %s",line_nr, args[0]);
						
						MV=LookupMesh (args[0],  Meshes, Nm);							
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);											
						nc=DeleteUnusedProperties(&(MV->M));
						Print(NORMAL,"            ---> %d areas removed, mesh now contains %d areas",nc, MV->M.Na);	
						FreeArgs (args, 1);
						break;					
					}
					case ADDEL:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
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
					case PURGERES:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Purging all simulated results in mesh %s",line_nr, word);
													
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr, word);		
						
						PurgeResults(&(MV->M));			
						break;		
					}
					case PURGERESI:
					{
						int Vai;
						double Va=0;
						char **args;
						meshvar *MV;
						
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
						MV=LookupMesh (args[0],  Meshes, Nm);							
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);	
								
						if (args[1][0]=='i')
						{
							Vai=atoi(args[1]+1);
							if ((Vai>=0)&&(Vai<MV->M.res.Nva))
								Va=MV->M.res.Va[Vai];								
						}
						else
						{	
							Va=atof(args[1]);
							Vai=FindVa(Va, MV->M.res.Va, MV->M.res.Nva);
							Va=MV->M.res.Va[Vai];	
						}
						if ((Vai<0)||(Vai>=MV->M.res.Nva))
						{
							Warning("* line %3d: Result index %d not available, nothing to purge\n",line_nr, Vai);
						}
						else
						{
							Print(NORMAL, "* line %3d: Purging result for %e V at index %d in mesh %s",line_nr, Va, Vai, args[0]);
							PurgeResultAtIndex(&(MV->M),Vai);
						}
						FreeArgs (args, 2);
						break;					
					}
					case SPLITX:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all elements in %s in x-direction",line_nr, word);
							fflush(stdout);
							SplitMeshX(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected elements in %s in x-direction",line_nr, word);
							fflush(stdout);
							SplitListX(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);				
						break;
					}
					case SPLITY:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all elements in %s in y-direction",line_nr, word);
							fflush(stdout);
							SplitMeshY(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected elements in %s in y-direction",line_nr, word);
							fflush(stdout);
							SplitListY(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);					
						break;
					}
					case SPLITXY:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all elements in %s in x- and y-direction",line_nr, word);
							fflush(stdout);
							SplitMeshXY(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected elements in %s in x- and y-direction",line_nr, word);
							fflush(stdout);
							SplitListXY(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);					
						break;
					}
					case SPLITLONG:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);		
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all elements in %s in the longest direction",line_nr, word);
							fflush(stdout);
							SplitMeshLong(&(MV->M));
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected elements in %s in the longest direction",line_nr, word);
							fflush(stdout);
							SplitListLong(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
						}
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);
						break;
					}
					case SPLITCOARSE:
					{
						meshvar *MV;
						double dx, dy;
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
							
						dx=atof(args[1]);
						dy=atof(args[2]);
											
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Splitting all elements in %s until all edges shorter than %e in x-direction and %e in y-direction",line_nr, MV->name, dx, dy);
							fflush(stdout);
							SplitMeshWhileCoarse(&(MV->M), dx, dy);
						}
						else
						{
							Print(NORMAL,"* line %3d: Splitting selected elements in %s until all edges shorter than %e in x-direction and %e in y-direction",line_nr, MV->name, dx, dy);
							fflush(stdout);
							SplitListWhileCoarse(&(MV->M), MV->nodes, dx, dy);
							MV->nodes[0]=0;
						}	
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);
						FreeArgs (args, 3);											
						break;
					}
					case MESH_SCALEMOVE:
					{	
						meshvar *MV;
						char **args;
						double fx,fy, dx, dy;
						args=GetArgs (&begin, 5, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
							
						fx=atof(args[1]);
						fy=atof(args[2]);
						dx=atof(args[3]);
						dy=atof(args[4]);
						Print(NORMAL,"* line %3d: Scaling and moving mesh %s, scaling: fx=%e, fy=%e moving: dx=%e dy=%e",line_nr, MV->name, fx,fy, dx,dy);
						Mesh_ScaleMove(&(MV->M), fx, fy, dx, dy);
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 5);										
						break;
					}
					case MESH_ROTATE:
					{	
						meshvar *MV;
						char **args;
						double x,y;
						int d;
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
							
						x=atof(args[1]);
						y=atof(args[2]);
						d=atoi(args[3]);
						if (d%90)
							Error("* line %3d: requested rotation over %d degrees is not a multiple of 90",line_nr,d);	
							
						Print(NORMAL,"* line %3d: Rotating mesh %s over %d degrees around point (%e, %e)",line_nr, MV->name, d, x,y);
						RotateMesh(&(MV->M), x, y, d);
						GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						FreeArgs (args, 4);												
						break;
					}
					case MESH_FLIPX:
					{		
						meshvar *MV;
						char **args;
						double xx1, yy1, xx2, yy2;
						
						args=GetArgs (&begin, 1, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
							
						Print(NORMAL,"* line %3d: Flip x-coordinates for mesh %s",line_nr, MV->name);
						GetMeshBB(&(MV->M), &xx1, &yy1, &xx2, &yy2);
						Mesh_ScaleMove(&(MV->M), -1.0, 1.0, (xx1+xx2), 0);
						FreeArgs (args, 1);																				
						break;
					}
					case MESH_FLIPY:
					{		
						meshvar *MV;
						char **args;
						double xx1, yy1, xx2, yy2;
						args=GetArgs (&begin, 1, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						Print(NORMAL,"* line %3d: Flip y-coordinates for mesh %s",line_nr, MV->name);
						GetMeshBB(&(MV->M), &xx1, &yy1, &xx2, &yy2);
						Mesh_ScaleMove(&(MV->M), 1.0, -1.0, 0, (yy1+yy2));
						FreeArgs (args, 1);												
						break;
					}
					case MESH_SETBB:
					{	
						meshvar *MV;
						char **args;
						double x1,x2,y1,y2;
						double fx, fy, dx, dy;
						int FixR;
						args=GetArgs (&begin, 6, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);		
						x1=atof(args[1]);			
						y1=atof(args[2]);		
						x2=atof(args[3]);		
						y2=atof(args[4]);
						Print(NORMAL,"* line %3d: Fitting mesh %s to bounding box (%e,%e) (%e,%e)",line_nr, MV->name, x1,y1,x2,y2);
						FixR=atoi(args[5]);
						if (FixR)	
							Print(NORMAL,"*           Preserving ascpect ratio");
						else
							Print(NORMAL,"*           Not preserving ascpect ratio");
						
						SetMeshBB(&(MV->M), x1, y1, x2, y2, FixR, &fx, &fy, &dx, &dy);
						MV->x1=x1;
						MV->y1=y1;
						MV->x2=x2;
						MV->y2=y2;
						Print(NORMAL,"* line %3d: Mesh scaled by fx=%e, fy=%e",line_nr, fx,fy);
						Print(NORMAL,"* line %3d: Mesh translated by dx=%e, dy=%e",line_nr, dx,dy);
						Print(NORMAL,"* line %3d: Variables Fx, Fy, Dx, and Dy set accordingly",line_nr);
						DefineVar("Fx", fx);
						DefineVar("Fy", fy);
						DefineVar("Dx", dx);
						DefineVar("Dy", dy);						
						FreeArgs (args, 6);												
						break;
					}
					case POLY_SCALEMOVE:
					{	
						char **args;
						double fx,fy, dx, dy;
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
											
						if (!P.N)
							Warning("* line %3d: No polygon defined\n",line_nr);		
							
						fx=atof(args[0]);
						fy=atof(args[1]);
						dx=atof(args[2]);
						dy=atof(args[3]);
						Print(NORMAL,"* line %3d: Scaling and moving polygon, scaling: fx=%e, fy=%e moving: dx=%e dy=%e",line_nr, fx,fy, dx,dy);
						Poly_ScaleMove(P, fx, fy, dx, dy);
						FreeArgs (args, 4);										
						break;
					}
					case POLY_ROTATE:
					{	
						char **args;
						double x,y;
						int d;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
											
						if (!P.N)
							Warning("* line %3d: No polygon defined\n",line_nr);	
							
						x=atof(args[0]);
						y=atof(args[1]);
						d=atoi(args[2]);
							
						Print(NORMAL,"* line %3d: Rotating polygon over %d degrees around point (%e, %e)",line_nr,d, x,y);
						Poly_Rotate(P, x, y, d);
						FreeArgs (args, 3);												
						break;
					}
					case POLY_FLIPX:
					{											
						if (!P.N)
							Warning("* line %3d: No polygon defined\n",line_nr);	
						else	
						{
							Print(NORMAL,"* line %3d: Flip polygon x-coordinates",line_nr);
							Poly_FlipX(P);
						}														
						break;
					}
					case POLY_FLIPY:
					{												
						if (!P.N)
							Warning("* line %3d: No polygon defined\n",line_nr);	
						else	
						{
							Print(NORMAL,"* line %3d: Flip polygon y-coordinates",line_nr);
							Poly_FlipY(P);
						}														
						break;	
					}
					case POLY_SETBB:
					{	
						char **args;
						double x1,x2,y1,y2;
						double fx, fy, dx, dy;
						int FixR;
						args=GetArgs (&begin, 5, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
											
						if (!P.N)
						{
							Warning("* line %3d: No polygon defined\n",line_nr);						
							FreeArgs (args, 5);	
							break;
						}
						x1=atof(args[0]);			
						y1=atof(args[1]);		
						x2=atof(args[2]);		
						y2=atof(args[3]);
						Print(NORMAL,"* line %3d: Fitting polygon to bounding box (%e,%e) (%e,%e)",line_nr, x1,y1,x2,y2);
						FixR=atoi(args[4]);
						if (FixR)	
							Print(NORMAL,"*           Preserving ascpect ratio");
						else
							Print(NORMAL,"*           Not preserving ascpect ratio");
						
						Poly_SetBoundingBox(P, x1, y1, x2, y2, FixR, &fx, &fy, &dx, &dy);
						Print(NORMAL,"* line %3d: Polygon scaled by fx=%e, fy=%e",line_nr, fx,fy);
						Print(NORMAL,"* line %3d: Polygon translated by dx=%e, dy=%e",line_nr, dx,dy);
						Print(NORMAL,"* line %3d: Variables Fx, Fy, Dx, and Dy set accordingly",line_nr);
						DefineVar("Fx", fx);
						DefineVar("Fy", fy);
						DefineVar("Dx", dx);
						DefineVar("Dy", dy);						
						FreeArgs (args, 5);												
						break;
					}
					case RESOLVPOLY:
					{
						meshvar *MV;
						double d;
						int loop;
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);		
										
						if (P.N)
						{
							d=atof(args[1]);
							loop=atoi(args[2]);
						
							Print(NORMAL,"* line %3d: Resolving polygon in mesh %s  down to %e",line_nr, MV->name, d);
							ResolvContour(P, &(MV->M), loop, d);
							Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);
						}
						else
							Print(NORMAL,"* line %3d: Cannot resolve polygon in mesh %s, no polygon defined",line_nr, MV->name);
						
						FreeArgs (args, 3);											
						break;
					}
					case RESOLVCIRC:
					{
						meshvar *MV;
						double d,x,y,r;
						char **args;
						args=GetArgs (&begin, 5, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);		
							
						d=atof(args[1]);
						x=atof(args[2]);
						y=atof(args[3]);
						r=atof(args[4]);
						
						Print(NORMAL,"* line %3d: Resolving circle in mesh %s  down to %e",line_nr, MV->name, d);
						ResolvCircle(x,y,r, &(MV->M), d);
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);
						FreeArgs (args, 5);											
						break;
					}
					case RESOLVRECT:
					{
						meshvar *MV;
						double d,x1,y1,x2, y2;
						char **args;
						args=GetArgs (&begin, 6, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);		
							
						d=atof(args[1]);
						x1=atof(args[2]);
						y1=atof(args[3]);
						x2=atof(args[4]);
						y2=atof(args[5]);
						
						Print(NORMAL,"* line %3d: Resolving rectangle in mesh %s  down to %e",line_nr, MV->name, d);
						ResolvRect(x1,y1,x2, y2, &(MV->M), d);
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);
						FreeArgs (args, 6);											
						break;
					}
					case SIMPLIFY_MESH:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);							
						Print(NORMAL,"* line %3d: Simplifying mesh of %d elements",line_nr, MV->M.Nn);
						fflush(stdout);		
						Chunkify(&(MV->M));
						MV->nodes[0]=0;
						Print(NORMAL,"            ---> Mesh %s consists of %d elements",MV->name, MV->M.Nn);	
						break;
					}
					case LOADMESH:
					{
						int len;
						char *name;
						mesh Mnew;
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
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
						meshvar *MV;						
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL, "* line %3d: Saving mesh %s to file %s",line_nr,args[0], args[1]);
														
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						WriteMesh(args[1],&(MV->M));
						FreeArgs (args, 2);						
						break;
					}
					/********************************* Secion Output */
					/********************************* SubSecion mesh data */
					case PRINTMESH:
					{
						meshvar *MV;						
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
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
						args=GetArgs (&begin, 2, Meshes, Nm);
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
						args=GetArgs (&begin, 2, Meshes, Nm);
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
						args=GetArgs (&begin, 2, Meshes, Nm);
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
					case PRINTSOLPAR:
					{
						meshvar *MV;			
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print solar cell parameters for mesh %s to file %s",line_nr,args[0], args[1]);
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						PrintSolPar(args[1],&(MV->M));
						FreeArgs (args, 2);	
						break;
					}
					case PRINTIV:
					{
						meshvar *MV;			
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print simulated I-V pairs of mesh %s to file %s",line_nr,args[0], args[1]);
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						PrintIV(args[1],&(MV->M));
						FreeArgs (args, 2);	
						break;
					}
					case PRINTPROBE:
					{
						meshvar *MV;			
						char **args;
						double x,y;
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print simulated probe voltages of mesh %s to file %s",line_nr,args[0], args[3]);
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						x=atof(args[1]);
						y=atof(args[2]);						
						PrintProbe(args[3],&(MV->M),x,y);
						FreeArgs (args, 4);	
						break;
					}
					case PRINTINIPIV:
					{
						meshvar *MV;			
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print simulated I-V pairs of mesh %s to file %s",line_nr,args[0], args[1]);
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						PrintInIp(args[1],&(MV->M), MV->nodes);
						FreeArgs (args, 2);	
						break;
					}
					case PRINTPARS:
					{
						meshvar *MV;				
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print parameters per area in mesh %s to file %s",line_nr,args[0], args[1]);
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
								
						PrintPars(args[1], &(MV->M));
						FreeArgs (args, 2);	
						break;
					}
					case PRINTLOCALJV:
					{
						
						meshvar *MV;				
						char **args;
						double x, y, Vstart, Vend;
						int Nstep, el;
						args=GetArgs (&begin, 8, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Print local JV characteristics in mesh %s to file %s",line_nr,args[0], args[7]);
							
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						
						x=atof(args[1]);
						y=atof(args[2]);
						el=atoi(args[3]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
						Vstart=atof(args[4]);
						Vend=atof(args[5]);
						Nstep=atoi(args[6]);
						Print(NORMAL, "* line %3d: Position (%e,%e), inter electrode index %d",line_nr,x,y,el);
						PrintLocalJV(args[7], MV->M, x, y, el, Vstart, Vend, Nstep);
						FreeArgs (args, 8);	
						break;
					}
					case SAVEVARS:
					{
						int j;					
						char **args;
						FILE *varf;
						args=GetArgs (&begin, 1, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						
						if ((varf=fopen(args[0],"w"))==NULL)
							Warning("Warning: Cannot open file %s on line %d, cannot save variables\n", args[0], line_nr);
						else
						{
							for (j=0;j<var_count;j++)
								fprintf(varf, "%s\t=\t%.12e\n", var_names[j], var_values[j]);
							fclose(varf);
						}						
						FreeArgs (args, 1);	
						break;
					}
					case LOADVARS:
					{
						char *vname, *vvalue, *vline; 			
						char **args;
						FILE *varf;
						args=GetArgs (&begin, 1, Meshes, Nm);
						if ((varf=fopen(args[0],"r"))==NULL)
							Warning("Warning: Cannot open file %s on line %d, cannot load variables\n", args[0], line_nr);
						else
						{
							while(feof(varf)==0)
							{
				
								vline=GetLine(varf);
								GetNameValue(vline, &(vname), &(vvalue));
								if (vname && vvalue)
								{
									Print(NORMAL,"* line %3d: defining \"%s=%s\"",line_nr, vname, vvalue);	
									DefineVar(vname, atof(vvalue));
								}
								free(vline);	
							}
							fclose(varf);		
						}
						FreeArgs (args, 1);
						break;	
					}
					case SURFCOLCUR:
					case SURFDCOLCUR:
					{
						meshvar *MV;
						int el, Nx, Ny, RI;
						double Va, x1, x2, y1, y2;					
						char **args;
						args=GetArgs (&begin, 11, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export differential collected Currtent Density in mesh %s to file %s",line_nr,args[0], args[1]);
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						el=atoi(args[2]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
						Va=atof(args[3]);
						
						x1=atof(args[4]);
						y1=atof(args[5]);
						x2=atof(args[6]);
						y2=atof(args[7]);
						Nx=atoi(args[8]);
						Ny=atoi(args[9]);
						RI=atoi(args[10]);

						PrintLocallyCollectedCurrent(args[1], &(MV->M), x1, y1, x2, y2, Nx, Ny, Va, el, (key==SURFDCOLCUR), RI);
						FreeArgs (args, 11);	
						break;
					}
					case SURFVPLOT:
					{
						meshvar *MV;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export potentials from mesh %s to file %s",line_nr,args[0], args[8]);
							
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);
																	
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						if (args[7][0]=='i')
							Vai=atoi(args[7]+1);
						else
						{	
							Va=atof(args[7]);
							Vai=FindVa(Va, MV->M.res.Va, MV->M.res.Nva);
						}
						if ((Vai>=0)&&(Vai<MV->M.res.Nva))
						{
							Print(NORMAL,"            -->  Using simulation at %e V", MV->M.res.Va[Vai]);
							SurfVPlot(args[8], &(MV->M), Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present.\n", line_nr);	
							
						FreeArgs (args, 9);	
						break;
					}
					case SURFPPLOT:
					{
						meshvar *MV;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export power density from mesh %s  to file %s",line_nr,args[0], args[8]);
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);
													
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						if (args[7][0]=='i')
							Vai=atoi(args[7]+1);
						else
						{	
							Va=atof(args[7]);
							Vai=FindVa(Va, MV->M.res.Va, MV->M.res.Nva);
						}
						
						
						if ((Vai>=0)&&(Vai<MV->M.res.Nva))
						{
							Print(NORMAL,"            -->  Using simulation at %e V", MV->M.res.Va[Vai]);
							SurfPPlot(args[8], &(MV->M), Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present\n", line_nr);	
						FreeArgs (args, 9);	
						break;
					}
					case SURFJPLOT:
					{
						meshvar *MV;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export current densities from mesh %s to file %s",line_nr,args[0], args[8]);
							
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);
																	
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						if (args[7][0]=='i')
							Vai=atoi(args[7]+1);
						else
						{	
							Va=atof(args[7]);
							Vai=FindVa(Va, MV->M.res.Va, MV->M.res.Nva);
						}
						
						if ((Vai>=0)&&(Vai<MV->M.res.Nva))
						{
							Print(NORMAL,"            -->  Using simulation at %e V", MV->M.res.Va[Vai]);
							SurfJPlot(args[8], &(MV->M), Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present.\n", line_nr);	
							
						FreeArgs (args, 9);	
						break;
					}
					case SURFVJPLOT:
					{
						meshvar *MV;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export junction voltages from mesh %s to file %s",line_nr,args[0], args[8]);
							
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);
																	
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						if (args[7][0]=='i')
							Vai=atoi(args[7]+1);
						else
						{	
							Va=atof(args[7]);
							Vai=FindVa(Va, MV->M.res.Va, MV->M.res.Nva);
						}
						
						if ((Vai>=0)&&(Vai<MV->M.res.Nva))
						{
							Print(NORMAL,"            -->  Using simulation at %e V", MV->M.res.Va[Vai]);
							SurfVjPlot(args[8], &(MV->M), Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present.\n", line_nr);	
							
						FreeArgs (args, 9);	
						break;
					}
					case SURFEPLOT:
					{
						meshvar *MV;
						double x1, y1, x2, y2, Va;
						int Nx, Ny, Vai;				
						char **args;
						args=GetArgs (&begin, 9, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export electric fields from mesh %s to file %s",line_nr,args[0], args[8]);
							
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);
																	
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						if (args[7][0]=='i')
							Vai=atoi(args[7]+1);
						else
						{	
							Va=atof(args[7]);
							Vai=FindVa(Va, MV->M.res.Va, MV->M.res.Nva);
						}
						if ((Vai>=0)&&(Vai<MV->M.res.Nva))
						{
							Print(NORMAL,"            -->  Using simulation at %e V", MV->M.res.Va[Vai]);
							SurfEPlot(args[8], &(MV->M), Vai, x1, y1, x2, y2, Nx, Ny);
						}
						else
							Warning("\n* line %3d: Warning: no data present.\n", line_nr);	
							
						FreeArgs (args, 9);	
						break;
					}
					case SURFDEFPLOT:
					{
						meshvar *MV;
						double x1, y1, x2, y2;
						int Nx, Ny;				
						char **args;
						args=GetArgs (&begin, 8, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL, "* line %3d: Export area definition from mesh %s to file %s",line_nr,args[0], args[7]);
							
								
						x1=atof(args[1]);
						y1=atof(args[2]);
						x2=atof(args[3]);
						y2=atof(args[4]);
						
						Nx=atoi(args[5]);	
						Ny=atoi(args[6]);
																	
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						SurfDefPlot(args[7], &(MV->M), x1, y1, x2, y2, Nx, Ny);
							
						FreeArgs (args, 8);	
						break;
					}
					/********************************* Secion Node Selection */
					case LOAD_POLY:
					{		
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL, "* line %3d: Load polygon from file %s",line_nr,word);
						if (P.N)
						{
							free(P.x);
							free(P.y);
							free(P.BR);
						}
						P=ReadPoly(word);
						Print(NORMAL, "            -->  Polygon with %d segments loaded",P.N);
						
						break;
					}	
					case DEF_POLY:
					{			
						char **args;
						int Na=50;
						/* check for end of table */
						Print(NORMAL, "* line %3d: Define polygon",line_nr);						
						if (P.N)
						{
							free(P.x);
							free(P.y);
							free(P.BR);
						}	
						P.N=0;
						P.x=malloc(Na*sizeof(double));
						P.y=malloc(Na*sizeof(double));
						P.BR=malloc(LISTBLOCK*sizeof(int));
						P.BR[0]=0;
    						fgets(line, MAXSTRLEN-1, f);
						line_nr++;
						key=NONE;
						while((feof(f)==0)&&(key!=DEF_POLY))
						{
							begin=Begin(line);
							if((begin)&&((*begin)!='#'))
							{
								GetWord (begin, word, Meshes, Nm);
								key=LookupKey (word,  KeyTable);
								if (key==DEF_POLY)
									break;
								args=GetArgs (&begin, 2, Meshes, Nm);
								P.x[P.N]=atof(args[0]);
								P.y[P.N]=atof(args[1]);
								P.N++;
								if (Na-1==P.N)
								{
									Na+=50;
									P.x=realloc(P.x, Na*sizeof(double));	
									P.y=realloc(P.y, Na*sizeof(double));					
								}
								FreeArgs (args, 2);
    							}
							else if (P.N) /* no use breaking an empty polygon. */
								P.BR=AddToList(P.BR, P.N);
							fgets(line, MAXSTRLEN-1, f);
							line_nr++;							
						}
						
						P.BR=AddToList(P.BR, P.N);
						
						if ((feof(f))||(!begin))
							goto premature_end;
												
						if (P.N)
						{
							P.x=realloc(P.x, (P.N+1)*sizeof(double));	
							P.y=realloc(P.y, (P.N+1)*sizeof(double));
						}
						else
						{
							free(P.x);
							free(P.y);
							free(P.BR);
						}
						begin=GetWord (begin, word, Meshes, Nm);
						Print(NORMAL, "            -->  Polygon with %d segments defined",P.N);
						break;
					}		
					case WHILE:
					{			
						char **args;
						int go=0;
						args=GetArgs (&begin, 3, Meshes, Nm);
						Print(NORMAL, "* line %3d: while-loop (%s %s %s)",line_nr, args[0], args[1], args[2]);
						if (strncmp(args[1], ">=",2)==0)
							go=(atof(args[0])>=atof(args[2]));
						else if (strncmp(args[1], "<=",2)==0)
							go=(atof(args[0])<=atof(args[2]));
						else if (strncmp(args[1], "==",2)==0)
						{
							double a,b;
							a=atof(args[0]);
							b=atof(args[2]);
							go=(fabs(a-b)/(fabs(a)+fabs(b)+1e-300)<1e-12);
						}
						else if (strncmp(args[1], "!=",2)==0)
						{
							double a,b;
							a=atof(args[0]);
							b=atof(args[2]);
							go=(fabs(a-b)/(fabs(a)+fabs(b)+1e-300)>1e-12);
						}
						else if (strncmp(args[1], ">",1)==0)
							go=(atof(args[0])>atof(args[2]));
						else if (strncmp(args[1], "<",1)==0)
							go=(atof(args[0])<atof(args[2]));
						else
							Warning("Warning at line %3d: ill formatted conditional (%s) for a while loop\n",line_nr, args[1]);
						
						FreeArgs (args, 3);	
						
						if (go)
						{
							Print(NORMAL, "            -->  Entering loop");							
							/* put filepointers and linenumbers on the stack */
							if (loop+1==LOOPALLOCBLOCK-1)
							{
								loop_allocated+=LOOPALLOCBLOCK;
								loop_stack=realloc(loop_stack, loop_allocated*sizeof(long int));
								loop_line_nr=realloc(loop_line_nr, loop_allocated*sizeof(int));
							}
							loop_line_nr[loop]=line_nr-1;
							loop_stack[loop]=last_pos;
							loop++;						
						}
						else
						{
							/* scan till end of while loop */
							fgets(line, MAXSTRLEN-1, f);
							line_nr++;
							go=1;		
							while((feof(f)==0)&&(go>0))
							{							
								begin=Begin(line);
								while((begin)&&((*begin)!='#'))
								{
									begin=GetWord_noexp (begin, word);
									key=LookupKey (word,  KeyTable);
									if (key==ENDWHILE)
										go--;
									if (key==WHILE)
										go++;
								}
								if (go)
								{
									fgets(line, MAXSTRLEN-1, f);
									line_nr++;
								}			
							}
							Print(NORMAL, "            -->  Exit loop at line %d", line_nr);
						}
						
						break;
					}		
					case ENDWHILE:
					{		
						if (!loop)
						{
							Print(NORMAL, "* line %3d: no while loop to end",line_nr);
							goto premature_end;						
						}
						loop--;		
						fseek(f, loop_stack[loop], SEEK_SET);
						line_nr=loop_line_nr[loop];
						break;
					}
					case IF:
					{			
						char **args;
						int go=0;
						args=GetArgs (&begin, 3, Meshes, Nm);
						Print(NORMAL, "* line %3d: If statement (%s %s %s)",line_nr, args[0], args[1], args[2]);
						if (strncmp(args[1], ">=",2)==0)
							go=(atof(args[0])>=atof(args[2]));
						else if (strncmp(args[1], "<=",2)==0)
							go=(atof(args[0])<=atof(args[2]));
						else if (strncmp(args[1], "==",2)==0)
						{
							double a,b;
							a=atof(args[0]);
							b=atof(args[2]);
							go=(fabs(a-b)/(fabs(a)+fabs(b)+1e-300)<1e-12);
						}
						else if (strncmp(args[1], "!=",2)==0)
						{
							double a,b;
							a=atof(args[0]);
							b=atof(args[2]);
							go=(fabs(a-b)/(fabs(a)+fabs(b)+1e-300)>1e-12);
						}
						else if (strncmp(args[1], ">",1)==0)
							go=(atof(args[0])>atof(args[2]));
						else if (strncmp(args[1], "<",1)==0)
							go=(atof(args[0])<atof(args[2]));
						else
							Warning("Warning at line %3d: ill formatted conditional (%s)\n",line_nr, args[1]);
						
						FreeArgs (args, 3);	
						
						if (!go)
						{
							/* scroll on to the next endif or else */
							fgets(line, MAXSTRLEN-1, f);
							line_nr++;
							go=1;		
							while((feof(f)==0)&&(go>0))
							{							
								begin=Begin(line);
								while((begin)&&((*begin)!='#'))
								{
									begin=GetWord_noexp (begin, word);
									key=LookupKey (word,  KeyTable);
									if (key==ENDIF)
										go--;
									else if (key==IF)
										go++;
									else if ((go==1)&&(key==ELSE))
										go=0;
								}
								if (go)
								{
									fgets(line, MAXSTRLEN-1, f);
									line_nr++;
								}			
							}
							line_nr++;
						}
						
						break;
					}	
					case ELSE:
					{
						int go=1;
						/* scroll on to the next endif */
						fgets(line, MAXSTRLEN-1, f);
						line_nr++;
						while((feof(f)==0)&&(go>0))
						{							
							begin=Begin(line);
							while((begin)&&((*begin)!='#'))
							{
								begin=GetWord_noexp (begin, word);
								key=LookupKey (word,  KeyTable);
								if (key==ENDIF)
									go--;
								else if (key==IF)
									go++;
							}
							if (go)
							{
								fgets(line, MAXSTRLEN-1, f);
								line_nr++;
							}			
						}
						line_nr++;
					}
					case ENDIF:
						break;
					case SELECT_RECT:
					{
						double x1,x2,y1,y2;			
						char **args;
						args=GetArgs (&begin, 5, Meshes, Nm);
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
						args=GetArgs (&begin, 6, Meshes, Nm);
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
						args=GetArgs (&begin, 4, Meshes, Nm);
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
						args=GetArgs (&begin, 5, Meshes, Nm);
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
						begin=GetWord (begin, word, Meshes, Nm);
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
						args=GetArgs (&begin, 3, Meshes, Nm);
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

						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Make area selection of elements in mesh %s", line_nr, word);
						SelectAreaNodes (word,  Meshes, Nm);
						break;
					}	
					case INVERTSELECT:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);	
										
						if (MV->nodes[0]==0)
							Error("* line %3d: No selection to invert\n",line_nr);	
						else
						{
							int *oldlist, i;
							oldlist=DuplicateList(MV->nodes);
							Print(NORMAL,"* line %3d: Inverting selection in mesh %s", line_nr, MV->name);	
							fflush(stdout);
							MV->nodes[0]=0;
							for (i=0;i<MV->M.Nn;i++)
								if (!IsInList(oldlist, MV->M.nodes[i].id))
									MV->nodes=AddToList(MV->nodes, MV->M.nodes[i].id);
						}				
						break;
					}
					case DESELECT:
					{
						meshvar *MV;
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMesh (word,  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);	
						Print(NORMAL,"* line %3d: Deselecting selection in mesh %s", line_nr, MV->name);	
						MV->nodes[0]=0;
						MV->setsel=0;
						break;
					}
					
					/***************** Section Local Properties */
					case ASSIGN_PROP:
					{
						meshvar *MV;
						int P;
						
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;								
						MV=LookupMeshArea (word,  Meshes, Nm, &P);
						
						if (P<0)
							Error("* line %3d: Area %s is not defined\n", line_nr, word);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Assigning all elements in the mesh to area %s",line_nr, word);
							fflush(stdout);
							AssignPropertiesMesh(&(MV->M), P);
						}
						else
						{
							Print(NORMAL,"* line %3d: Assigning selected elements to area %s",line_nr, word);
							fflush(stdout);
							AssignProperties(&(MV->M), MV->nodes, P);
						}					
						break;
					}
					case DELETE_ELEMENTS:
					{
						meshvar *MV;
						
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;				
						MV=LookupMesh (word,  Meshes, Nm);
							
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,word);	
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							Print(NORMAL,"* line %3d: Cannot delete elements, no elements selected in mesh %s",line_nr, word);
						}
						else
						{
							Print(NORMAL,"* line %3d: deleting selected elements in mesh %s",line_nr, word);
							fflush(stdout);
							CleanUpMesh(&(MV->M), MV->nodes);
							MV->nodes[0]=0;
							GetMeshBB(&(MV->M), &(MV->x1), &(MV->y1), &(MV->x2), &(MV->y2));
						}					
						break;
					}
					case SET_REL:
					{
						meshvar *MV;
						int P, el;
						double R;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
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
					case SET_SEL_REL:
					{
						
						meshvar *MV;
						int el;
						double R;
						char *area;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	

							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);					
						R=atof(args[2]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting R for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[a].name, args[2]);
								MV->M.P[a].Rel[el]=R;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting R for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[AreaList[a]].name, args[2]);
								MV->M.P[AreaList[a]].Rel[el]=R;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 3);					
						break;
					}
					case SET_RVP:
					{
						meshvar *MV;
						int P, el;
						double R;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
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
					case SET_SEL_RVP:
					{
						
						meshvar *MV;
						int el;
						double R;
						char *area;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);					
						R=atof(args[2]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting Rvp for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[a].name, args[2]);
								MV->M.P[a].Rvp[el]=R;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting Rvp for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[AreaList[a]].name, args[2]);
								MV->M.P[AreaList[a]].Rvp[el]=R;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 3);					
						break;
					}
					case SET_RVN:
					{
						meshvar *MV;
						int P, el;
						double R;	
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
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
					case SET_SEL_RVN:
					{
						
						meshvar *MV;
						int el;
						double R;
						char *area;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);					
						R=atof(args[2]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting Rvn for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[a].name, args[2]);
								MV->M.P[a].Rvn[el]=R;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting Rvn for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[AreaList[a]].name, args[2]);
								MV->M.P[AreaList[a]].Rvn[el]=R;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 3);					
						break;
					}
					case SET_RVG:
					{
						meshvar *MV;
						int P, el;
						double R;	
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
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
							
						Print(NORMAL,"* line %3d: Setting Rvg for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[P].name, args[2]);
														
						R=atof(args[2]);
						MV->M.P[P].Rvg[el]=R;
						FreeArgs (args, 3);	
						break;
					}
					case SET_SEL_RVG:
					{
						
						meshvar *MV;
						int el;
						double R;
						char *area;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);					
						R=atof(args[2]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting Rvg for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[a].name, args[2]);
								MV->M.P[a].Rvg[el]=R;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting Rvg for electrode %d in %s.%s to %s",line_nr, el, MV->name, MV->M.P[AreaList[a]].name, args[2]);
								MV->M.P[AreaList[a]].Rvg[el]=R;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 3);					
						break;
					}
					case SET_JV:
					{
						meshvar *MV;
						int P, el;
						polygon JV;							
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
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
						if (MV->M.P[P].conn[el].ParStruct)
							free(MV->M.P[P].conn[el].ParStruct);
						MV->M.P[P].conn[el].ParSize=0;
						MV->M.P[P].conn[el].ParStruct=NULL;
						free(JV.BR); /* note that the allocated x and y arrays are now in use and should not be freed*/
						FreeArgs (args, 3);						
						break;
					}
					case SET_SEL_JV:
					{
						
						meshvar *MV;
						int el;
						polygon JV;
						char *area;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
						
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter-electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
							
						JV=ReadPoly(args[2]);
						if (JV.N<2)
							Error("JV characteristic from file %s is too short\nAt least two voltage current pairs are requires\n", word);
						/* ensure monotoneous increasing voltages */
						BubbleSortJV(JV.N, JV.x, JV.y);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								int ii;	
								Print(NORMAL,"* line %3d: Setting tabular JV between electrodes %d and %d in %s.%s to data from file %s",line_nr, el, el+1, MV->name, MV->M.P[a].name, args[2]);
								free(MV->M.P[a].conn[el].V);
								free(MV->M.P[a].conn[el].J);
								
								MV->M.P[a].conn[el].V=malloc((JV.N+1)*sizeof(double));
								MV->M.P[a].conn[el].J=malloc((JV.N+1)*sizeof(double));
								for (ii=0;ii<JV.N;ii++)
								{
									MV->M.P[a].conn[el].V[ii]=JV.x[ii];
									MV->M.P[a].conn[el].J[ii]=JV.y[ii];
								}
								MV->M.P[a].conn[el].N=JV.N;
								MV->M.P[a].conn[el].model=JVD;	
								if (MV->M.P[a].conn[el].ParStruct)
									free(MV->M.P[a].conn[el].ParStruct);
								MV->M.P[a].conn[el].ParSize=0;
								MV->M.P[a].conn[el].ParStruct=NULL;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{		
								int ii;					
								Print(NORMAL,"* line %3d: Setting tabular JV between electrodes %d and %d in %s.%s to data from file %s",line_nr, el, el+1, MV->name, MV->M.P[AreaList[a]].name, args[2]);
								free(MV->M.P[AreaList[a]].conn[el].V);
								free(MV->M.P[AreaList[a]].conn[el].J);
								MV->M.P[AreaList[a]].conn[el].V=malloc((JV.N+1)*sizeof(double));
								MV->M.P[AreaList[a]].conn[el].J=malloc((JV.N+1)*sizeof(double));
								for (ii=0;ii<JV.N;ii++)
								{
									MV->M.P[AreaList[a]].conn[el].V[ii]=JV.x[ii];
									MV->M.P[AreaList[a]].conn[el].J[ii]=JV.y[ii];
								}
								MV->M.P[AreaList[a]].conn[el].N=JV.N;
								MV->M.P[AreaList[a]].conn[el].model=JVD;
								if (MV->M.P[AreaList[a]].conn[el].ParStruct)
									free(MV->M.P[AreaList[a]].conn[el].ParStruct);
								MV->M.P[AreaList[a]].conn[el].ParSize=0;	
								MV->M.P[AreaList[a]].conn[el].ParStruct=NULL;
							}							
							free(newarea);
							free(AreaList);
						}
						free(JV.BR);	
						FreeArgs (args, 3);					
						break;
					}
					case SET_2DJV:
					{
						meshvar *MV;
						int P, el;						
						char **args;
						args=GetArgs (&begin, 8, Meshes, Nm);
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
						
						if (MV->M.P[P].conn[el].ParStruct)
							free(MV->M.P[P].conn[el].ParStruct);
						MV->M.P[P].conn[el].ParStruct=InitOneTwoDiodeStruct(&(MV->M.P[P].conn[el].ParSize));
						
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->J01=atof(args[2]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->J02=atof(args[3]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Jph=atof(args[4]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Rs=atof(args[5]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Rsh=atof(args[6]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Eg=atof(args[7]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->nid1=1.0;
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->nid2=2.0;							
						MV->M.P[P].conn[el].model=TWOD;	
						FreeArgs (args, 8);		
						break;
					}
					case SET_SEL_2DJV:
					{
						meshvar *MV;
						int el;
						char *area;		
						char **args;
						double J01,J02,Jph,Rs,Rsh,Eg;
						args=GetArgs (&begin, 8, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);
						J01=atof(args[2]);
						J02=atof(args[3]);
						Jph=atof(args[4]);
						Rs=atof(args[5]);
						Rsh=atof(args[6]);
						Eg=atof(args[7]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting 2 diode model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[a].name);								
								if (MV->M.P[a].conn[el].ParStruct)
									free(MV->M.P[a].conn[el].ParStruct);
								MV->M.P[a].conn[el].ParStruct=InitOneTwoDiodeStruct(&(MV->M.P[a].conn[el].ParSize));
						
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->J01=J01;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->J02=J02;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Jph=Jph;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Rs=Rs;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Rsh=Rsh;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Eg=Eg;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->nid1=1.0;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->nid2=2.0;						
								MV->M.P[a].conn[el].model=TWOD;						
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
									MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting 2 diode model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[AreaList[a]].name);
								if (MV->M.P[AreaList[a]].conn[el].ParStruct)
									free(MV->M.P[AreaList[a]].conn[el].ParStruct);
								MV->M.P[AreaList[a]].conn[el].ParStruct=InitOneTwoDiodeStruct(&(MV->M.P[AreaList[a]].conn[el].ParSize));
						
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->J01=J01;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->J02=J02;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jph=Jph;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Rs=Rs;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Rsh=Rsh;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Eg=Eg;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->nid1=1.0;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->nid2=2.0;			
								MV->M.P[AreaList[a]].conn[el].model=TWOD;	
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 8);					
						break;
					}
					case SET_1DJV:
					{
						meshvar *MV;
						int P, el;						
						char **args;
						args=GetArgs (&begin, 8, Meshes, Nm);
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
						
						if (MV->M.P[P].conn[el].ParStruct)
							free(MV->M.P[P].conn[el].ParStruct);
						MV->M.P[P].conn[el].ParStruct=InitOneTwoDiodeStruct(&(MV->M.P[P].conn[el].ParSize));
						
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->J01=atof(args[2]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->nid1=atof(args[3]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Jph=atof(args[4]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Rs=atof(args[5]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Rsh=atof(args[6]);
						((OneTwoDiode *)(MV->M.P[P].conn[el].ParStruct))->Eg=atof(args[7]);
						MV->M.P[P].conn[el].model=ONED;	
						FreeArgs (args, 8);								
						break;
					}
					case SET_SEL_1DJV:
					{
						meshvar *MV;
						int el;
						char *area;		
						char **args;
						double J01,nid1,Jph,Rs,Rsh,Eg;
						args=GetArgs (&begin, 8, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);	
						J01=atof(args[2]);
						nid1=atof(args[3]);
						Jph=atof(args[4]);
						Rs=atof(args[5]);
						Rsh=atof(args[6]);
						Eg=atof(args[7]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting 1 diode model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[a].name);
								if (MV->M.P[a].conn[el].ParStruct)
									free(MV->M.P[a].conn[el].ParStruct);
								MV->M.P[a].conn[el].ParStruct=InitOneTwoDiodeStruct(&(MV->M.P[a].conn[el].ParSize));
						
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->J01=J01;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->nid1=nid1;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Jph=Jph;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Rs=Rs;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Rsh=Rsh;
								((OneTwoDiode *)(MV->M.P[a].conn[el].ParStruct))->Eg=Eg;	
								MV->M.P[a].conn[el].model=ONED;	
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting 1 diode model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[AreaList[a]].name);
								if (MV->M.P[AreaList[a]].conn[el].ParStruct)
									free(MV->M.P[AreaList[a]].conn[el].ParStruct);
								MV->M.P[AreaList[a]].conn[el].ParStruct=InitOneTwoDiodeStruct(&(MV->M.P[AreaList[a]].conn[el].ParSize));
						
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->J01=J01;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->nid1=nid1;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jph=Jph;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Rs=Rs;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Rsh=Rsh;
								((OneTwoDiode *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Eg=Eg;		
								MV->M.P[AreaList[a]].conn[el].model=ONED;	
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 8);					
						break;
					}
					case SET_PTJV:
					{
						meshvar *MV;
						int P, el;						
						char **args;
						args=GetArgs (&begin, 11, Meshes, Nm);
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
							
						Print(NORMAL,"* line %3d: Setting photo-transistor model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[P].name);
						
						if (MV->M.P[P].conn[el].ParStruct)
							free(MV->M.P[P].conn[el].ParStruct);
						MV->M.P[P].conn[el].ParStruct=InitPhotoTransistorStruct(&(MV->M.P[P].conn[el].ParSize));
						
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Jsbe=atof(args[2]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Jsbc=atof(args[3]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Jph=atof(args[4]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Rs=atof(args[5]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Rsh=atof(args[6]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->EgBE=atof(args[7]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->PhiBC=atof(args[8]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Bf=atof(args[9]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Vaf=atof(args[10]);
						MV->M.P[P].conn[el].model=PHOTOT;	
						FreeArgs (args, 11);								
						break;
					}
					case SET_PTEXTJV:
					{
						meshvar *MV;
						int P, el;						
						char **args;
						args=GetArgs (&begin, 12, Meshes, Nm);
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
						if (MV->M.P[P].conn[el].model!=PHOTOT)							
							Error("* line %3d: Setting extra paremeters for the phototransistor model requires the model to be set first\n", line_nr);
							
						Print(NORMAL,"* line %3d: Setting photo-transistor model extra parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[P].name);
						
						
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Jse=atof(args[2]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Jsc=atof(args[3]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Nf=atof(args[4]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Nr=atof(args[5]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Ne=atof(args[6]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Nc=atof(args[7]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->Var=atof(args[8]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->XTIBE=atof(args[9]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->XTIBC=atof(args[10]);
						((PhotoTransistor *)(MV->M.P[P].conn[el].ParStruct))->XTB=atof(args[11]);
						FreeArgs (args, 12);								
						break;
					}
					case SET_SEL_PTJV:
					{
						meshvar *MV;
						int el;
						char *area;		
						char **args;
						args=GetArgs (&begin, 11, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);	
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting photo-transistor model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[a].name);
								if (MV->M.P[a].conn[el].ParStruct)
									free(MV->M.P[a].conn[el].ParStruct);
								MV->M.P[a].conn[el].ParStruct=InitPhotoTransistorStruct(&(MV->M.P[a].conn[el].ParSize));
						
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Jsbe=atof(args[2]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Jsbc=atof(args[3]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Jph=atof(args[4]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Rs=atof(args[5]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Rsh=atof(args[6]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->EgBE=atof(args[7]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->PhiBC=atof(args[8]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Bf=atof(args[9]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Vaf=atof(args[10]);
								MV->M.P[a].conn[el].model=PHOTOT;								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting photo-transistor model parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[AreaList[a]].name);
								if (MV->M.P[AreaList[a]].conn[el].ParStruct)
									free(MV->M.P[AreaList[a]].conn[el].ParStruct);
								MV->M.P[AreaList[a]].conn[el].ParStruct=InitPhotoTransistorStruct(&(MV->M.P[AreaList[a]].conn[el].ParSize));
						
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jsbe=atof(args[2]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jsbc=atof(args[3]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jph=atof(args[4]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Rs=atof(args[5]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Rsh=atof(args[6]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->EgBE=atof(args[7]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->PhiBC=atof(args[8]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Bf=atof(args[9]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Vaf=atof(args[10]);
								MV->M.P[AreaList[a]].conn[el].model=PHOTOT;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 11);					
						break;
					}
					case SET_SEL_PTEXTJV:
					{
						meshvar *MV;
						int el;
						char *area;		
						char **args;
						args=GetArgs (&begin, 12, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel))
							Error("* line %3d: Invalid electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-1);	
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								if (MV->M.P[a].conn[el].model!=PHOTOT)							
									Error("* line %3d: Setting extra paremeters for the photo-transistor model requires the model to be set first in %s.%s\n", line_nr, MV->name, MV->M.P[a].name);
									
								Print(NORMAL,"* line %3d: Setting photo-transistor model extra parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[a].name);
						
						
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Jse=atof(args[2]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Jsc=atof(args[3]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Nf=atof(args[4]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Nr=atof(args[5]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Ne=atof(args[6]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Nc=atof(args[7]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->Var=atof(args[8]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->XTIBE=atof(args[9]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->XTIBC=atof(args[10]);
								((PhotoTransistor *)(MV->M.P[a].conn[el].ParStruct))->XTB=atof(args[11]);		
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								if (MV->M.P[AreaList[a]].conn[el].model!=PHOTOT)							
									Error("* line %3d: Setting extra paremeters for the photo-transistor model requires the model to be set first in %s.%s\n", line_nr, MV->name, MV->M.P[AreaList[a]].name);
								Print(NORMAL,"* line %3d: Setting photo-transistor model extra parameters between electrodes %d and %d in %s.%s",line_nr, el, el+1, MV->name, MV->M.P[AreaList[a]].name);
						
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jse=atof(args[2]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Jsc=atof(args[3]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Nf=atof(args[4]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Nr=atof(args[5]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Ne=atof(args[6]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Nc=atof(args[7]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->Var=atof(args[8]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->XTIBE=atof(args[9]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->XTIBC=atof(args[10]);
								((PhotoTransistor *)(MV->M.P[AreaList[a]].conn[el].ParStruct))->XTB=atof(args[11]);
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 12);					
						break;
					}
					case SET_R:
					{
						meshvar *MV;
						int P, el;
						double R;						
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
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
						if (MV->M.P[P].conn[el].ParStruct)
							free(MV->M.P[P].conn[el].ParStruct);
						MV->M.P[P].conn[el].ParSize=0;
						MV->M.P[P].conn[el].ParStruct=NULL;	
						FreeArgs (args, 3);		
						break;
					}
					case SET_SEL_R:
					{
						
						meshvar *MV;
						int el;
						double R;
						char *area;		
						char **args;
						args=GetArgs (&begin, 3, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
							
						el=atoi(args[1]);
						if ((el<0)||(el>=MV->M.Nel-1))
							Error("* line %3d: Invalid inter-electrode index: Index %i is not in the valid range of 0 %i\n", line_nr, el, MV->M.Nel-2);
										
						R=atof(args[2]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting R between electrodes %d and %d in %s.%s to %s",line_nr, el, el+1, MV->name, MV->M.P[a].name, args[2]);
								if (MV->M.P[a].conn[el].V)
									free(MV->M.P[a].conn[el].V);
								if (MV->M.P[a].conn[el].J)
									free(MV->M.P[a].conn[el].J);
								MV->M.P[a].conn[el].V=malloc(3*sizeof(double));
								MV->M.P[a].conn[el].J=malloc(3*sizeof(double));
								MV->M.P[a].conn[el].N=2;	
								MV->M.P[a].conn[el].V[0]=-1.0;
								MV->M.P[a].conn[el].J[0]=-1.0/R;
								MV->M.P[a].conn[el].V[1]=1.0;
								MV->M.P[a].conn[el].J[1]=1.0/R;
								MV->M.P[a].conn[el].model=JVD;	
								if (MV->M.P[a].conn[el].ParStruct)
									free(MV->M.P[a].conn[el].ParStruct);
								MV->M.P[a].conn[el].ParSize=0;
								MV->M.P[a].conn[el].ParStruct=NULL;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting R between electrodes %d and %d in %s.%s to %s",line_nr, el, el+1, MV->name, MV->M.P[AreaList[a]].name, args[2]);
								if (MV->M.P[AreaList[a]].conn[el].V)
									free(MV->M.P[AreaList[a]].conn[el].V);
								if (MV->M.P[AreaList[a]].conn[el].J)
									free(MV->M.P[AreaList[a]].conn[el].J);
								MV->M.P[AreaList[a]].conn[el].V=malloc(3*sizeof(double));
								MV->M.P[AreaList[a]].conn[el].J=malloc(3*sizeof(double));
								MV->M.P[AreaList[a]].conn[el].N=2;	
								MV->M.P[AreaList[a]].conn[el].V[0]=-1.0;
								MV->M.P[AreaList[a]].conn[el].J[0]=-1.0/R;
								MV->M.P[AreaList[a]].conn[el].V[1]=1.0;
								MV->M.P[AreaList[a]].conn[el].J[1]=1.0/R;
								MV->M.P[AreaList[a]].conn[el].model=JVD;
								if (MV->M.P[AreaList[a]].conn[el].ParStruct)
									free(MV->M.P[AreaList[a]].conn[el].ParStruct);
								MV->M.P[AreaList[a]].conn[el].ParSize=0;	
								MV->M.P[AreaList[a]].conn[el].ParStruct=NULL;	
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 3);					
						break;
					}
					case SET_T:
					{
						meshvar *MV;
						int P;
						double T;						
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
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
					case SET_SEL_T:
					{
						
						meshvar *MV;
						double T;
						char *area;		
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
											
						T=atof(args[1]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting initial temperature in %s.%s to %s",line_nr, MV->name, MV->M.P[a].name, args[1]);
								MV->M.P[a].T=T;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting initial temperature in %s.%s to %s",line_nr, MV->name, MV->M.P[AreaList[a]].name, args[1]);
								MV->M.P[AreaList[a]].T=T;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 2);					
						break;
					}
					case SET_VG:
					{
						meshvar *MV;
						int P;
						double Vg;						
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
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
														
						Print(NORMAL,"* line %3d: Setting ground-level potential in %s to %s",line_nr, args[0], args[1]);
						Vg=atof(args[1]);
						
						MV->M.P[P].Vg=Vg;
						FreeArgs (args, 2);		
						break;
					}
					case SET_SEL_VG:
					{
						
						meshvar *MV;
						double Vg;
						char *area;		
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
											
						Vg=atof(args[1]);
						
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							/* change R_el for all areas */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Setting ground-level potential in %s.%s to %s",line_nr, MV->name, MV->M.P[a].name, args[1]);
								MV->M.P[a].Vg=Vg;
								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Setting ground-level potential in %s.%s to %s",line_nr, MV->name, MV->M.P[AreaList[a]].name, args[1]);
								MV->M.P[AreaList[a]].Vg=Vg;
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 2);					
						break;
					}
					case SET_SPLITY:
					{
						meshvar *MV;
						int P;						
						char **args;
						args=GetArgs (&begin, 1, Meshes, Nm);
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
					case SET_SEL_SPLITY:
					{
						
						meshvar *MV;
						char *area;		
						char **args;
						args=GetArgs (&begin, 1, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
													
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Toggle split-y parameter of %s.%s",line_nr, MV->name, MV->M.P[a].name);
								if (MV->M.P[a].SplitY)
									MV->M.P[a].SplitY=0;
								else
									MV->M.P[a].SplitY=1;								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								/* assign element to area */
								AreaList=AddToList(AreaList, a);
								MV->M.nodes[MV->nodes[E]].P=a;
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Toggle split-y parameter of %s.%s",line_nr,MV->name, MV->M.P[AreaList[a]].name);
								if (MV->M.P[AreaList[a]].SplitY)
									MV->M.P[AreaList[a]].SplitY=0;
								else
									MV->M.P[AreaList[a]].SplitY=1;	
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 1);					
						break;
					}					
					case SET_SPLITX:
					{
						meshvar *MV;
						int P;					
						char **args;
						args=GetArgs (&begin, 1, Meshes, Nm);
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
					case SET_SEL_SPLITX:
					{
						
						meshvar *MV;
						char *area;		
						char **args;
						args=GetArgs (&begin, 1, Meshes, Nm);
						if (args==NULL)
							goto premature_end;

						area=args[0];
						while (((*area)!='.')&&(*(area+1)))
							area++;						
						area++;	
						if (!(*area))
							Error("* line %3d: Expecting area indication in the form <mesh>.<area_modifier>\n", line_nr);
						(*(area-1))='\0';	
											
						MV=LookupMesh (args[0],  Meshes, Nm);
						if (!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);	
													
												
						if (MV->nodes[0]==0)
						{
							/* all nodes */
							int a;
							for (a=0;a<MV->M.Na;a++)
							{
								Print(NORMAL,"* line %3d: Toggle split-y parameter of %s.%s",line_nr, MV->name, MV->M.P[a].name);
								if (MV->M.P[a].SplitX)
									MV->M.P[a].SplitX=0;
								else
									MV->M.P[a].SplitX=1;								
							}
						}
						else
						{
							/* create new area for all areas in selection and assign elements to it */
							/* add_area to the area names */
							int E, a, l1, l2;
							int *AreaList;
							char *newarea;
							newarea=malloc(MAXSTRLEN*sizeof(char));
							AreaList=malloc(LISTBLOCK*sizeof(int));
							AreaList[0]=0;
							for (E=1;E<=MV->nodes[0];E++)
							{
								a=MV->M.nodes[MV->nodes[E]].P;
								if (!MV->setsel)
								{
									/* create new area string */
									l1=strlen(MV->M.P[a].name);
									l2=strlen(area);
									if (l1+l2>MAXSTRLEN)
										Error("* line %3d: New area-name exceeds %d characters\n", line_nr, MAXSTRLEN);	
										
									strncpy(newarea, MV->M.P[a].name,l1);
									strncpy(newarea+l1, area,l2+1);								
									a=FindProperties(MV->M, newarea);
									if (a<0)
									{
										/* create new area */
										Print(NORMAL,"* line %3d: Creating new area %s.%s",line_nr, MV->name, newarea);
										a=MV->M.Na;
										MV->M.Na++;
										MV->M.P=realloc(MV->M.P, (MV->M.Na+1)*sizeof(local_prop));
										DuplicateProperties(&(MV->M), MV->M.P+a, MV->M.P+MV->M.nodes[MV->nodes[E]].P);
										MV->M.P[a].name=realloc(MV->M.P[a].name, (strlen(newarea)+2)*sizeof(char));
										strncpy(MV->M.P[a].name, newarea,strlen(newarea)+1);
									}
									/* assign element to area */
									MV->M.nodes[MV->nodes[E]].P=a;
								}
								AreaList=AddToList(AreaList, a);
								
							}
							MV->setsel=1;
							for (a=1;a<=AreaList[0];a++)
							{
								Print(NORMAL,"* line %3d: Toggle split-y parameter of %s.%s",line_nr, MV->name, MV->M.P[AreaList[a]].name);
								if (MV->M.P[AreaList[a]].SplitX)
									MV->M.P[AreaList[a]].SplitX=0;
								else
									MV->M.P[AreaList[a]].SplitX=1;	
							}							
							free(newarea);
							free(AreaList);
						}	
						FreeArgs (args, 1);					
						break;
					}				
					/********************************* Solving*/		
					case MAXITER:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting maximum number of iterations to %s",line_nr, word);							
						Numeric_Settings.max_iter=atoi(word);
						if (Numeric_Settings.max_iter<=0)
						{
							Error("* line %3d: Invalid maxiter value %d\n", line_nr,Numeric_Settings.max_iter);
							exit(1);
						}
						break;		
					case TOLV:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting absolute voltage tolerance to %s",line_nr, word);							
						Numeric_Settings.tol_v_abs=atof(word);
						if (Numeric_Settings.tol_v_abs<=1e-16)
						{
							Error("* line %3d: Invalid absolute voltage tolerance value %e\n", line_nr,Numeric_Settings.tol_v_abs);
							exit(1);
						}
						break;		
					case RELTOLV:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting relative voltage tolerance to %s",line_nr, word);							
						Numeric_Settings.tol_v_rel=atof(word);
						if (Numeric_Settings.tol_v_rel<=1e-16)
						{
							Error("* line %3d: Invalid relative voltage tolerance value %e\n", line_nr,Numeric_Settings.tol_v_rel);
							exit(1);
						}
						break;		
					case TOLKCL:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;	
						Print(NORMAL,"* line %3d: Setting absolute KCL tolerance to %s",line_nr, word);							
						Numeric_Settings.tol_kcl_abs=atof(word);
						if (Numeric_Settings.tol_kcl_abs<=1e-16)
						{
							Error("* line %3d: Invalid absolute KCL tolerance value %e\n", line_nr,Numeric_Settings.tol_kcl_abs);
							exit(1);
						}
						break;		
					case RELTOLKCL:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting relative KCL tolerance to %s",line_nr, word);								
						Numeric_Settings.tol_kcl_rel=atof(word);
						if (Numeric_Settings.tol_kcl_rel<=1e-16)
						{
							Error("* line %3d: Invalid relative KCL tolerance value %e\n", line_nr,Numeric_Settings.tol_kcl_rel);
							exit(1);
						}
						break;
					case NLINSEARCH:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting number of steps for a linear search in the Newton direction to %s",line_nr, word);								
						Numeric_Settings.N_lin_search=atoi(word);
						if (Numeric_Settings.N_lin_search<0)
						{
							Error("* line %3d: Invalid N_lin_search value %d\n", line_nr,Numeric_Settings.N_lin_search);
							exit(1);
						}
						break;	
					case GMINSTEP:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting number of Gmin steps to %s (Gmin stepping is triggered at convergence problems)",line_nr, word);								
						Numeric_Settings.GminStep=atoi(word);
						if (Numeric_Settings.GminStep<0)
						{
							Error("* line %3d: Invalid GminStep value %d\n", line_nr,Numeric_Settings.GminStep);
							exit(1);
						}
						break;		
					case GMINSTART:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting number of Gmin start iterations to %s",line_nr, word);								
						Numeric_Settings.GminStart=atoi(word);
						if (Numeric_Settings.GminStart<0)
						{
							Error("* line %3d: Invalid GminStart value %d\n", line_nr,Numeric_Settings.GminStart);
							exit(1);
						}
						break;	
					case GMINMAX:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting maximum Gmin value to %s",line_nr, word);								
						Numeric_Settings.GminMax=atof(word);
						if (Numeric_Settings.GminMax<=1e-16)
						{
							Error("* line %3d: Invalid maximum Gmin value %e\n", line_nr,Numeric_Settings.GminMax);
							exit(1);
						}
						break;	
					case GMINFAC:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting Gmin stepping factor to %s",line_nr, word);								
						Numeric_Settings.GminFac=atof(word);
						if (Numeric_Settings.GminFac<=1)
						{
							Error("* line %3d: Invalid maximum GminFac value %e\n", line_nr,Numeric_Settings.GminFac);
							exit(1);
						}
						break;	
					case GMIN:
						begin=GetWord (begin, word, Meshes, Nm);
						if(word[0]=='\0')
							goto premature_end;
						Print(NORMAL,"* line %3d: Setting Gmin to %s",line_nr, word);								
						Numeric_Settings.Gmin=atof(word);
						if (Numeric_Settings.Gmin<=0.0)
						{
							Error("* line %3d: Invalid maximum Gmin value %e\n", line_nr,Numeric_Settings.Gmin);
							exit(1);
						}
						break;		
					case SOLVE:
					{
						meshvar *MV;
						double Vstart, Vend;
						int Nstep;				
						char **args;
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL,"* line %3d: Solving potentials in mesh %s",line_nr, args[0]);						
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						
						Vstart=atof(args[1]);
						Vend=atof(args[2]);
						Nstep=atoi(args[3]);
						
						SolveVa(&(MV->M), Vstart, Vend, Nstep);
						FreeArgs (args, 4);
						break;
					}	
					case REFINEOC:
					{
						meshvar *MV;
						double tol_v, tol_i;
						int Niter;				
						char **args;
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL,"* line %3d: Refining open circuit point for mesh %s",line_nr, args[0]);						
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						tol_i=atof(args[1]);
						tol_v=atof(args[2]);
						Niter=atoi(args[3]);
						
						RefineOC(&(MV->M), tol_i, tol_v, Niter);
						FreeArgs (args, 4);
						
						break;
					}		
					case REFINEMPP:
					{
						meshvar *MV;
						double tol_v, tol_i;
						int Niter;			
						char **args;
						args=GetArgs (&begin, 4, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
							
						Print(NORMAL,"* line %3d: Refining maximum power-point for mesh %s",line_nr, args[0]);						
						MV=LookupMesh (args[0],  Meshes, Nm);
						if(!MV)
							Error("* line %3d: Mesh \"%s\" does not exist\n",line_nr,args[0]);
						
						
						tol_i=atof(args[1]);
						tol_v=atof(args[2]);
						Niter=atoi(args[3]);
						RefineMPP(&(MV->M), tol_i, tol_v, Niter);
						FreeArgs (args, 4);
						break;
					}
					case ADAPTIVE_SOLVE:
					{
						meshvar *MV;
						double Va, rel_th;
						int Na;				
						char **args;
						args=GetArgs (&begin, 4, Meshes, Nm);
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
						AdaptiveSolveVa(&(MV->M), Va, rel_th, Na);
						FreeArgs (args, 4);
						break;
					
					}
					case TIC:
						Print(NORMAL,"* line %3d: Starting timer",line_nr);	
						Print(NORMAL, "------------------------timer start-----------------------------");
						tic = clock();
						break;
					case TOC:
						toc = clock();
						Print(NORMAL, "-----------------------timer readout----------------------------");
						Print(NORMAL,"* line %3d: %e seconds since last tic",line_nr, ((double) (toc - tic)) / CLOCKS_PER_SEC);
						break;
					case EXPR_DEF:	
					{	
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						Print(NORMAL,"* line %3d: defining \"%s=%s\"",line_nr, args[0], args[1]);	
						DefineVar(args[0], atof(args[1]));
						FreeArgs (args, 2);
						break;
					}
					case EXPR_IFDEF:	
					{	
						char **args;
						args=GetArgs (&begin, 2, Meshes, Nm);
						if (args==NULL)
							goto premature_end;
						if (!IsDefined(args[0]))
						{
							Print(NORMAL,"* line %3d: defining \"%s=%s\"",line_nr, args[0], args[1]);	
							DefineVar(args[0], atof(args[1]));
						}
						else
							Print(NORMAL,"* line %3d: variable %s is already defined",line_nr, args[0]);
							
						FreeArgs (args, 2);
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
					case EXIT:
						Print(NORMAL,"* line %3d: Exit command encountered",line_nr);					
						ExitPVMOS=0;
						break;
					case DEBUGCOMMAND:
						/* I sometimes insert things here for testing purposes,  it is not intended for normal use */
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
		last_pos=ftell(f);
    		fgets(line, MAXSTRLEN-1, f);
		line_nr++;
	}
	fclose(f);
	free(loop_stack);
	free(loop_line_nr);
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
