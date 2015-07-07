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
 * lists are sorted arrays of unique (i.e. a certain integer     *
 * value can occur only once in a given list) integers           * 
 * Here are routines to create list, add/remove elements to/from *
 * a list, and find elements in a list                           * 
 *                                                               *            
 *****************************************************************/     
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "list.h"
#include "main.h"
#include "utils.h"
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))


/* I use many lists for the bookkeeping. The format of a list is:
   A list is an array of type int
   The first element is the number of elements in the list, the resit is the list (obviously)
*/
int FindInList(int id, int *list, int *index)
/* finds an element in the list, returns the index of this element. If the element is not present
   the index points to where the element *should* be it it were in the list */
{
	int max, min, i;
	if ((list[0]==0)||(list[1]>id))
	{
		(*index)=1;
		return 0;
	}
	if (list[list[0]]<id)
	{
		(*index)=list[0]+1;
		return 0;
	}		
	min=1;
	max=list[0];
	while(max-min>1)
	{
		i=(min+max)/2;
		if (list[i]==id)
		{
			(*index)=i;
			return 1;		
		}
		if (list[i]<id)
			min=i;
		else
			max=i;
	}
	if (list[min]<id)
	{
		(*index)=max;
		return (list[max]==id);
	}
	(*index)=min;		
	return (list[min]==id);
}

int *AddToList(int *list, int id)
/*  add a node id to a node list */
{
	int i;
	
	if (!FindInList(id, list, &i))
	{
		int j;
		if ((list[0]+3)%LISTBLOCK==0)
			list=realloc(list,(list[0]+3+LISTBLOCK)*sizeof(int));
		for (j=list[0];j>=i;j--)
			list[j+1]=list[j];
		list[i]=id;
		list[0]++;
		
	}
	return list;
}

int *RemoveFromList(int *list, int id)
/* remove a node id from a node list */
{
	int i;
	
	if (!FindInList(id, list, &i))
		return list;
	for (;i<list[0];i++)
		list[i]=list[i+1];
	list[0]--;
	if ((list[0]+2)%LISTBLOCK==0)
		list=realloc(list,(list[0]+2+LISTBLOCK)*sizeof(int));
	
	return list;
}

int *AddListToList(int *dest, int *source)
/* add a node id to a node list */
{
	int i;
	
	for (i=1;i<=source[0];i++)
		dest=AddToList(dest, source[i]);
	return dest;
}
int *RemoveListFromList(int *dest, int *source)
/* add a node id to a node list */
{
	int i;
	
	for (i=1;i<=source[0];i++)
		dest=RemoveFromList(dest, source[i]);
	return dest;
}

int IsInList(int *list, int id)
{
	int max, min, i;
	if ((list[0]==0)||(list[1]>id))
		return 0;
	if (list[list[0]]<id)
		return 0;
	min=1;
	max=list[0];
	while(max-min>1)
	{
		i=(min+max)/2;
		if (list[i]==id)
			return 1;
		if (list[i]<id)
			min=i;
		else
			max=i;
	}
	return ((list[min]==id)||(list[max]==id));
}

void SortList(int * list)
{
	int a, b, c, s;
	for(a=1;a<list[0];a++)
	{
		s=1;
		for(b=1;b<list[0]-a;b++)
		{
			if (list[b]>list[b+1])
			{
				c=list[b];
				list[b]=list[b+1];
				list[b+1]=c;
				s=0;
			}
		}
		if (s)
			break;
	}
}

int *DuplicateList(int *list)
/* duplicate a node list */
{
	int *res;
	res=malloc(((list[0]+2)/LISTBLOCK+1)*LISTBLOCK*sizeof(int));
	res=memcpy(res, list, (list[0]+1)*sizeof(int));
	return res;
}
