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
 * some small utilitiy routines, mostly regarding printing       *
 * output with respect to the set verbosity level                *         
 *                                                               *            
 *****************************************************************/     
#ifdef __MINGW32__ 
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include "main.h"
#include "mesh2d.h"
#include "utils.h"
#include "config.h"
#define LINE_WIDTH 65
#define LINE_SKIP 12

void Print(VERB_LEVEL v_level, const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	if (v_level<=verbose)
	{
		int numc;
		char * s;
		s=malloc((LINE_WIDTH)*sizeof(char));
		numc=vsnprintf(s, LINE_WIDTH*sizeof(char), format_str, ap);
      		/* on one particular system ap changed upon execustion of vsnprintf! 
		   I do not know whether it is a bug or that I was wrong to assume
		   it does not change. Anyway, my fix is multiple calls to va_start 
		*/
      		va_start (ap, format_str); 
		if (numc>=LINE_WIDTH)
		{
			int br, i, minl=0;
			char *newline;
			
			s=realloc(s,(numc+1)*sizeof(char));
			numc=vsnprintf(s, (numc+1)*sizeof(char), format_str, ap);
      			va_start (ap, format_str);
			newline=s;
			while (strlen(newline)>LINE_WIDTH)
			{
				br=LINE_WIDTH-1;
				while ((!isblank(newline[br]))&&(br>minl))
					br--;
				if (br>minl)
				{
					newline[br]='\0';
					printf("%s\n",newline);
					for (i=0;i<LINE_SKIP;i++)
						newline[br-i]=' ';
					newline=newline+br+1-i;
				}
				else
				{
					char c;
					br=LINE_WIDTH-1;
					c=newline[br];
					newline[br]='\0';
					printf("%s\n",newline);					
					newline[br]=c;
					for (i=1;i<LINE_SKIP;i++)
						newline[br-i]=' ';
					newline=newline+br+1-i;
				}
				minl=LINE_SKIP;
			}
			printf("%s\n",newline);
		}
		else
			printf("%s\n",s);	
		free(s);
		/* vprintf(format_str, ap);*/
	}
}

void Print_nlb(VERB_LEVEL v_level, const char *format_str, ...)
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
	printf("PVMOS Version %s   %s   B.E. Pieters \n", VERSION, __DATE__); 
	printf("IEK-5 Photovoltaik Forschungszentrum Juelich, Germany\n"); 
	printf("\n");
}

char *TimeString(void)
{
	time_t current_time;
	char* c_time_string;
	int len;
	 
	current_time = time(NULL);
	if (current_time == ((time_t)-1))
	{
		Warning("Warning: time is unknown. I have no idea when we are!\n");
		c_time_string=calloc(1, sizeof(char));
 		return c_time_string;
	}
	 
	/* Convert to local time format. */
	c_time_string = ctime(&current_time);
	 
	if (c_time_string == NULL)
	{
		Warning("Warning: time is unknown. I have no idea when we are!\n");
		c_time_string=calloc(1, sizeof(char));
 		return c_time_string;
	}
	
	len=strlen(c_time_string);
	c_time_string=malloc((len+1)*sizeof(char));
	strncpy(c_time_string, ctime(&current_time), len);
	c_time_string[len-1]='\0';
	
 	return c_time_string;
}

void PrintFileHeader(FILE *f)
{
	char* c_time_string;
	fprintf(f, "#   _____   __   __  __  ___  ___                                  \n");
	fprintf(f, "#  | _ \\ \\ / /__|  \\/  |/ _ \\/ __|                                 \n");
	fprintf(f, "#  |  _/\\ V /___| |\\/| | (_) \\__ \\                                 \n");
	fprintf(f, "#  |_|   \\_/    |_|  |_|\\___/|___/    PV-Module Simulator          \n");
	fprintf(f, "# ________________________________________________________________\n");
	fprintf(f, "# PV-MOS Version %s   %s   B.E. Pieters \n", VERSION, __DATE__); 
	fprintf(f, "# IEK-5 Photovoltaik Forschungszentrum Juelich, Germany\n"); 
	c_time_string=TimeString();
	fprintf(f, "# Date: %s\n#\n", c_time_string); 
	free(c_time_string);
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

char *GetLine(FILE *f)
{
	/* reads a line from a file and returns it as an allocated and null terminated string */
	int ll=0, al=20; /* line length, allocated length */
	char *line, c;
	line=malloc(al*sizeof(char));
	
	c=fgetc(f);
	while (feof(f)==0)
	{
		if (c=='\n')
			break;
		ll++;
		if(ll==al)
		{
			al+=20;
			line=realloc(line, al*sizeof(char));
		}
		line[ll-1]=c;
		c=fgetc(f);	
	}
	line[ll]='\0';
	line=realloc(line, (ll+1)*sizeof(char));
	return line;
}

void GetNameValue(char *string, char **name, char **value)
{
	/* Separates the name from the value in a string containing a variable definition */
	/* name and value are separated by a = sign */
	/* blanks at the start of the line are ignored */
	/* blanks at the end of the line are ignored */
	/* blanks before = are ignored */
	/* blanks after = are ignored */
	/* names and values can be quoted (e.g. to add blanks which would otherwise be ignored) */
	/* characters can be escaped (e.g. to add quotation marks to a name or value) */
	
	int quoted=0, esc=0;
	char *name_end;
	
	/* strip blanks at the start */
	while(((*string)==' ') || ((*string)=='\t') )
		string++;
	
	(*name)=string;
	name_end=string;
	
	while (*name_end)
	{
		if (esc)
			esc=0;
		else if (quoted)
		{
			if(*name_end==quoted)
				quoted=0;
		}
		else
		{
			if ( ((*name_end)=='\\'))
			{
				esc=1;
				while (name_end[esc-1])
				{
					name_end[esc-1]=name_end[esc];
					esc++;
					printf("%s\n",name_end);
				}				
				
			}
			else if ( ((*name_end)=='\"') || ((*name_end)=='\'') )
				quoted=(*name_end);
			else if ((*name_end) == '=' )
			{
				char *value_start;
				int value_len;
				
				/* walk back over whitespaces */
				value_start=name_end;
				do
					name_end--;
				while (( ((*name_end)==' ') || ((*name_end)=='\t') ) && name_end>string);
				
				if( (*(*name)==*name_end) && (((*name_end)=='\"') || ((*name_end)=='\'')))
				{
					/* skip quotes */
					(*name)++;
					name_end--;
				}
				*(name_end+1)='\0';
				
				/* walk forward over whitespaces */				
				do
					value_start++;
				while ( ((*value_start)==' ') || ((*value_start)=='\t') );
				
				/* deal with escapes */
				while(value_start[esc])
				{
					if (value_start[esc]=='\\')
					{
						int j=esc+1;
						while (value_start[j-1])
						{
							value_start[j-1]=value_start[j];
							j++;
						}						
					}
					esc++;
				}
				
				
				value_len=strlen(value_start);
				
				value_len--;
				while ((value_start[value_len]==' ') || (value_start[value_len]=='\t'))
					value_len--;
				
				if(( value_start[0]==value_start[value_len]) && ((value_start[0]=='\"') || (value_start[0]=='\'')))
				{
					/* skip quotes */
					value_start++;
					value_start[value_len-1]='\0';
				}							
				(*value)=value_start;
				return;
			}
		}
		name_end++;
	}
	(*name)=NULL;
	(*value)=NULL;
}
