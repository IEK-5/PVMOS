#include <octave/oct.h>
#include <octave/variables.h>
#include <octave/ov-struct.h>
extern "C"
{
  #include "mesh2d.h"
  #include "main.h"
  #include "utils.h"
  #include "diode.h"
  #include "phototransistor.h"
}


VERB_LEVEL verbose=NORMAL;
int fixverb=0;
typedef struct {
	const char *name;
	diode_model D;
} ElModel;

const ElModel ElModelTable[] =
{
      	{"JVD", JVD},
	{"ONED", ONED},
	{"TWOD", TWOD},
	{"PHOTOT", PHOTOT},
      	{NULL, JVD}
};

octave_value GetValue(octave_map Area_Def, int i,  char const * name)
{
	if (! error_state && Area_Def.contains (name))
	{
		octave_map::const_iterator index = Area_Def.seek (name);
		if (i>=Area_Def.contents(index).numel ())
		{
			fprintf(stderr, "Error: index %i out of bound %i in GetValue\n", i+1, Area_Def.contents(index).numel ());
			error_state=1;
			return octave_value();
		}
		octave_value tmp = Area_Def.contents(index)(i);
		return tmp;
	}
	error("Struct does not contain field %s\n", name);
	error_state=1;
	return octave_value();
}

char * Getstring(charMatrix string)
/* There is probably better ways of dong this but I simply want simple plain null terminated pointers to characters */
/* This fetached the data from a charMatrix and coppies it to a freshly allocated string and null-terminetes it */
/* do not forget to free the allocated string after use. */
{
	char * res;
	res=(char *) malloc((string.cols()+1)*sizeof(char));
	memcpy(res,string.fortran_vec(),string.cols()*sizeof(char) );
	res[string.cols()]='\0';
	return res;
}
octave_value VectorIndex(Matrix v, int i)
{
	octave_value tmp;
	if (v.cols()==1)
	{
		if ((i<0)||(i>=v.rows()))
			i=v.rows()-1;
		tmp=v(i,0);
	}
	else
	{
		if ((i<0)||(i>=v.cols()))
			i=v.cols()-1;
		tmp=v(0,i);
	}
	return tmp;
}


int Min(Matrix m)
{
	int i,j;
	double res;
	res=(double)m(0,0);
	for(i=0;i<m.cols();i++)
		for (j=0;j<m.rows();j++)
			if (res>(double)m(j,i))
				res=(double)m(j,i);
				
	return (int)res;
}

DEFUN_DLD (mkpvmosmesh, args, nargout,"\
mkpvmosmesh(Na, Nel, Area_Index, Area_Def, x, y, filename)\n\
Generate a regular mesh for PVMOS from octave data structures and write it to a file.\n\
Na          - The number or area's in the mesh\n\
Nel         - The number or electrodes\n\
Area_Index  - A matrix with index number reffering to areas\n\
Area_Def    - The area definition struct array describing the properties of each area\n\
x           - vector with x coordinates of element boundaries\n\
y           - vector with y coordinates of element boundaries\n\
filename    - filename of the file to dump the resulting PVMOS mesh in")
{
	mesh M;
	int i, j, k, N;
  	int nargin = args.length ();
  	if (nargin != 7)
	{
    		print_usage ();
		return octave_value_list ();
	}
	int Na = args(0).int_value();
	int Nel = args(1).int_value();
	Matrix Area_Index = args(2).matrix_value();
        octave_map Area_Def = args(3).map_value ();
	Matrix x=args(4).matrix_value();
	Matrix y=args(5).matrix_value();
	char * string;
 	charMatrix fn = args(6).char_matrix_value ();
	octave_value res;
	
	if (Nel<2)
	{
    		fprintf(stderr, "Error: PVMOS meshes have a minimum of 2 electrodes\n");	
		FreeMesh(&M);
		return octave_value_list ();
	}
	
    	printf("Init mesh\n");
	
	
	/* the default area for the newly created mesh is equal to the name of the first defined area */
	res=GetValue(Area_Def, 0,  "name");
	if (error_state)
	{
		FreeMesh(&M);
		return octave_value();
	}
	string=Getstring(res.char_matrix_value());
	M=InitMesh(string, 0,  1, 0, 1, Area_Index.cols(), Area_Index.rows());

	
	
    	printf("Setting the number of electrodes to %d\n", Nel);
	/* set the right number of electrodes */
	for (i=0;i<Nel-2;i++)
		AddElectrode(&M);
	
    	printf("Creating %d areas\n", Na);
	
	/* we loop through the struct array to define the areas */
	/* create the new areas */
	for (i=0;i<Na;i++)
	{
		/* allow user intterupts to shut things down */
    		OCTAVE_QUIT;
    		printf("Area %i of %i: ", i, Na);
		fflush(stdout);
		
		res=GetValue(Area_Def, i,  "name");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		string=Getstring(res.char_matrix_value());
    		printf("creating area %s\n", string);
		NewProperties(&M, string);
		free(string);
		
		res=GetValue(Area_Def, i,  "Rel");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		Matrix Rel=res.matrix_value();
		if (Nel!=Rel.numel ())
		{
			fprintf(stderr, "Error: length of Rel is unequal to the number of electrodes\nNel=%i and Rel is %i long\n", Nel, Rel.numel ());
			exit(1);
		}
		for (j=0;j<Nel;j++)
			M.P[i].Rel[j]=Rel(j);
			
		res=GetValue(Area_Def, i,  "Rvp");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		Matrix Rvp=res.matrix_value();
		if (Nel!=Rvp.numel ())
		{
			fprintf(stderr, "Error: length of Rvp is unequal to the number of electrodes\nNel=%i and Rel is %i long\n", Nel, Rvp.numel ());
			exit(1);
		}
		for (j=0;j<Nel;j++)
			M.P[i].Rvp[j]=Rvp(j);
			
		res=GetValue(Area_Def, i,  "Rvn");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		Matrix Rvn=res.matrix_value();
		if (Nel!=Rvn.numel ())
		{
			fprintf(stderr, "Error: length of Rvp is unequal to the number of electrodes\nNel=%i and Rel is %i long\n", Nel, Rvn.numel ());
			exit(1);
		}
		for (j=0;j<Nel;j++)
			M.P[i].Rvn[j]=Rvn(j);
			
		res=GetValue(Area_Def, i,  "conn");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		octave_map elcon=res.map_value();
		if (Nel-1!=elcon.numel ())
		{
			fprintf(stderr, "Error: length of the connection struct is unequal to the number of electrodes - 1\nNel=%i and connection struct is %i long\n", Nel, elcon.numel ());
			exit(1);
		}
		for (j=0;j<Nel-1;j++)
		{
			ElModel *model;
			res=GetValue(elcon, j,  "model");
			if (error_state)
			{
				FreeMesh(&M);
				return octave_value();
			}
			model=(ElModel *) ElModelTable;
			
      			while(model->name)
			{
				fflush(stdout);
				if (strncmp (res.char_matrix_value().fortran_vec(), model->name, 4) == 0)
					M.P[i].conn[j].model=model->D;
				model++;
			}
			
			switch (M.P[i].conn[j].model)
			{
				case JVD:
				{
					res=GetValue(elcon, j,  "VJ");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					Matrix JV=res.matrix_value();
					M.P[i].conn[j].N=JV.rows();
					
					N=JV.rows();
					if ((N<2)||(JV.cols()!=2))
					{
						fprintf(stderr, "JV table should be a 2 columns and N rows with N>1\n");
						exit(1);
					}
					M.P[i].conn[j].V=(double *)realloc(M.P[i].conn[j].V, (N+1)*sizeof(double));
					M.P[i].conn[j].J=(double *)realloc(M.P[i].conn[j].J, (N+1)*sizeof(double));
		
					M.P[i].conn[j].N=N;
					for (k=0;k<N;k++)
					{
						M.P[i].conn[j].V[k]=JV(k,0);
						M.P[i].conn[j].J[k]=JV(k,1);
					}
					break;
				}
				case ONED:
				case TWOD:
					M.P[i].conn[j].ParStruct=InitOneTwoDiodeStruct(&(M.P[i].conn[j].ParSize));
					res=GetValue(elcon, j,  "J01");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->J01=res.double_value();
					res=GetValue(elcon, j,  "J02");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->J02=res.double_value();
					res=GetValue(elcon, j,  "Jph");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->Jph=res.double_value();
					res=GetValue(elcon, j,  "nid1");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->nid1=res.double_value();
					res=GetValue(elcon, j,  "nid2");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->nid2=res.double_value();
					res=GetValue(elcon, j,  "Eg");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->Eg=res.double_value();
					res=GetValue(elcon, j,  "Rs");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->Rs=res.double_value();
					res=GetValue(elcon, j,  "Rsh");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((OneTwoDiode *) M.P[i].conn[j].ParStruct)->Rsh=res.double_value();
					break;
				case PHOTOT:
					M.P[i].conn[j].ParStruct=InitPhotoTransistorStruct(&(M.P[i].conn[j].ParSize));
					res=GetValue(elcon, j,  "Jsbe");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Jsbe=res.double_value();
					
					res=GetValue(elcon, j,  "Jsbc");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Jsbc=res.double_value();
					
					res=GetValue(elcon, j,  "Jph");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Jph=res.double_value();
					
					res=GetValue(elcon, j,  "Rs");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Rs=res.double_value();

					res=GetValue(elcon, j,  "Rsh");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Rsh=res.double_value();

					res=GetValue(elcon, j,  "EgBE");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->EgBE=res.double_value();

					res=GetValue(elcon, j,  "PhiBC");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->PhiBC=res.double_value();

					res=GetValue(elcon, j,  "Bf");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Bf=res.double_value();

					res=GetValue(elcon, j,  "Vaf");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Vaf=res.double_value();

					res=GetValue(elcon, j,  "Jse");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Jse=res.double_value();

					res=GetValue(elcon, j,  "Jsc");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Jsc=res.double_value();

					res=GetValue(elcon, j,  "Nf");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Nf=res.double_value();

					res=GetValue(elcon, j,  "Nr");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Nr=res.double_value();

					res=GetValue(elcon, j,  "Ne");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Ne=res.double_value();

					res=GetValue(elcon, j,  "Nc");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Nc=res.double_value();

					res=GetValue(elcon, j,  "Var");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->Var=res.double_value();

					res=GetValue(elcon, j,  "XTIBE");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->XTIBE=res.double_value();

					res=GetValue(elcon, j,  "XTIBC");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->XTIBC=res.double_value();

					res=GetValue(elcon, j,  "XTB");
					if (error_state)
					{
						FreeMesh(&M);
						return octave_value();
					}
					((PhotoTransistor *)(M.P[i].conn[j].ParStruct))->XTB=res.double_value();
					break;
					
			}
			
		}
			
		res=GetValue(Area_Def, i,  "T");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		M.P[i].T=res.double_value();
		
		res=GetValue(Area_Def, i,  "SplitX");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		M.P[i].SplitX=res.int_value();
		res=GetValue(Area_Def, i,  "SplitY");
		if (error_state)
		{
			FreeMesh(&M);
			return octave_value();
		}
		M.P[i].SplitY=res.int_value();
	}
    	printf("Setting element properties (assigning to areas and coordinates)\n");
	int Nn=0;
	double x1=0;
	double x2=0;
	double y1=0;
	double y2=0;
	int min=Min(Area_Index);
	for (i=0;i<Area_Index.cols();i++)
	{
		x1=VectorIndex(x, i).double_value();
		x2=VectorIndex(x, i+1).double_value();
		for (j=0;j<Area_Index.rows();j++)
		{
			y1=VectorIndex(y, j).double_value();
			y2=VectorIndex(y, j+1).double_value();
			M.nodes[Nn].P=(int)Area_Index(j,i)-min;
			
			M.nodes[Nn].x1=x1;
			M.nodes[Nn].x2=x2;
			M.nodes[Nn].y1=y1;
			M.nodes[Nn].y2=y2;
			Nn++;
		}
	}
	string=Getstring(fn);
    	printf("Writing mesh to file %s\n", string);
	WriteMesh(string, &M);	
	free(string);
	FreeMesh(&M);
	return octave_value();
}
