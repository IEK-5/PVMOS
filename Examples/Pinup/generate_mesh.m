#!/usr/bin/octave -q
pkg load mkpvmosmesh

function img=ReadImageIntensity(fn)
	img=imread(fn);
	if (ndims(img)>2)
		img = rgb2ind (img);
	endif
	img=double(img);
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%
function AreaDef=DefaultAreaDef(Nel)
% generate the default material properties struct. Input: Nel, the number of electrodes.
	AreaDef.name="default";				% name of the area definition
	AreaDef.Rel=ones(1,Nel);			% electrode sheet resistance [Ohm]
	AreaDef.Rvp=-ones(1,Nel);			% Resistance to the positive node in [Ohm cm^2], if the value is negative there is no direct connection to the positive nde
	AreaDef.Rvn=-ones(1,Nel);			% Resistance to the negative node in [Ohm cm^2], if the value is negative there is no direct connection to the negative nde
	for i=1:Nel-1					% Define the connections between the electrodes
		AreaDef.conn(i).model="ONED";		% Connection model, can be ONED (one diode model) TWOD (two diode model) or JVD (Tabular JV data)
		AreaDef.conn(i).J01=1e-12;		% Diode 1 saturation current density (used in ONED and TWOD) [A/cm^2]
		AreaDef.conn(i).J02=1e-8;		% Diode 2 saturation current density (used only in TWOD) [A/cm^2]
		AreaDef.conn(i).Jph=0; 			% Photo-current density current (used in ONED and TWOD) [A/cm^2]
		AreaDef.conn(i).nid1=1; 		% Ideality factor for diode 1 (used in ONED and TWOD)
		AreaDef.conn(i).nid2=2;  		% Ideality factor for diode 2 (used in TWOD)
		AreaDef.conn(i).Eg=1.12;  		% Bandgap of the diode (used for temperature dependency, used in ONED and TWOD) [eV]
		AreaDef.conn(i).Rs=1e-5;   		% Series resistance (used in ONED and TWOD) [Ohm cm^2]
		AreaDef.conn(i).Rsh=1e3;    		% Shunt resistance (used in ONED and TWOD) [Ohm cm^2]
		AreaDef.conn(i).N=2;     		% Number of JV datapoints (for JVD model)
		AreaDef.conn(i).VJ=[-1,-1;1,1];     	% Tabular JV data (for JVD model) Column 1: [V] Column 2: [A/cm^2]
	endfor
	AreaDef.T=300;   				% Initial temperature in the area
	AreaDef.SplitX=1;    				% Allow element to be split in x-direction during adaptive meshing (0 false, 1: true)
	AreaDef.SplitY=1;    				% Allow element to be split in y-direction during adaptive meshing (0 false, 1: true)
endfunction

function [Area_Index, AreaDef]=CleanupAreas(Area_Index,AreaDef)
% cleanup any unused area definitions
	set=unique(Area_Index);
	a=ones(size(Area_Index));
	m=max(set);
	for i=1:length(set)
		ii=(Area_Index==set(i));
		Area_Index(ii)=a(ii)*(m+i);
		NewDef(i)=AreaDef(set(i));
	endfor
	Area_Index.-=m;
	AreaDef=NewDef;		
endfunction

function [Area_Index, AreaDef]=AddImage(Area_Index,AreaDef, img, modlist, name)
% this function helps to define regular PVMOS meshes.
	printf("===========================================================================\n");
	printf("Adding image %s\n",name);
	printf("---------------------------------------------------------------------------\n");
	for j=1:length(modlist)
		printf("%s;\n",modlist{j})
	endfor
	% make sure img is boolean type
	img=(img>0);
	if (length(Area_Index)==0)
		Area_Index=img.+1;
		m=1;
	else
		m=max(max(Area_Index));
		img.*=m;
		Area_Index.+=img;
		
	endif
	Na=max(max(Area_Index))
	for i=m+1:Na
		AreaDef(i)=AreaDef(i-m);
		AreaDef(i).name=[AreaDef(i-m).name,name];
		for j=1:length(modlist)
			c=sprintf("AreaDef(%d).%s;",i,modlist{j});
			eval(c);
		endfor
	endfor
	[Area_Index, AreaDef]=CleanupAreas(Area_Index,AreaDef);
	for i=m+1:max(max(Area_Index))
		printf("New Area: %s\n",AreaDef(i).name);
	endfor
	printf("Created %i new areas\n",max(max(Area_Index))-m);
	printf("===========================================================================\n\n");
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parametrization %%%%
Nel=3;						% this example uses three electrodes. Electrode 1 is the 
						% front contact, 2: back contact, 3 contacting foil
AreaDef=DefaultAreaDef(Nel);			% Init the default parameters
AreaDef.Rel=[0.01, 0.005, 40]; 			% Set the electrode sheet resistances, foil, back contact, front (emitter)
AreaDef.conn(1).model="JVD";			% Set connection 2 (between electrodes 2 and 3) to tabular JV data
AreaDef.conn(1).N=2; 				% Set the number of data points for connection 1
AreaDef.conn(1).VJ=[-1,-1e-10;1,1e-10]; 	% Set connection 1 is ohmic with a high resistance (1e10 Ohm cm^2)
AreaDef.conn(2).model="ONED";			% Set connection 2 (between electrodes 2 and 3) to a 1 diode model
AreaDef.conn(2).J01=1e-12;			% Set the satuartion current density for connection 2
AreaDef.conn(2).Jph=0.037;			% Set the photo-current density for connection 2
AreaDef.conn(2).nid1=1;				% Set the ideality factor for connection 2
AreaDef.conn(2).Rs=1e-5;			% Set the series resistance for connection 2
AreaDef.conn(2).Rsh=1e5;			% Set the shunt resistance for connection 2


% the default parameters we just defined apply everywhere, we will however locally change the properties
% To this end we load a series of images which we use to select elements in the mesh.
% We load the first image, which represents the front melatization pattern:
img=ReadImageIntensity("pinup_frontmetal.png");
% The parts in the image that are black are the areas where we have a front metal contact. These elemnts 
% in the matrix are 0. This command will change the matrix to a matrix with 0 where there is no metal and
% 1 where there is metal:
img=(img==0);

% this variable is a list of definitions that will be applied to those elements which have a front metal
% In this case that is a reduction in the sheet resistance of electrode 1 and setting the photo-current to 0
modlist={"Rel(3)=0.01"; "conn(2).Jph=0"};
% For this first image we pass an empty matrix for the Area_Index matrix, this will signal to initialize the 
% mesh to a mesh with the dimensions of the img matrix.
[Area_Index, AreaDef]=AddImage([],AreaDef, img, modlist, "front_metal");

% the next image marks the areas where the back metalization is removed around the metal-wrap-through-holes
img=ReadImageIntensity("pinup_backisolation.png");
img=(img==0);
modlist={"Rel(2)=1e10"};	% make electrode 2 isolating
[Area_Index, AreaDef]=AddImage(Area_Index, AreaDef, img, modlist, "no_backmetal");

% Marks the vias
img=ReadImageIntensity("pinup_vias.png");
img=(img==0);
% At the vias we do not only connect the front to the back contact but also connect to the contacting foil (electrode 1)
modlist={"Rel(2)=0.01"; "conn(2).model=\"JVD\""; "conn(2).VJ=[-1,-1e8;1,1e8]"; "conn(1).model=\"JVD\""; "conn(1).VJ=[-1,-1e8;1,1e8]"};
[Area_Index, AreaDef]=AddImage(Area_Index, AreaDef, img, modlist, "vias_p-contact"); 

% This marks the spots where the back contact is connected to the contact foil
img=ReadImageIntensity("pinup_back2contactfoil.png");
img=(img==0);
modlist={"conn(1).model=\"JVD\""; "conn(1).VJ=[-1,-1e8;1,1e8]"};
[Area_Index, AreaDef]=AddImage(Area_Index, AreaDef, img, modlist, "p-contact");


% This marks areas where the contact foil is isolating to separate the positive and negative contacts
img=ReadImageIntensity("pinup_foilisolation.png");
img=(img==0);
modlist={"Rel(1)=1e15"};
[Area_Index, AreaDef]=AddImage(Area_Index, AreaDef, img, modlist, "contact_foil_isolation");

% This marks areas where the contact foil contacted to the positive electrode
img=ReadImageIntensity("pinup_contactfoil_vp.png");
img=(img==0);
modlist={"Rvn(1)=1e-8"};
[Area_Index, AreaDef]=AddImage(Area_Index, AreaDef, img, modlist, "contact_foil_vp");

% This marks areas where the contact foil contacted to the negative electrode
img=ReadImageIntensity("pinup_contactfoil_vn.png");
img=(img==0);
modlist={"Rvp(1)=1e-8"};
[Area_Index, AreaDef]=AddImage(Area_Index, AreaDef, img, modlist, "contact_foil_vn");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate the mesh %%%%
Na=max(max(Area_Index));
Nx=length(Area_Index(1,:))
Ny=length(Area_Index(:,1))
x=0:10/Nx:10;
y=0:10/Ny:10;
% make a 10x10 cm^2 mesh
mkpvmosmesh(Na, Nel, Area_Index, AreaDef, x, y, "pinupmesh2.bin")

