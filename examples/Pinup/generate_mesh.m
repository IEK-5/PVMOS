#!/usr/bin/octave -q
pkg load pvmos-mesh


% reads an image and retuns it as a 2D double array
function img=ReadImageIntensity(fn)
	img=imread(fn);
	if (ndims(img)>2)
		img = rgb2ind (img);
	endif
	img=double(img);
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parametrization %%%%
Nel=3;						% this example uses three electrodes. Electrode 1 is the 
						% front contact, 2: back contact, 3 contacting foil
AreaDef=PVMOS_AreaProperties(Nel);		% Init the default parameters
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


img=ReadImageIntensity("pinup_frontmetal.png");	% loads a bitmap image with the front metalization pattern
mask=(img==0);					% create a boolean mask with 1 where there is metal
						% (was black in the image)

modlist={"Rel(3)=0.01"; "conn(2).Jph=0"};	% these are the *changes* for the front matalization
						% The resistance of the front electrode is reduced
						% and the photocurrent density is set to 0
						
% create the area index, and apply the front metalization mask
[Area_Index, AreaDef]=PVMOS_ModifyMask([],AreaDef, mask, modlist, "front_metal");



img=ReadImageIntensity("pinup_backisolation.png");% loads a bitmap image with where the back metalization is 
						% removed around the metal-wrap-through-holes
mask=(img==0);
modlist={"Rel(2)=1e10"};			% make electrode 2 isolating
[Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index, AreaDef, mask, modlist, "no_backmetal");



img=ReadImageIntensity("pinup_vias.png");	% loads a bitmap image with the vias
mask=(img==0);

% At the vias we do not only connect the front to the back contact but also connect to the contacting foil (electrode 1)
modlist={"Rel(2)=0.01"; "conn(2).model=\"JVD\""; "conn(2).VJ=[-1,-1e8;1,1e8]"; "conn(1).model=\"JVD\""; "conn(1).VJ=[-1,-1e8;1,1e8]"};
[Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index, AreaDef, mask, modlist, "vias_p-contact"); 



% This marks the spots where the back contact is connected to the contact foil
img=ReadImageIntensity("pinup_back2contactfoil.png");	% loads a bitmap image with the connections of the back to the contacting foil
mask=(img==0);
modlist={"conn(1).model=\"JVD\""; "conn(1).VJ=[-1,-1e8;1,1e8]"};
[Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index, AreaDef, mask, modlist, "p-contact");


img=ReadImageIntensity("pinup_foilisolation.png");	% loads a bitmap image with the foil isolation
mask=(img==0);
modlist={"Rel(1)=1e15"};
[Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index, AreaDef, mask, modlist, "contact_foil_isolation");


% This marks areas where the contact foil contacted to the positive electrode
img=ReadImageIntensity("pinup_contactfoil_vp.png");	% loads a bitmap image with the foil contact to the positive node
mask=(img==0);
modlist={"Rvn(1)=1e-8"};
[Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index, AreaDef, mask, modlist, "contact_foil_vp");


img=ReadImageIntensity("pinup_contactfoil_vn.png");	% loads a bitmap image with the foil contact to the ground node
mask=(img==0);
modlist={"Rvp(1)=1e-8"};
[Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index, AreaDef, mask, modlist, "contact_foil_vn");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate the mesh %%%%
Na=max(max(Area_Index));
Nx=length(Area_Index(1,:))
Ny=length(Area_Index(:,1))
% make a 10x10 cm^2 mesh
x=0:10/Nx:10;
y=0:10/Ny:10;
mkpvmosmesh(Na, Nel, Area_Index, AreaDef, x, y, "pinupmesh.bin")

