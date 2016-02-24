function AreaDef=PVMOS_AreaProperties(Nel)
% generate the default material properties struct. Input: Nel, the number of electrodes.
	AreaDef.name="default";				% name of the area definition
	AreaDef.Rel=ones(1,Nel);			% electrode sheet resistance [Ohm]
	AreaDef.Rvp=-ones(1,Nel);			% Resistance to the positive node in [Ohm cm^2], if the value is negative there is no direct connection to the positive nde
	AreaDef.Rvn=-ones(1,Nel);			% Resistance to the negative node in [Ohm cm^2], if the value is negative there is no direct connection to the negative nde
	for i=1:Nel-1					% Define the connections between the electrodes
		AreaDef.conn(i).model="ONED";		% Connection model, can be ONED (one diode model) TWOD (two diode model) or JVD (Tabular JV data)
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for 1D, 2D and Phototransistor models
		AreaDef.conn(i).Jph=0; 			% Photo-current density current  [A/cm^2]
		AreaDef.conn(i).Rs=1e-5;   		% Series resistance              [Ohm cm^2]
		AreaDef.conn(i).Rsh=1e3;    		% Shunt resistance               [Ohm cm^2]		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for the 1D and 2D Models
		AreaDef.conn(i).J01=1e-12;		% Diode 1 saturation current density (used in ONED and TWOD) [A/cm^2]
		AreaDef.conn(i).J02=1e-8;		% Diode 2 saturation current density (used only in TWOD) [A/cm^2]
		AreaDef.conn(i).nid1=1; 		% Ideality factor for diode 1 (used in ONED and TWOD)
		AreaDef.conn(i).nid2=2;  		% Ideality factor for diode 2 (used in TWOD)
		AreaDef.conn(i).Eg=1.12;  		% Bandgap of the diode (used for temperature dependency, used in ONED and TWOD) [eV]
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for the Photo-Transistor model
		AreaDef.conn(i).Jsbc=5e-4;
		AreaDef.conn(i).Jsbe=5e-15;
		AreaDef.conn(i).Jse=0.0;
		AreaDef.conn(i).Jsc=0.0;
		AreaDef.conn(i).Nf=1.0;
		AreaDef.conn(i).Nr=1.0;
		AreaDef.conn(i).Ne=1.5;
		AreaDef.conn(i).Nc=2.0;
		AreaDef.conn(i).Vaf=inf;
		AreaDef.conn(i).Var=inf;
		AreaDef.conn(i).Bf=3.0;
		AreaDef.conn(i).EgBE=1.2;
		AreaDef.conn(i).PhiBC=0.3;
		AreaDef.conn(i).XTIBE=3.0;
		AreaDef.conn(i).XTIBC=3.0;
		AreaDef.conn(i).XTB=0.0;

		AreaDef.conn(i).N=2;     		% Number of JV datapoints (for JVD model)
		AreaDef.conn(i).VJ=[-1,-1;1,1];     	% Tabular JV data (for JVD model) Column 1: [V] Column 2: [A/cm^2]
	endfor
	AreaDef.T=300;   				% Initial temperature in the area
	AreaDef.SplitX=1;    				% Allow element to be split in x-direction during adaptive meshing (0 false, 1: true)
	AreaDef.SplitY=1;    				% Allow element to be split in y-direction during adaptive meshing (0 false, 1: true)
endfunction
