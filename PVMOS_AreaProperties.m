% The Area Properties Structure:
% Member		Default Value	Unit		Description:	
% ---------------------------------------------------------------------------
% name		"default"	(string)	Area Name 
% Rel  		1 		[Ohm]		Array of length Nel with 
% 						the electrode sheet 
% 						resistances
% Rvp		-1		[Ohm cm^2]	Array of length Nel with 
% 						the electrode resistances 
% 						to the positive node. 
% 						No connection is made when 
% 						this value is negative
% Rvn		-1		[Ohm cm^2]	Array of length Nel with 
% 						the electrode resistances 
% 						to the ground node. No
% 						connection is made when 
% 						this value is negative
% T		300		[K]		Temperature in the area
% SplitX		1		[-]		Boolean value whether 
% 						the element may be split 
% 						in the x-direction when 
% 						adaptive the mesh 
% 						(0 false, 1: true)
% SplitY		1		[-]		Boolean value whether 
% 						the element may be split 
% 						in the y-direction when 
% 						adaptive the mesh 
% 						(0 false, 1: true)
% conn		-		(structure)	Structure array of length
% 						Nel-1 describing the 
% 						inter-electrode connections
% 
% The \"conn\" Structure (inter-electrode connection): 
% Member		Default Value	Unit		Description:	
% ---------------------------------------------------------------------------
% model		"ONED"		(string)	String describing the 
% 						connection model. The 
% 						following models are 
% 						available:
% 						"ONED"	 One Diode Model
% 						"TWOD"	 Two Diode Model
% 						"PHOTOT" Photo-Transistor
% 							 Model
% 						"JVD"	 Tabular JV data
% Jph		0		[A/cm^2]	Photo-current density. 
% 						Used in models:
% 						ONED, TWOD, PHOTOT
% Rs		1e-5		[Ohm cm^2]	Series resistance.
% 						Used in models:
% 						ONED, TWOD, PHOTOT
% Rsh		1e3		[Ohm cm^2]	Shunt resistance
% 						Used in models:
% 						ONED, TWOD, PHOTOT		
% J01		1e-12		[A/cm^2]	Diode 1 saturation current
% 						density. Used in models:
% 						ONED, TWOD	
% J02		1e-8		[A/cm^2]	Diode 2 saturation current
% 						density. Used in models:
% 						TWOD
% nid1		1		[-]		Diode 1 ideality factor.
% 						Used in models:
% 						ONED, TWOD
% nid2		2		[-]		Diode 2 ideality factor.
% 						Used in models:
% 						TWOD
% Eg		1.12		[eV]		Bandgap of the diode 
% 						(temperature dependencies).
% 						Used in models:
% 						ONED, TWOD
% Jsbc		5e-4		[A/cm^2]	BC junction saturation 
% 						current density. Used in 
% 						models:
% 						PHOTOT
% Jsbe		5e-15		[A/cm^2]	BE junction saturation 
% 						current density. Used in 
% 						models:
% 						PHOTOT
% Jse		0.0		[A/cm^2]	BE leakage junction 
% 						saturation current density.
% 						Used in models:
% 						PHOTOT
% Jsc		0.0		[A/cm^2]	BC leakage junction 
% 						saturation current density.
% 						Used in models:
% 						PHOTOT
% Nf		1.0		[-]		BE junction ideality factor
% 						Used in models:
% 						PHOTOT
% Nr		1.0		[-]		BC junction ideality factor
% 						Used in models:
% 						PHOTOT
% Ne		1.0		[-]		BE leakage junction 
% 						ideality factor. Used in 
% 						models:
% 						PHOTOT
% Nc		1.0		[-]		BC leakage junction 
% 						ideality factor. Used in 
% 						models:
% 						PHOTOT
% Vaf		inf		[V]		Forward operation Early 
% 						voltage. Used in models:
% 						PHOTOT
% Var		inf		[V]		Reverse operation Early 
% 						voltage. Used in models:
% 						PHOTOT
% Bf		3.0		[-]		Forward operation gain
% EgBE		1.2		[eV]		Band gap BE junction. 
% 						Used in models:
% 						PHOTOT
% PhiBC		0.3		[eV]		BC junction barrier height
% XTIBE		3.0		[-]		Temperature dependency of 
% 						the BE junction saturation 
% 						current density (T^XTIBE)
% 						Used in models:
% 						PHOTOT
% XTIBC		3.0		[-]		Temperature dependency of 
% 						the BC junction saturation 
% 						current density (T^XTIBC)
% 						Used in models:
% 						PHOTOT
% XTB		0.0		[-]		Temperature dependency of 
% 						the forward gain (T^XTB)
% 						Used in models:
% 						PHOTOT
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
