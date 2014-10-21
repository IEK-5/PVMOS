#!/usr/bin/octave -q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example ThinFilm Module:
% In this example we simulate a 40x40 cm^2 thin-film a-Si:H module with 40 cell stripes and  
% a random distribution of defect. 
%
% The PVMOS script thin_film40x40.mos creates a mesh for a thin-film module *without* defects 
% The script in this file serves to take this mesh and shoot defects in it.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our module has two electrodes, TCO with 5 Ohm and a metal back contact with 0.5 Ohm
% Note that these values are, somewhat inconveniently, also  defined in thin_film40x40.mos, i.e.
% if you want to change these values you have to do so at two locations.
% The sheet resistances of these two electrodes are stored in the R array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global R=[5,0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                            %%
%%                                          Functions:                                        %%
%%                                                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes a PVMOS input file including shunts
% input:
% fn: filename
% shunts matrix, each row specifies a local defect
% column1 x coordinate
% column2 y coordinate
% column3 radius
% column4 resistance [Ohm cm^2]
% output: a PVMOS input file in "fn"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteMOSFile(fn, shunts)
	global R;

	f=fopen(fn, "w");
	
	% load the defect free mesh generated with the thin_film.mos PVMOS script
	fprintf(f,"loadmesh module_defectfree.bin module\n");	
	
	% we loop through the shunts to add em, if shunts is empty we do a defect free simulation
	if (length(shunts))	
		% shunts are defined
		for i=1:length(shunts(:,1))
			% create a new area for the i-th shunt and set the parameters accordingly
			fprintf(f,"set_R\t\tmodule.shunt%i 	0 %e\n", i, shunts(i,4));
			for j=1:length(R)
				fprintf(f,"set_Rel\t\tmodule.shunt%i %i %e\n", i, j-1, R(j));
			endfor
			
			% we refine the mesh locally to add the shunt
			% to this end we select several circular areas and specify the minimal resolution in these circles
			% radius of the area around the shunt
			rr=2;
			% maximum element length or width, all elemnts which exceed this length or width are split in smaller 
			% elements accordingly
			d=0.2;
			fprintf(f,"select_circ %e %e %e module\n", shunts(i,1), shunts(i,2), rr);
			fprintf(f,"split_coarse module %e\n", d);
			% make sure the radius of the shunt is at least 10 elements long
			% to make a somewhat reasonable definition of a circle
			while (shunts(i,3)/d<10)
				rr/=2;
				d/=2;
				fprintf(f,"select_circ %e %e %e module\n", shunts(i,1), shunts(i,2), rr);
				fprintf(f,"split_coarse module %e\n", d);
			endwhile
			% select and assign the area that is to be shunted
			fprintf(f,"select_circ %e %e %e module\n", shunts(i,1), shunts(i,2), shunts(i,3));
			fprintf(f,"assign_properties module.shunt%i\n", i);
		endfor
		fprintf(f,"solve 	module 	0 25 	3\n");
		% adapt the mesh at 25V
		fprintf(f,"adaptive_solve 	module 	25 	0.3 	2\n");
	else
		i=0;
	endif
	% solve IV from 0-35 V in 36 steps (1V steps including both 0 and 35 V)
	fprintf(f,"solve module 0 35 36\n");
	fprintf(f,"surfVplot module 0 0 40 40 1000 1000 25 module%i_25.dat\n", i);
	fprintf(f,"surfVplot module 0 0 40 40 1000 1000 0 module%i_00.dat\n", i);
	fprintf(f,"surfVplot module 0 0 40 40 1000 1000 30 module%i_30.dat\n", i);
	fprintf(f,"printIV module IV%i.dat\n",i);
	fclose(f);
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random shunts to feed to WriteMOSFile
% input:
% N: Number of shunts to be generated
% w: width of the module
% l: length of the module
% r1-r2: range of defect radii
% R1-R2: range of defect resistances in Ohm cm2
% output: shunt matrix with N rows and where
% column1 x coordinate
% column2 y coordinate
% column3 radius
% column4 resistance [Ohm cm^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shunts=RandomShunt(N, w, l, r1, r2, R1, R2)
	shunts=[];
	for i=1:N
		x=unifrnd (0, w);
		y=unifrnd (0, l);
		r=unifrnd (r1, r2);
		R=unifrnd (R1, R2);
		shunts=[shunts;x,y,r,R];
	endfor
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                            %%
%%                                          Doin' Stuff:                                      %%
%%                                                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run PVMOS to generate the defect free mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system("../../pvmos thin_film40x40.mos");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a 100 shunts: 40x40, radii between 0.01-0.3 cm, resistances of 0.1-100 Ohm cm^2
% Note that for a 100 defects you should expect some calculation time and memory usage! 
% Due to the random nature of the defects unexpected things may happen...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsh=100;
shunts=RandomShunt(Nsh, 40, 40, 0.01, 0.3, 1e-1, 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save shunt distribution for future reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save shunts.dat shunts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation, throw out some defects, re-run till no more defects left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (Nsh>1)
	WriteMOSFile("thin_film_shunted.mos", shunts);
	system("./pvmos thin_film_shunted.mos");
	Nsh=floor(Nsh/2+0.5)
	shunts=shunts(1:Nsh,:);
endwhile

WriteMOSFile("thin_film_shunted.mos", shunts);
system("./pvmos thin_film_shunted.mos");	
WriteMOSFile("thin_film_shunted.mos", []);
system("./pvmos thin_film_shunted.mos");	
