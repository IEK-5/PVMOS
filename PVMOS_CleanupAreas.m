function [Area_Index, AreaDef]=PVMOS_CleanupAreas(Area_Index,AreaDef)
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
