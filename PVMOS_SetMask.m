function [Area_Index, AreaDef]=PVMOS_SetMask(Area_Index,AreaDef, mask, MaskArea)
% this function helps to define regular PVMOS meshes.
	printf("===========================================================================\n");
	printf("Set mask %s\n",MaskArea.name);
	printf("---------------------------------------------------------------------------\n");
	% make sure mask is boolean type
	mask=(mask>0);
	if (length(Area_Index)==0)
		Area_Index=mask.+1;
		m=1;
	else
		m=max(max(Area_Index));
		Area_Index(mask)=m+1;		
	endif
	
	
	Na=m+1
	AreaDef(m+1)=MaskArea;
	[Area_Index, AreaDef]=PVMOS_CleanupAreas(Area_Index,AreaDef);
	printf("===========================================================================\n\n");
endfunction
