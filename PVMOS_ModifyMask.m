function [Area_Index, AreaDef]=PVMOS_ModifyMask(Area_Index,AreaDef, mask, modlist, name)
% this function helps to define regular PVMOS meshes.
	printf("===========================================================================\n");
	printf("Modify mask %s\n",name);
	printf("---------------------------------------------------------------------------\n");
	for j=1:length(modlist)
		printf("%s;\n",modlist{j})
	endfor
	% make sure mask is boolean type
	mask=(mask>0);
	if (length(Area_Index)==0)
		Area_Index=mask.+1;
		m=1;
	else
		m=max(max(Area_Index));
		mask.*=m;
		Area_Index.+=mask;
		
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
	[Area_Index, AreaDef]=PVMOS_CleanupAreas(Area_Index,AreaDef);
	for i=m+1:max(max(Area_Index))
		printf("New Area: %s\n",AreaDef(i).name);
	endfor
	printf("Created %i new areas\n",max(max(Area_Index))-m);
	printf("===========================================================================\n\n");
endfunction
