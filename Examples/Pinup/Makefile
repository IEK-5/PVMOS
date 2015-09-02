all:
	make -C ../../ mkpvmosmesh
	./generate_mesh.m
	../../pvmos pinup.mos
clean:
	rm *.dat mkpvmosmesh.oct *.bin
