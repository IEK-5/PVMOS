all:
	make -C ../../ mkpvmosmesh
	cp ../../mkpvmosmesh.oct .
	./generate_mesh.m
	../../pvmos pinup.mos
clean:
	rm *.dat mkpvmosmesh.oct *.bin
