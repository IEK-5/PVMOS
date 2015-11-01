# Makefile for PhotoVoltaic MOdule Simulator (PVMOS)
src=main.c parse.c utils.c mesh2d.c solve.c list.c select_nodes.c dataexport.c diode.c phototransistor.c expr.c mkpvmosmesh.cc
hdr=parse.h mesh2d.h parsedef.h utils.h main.h solve.h list.h select_nodes.h dataexport.h diode.h phototransistor.h expr.h
obj=main.o parse.o utils.o mesh2d.o solve.o list.o select_nodes.o dataexport.o diode.o phototransistor.o expr.o

src_util=MeshHasher.c md5.c
hdr_util=md5.h
obj_util=md5.o

CC=gcc
target=pvmos
VERSION=0.78


CFLAGS=-Ofast -flto -Wall -fPIC
# CFLAGS=-Og -g -Wall -fPIC


LFLAGS= -flto -lcholmod -lopenblas-r0.2.14 -lm -lmatheval
# if linking to OpenBLAS we need to set the numthreads to 1
# Comment this if you are not using OpenBLAS or believe OpenBLAS should run on more threads (at least for OpenBLAS version 0.2.13 this is a bad idea)
# Seems that perhaps newer versions fix this problem, will test sometime...
OPENBLAS=""
WITH_LIBMATHEVAL=""

pvmos: newversion $(obj)
	$(CC) -o $(target)  $(obj) $(LFLAGS)
install: pvmos doc
	cp $(target) /usr/bin/
utils.o: utils.c mesh2d.h utils.h
	$(CC) $(CFLAGS)   -c -o utils.o -DVERSION=\"$(VERSION)\" utils.c
main.o: parse.h main.c main.h
ifdef OPENBLAS
	$(CC) $(CFLAGS)   -c -o main.o -DOPENBLAS main.c
else
	$(CC) $(CFLAGS)   -c main.c
endif	
parse.o: utils.h parsedef.h parse.h parse.c mesh2d.h solve.h list.h select_nodes.h dataexport.h expr.h diode.h phototransistor.h
mesh2d.o: mesh2d.h mesh2d.c utils.h meshhash.h
solve.o: mesh2d.h utils.h solve.c diode.h phototransistor.h
list.o: list.c utils.h
select_nodes.o: select_nodes.c list.h utils.h mesh2d.h
dataexport.o: dataexport.c list.h utils.h mesh2d.h diode.h phototransistor.h
diode.o: diode.c diode.h
expr.o: expr.c expr.h	
ifdef WITH_LIBMATHEVAL
	$(CC) $(CFLAGS)   -c -o expr.o -DWITH_LIBMATHEVAL expr.c
else
	$(CC) $(CFLAGS)   -c -o expr.o expr.c
endif	
meshhash.h: mesh2d.h mesh2d.c MeshHasher.c list.o utils.o md5.o diode.o phototransistor.o diode.h phototransistor.h
	# Make a hash (md5 sum) of a standard small mesh as a signature of the current mesh data structure
	# This hash is used to test compatibility of binary mesh files
	# we first make a dummy meshhash.h to compile a dummy mesh2d.o
	echo "#define _HAS_MESHHASH" > meshhash.h
	echo "int NMESHHASH=1;" >> meshhash.h
	echo "unsigned char MESHHASH[] = { 0 };" >> meshhash.h 
	$(CC) -Og -g -Wall -fPIC -lm   -c -o mesh2d.o mesh2d.c
	# Build the mesh hasher
	$(CC) -Og -g -Wall -fPIC -lm -g -Wall  -o MeshHasher mesh2d.o  utils.o list.o md5.o diode.o phototransistor.o MeshHasher.c
	# generate the hash
	./MeshHasher meshhash.h
	# mesh2d.o needs to be recompiled, with the newly created meshhash.h
	rm mesh2d.o
	# If this recipe fails the meshhash.h file is still created but with a bogus hash
	# Thus we need to make sure the meshhash.h file is deleted if this recipe fails. Turns out there
	# is a special target for that, .DELETE_ON_ERROR. Include it to make sure mesh2d.c is never compiled with
	# the dummy meshhash.h 
.DELETE_ON_ERROR:
mkpvmosmesh: mkpvmosmesh.pkg
	octave --eval "pkg install pvmos-mesh-$(VERSION).tar.gz"
mkpvmosmesh.pkg: mesh2d.c utils.c list.c main.h mesh2d.h utils.h list.h meshhash.h phototransistor.h diode.h  phototransistor.c diode.c 
	mkdir -p pvmos-mesh-$(VERSION)/src
	mkdir -p pvmos-mesh-$(VERSION)/inst
	cp PVMOS_*.m pvmos-mesh-$(VERSION)/inst/
	sed -i 's/^VERSION[^\n]\+/VERSION=$(VERSION)/g' Makefile_mkpvmosmesh
	cp mkpvmosmesh.cc diode.c diode.h phototransistor.c phototransistor.h mesh2d.c utils.c list.c main.h mesh2d.h utils.h list.h meshhash.h pvmos-mesh-$(VERSION)/src/
	cp Makefile_mkpvmosmesh pvmos-mesh-$(VERSION)/src/Makefile
	echo "Name: pvmos-mesh" > pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Version: $(VERSION)" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Date: $(shell date)" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Author: B.E.Pieters" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Maintainer: B.E.Pieters" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Categories: Create PVMOS meshes in Octave >>mkpvmosxmesh" >>  pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Title: pvmos-mesh" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "Description: A mesh generator for PVMOS" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "License: GPLv3+" >> pvmos-mesh-$(VERSION)/DESCRIPTION
	echo "The GPLv3+ applies to all files in thic package" > pvmos-mesh-$(VERSION)/COPYING
	tar -zcvf pvmos-mesh-$(VERSION).tar.gz pvmos-mesh-$(VERSION)
	rm -rf  pvmos-mesh-$(VERSION)
newversion:
	sed -i 's/version [0-9\.]\+/version $(VERSION)/g' README.md
doc:
	cd Doc/;echo "\newcommand{\version}{$(VERSION)}" > version.tex
	cd Doc/;pdflatex PVMOS_manual.tex
	cd Doc/;pdflatex PVMOS_manual.tex	
cleancopy:
	mkdir -p CleanCopy
	cp $(src) CleanCopy
	cp $(hdr) CleanCopy
	cp $(src_util) CleanCopy
	cp $(hdr_util) CleanCopy
	cp README.md CleanCopy
	cp Makefile CleanCopy
clean:
	-rm *.o *.oct $(target) pvmos-mesh-$(VERSION).tar.gz MeshHasher
