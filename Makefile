# PhotoVoltaic MOdule Simulator (PVMOS)
src=main.c parse.c utils.c mesh2d.c solve.c list.c select_nodes.c dataexport.c diode.c expr.c mkpvmosmesh.cc
hdr=parse.h mesh2d.h parsedef.h utils.h main.h solve.h list.h select_nodes.h dataexport.h diode.h expr.h
obj=main.o parse.o utils.o mesh2d.o solve.o list.o select_nodes.o dataexport.o diode.o expr.o

CC=gcc
target=pvmos

# CFLAGS=-O3
CFLAGS=-Ofast -flto -Wall -fPIC
# CFLAGS=-Og -g -Wall -fPIC
LFLAGS= -flto -lcholmod -lopenblas-r0.2.14 -lm -lmatheval
# if linking to OpenBLAS we need to set the numthreads to 1
# Uncomment this if you are not using OpenBLAS or believe OpenBLAS should run on more threads (at least for OpenBLAS version 0.2.13 this is a bad idea)
OPENBLAS=""
WITH_LIBMATHEVAL=""

# CFLAGS=-Og -g -Wall -fPIC
# LFLAGS=-lcholmod -lopenblas-r0.2.13 -lm
# LFLAGS= -lcholmod -L/usr/lib64/libblas.so.3 -lm
#LFLAGS= -lcholmod -L"/usr/local/cuda-5.5/targets/x86_64-linux/lib/" -L"/usr/lib64/nvidia-bumblebee/" -lcuda -lcudart -lcublas -lcufft -lm
VERSION=0.73
Y=$(shell date +%Y)
OY=$(shell cat README.md |grep "Dr. Bart E. Pieters"|egrep -o '201[0-9]')
NEWYEAR=$(shell [ "$(OY)" != "$(Y)" ] && echo true)
pvmos: newversion newyear $(obj)
	$(CC) -o $(target)  $(obj) $(LFLAGS)
install: pvmos
	cp $(target) /usr/bin/
utils.o: utils.c mesh2d.h utils.h
	$(CC) $(CFLAGS)   -c -o utils.o -DVERSION=\"$(VERSION)\" utils.c
main.o: parse.h main.c main.h
ifdef OPENBLAS
	$(CC) $(CFLAGS)   -c -o main.o -DOPENBLAS main.c
else
	$(CC) $(CFLAGS)   -c main.c
endif	
parse.o: utils.h parsedef.h parse.h parse.c mesh2d.h solve.h list.h select_nodes.h dataexport.h expr.h
mesh2d.o: mesh2d.h mesh2d.c utils.h meshhash.h
meshhash.h: mesh2d.h MeshHasher.c list.o utils.o 
	# this part is only executed if the mesh data structures are changed
	# First create a dummy meshhash.h file to make a dummy mesh2d.o file
	echo "#define _HAS_MESHHASH" > meshhash.h
	echo "int NMESHHASH=1;" >> meshhash.h
	echo "unsigned char MESHHASH[] = { 0 };" >> meshhash.h 
	$(CC) $(CFLAGS)   -c -o mesh2d.o mesh2d.c
	# Build the mesh hasher
	$(CC) $(CFLAGS) -lssl -lcrypto -g -Wall  -o MeshHasher mesh2d.o  utils.o list.o MeshHasher.c
	# generate the hash
	./MeshHasher meshhash.h
	# remove the dummy mesh2d.o so it is properly compiled after this.
	rm mesh2d.o		
solve.o: mesh2d.h utils.h solve.c diode.h
list.o: list.c utils.h
select_nodes.o: select_nodes.c list.h utils.h mesh2d.h
dataexport.o: dataexport.c list.h utils.h mesh2d.h
diode.o: diode.c diode.h
expr.o: expr.c expr.h
ifdef WITH_LIBMATHEVAL
	$(CC) $(CFLAGS)   -c -o expr.o -DWITH_LIBMATHEVAL expr.c
else
	$(CC) $(CFLAGS)   -c -o expr.o expr.c
endif	
mkpvmosmesh: mesh2d.o utils.o list.o main.h
	mkoctfile -v  mkpvmosmesh.cc mesh2d.o utils.o list.o
mkpvmosmesh.pkg: mesh2d.c utils.c list.c main.h mesh2d.h utils.h list.h 
	mkdir -p mkpvmosmesh-$(VERSION)/src
	sed -i 's/^VERSION[^\n]\+/VERSION=$(VERSION)/g' Makefile_mkpvmosmesh
	cp mkpvmosmesh.cc mesh2d.c utils.c list.c main.h mesh2d.h utils.h list.h mkpvmosmesh-$(VERSION)/src/
	cp Makefile_mkpvmosmesh mkpvmosmesh-$(VERSION)/src/Makefile
	echo "Name: mkpvmosmesh" > mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Version: $(VERSION)" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Date: $(shell date)" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Author: B.E.Pieters" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Maintainer: B.E.Pieters" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Categories: Create PVMOS meshes in Octave >>mkpvmosxmesh" >>  mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Title: mkpvmosmesh" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "Description: A mesh generator for PVMOS" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "License: GPLv3+" >> mkpvmosmesh-$(VERSION)/DESCRIPTION
	echo "The GPLv3+ applies to all files in thic package" > mkpvmosmesh-$(VERSION)/COPYING
	tar -zcvf mkpvmosmesh-$(VERSION).tar.gz mkpvmosmesh-$(VERSION)
	rm -rf  mkpvmosmesh-$(VERSION)

newversion:
	sed -i 's/version [^\n]\+/version $(VERSION)/g' README.md
newyear:
ifeq ($(NEWYEAR),true)
	sed -i 's/Dr. Bart E. Pieters 201[0-9]\+/Dr. Bart E. Pieters $(Y)/g' README.md
	find . -maxdepth 1 -name '*.[ch]' -exec sed -i 's/Dr. Bart E. Pieters 201[0-9]\+/Dr. Bart E. Pieters $(Y)/g' {} \;
endif
cleancopy:
	mkdir -p CleanCopy
	cp $(src) CleanCopy
	cp $(hdr) CleanCopy
	cp README.md CleanCopy
	cp Makefile CleanCopy
clean:

	-rm *.o *.oct $(target) mkpvmosmesh-$(VERSION).tar.gz MeshHasher
