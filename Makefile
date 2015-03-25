# PhotoVoltaic MOdule Simulator (PVMOS)
src=main.c parse.c utils.c mesh2d.c solve.c list.c select_nodes.c dataexport.c diode.c expr.c
hdr=parse.h mesh2d.h parsedef.h utils.h main.h solve.h list.h select_nodes.h dataexport.h diode.h expr.h
obj=main.o parse.o utils.o mesh2d.o solve.o list.o select_nodes.o dataexport.o diode.o expr.o

CC=gcc
target=pvmos

# CFLAGS=-O3
# CFLAGS=-Ofast -flto -Wall -fPIC
CFLAGS=-Og -g -Wall -fPIC
LFLAGS= -flto -lcholmod -lopenblas-r0.2.13 -lm -lmatheval
# if linking to OpenBLAS we need to set the numthreads to 1
# Uncomment this if you are not using OpenBLAS or believe OpenBLAS should run on more threads (at least for OpenBLAS version 0.2.13 this is a bad idea)
OPENBLAS=""

# CFLAGS=-Og -g -Wall -fPIC
# LFLAGS=-lcholmod -lopenblas-r0.2.13 -lm
# LFLAGS= -lcholmod -L/usr/lib64/libblas.so.3 -lm
#LFLAGS= -lcholmod -L"/usr/local/cuda-5.5/targets/x86_64-linux/lib/" -L"/usr/lib64/nvidia-bumblebee/" -lcuda -lcudart -lcublas -lcufft -lm
VERSION=0.63

pvmos: $(obj)
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
mesh2d.o: mesh2d.h mesh2d.c utils.h
solve.o: mesh2d.h utils.h solve.c diode.h
list.o: list.c utils.h
select_nodes.o: select_nodes.c list.h utils.h mesh2d.h
dataexport.o: dataexport.c list.h utils.h mesh2d.h
diode.o: diode.c diode.h
expr.o: expr.c expr.h
mkpvmosmesh: mesh2d.o utils.o list.o
	mkoctfile -v  mkpvmosmesh.cc mesh2d.o utils.o list.o
clean:
	-rm *.o *.oct $(target)
