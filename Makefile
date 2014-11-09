# PhotoVoltaic MOdule Simulator (PVMOS)
src=main.c parse.c utils.c mesh2d.c solve.c list.c select_nodes.c dataexport.c diode.c
hdr=parse.h mesh2d.h parsedef.h utils.h main.h solve.h list.h select_nodes.h dataexport.h diode.h
obj=main.o parse.o utils.o mesh2d.o solve.o list.o select_nodes.o dataexport.o diode.o

CC=gcc
target=pvmos

# CFLAGS=-O3
CFLAGS=-Og -Wall -g -fPIC
LFLAGS= -lcholmod -lopenblas -lm
#LFLAGS= -lcholmod -L"/usr/local/cuda-5.5/targets/x86_64-linux/lib/" -L"/usr/lib64/nvidia-bumblebee/" -lcuda -lcudart -lcublas -lcufft -lm
VERSION=0.52

all: $(obj)
	$(CC) -o $(target)  $(obj) $(LFLAGS)
utils.o: utils.c mesh2d.h utils.h
	$(CC) $(CFLAGS)   -c -DVERSION=\"$(VERSION)\" utils.c
main.o: parse.h main.c main.h
parse.o: utils.h parsedef.h parse.h parse.c mesh2d.h solve.h list.h select_nodes.h dataexport.h
mesh2d.o: mesh2d.h mesh2d.c utils.h
solve.o: mesh2d.h utils.h solve.c diode.h
list.o: list.c utils.h
select_nodes.o: select_nodes.c list.h utils.h mesh2d.h
dataexport.o: dataexport.c list.h utils.h mesh2d.h
diode.o: diode.c diode.h
mkpvmosmesh: mesh2d.o utils.o list.o
	mkoctfile -v  mkpvmosmesh.cc mesh2d.o utils.o list.o
clean:
	-rm *.o *.oct $(target)
