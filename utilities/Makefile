# PhotoVoltaic MOdule Simulator (PVMOS)
CC=gcc
target=polytransform

# CFLAGS=-O3
CFLAGS=-Og -Wall -g
LFLAGS= -lm
VERSION=0.1

polytransform: polytransform.c
	$(CC) -o $(target)  polytransform.c $(LFLAGS)
install: polytransform
	cp $(target) /usr/bin/
clean:
	-rm $(target)
