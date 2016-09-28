# what flags you want to pass to the C compiler & linker
# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile

bin_PROGRAMS = pvmos polytransform
pvmos_SOURCES = main.c parse.c utils.c mesh2d.c solve.c list.c polygon.c select_nodes.c dataexport.c diode.c phototransistor.c expr.c parse.h mesh2d.h parsedef.h utils.h main.h solve.h list.h polygon.h select_nodes.h dataexport.h diode.h phototransistor.h expr.h meshhash.h
polytransform_SOURCES=polytransform.c
EXTRA_pvmos_SOURCES = malloc-2.8.6.c
EXTRA_DIST = PVMOS_AreaProperties.m PVMOS_CleanupAreas.m PVMOS_ModifyMask.m PVMOS_SetMask.m mkpvmosmesh.cc Makefile_pvmosmesh.in

BUILT_SOURCES = meshhash.h
noinst_PROGRAMS = MeshHasher
MeshHasher_SOURCES = md5.c md5.h mesh2d.c mesh2d.h  utils.c utils.h list.c list.h diode.c diode.h phototransistor.c phototransistor.h MeshHasher.c	

meshhash.h: MeshHasher$(EXEEXT) $(MeshHasher_SOURCES)
	./MeshHasher $@

MeshHasher$(EXEEXT): $(MeshHasher_SOURCES) utils.o list.o md5.o diode.o phototransistor.o
	echo '#define _HAS_MESHHASH' > meshhash.h
	echo 'int NMESHHASH=1;' >> meshhash.h
	echo 'unsigned char MESHHASH[] = { 0 };' >> meshhash.h
	$(CC) -Og -g -Wall -fPIC -c -o mesh2d.o $(srcdir)/mesh2d.c
	$(CC) -Og -g -Wall -fPIC -g -Wall  -o MeshHasher mesh2d.o  utils.o list.o md5.o diode.o phototransistor.o $(srcdir)/MeshHasher.c -lm
	rm mesh2d.o


if DMALLOC
# my somewhat hackish solution to getting dmalloc compiled and linked to
pvmos_LDADD = malloc-2.8.6.o
malloc-2.8.6.o:
	$(CC) -O  -finline-functions -frename-registers -foptimize-sibling-calls  -fstrength-reduce  -fcse-follow-jumps  -fcse-skip-blocks  -frerun-cse-after-loop  -frerun-loop-opt  -fgcse   -fgcse-lm   -fgcse-sm  -fdelete-null-pointer-checks  -fexpensive-optimizations  -fregmove  -fschedule-insns  -fschedule-insns2  -fsched-interblock -fsched-spec  -fcaller-saves  -fpeephole2  -freorder-blocks  -freorder-functions  -fstrict-aliasing  -falign-functions  -falign-jumps  -falign-loops  -falign-labels    -Wall -fPIC  -c -o malloc-2.8.6.o malloc-2.8.6.c
endif

if PVMOSMESH
dist_pkgdata_DATA = octave-pvmos-mesh-$(VERSION).0.tar.gz
Makefile_pvmosmesh: $(srcdir)/Makefile_pvmosmesh.in
	sed 's|SRCDIR|'$(srcdir)'|g' $(srcdir)/Makefile_pvmosmesh.in > Makefile_pvmosmesh

octave-pvmos-mesh-$(VERSION).0.tar.gz: meshhash.h Makefile_pvmosmesh
	sed -i 's/^VERSION[^\n]\+/VERSION=$(VERSION)/g' Makefile_pvmosmesh; \
	make -f Makefile_pvmosmesh
endif

clean-local:
	-rm -f meshhash.h pvmos-mesh-$(VERSION).tar.gz Makefile_pvmosmesh
