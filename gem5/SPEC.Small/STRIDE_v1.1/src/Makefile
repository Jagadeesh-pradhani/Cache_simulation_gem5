#
#          @@@@@      @@@      @@@@@   @       @ @@@@@@@@@
#         @     @    @   @    @     @  @       @ @
#       @       @  @     @  @       @ @       @ @
#       @         @       @ @         @@@@@@@@@ @@@@@@@
#       @       @ @@@@@@@@@ @       @ @       @ @
#         @     @  @       @  @     @  @       @ @
#          @@@@@   @       @   @@@@@   @       @ @@@@@@@@@
#
#
#  This  makefile makes the Livermore Computing cache trashing benchmarks.
#       When compiling these benchmarks please generate  a  compiler
#       listing (if it   is  possible  :-) and  include   these when
#       returning the results.
#
#               Send results to:
#                       Dr. Mark K. Seager
#                       LLNL/Integrated Computing and Communications Dept
#                       PO Box 808, L-554
#                       Livermore, CA 94551-9900
#                       seager@llnl.gov
#                       925-423-3141-Voice
#                       925-423-8911-Fax
#
# What shell to use...
SHELL=/bin/bash
#
# What are your Fortran 77 and C compilers called?
# INTEL COMPILERS
#FC=ifort
#CC=icc
# GNU COMPILERS
#FC=f77
#CC=gcc

SYSROOT=$(ARMTOOLCHAINDIR)/gcc-linaro-5.4.1-2017.05-x86_64_aarch64-linux-gnu/aarch64-linux-gnu/libc
CC=$(ARMTOOLCHAINDIR)/gcc-linaro-5.4.1-2017.05-x86_64_aarch64-linux-gnu/bin/aarch64-linux-gnu-gcc

# IBM COMPILERS
#FC=mpxlf -qlist -qsource
#FC=xlf -qlist -qsource
#CC=mpxlc
#CC=xlc
# PathScale compilers
#FC=mpipathf90
#CC=pathcc
#
# Use the highest optimization flags your compiler can stand...
#GFLAGS=-g -ffast-math -funroll-loops
#GCFLAGS=-g -ffast-math -DLINUX
# Itanium2 g
#GFLAGS=-g -tpp2 -pad
#GCFLAGS=-g -tpp2 -pad -DLINUX
# Pentium4 with SSE2
#GFLAGS=-g -tpp7 -xW -pad -Zp16
#GCFLAGS=-g -tpp7 -xW -align -Zp16 -DLINUX
# f77 and GCC for Pentium 4 with SSE2
#FFLAGS=-O3 -ffast-math -funroll-loops -fexpensive-optimizations -march=i686 -mcpu=i686 -msse2 -malign-double
#CFLAGS=-O3 -ffast-math -DLINUX -march=i686 -mcpu=i686 -msse2 -malign-double
# Directives for INTEL 8.1 Compiler on IA-64 with Linux
#FFLAGS=-O3 -tpp2 -IPF_fma -unroll -pad -Zp16 -opt_report -opt_report_phaseall -vec_report5
#CFLAGS=-O3 -tpp2 -IPF_fma -unroll -Zp16 -DLINUX
# Directives for INTEL Compiler on P4 USE SSE2 with Linux
#FFLAGS=-O3 -tpp7 -xW -pad -Zp16 -opt_report
#CFLAGS=-O3 -tpp7 -xW -align -march=pentium4 -mcpu=pentium4 -Zp16 -DLINUX -opt_report -opt_report_phaseall -vec_report3
# Directives for INTEL Compiler on P4 DONT USE SSE2 with Linux
#FFLAGS=-O3 -tpp7 -pad -Zp16
#CFLAGS=-O3 -tpp7 -Zp16 -DLINUX
# Directives for INTEL Compiler on P3 with Linux
#FFLAGS=-O3 -tpp6 -xK -pad -Zp16
#CFLAGS=-O3 -tpp6 -xK -Zp16
#
# Directives for BGL systems
#FFLAGS=-O4 -qhot -qtune=440 -qarch=440d
#CFLAGS=-O3 -D_BGL -qtune=440 -qarch=440d
#GFLAGS=-g -qtune=440 -qarch=440d -D
#GCFLAGS=-g -D_BGL -qtune=440 -qarch=440d
# If you need any special load flags put them here...
#LDFLAGS=-L/bgl/local/lib -lmassv -lm
#LDFFLAGS=-L/bgl/local/lib -lmassv
#LDFFLAGS=-L/bglhome/bwallen/libs/lib -lmassv
LDCFLAGS=-lm
#
# Directives for IBM Aix Power5 SMP systems
#FFLAGS=-O4 -qhot -qtune=pwr5 -qarch=pwr5 -DAIX
#FFLAGS=-O4 -qhot -qtune=pwr5 -qarch=pwr5 -WF,-DAIX
#FFLAGS=-O2  -qtune=pwr5 -qarch=pwr5 -WF,-DAIX -qmaxmem=-1
#CFLAGS=-O3 -DAIX -qtune=pwr5 -qarch=pwr5
#CFLAGS=-O2 -DAIX -qtune=pwr5 -qarch=pwr5
#GFLAGS=-g -qtune=pwr5 -qarch=pwr5 -DAIX
#GFLAGS=-g -qtune=pwr5 -qarch=pwr5 -WF,-DAIX
#GCFLAGS=-g -DAIX -qtune=pwr5 -qarch=pwr5
# Directives for IBM Aix Power4 SMP systems
#FFLAGS=-O4 -qhot -qsmp=omp -qreport=smplist -qtune=pwr4 -qarch=pwr4 -DAIX
#CFLAGS=-O3 -DAIX -qsmp  -qtune=pwr4 -qarch=pwr4 -DDGEMM -DMXMD16 -D_OPENMP
#FFLAGS=-O4 -qhot
#FFLAGS=-O5 -fast -dalign -fns -xprofile=collect
#FFLAGS=-O -qhsflt -Pv
#FFLAGS=-O2 -dalign -Mvect -Minfo=loop -Mx,19,1
# If you need any special load flags put them here...
#LDFLAGS=-lm
#
# Directives for PathScale Compilers
FFLAGS=-O3 -DLINUX
CFLAGS=-O3 -DLINUX
GFFLAGS=-g -DLINUX
GCFLAGS=-g -DLINUX

# Things to make and run...
STRID3=strid3c.Opt strid3c.g
VECOP=vecopc.Opt vecopc.g
CACHE=cachec.Opt cachec.g
STRIDDOT=striddot.Opt striddot.g
CACHEDOT=cachedot.Opt cachedot.g

TARGETS=$(STRID3) $(VECOP) $(CACHE)

.f.o:
	$(FC) ${FFLAGS} -c $*.f

all : $(TARGETS)

clean:
	rm -f *.o $(TARGETS) *.out core *.list *.l
#
#    Optimized Codes...
#
#
strid3.Opt : strid3.f secf.o
	$(FC) $(FFLAGS) strid3.f secf.o -o strid3.Opt $(LDFFLAGS)
strid3c.Opt : strid3.c secc.o
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static strid3.c secc.o -o strid3c.Opt $(LDCFLAGS)
vecop.Opt : vecop.f secf.o
	$(FC) $(FFLAGS) vecop.f secf.o -o vecop.Opt $(LDFFLAGS)
vecopc.Opt : vecop.c secc.o
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static vecop.c secc.o -o vecopc.Opt $(LDCFLAGS)
cache.Opt : cache.f secf.o
	$(FC) $(FFLAGS) cache.f secf.o -o cache.Opt $(LDFFLAGS)
cachec.Opt : cache.c secc.o
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static cache.c secc.o -o cachec.Opt $(LDCFLAGS)
striddot.Opt : striddot.f secf.o
	$(FC) $(FFLAGS) striddot.f secf.o -o striddot.Opt $(LDFFLAGS)
cachedot.Opt : cachedot.f secf.o
	$(FC) $(FFLAGS) cachedot.f secf.o -o cachedot.Opt $(LDFFLAGS)
#
#    Unoptimized codes...
#
strid3.g : strid3.f secf.o
	$(FC) $(GFLAGS) strid3.f secf.o -o strid3.g $(LDFFLAGS)
strid3c.g : strid3.c secc.o
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static strid3.c secc.o -o strid3c.g $(LDCFLAGS)
vecop.g : vecop.f secf.o
	$(FC) $(GFLAGS) vecop.f secf.o -o vecop.g $(LDFFLAGS)
vecopc.g : vecop.c secc.o
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static vecop.c secc.o -o vecopc.g $(LDCFLAGS)
cache.g : cache.f secf.o
	$(FC) $(GFLAGS) cache.f secf.o -o cache.g $(LDFFLAGS)
cachec.g : cache.c secc.o
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static cache.c secc.o -o cachec.g $(LDCFLAGS)
striddot.g : striddot.f secf.o
	$(FC) $(GFLAGS) striddot.f secf.o -o striddot.g $(LDFFLAGS)
cachedot.g : cachedot.f secf.o
	$(FC) $(GFLAGS) cachedot.f secf.o -o cachedot.g $(LDFFLAGS)

cache.sh : READ.ME makefile cache.f cache.c strid3.f strid3.c \
       vecop.f vecop.c secc.c secf.c run.cache run.stride
	shar READ.ME makefile cache.f cache.c strid3.f strid3.c \
               vecop.f vecop.c striddot.f cachedot.f \
               secc.c secf.c run.cache run.stride > cache.sh
