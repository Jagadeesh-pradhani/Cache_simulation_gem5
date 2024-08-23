###############################################################################
# SPECtools2000
#   Version 0.3
#   This makefile should NOT be changed other than to customize for a 
#   particular benchmark.  Users may overide variables by placing them
#   in the Makefile.spec file.  Or by using the "runspec" program to 
#   build the Makefile.spec file for them.
###############################################################################
all:
	#gcc -fopenmp -o quake quake.c -lm
	gcc -m32 -o quake quake.c -lm
run:
	./quake < inp.in
