###############################################################################
# SPECtools2000
#   Version 0.3
#   This makefile should NOT be changed other than to customize for a 
#   particular benchmark.  Users may overide variables by placing them
#   in the Makefile.spec file.  Or by using the "runspec" program to 
#   build the Makefile.spec file for them.
###############################################################################
#	BENCHMARK	-- standard definitions for this benchmark

all:

	f95 -m32 -c -o blts.o                blts.f
	f95 -m32 -c -o buts.o                buts.f
	f95 -m32 -c -o domain.o                domain.f
	f95 -m32 -c -o erhs.o                erhs.f
	f95 -m32 -c -o error.o                error.f
	f95 -m32 -c -o exact.o                exact.f
	f95 -m32 -c -o jacld.o                jacld.f
	f95 -m32 -c -o jacu.o                jacu.f
	f95 -m32 -c -o l2norm.o                l2norm.f
	f95 -m32 -c -o applu.o                applu.f
	f95 -m32 -c -o pintgr.o                pintgr.f
	f95 -m32 -c -o print_results.o                print_results.f
	f95 -m32 -c -o read_input.o                read_input.f
	f95 -m32 -c -o rhs.o                rhs.f
	f95 -m32 -c -o setbv.o                setbv.f
	f95 -m32 -c -o setcoeff.o                setcoeff.f
	f95 -m32 -c -o setiv.o                setiv.f
	f95 -m32 -c -o ssor.o                ssor.f
	f95 -m32 -c -o verify.o                verify.f
	f95 -m32   blts.o buts.o domain.o erhs.o error.o exact.o jacld.o jacu.o l2norm.o applu.o pintgr.o print_results.o read_input.o rhs.o setbv.o setcoeff.o setiv.o ssor.o verify.o     -o applu
run :
	./applu < applu.in

