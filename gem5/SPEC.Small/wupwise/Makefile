###############################################################################
# SPECtools2000
#   Version 0.3
#   This makefile should NOT be changed other than to customize for a 
#   particular benchmark.  Users may overide variables by placing them
#   in the Makefile.spec file.  Or by using the "runspec" program to 
#   build the Makefile.spec file for them.
###############################################################################
#	BENCHMARK	-- standard definitions for this benchmark

NUMBER    = 310
NAME      = wupwise_m
EXEBASE   = wupwise
BENCHLANG = F
SOURCES   = dcabs1.f dlaran.f dlarnd.f dznrm2.f gammul.f init.f lsame.f \
		matmul.f muldeo.f muldoe.f phinit.f rndcnf.f rndphi.f \
		su3mul.f uinith.f wupwise.f xerbla.f \
		zaxpy.f zcopy.f zdotc.f zgemm.f zscal.f

all: 
	f95 -c *.f
	f95 -o wupwise *.o 
run : 
	./wupwise < wupwise.in

# Include standard template for this suite
-include ../../../../Makefile.defaults
-include ../../../Makefile.defaults
-include ../../Makefile.defaults
-include ../Makefile.defaults


nomakefiles:
	@echo "*** Whoops!  I didn't find a Makefile.defaults file."
	@echo "*** Your distribution is probably corrupt."
	@echo
