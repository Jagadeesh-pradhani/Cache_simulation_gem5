# Include standard template for this suite
-include $(SPEC)/benchspec/Makefile.defaults

CC=gcc 

all: 
	make linux
nomakefiles:
	@echo "*** Whoops!  I didn't find a Makefile.defaults file."
	@echo "*** Your distribution is probably corrupt."
	@echo
linux:
	make clean
	make alt-compile CC=gcc
ssmips:
	make clean
	make alt-compile CC=sslittle-na-sstrix-gcc 
sparc:
	make clean
	make alt-compile CC=sparc-elf-gcc
cradle:
	make clean
	make alt-compile CC=umsgcc
	
alt-compile:
	echo "CC =$(CC)"
	$(CC) -c -I./ *.c
	$(CC) -o mcf  *.o
clean:
	rm -rf *.o
run :
	./mcf inp.in	
	
