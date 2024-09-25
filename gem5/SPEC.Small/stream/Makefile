
SYSROOT=$(ARMTOOLCHAINDIR)/gcc-linaro-5.4.1-2017.05-x86_64_aarch64-linux-gnu/aarch64-linux-gnu/libc
CC=$(ARMTOOLCHAINDIR)/gcc-linaro-5.4.1-2017.05-x86_64_aarch64-linux-gnu/bin/aarch64-linux-gnu-gcc
CFLAGS = -O2

FF = f77
FFLAGS = -O2

all: stream_c.exe

stream_f.exe: stream.f mysecond.o
	$(CC) $(CFLAGS) -c mysecond.c
	$(FF) $(FFLAGS) -c stream.f
	$(FF) $(FFLAGS) stream.o mysecond.o -o stream_f.exe

stream_c.exe: stream.c
	$(CC) --sysroot=$(SYSROOT) --no-sysroot-suffix --static stream.c -o stream_c.exe

clean:
	rm -f stream_f.exe stream_c.exe *.o

# an example of a more complex build line for the Intel icc compiler
stream.icc: stream.c
	icc -O3 -xCORE-AVX2 -ffreestanding -qopenmp -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 stream.c -o stream.omp.AVX2.80M.20x.icc
