
export OMP_NUM_THREADS=30

#echo "Running matmul"
#../src/m2s -cores 36 ./matmul/matmul 
#mv www.txt  matmul.cctrace


echo "Running fma3d"
../../src/m2s -cores 30 ./fma3d < fma3d.in
mv www.txt fma3d.cctrace

