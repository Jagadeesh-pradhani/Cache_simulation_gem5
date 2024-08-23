
export OMP_NUM_THREADS=36

#echo "Running matmul"
#../src/m2s -cores 36 ./matmul/matmul 
#mv www.txt  matmul.cctrace


echo "Running Grid"
../../src/m2s -cores 36 ./applu < applu.in
mv www.txt  mgrid.cctrace

