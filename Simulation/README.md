## Debugging in GEM5
1) All of the available debug flags, by running gem5 with the `--debug-help` parameter.
 ```
 build/X86/gem5.opt --debug-help
 ```

2) Using a debug parameter
   Add the debug flags with simulaiton commands.
   ```
   gem5/build/X86/gem5.opt --debug-flags=RubyCache --debug-file=output.out.gz gem5/configs/deprecated/example/se.py
   ```
   `--debug-flags=RubyCache` : This will trace the address allocated to all the ruby components
   `--debug-file=output.out.gz` : Store the trace file in `m5out` folder. this file may be sometimes very large.

## Simulation parameters
   Script             : `se.py`             <br>
   Cores              : `16`                <br>
   CPU                : `O3CPU`             <br>
   Caches             : `PrivateL1SharedL2` <br>
   Cache_type         : `ruby`              <br>
   L1I_Cache          : `64kB`              <br>
   L1D_Cache          : `64kB`              <br>
   L2_cache           : `512kB`             <br>
   Memory             : `512MB`             <br>
   L1_assoc           : `2-way`             <br>
   L2_assoc           : `8-way`             <br>
   Network            : `garnet`            <br>
   Topology           : `Mesh_XY`           <br>
   Gtid               : `4x4`               <br>
   Link_latency       : `1`                 <br>
   Router_latency     : `1`                 <br>
   Debug              : `RubyCache`         <br>
   Benchmark          : `SPEC-CPU2006`      <br>
    Workloads accross Cores:
    CPU0 : `dct`         CPU1 : `fft`         CPU2 : `heat`         CPU3 : `dct`  <br>
    CPU4 : `dct`         CPU5 : `fft`         CPU6 : `heat`         CPU7 : `dct`  <br>
    CPU8 : `dct`         CPU9 : `fft`         CPU10 : `heat`        CPU11 : `dct`  <br>
    CPU12 : `dct`        CPU13 : `fft`        CPU14 : `heat`        CPU15 : `dct`  <br>
   
   














Reference: <br>
https://www.gem5.org/documentation/learning_gem5/part2/debugging/
