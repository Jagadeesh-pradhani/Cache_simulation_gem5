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

## Simulation Parameters

- **Script**: `se.py`
- **Cores**: `16`
- **CPU**: `O3CPU`
- **Caches**: `PrivateL1SharedL2`
- **Cache Type**: `ruby`
- **L1I Cache**: `64kB`
- **L1D Cache**: `64kB`
- **L2 Cache**: `512kB`
- **Memory**: `512MB`
- **L1 Associativity**: `2-way`
- **L2 Associativity**: `8-way`
- **Network**: `garnet`
- **Topology**: `Mesh_XY`
- **Gtid**: `4x4`
- **Link Latency**: `1`
- **Router Latency**: `1`
- **Debug**: `RubyCache`
- **Benchmark**: `SPEC-CPU2006`

## Workloads Across Cores

| CPU Core | Workload |
|----------|----------|
| CPU0     | dct      |
| CPU1     | fft      |
| CPU2     | heat     |
| CPU3     | dct      |
| CPU4     | dct      |
| CPU5     | fft      |
| CPU6     | heat     |
| CPU7     | dct      |
| CPU8     | dct      |
| CPU9     | fft      |
| CPU10    | heat     |
| CPU11    | dct      |
| CPU12    | dct      |
| CPU13    | fft      |
| CPU14    | heat     |
| CPU15    | dct      |


   
   














Reference: <br>
https://www.gem5.org/documentation/learning_gem5/part2/debugging/
