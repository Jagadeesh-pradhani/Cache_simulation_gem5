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
   `--debug-file=output.out.gz` : Store the trace file in a zip. this file may be sometimes very large.














Reference: <br>
https://www.gem5.org/documentation/learning_gem5/part2/debugging/
