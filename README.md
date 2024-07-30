# Cache_prefetcher_simulation

## ruby_test.py

1) Navigate to the Build Directory:
   ```
   cd ~/gem5
   ```
2) Write the Configuration Script
   ```
   gedit ~/gem5/configs/example/ruby_test.py
   ```

   Paste the ruby_test.py code from this repo, and save the file.

3) Run the simulation
   ```
   ~/gem5/build/X86/gem5.opt ~/gem5/configs/example/ruby_test.py 
   ```
4) Output
   ![image](https://github.com/user-attachments/assets/f27eb63c-602a-4f30-9a9e-2c702d23a0f3)

5) Analyze the Results
   GEM5 will generate various output files and statistics that you can analyze to understand the performance of your system.
   Open file named 'm5out' for results.



## Extra

https://gem5bootcamp.github.io/gem5-bootcamp-env/modules/using%20gem5/stdlib/ <br>

https://pages.cs.wisc.edu/~swilson/gem5-docs/Prefetcher_8cc_source.html  <br>

https://github.com/dhschall/gem5-fdp/tree/fdp
