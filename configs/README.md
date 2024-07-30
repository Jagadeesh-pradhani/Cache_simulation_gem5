# RUBY

## 1) ruby_test.py

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


## 2) ruby_direct_test.py

1) Replace the Mesh_XY.py file from this repo to the file in the location "gem5/configs/topologies/Mesh_XY.py".
2) Navigate to the Build Directory:
   ```
   cd ~/gem5
   ```
3) Run the ruby script with Garnet network, <br>
   cpu's = 4 <br>
   dir's = 4 <br>
   rows = 2 <br>
   latency = 1 <br>

   ```
   ~/gem5/build/X86_MESI_Two_Level/gem5.opt ~/gem5/configs/example/ruby_direct_test.py --num-cpus=4 --num-dirs=4 --network=garnet --topology=Mesh_XY --mesh-rows=2 --link-latency=1 --router-latency=1
   ```
4) Output
   ![image](https://github.com/user-attachments/assets/4af03ef4-4ec9-4866-88f5-ac24018f99da)

5) Result
   Results can be viewd in m5out directory <br>
   config.pdf <br>
   ![image](https://github.com/user-attachments/assets/e6f9c5f6-8c93-4870-8bb8-e4d6d1a7e899)
   <br>
   ruby.pdf<br>
   ![image](https://github.com/user-attachments/assets/5aa3ae8b-8880-4ad7-b355-65be9b9e5bcc)
   <br>
   <br>
   Stats.txt
   ```
   system.ruby.network.average_flit_latency
   system.ruby.network.average_flit_network_latency
   system.ruby.network.flits_received::total
   system.ruby.network.flits_received 
   ```
   Flits ststs used for Calculations.



   
