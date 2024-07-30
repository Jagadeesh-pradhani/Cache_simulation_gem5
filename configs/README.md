# RUBY

## 1) ruby_test.py : normal run

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


## 2) ruby_direct_test.py :Mesh_XY topology

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
   <br>
   Prefetcher = Rubyprefetcher (default-strideprefetcher)

   ```
   ~/gem5/build/X86/gem5.opt ~/gem5/configs/example/ruby_direct_test.py --num-cpus=4 --num-dirs=4 --network=garnet --topology=Mesh_XY --mesh-rows=2 --link-latency=1 --router-latency=1
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
   Flits ststs used for Calculations. <br>

### Interconnected Parameters of Flits
1. *Flit Latency*: Look for the average, minimum, and maximum latencies of flits.
2. *Flit Throughput*: Metrics indicating the number of flits sent/received over time.
3. *Network Traffic*: Statistics related to the volume and distribution of flit traffic across the network.


### Performance Calculation
1. *Execution Time*: simSeconds and simTicks can provide the simulated execution time.
2. *Instructions Per Cycle (IPC)*: This can be derived if you have the total number of instructions executed and the number of cycles.
3. *Bandwidth Utilization*: Metrics like avgRdBWSys and avgWrBWSys give insights into read and write bandwidth utilization.


### Sample Data from the Stats File
Here are some of the relevant entries extracted from your stats-2.txt:

- system.mem_ctrls0.readReqs: Number of read requests accepted (Count)
- system.mem_ctrls0.writeReqs: Number of write requests accepted (Count)
- system.mem_ctrls0.avgRdBWSys: Average system read bandwidth in Byte/s
- system.mem_ctrls0.avgWrBWSys: Average system write bandwidth in Byte/s
- simSeconds: Number of seconds simulated (Second)
- simTicks: Number of ticks simulated (Tick)
- finalTick: Number of ticks from the beginning of the simulation (Tick)
- system.clk_domain.clock: Clock period in ticks (Tick)
- system.cpu.power_state.pwrStateResidencyTicks::UNDEFINED: Cumulative time in ticks in various power states (Tick)
   
