

![image](https://github.com/user-attachments/assets/fa66786e-3648-41aa-be6a-978273e8648f)


<p align="center">
    <img width="100px" height="20px" src="https://img.shields.io/badge/Ubuntu-22.04-orange?logo=Ubuntu&Ubuntu-22.04"
        alt="ubuntu" />
    <img width="100px" height="20px" src="https://img.shields.io/badge/gem5-23.1.0.0-blue?logo=ROS&ROS=noetic" alt="ROS" />
</p>

# Cache_prefetcher_simulation

This project involves modifying the gem5 simulator to implement a custom version of the MESI_Two_Level cache coherence protocol. The goal of this modification is to enhance cache efficiency in multicore processors by incorporating an energy-efficient prefetching mechanism. Specifically, the protocol is designed to search neighboring cores' L1 caches for data upon an L1 cache miss before accessing the L2 cache, thereby reducing average memory access time (AMAT) and improving overall system performance. The modified protocol has been tested on various SPEC CPU2006 workloads using a 16-core multicore system.

# Quick installation guide

1. Clone this repo
   ```bash
   git clone https://github.com/Jagadeesh-pradhani/Cache_simulation_gem5.git -b dev
   ```
2. Install dependencies
   ```bash
   cd Cache_simulation_gem5/
   ```
   ```bash
   sudo apt install build-essential git m4 scons zlib1g zlib1g-dev \
    libprotobuf-dev protobuf-compiler libprotoc-dev libgoogle-perftools-dev \
    python3-dev libboost-all-dev pkg-config python3-tk
   ```
   ```bash
   pip3 install requirements.txt
   ```
4. Build for **X86_MESI_Two_Level**
   ```bash
   scons build/X86_MESI_Two_Level/gem5.opt -j$(nproc)
   ```

# Simulation

1. Run the **run_simulation.py** file
   ```bash
   cd gem5/
   ```
   ```bash
   python3 run_simulation.py
   ```



Open each folders for details.

## References

https://gem5bootcamp.github.io/gem5-bootcamp-env/modules/using%20gem5/stdlib/ <br>

https://www.gem5.org/documentation/general_docs/ruby/interconnection-network/ <br>

https://pages.cs.wisc.edu/~swilson/gem5-docs/Prefetcher_8cc_source.html  <br>

https://github.com/dhschall/gem5-fdp/tree/fdp

