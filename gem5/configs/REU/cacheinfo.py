import argparse
import os
import sys


import m5
from m5.objects import *
from m5.util import addToPath, fatal

# Set the path to where your Ruby system is defined (adjust as necessary)
addToPath('../')
# addToPath('../ruby/')

# Import necessary Ruby components
from common import Options
from common.Caches import *
from common.cpu2000 import *
from ruby import Ruby

# Parse command line options

parser = argparse.ArgumentParser()
Options.addCommonOptions(parser)
Options.addSEOptions(parser)
Ruby.define_options(parser)

args = parser.parse_args()


# Create the system
system = System()

# Set the clock domain
system.clk_domain = SrcClockDomain()
system.clk_domain.clock = '1GHz'
system.clk_domain.voltage_domain = VoltageDomain()

# Set up the memory
system.mem_mode = 'timing'
system.mem_ranges = [AddrRange('512MB')]

# Create a simple CPU and connect it to the system
system.cpu = TimingSimpleCPU()

# Set up the Ruby system
Ruby.create_system(args, False, system)

# Assign L1 and L2 caches to the Ruby system (example setup)
system.ruby_system = system.ruby

# Create the L1 cache
system.cpu.icache = L1ICache(size='32kB', assoc=2)
system.cpu.dcache = L1DCache(size='32kB', assoc=2)

# Create the L2 cache
system.l2cache = L2Cache(size='256kB', assoc=8)

# Connect the CPU and the L1 caches
system.cpu.icache_port = system.cpu.icache.cpu_side
system.cpu.dcache_port = system.cpu.dcache.cpu_side

# Connect L1 to L2
system.cpu.icache.mem_side = system.l2cache.cpu_side
system.cpu.dcache.mem_side = system.l2cache.cpu_side

# Create a memory controller
system.mem_ctrl = MemCtrl()
system.mem_ctrl.dram = DDR3_1600_8x8()
system.mem_ctrl.dram.range = system.mem_ranges[0]

# Connect the L2 cache to the memory controller
system.l2cache.mem_side = system.mem_ctrl.port

# Create and connect a system port (for debugging)
system.system_port = system.cpu.icache_port

# Create the workload (simple "Hello World" application)
binary = 'tests/test-progs/hello/bin/x86/linux/hello'
system.cpu.workload = Process(executable=binary)
system.cpu.createThreads()

# Instantiate the simulation (this actually starts the simulation setup)
m5.instantiate()

# Print initial cache occupancy
print("Initial Cache Occupancy (Occupied Blocks):")
print(f"L1 ICache Occupied Size: {system.ruby_system.cacheMemory.getOccupiedSize()} blocks")
print(f"L1 DCache Occupied Size: {system.ruby_system.cacheMemory.getOccupiedSize()} blocks")
print(f"L2 Cache Occupied Size: {system.ruby_system.cacheMemory.getOccupiedSize()} blocks")

# Run the simulation
print("Running simulation...")
exit_event = m5.simulate()

# Print cache occupancy after simulation
print("Final Cache Occupancy (Occupied Blocks):")
print(f"L1 ICache Occupied Size: {system.ruby_system.cacheMemory.getOccupiedSize()} blocks")
print(f"L1 DCache Occupied Size: {system.ruby_system.cacheMemory.getOccupiedSize()} blocks")
print(f"L2 Cache Occupied Size: {system.ruby_system.cacheMemory.getOccupiedSize()} blocks")

# Check the simulation status
print('Exiting @ tick {} because {}'.format(m5.curTick(), exit_event.getCause()))
