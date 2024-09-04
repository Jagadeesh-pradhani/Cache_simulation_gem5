import m5
from m5.defines import buildEnv
from m5.objects import *
from m5.util import (
    addToPath,
    fatal,
    panic,
)

addToPath("../")
from common.FileSystemConfig import config_filesystem


from cache.caches import MyCacheSystem

system = System()

# Set the clock frequency of the system (and all of its children)
system.clk_domain = SrcClockDomain()
system.clk_domain.clock = "1GHz"
system.clk_domain.voltage_domain = VoltageDomain()

# Set up the system
system.mem_mode = "timing"  # Use timing accesses
system.mem_ranges = [AddrRange("512MB")]  # Create an address range

# Create a pair of simple CPUs
system.cpu = [X86TimingSimpleCPU() for i in range(4)]

# Create a DDR3 memory controller and connect it to the membus
system.mem_ctrl = MemCtrl()
system.mem_ctrl.dram = DDR3_1600_8x8()
system.mem_ctrl.dram.range = system.mem_ranges[0]

for cpu in system.cpu:
    cpu.createInterruptController()

system.caches = MyCacheSystem()
system.caches.setup(system, system.cpu, [system.mem_ctrl])


print("hello")

