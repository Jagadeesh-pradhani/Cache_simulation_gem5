import m5
from m5.objects import *

# System configuration
system = System()
system.clk_domain = SrcClockDomain()
system.clk_domain.clock = '1GHz'
system.clk_domain.voltage_domain = VoltageDomain()

# Memory configuration
system.mem_mode = 'timing'
system.mem_ranges = [AddrRange('512MB')]

# CPU configuration
system.cpu = TimingSimpleCPU()

# Cache configuration
system.cpu.icache = L1ICache(size='32kB')
system.cpu.dcache = L1DCache(size='32kB')
system.cpu.icache.connectCPU(system.cpu)
system.cpu.dcache.connectCPU(system.cpu)

# Prefetcher configuration
system.cpu.dcache.prefetcher = StridePrefetcher()

# L2 Cache configuration
system.l2cache = L2Cache(size='256kB')
system.l2bus = L2XBar()
system.cpu.dcache.connectBus(system.l2bus)
system.l2cache.connectCPUSideBus(system.l2bus)
system.l2cache.connectMemSideBus(system.membus)

# Memory bus
system.membus = SystemXBar()

# Connect CPU and Memory
system.cpu.icache_port = system.membus.slave
system.cpu.dcache_port = system.membus.slave
system.system_port = system.membus.slave

# Create root and instantiate system
root = Root(full_system=False, system=system)
m5.instantiate()

print("Beginning simulation!")
exit_event = m5.simulate()
print("Exiting @ tick {} because {}"
      .format(m5.curTick(), exit_event.getCause()))
