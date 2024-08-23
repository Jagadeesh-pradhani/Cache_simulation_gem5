import m5
from m5.objects import *
from caches import *
from m5.objects import GarnetNetwork

system = System()

#n-> number of tiles or cores
n = 4

system.clk_domain = SrcClockDomain()
system.clk_domain.clock = '1GHz'
system.clk_domain.voltage_domain = VoltageDomain()

system.mem_mode = 'timing'
system.mem_ranges = [AddrRange('512MB')]

# Create the CPU cores
system.cpu = [X86TimingSimpleCPU() for i in range(n)]

for i in range(n):
    system.cpu[i].icache = L1ICache(prefetcher = BOPPrefetcher())
    system.cpu[i].dcache = L1DCache(prefetcher = BOPPrefetcher())
    system.cpu[i].icache.connectCPU(system.cpu[i])
    system.cpu[i].dcache.connectCPU(system.cpu[i])

system.mainNOC = SystemXBar()
system.noc = [L2XBar() for _ in range(n)]

for i in range(n):
    system.cpu[i].icache.connectBus(system.noc[i])
    system.cpu[i].dcache.connectBus(system.noc[i])
    system.mainNOC.cpu_side_ports = system.noc[i].mem_side_ports

# Create an L2Cache with StridePrefetcher
system.l2cache = L2Cache(prefetcher = BOPPrefetcher())
system.l2cache.connectCPUSideBus(system.mainNOC)

#create a root mem bus
system.membus = SystemXBar()
system.l2cache.connectMemSideBus(system.membus)

for i in range(n):
    system.cpu[i].createInterruptController()
    system.cpu[i].interrupts[0].pio = system.membus.mem_side_ports
    system.cpu[i].interrupts[0].int_requestor = system.membus.cpu_side_ports
    system.cpu[i].interrupts[0].int_responder = system.membus.mem_side_ports

system.system_port = system.membus.cpu_side_ports

system.mem_ctrl = MemCtrl()
system.mem_ctrl.dram = DDR3_1600_8x8()
system.mem_ctrl.dram.range = system.mem_ranges[0]
system.mem_ctrl.port = system.membus.mem_side_ports

# binary = 'tests/test-progs/hello/bin/x86/linux/hello'
binary = [
    'tests/test-progs/hello/bin/x86/linux/hello',
    'tests/test-progs/hello/bin/x86/linux/hello',
    'tests/test-progs/hello/bin/x86/linux/hello',
    'tests/test-progs/hello/bin/x86/linux/hello'
    ]

# For gem5 V21 and beyond
system.workload = SEWorkload.init_compatible(binary[0])

# Assign a separate process to each CPU
for i in range(n):
    process = Process(pid = 100+i)
    process.cmd = [binary[i]]
    system.cpu[i].workload = process
    system.cpu[i].createThreads()

root = Root(full_system=False, system=system)
m5.instantiate()

print("Beginning simulation!")
exit_event = m5.simulate()

print('Exiting @ tick {} because {}'.format(m5.curTick(), exit_event.getCause()))