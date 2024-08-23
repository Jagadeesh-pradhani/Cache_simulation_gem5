from m5.objects import *
from m5.util import addToPath
import m5

# Add the common scripts to our path
addToPath('../')
from ruby import Ruby
addToPath('../common')
from common import Options

# Create a system
system = System()
system.clk_domain = SrcClockDomain()
system.clk_domain.clock = '1GHz'
system.clk_domain.voltage_domain = VoltageDomain()

# Create Ruby system
system.ruby = RubySystem()
system.mem_ranges = [AddrRange('512MB')]

# Define number of CPUs
num_cpus = 2

# Create CPU objects
system.cpu = [TimingSimpleCPU(cpu_id=i) for i in range(num_cpus)]

# Create L1 instruction and data caches
for i in range(num_cpus):
    system.cpu[i].icache = RubyCache(size='32kB', assoc=4)
    system.cpu[i].dcache = RubyCache(size='32kB', assoc=4)

# Create L2 cache
l2_cache = RubyCache(size='256kB', assoc=8)

# Create Ruby network
system.ruby.network = GarnetNetwork(ruby_system=system.ruby)
system.ruby.network.number_of_virtual_networks = 1

# Create sequencers and connect them to the caches
for i in range(num_cpus):
    # Create L1 cache controllers
    system.cpu[i].icache_controller = L1Cache_Controller(
        version=i, L1Icache=system.cpu[i].icache, L1Dcache=system.cpu[i].dcache,
        send_evictions=False, ruby_system=system.ruby)
    system.cpu[i].dcache_controller = L1Cache_Controller(
        version=i, L1Icache=system.cpu[i].icache, L1Dcache=system.cpu[i].dcache,
        send_evictions=False, ruby_system=system.ruby)
    
    # Create sequencer
    system.cpu[i].sequencer = RubySequencer(
        version=i, ruby_system=system.ruby)
    
    # Connect controllers to network
    system.cpu[i].icache_controller.mandatoryQueue = system.ruby.network.in_port
    system.cpu[i].dcache_controller.mandatoryQueue = system.ruby.network.in_port
    system.cpu[i].sequencer.in_port = system.ruby.network.in_port

    # Connect sequencer to cache controllers
    system.cpu[i].icache_controller.sequencer = system.cpu[i].sequencer
    system.cpu[i].dcache_controller.sequencer = system.cpu[i].sequencer

    # Append controllers to Ruby system
    system.ruby._cpu_ports.append(system.cpu[i].icache_controller)
    system.ruby._cpu_ports.append(system.cpu[i].dcache_controller)

# Create L2 cache controller
l2_cntrl = L2Cache_Controller(
    version=0, L2cache=l2_cache, send_evictions=False, ruby_system=system.ruby)
l2_cntrl.mandatoryQueue = system.ruby.network.in_port

# Append L2 controller to Ruby system
system.ruby._cpu_ports.append(l2_cntrl)

# Create memory controller
system.mem_ctrl = MemCtrl()
system.mem_ctrl.dram = DDR3_1600_8x8()
system.mem_ctrl.port = system.ruby._dma_ports[0]

# Instantiate the system
root = Root(full_system=False, system=system)
m5.instantiate()

# Run simulation
print("Beginning simulation!")
exit_event = m5.simulate()
print('Exiting @ tick %i because %s' % (m5.curTick(), exit_event.getCause()))
