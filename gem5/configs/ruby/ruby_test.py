import math

import m5
from m5.defines import buildEnv
from m5.objects import *
from m5.util import (
    addToPath,
    fatal,
)

from gem5.isas import ISA
from gem5.runtime import get_supported_isas

addToPath("../")

from common import (
    FileSystemConfig,
    MemConfig,
    ObjectList,
)
from network import Network
from topologies import *

class MyCache(RubyCache):
    def __init__(self, size, assoc):
        super().__init__()
        self.size = size
        self.assoc = assoc

# System configuration
system = System()
system.clk_domain = SrcClockDomain()
system.clk_domain.clock = '1GHz'
system.clk_domain.voltage_domain = VoltageDomain()

# Ruby configuration
system.ruby = RubySystem()
ruby = system.ruby

# L1 Cache
l1i_cache = MyCache(size='32kB', assoc=4)
l1d_cache = MyCache(size='32kB', assoc=4)

# L2 Cache
l2_cache = MyCache(size='256kB', assoc=8)

# Network configuration
network = Network()
ruby.network = network

# Create the cores
system.cpu = [TimingSimpleCPU(cpu_id=i) for i in range(2)]

# Connect the caches and the network
for i in range(2):
    l1_cntrl = L1Cache_Controller(version=i,
                                  L1Icache=l1i_cache,
                                  L1Dcache=l1d_cache,
                                  send_evictions=False,
                                  cache_response_latency=10,
                                  cache_request_latency=5)
    ruby.l1_cntrl = l1_cntrl

    l2_cntrl = L2Cache_Controller(version=i,
                                  L2cache=l2_cache,
                                  send_evictions=False,
                                  transitions_per_cycle=4,
                                  cache_request_latency=20,
                                  cache_response_latency=10)
    ruby.l2_cntrl = l2_cntrl

    cpu_sequencer = RubySequencer(version=i,
                                  icache=l1i_cache,
                                  dcache=l1d_cache)
    ruby.sequencer = cpu_sequencer

    # Connect the L1 and L2 controllers to the network
    l1_cntrl.mandatoryQueue = network.getNetworkInterface()
    l2_cntrl.mandatoryQueue = network.getNetworkInterface()

# Memory controller configuration
mem_ctrl = MemCtrl()
mem_ctrl.port = ruby._dma_ports[0]

# Simulation run
root = Root(full_system=False, system=system)
m5.instantiate()
print("Beginning simulation!")
exit_event = m5.simulate()
print('Exiting @ tick %i because %s' % (m5.curTick(), exit_event.getCause()))
