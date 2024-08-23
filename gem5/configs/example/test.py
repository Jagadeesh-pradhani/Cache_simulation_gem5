from m5.objects import *
from m5.util import addToPath
addToPath("../../")
from src.python.gem5.components.cachehierarchies.ruby.mesi_two_level_cache_hierarchy import MESITwoLevelCacheHierarchy
from src.python.gem5.components.memory import SingleChannelDDR3_1600
from src.python.gem5.components.processors.simple_processor import SimpleProcessor

addToPath("../")
from topologies import *

# Ensure you import the necessary components for your simulation
# For example, import SimpleBoard, MESITwoLevelCacheHierarchy, etc.

def create_system(options):
    # Create the system
    system = System()

    # Create the cache hierarchy
    cache_hierarchy = MESITwoLevelCacheHierarchy(
        l1d_size="32KiB",
        l1i_size="32KiB",
        l1d_assoc=8,
        l1i_assoc=8,
        l2_size="256KiB",
        l2_assoc=16,
        num_l2_banks=1,
    )

    # Create the memory and processor
    memory = SingleChannelDDR3_1600("1GiB")
    processor = SimpleProcessor(cpu_type=CPUTypes.TIMING, num_cores=4, isa=ISA.X86)

    # Create the board
    board = SimpleBoard(
        clk_freq="3GHz",
        processor=processor,
        memory=memory,
        cache_hierarchy=cache_hierarchy,
    )

    # Set up the network and topology
    network = Network()
    topology = Mesh_XY([system.cpu[i] for i in range(options.num_cpus)] +
                       [system.l2_cntrl[i] for i in range(options.num_cpus)])

    # Set network and topology parameters
    network.num_routers = options.num_cpus
    network.mesh_rows = options.mesh_rows
    network.link_latency = options.link_latency
    network.router_latency = options.router_latency

    # Create the topology
    topology.makeTopology(options, network, IntLink, ExtLink, Router)

    # Register topology with the system
    topology.registerTopology(options)

    # Set up the binary
    binary = obtain_resource("x86-hello64-static")
    board.set_se_binary_workload(binary)

    # Return the configured system
    return system

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--num-cpus", type=int, default=4, help="Number of CPUs")
parser.add_argument("--num-dirs", type=int, default=4, help="Number of directories")
parser.add_argument("--network", type=str, default="garnet", help="Network type")
parser.add_argument("--topology", type=str, default="Mesh_XY", help="Topology type")
parser.add_argument("--mesh-rows", type=int, default=2, help="Number of rows in the mesh")
parser.add_argument("--link-latency", type=int, default=1, help="Link latency")
parser.add_argument("--router-latency", type=int, default=1, help="Router latency")
options = parser.parse_args()

system = create_system(options)

# Set up the simulator and run the simulation
simulator = Simulator(board=system)
simulator.run()
