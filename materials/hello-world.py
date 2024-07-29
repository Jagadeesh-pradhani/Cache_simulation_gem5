from gem5.components.boards.simple_board import SimpleBoard
from gem5.components.cachehierarchies.classic.no_cache import NoCache
from gem5.components.cachehierarchies.classic.private_l1_private_l2_cache_hierarchy import PrivateL1PrivateL2CacheHierarchy
from gem5.components.cachehierarchies.ruby.mesi_two_level_cache_hierarchy import (
    MESITwoLevelCacheHierarchy,
)
from gem5.components.memory import SingleChannelDDR3_1600
from gem5.components.processors.simple_processor import SimpleProcessor
from gem5.components.processors.cpu_types import CPUTypes
from gem5.resources.resource import obtain_resource
from gem5.simulate.simulator import Simulator
from gem5.isas import ISA

from unique_cache_hierarchy.unique_cache_hierarchy import UniqueCacheHierarchy
from private_l1_shared_l2_cache.PrivateL1SharedL2 import PrivateL1SharedL2CacheHierarchy

# Obtain the components
cache_hierarchy = NoCache()
# cache_hierarchy = PrivateL1PrivateL2CacheHierarchy(
#     l1d_size = "32KiB",
#     l1i_size = "32KiB",
#     l2_size = "64KiB",
# )

# cache_hierarchy = MESITwoLevelCacheHierarchy(
#     l1d_size = "32KiB",
#     l1i_size = "32KiB",
#     l1d_assoc = 8,
#     l1i_assoc = 8,
#     l2_size = "256KiB",
#     l2_assoc = 16,
#     num_l2_banks = 1,
# )

# cache_hierarchy = UniqueCacheHierarchy()

cache_hierarchy = PrivateL1PrivateL2CacheHierarchy(
    l1d_size="32KiB", l1i_size="32KiB", l2_size="64KiB"
)

memory = SingleChannelDDR3_1600("1GiB")
processor = SimpleProcessor(cpu_type=CPUTypes.TIMING, num_cores=1, isa=ISA.X86)

board = SimpleBoard(
    clk_freq="3GHz",
    processor=processor,
    memory=memory,
    cache_hierarchy=cache_hierarchy,
)

binary = obtain_resource("x86-hello64-static")
board.set_se_binary_workload(binary)

# Setup the simulator and run the simulation.
simulator = Simulator(board=board)
simulator.run()
