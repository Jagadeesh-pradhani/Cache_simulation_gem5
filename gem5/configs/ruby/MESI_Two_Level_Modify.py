import math

import m5
from m5.defines import buildEnv
from m5.objects import *

from .Ruby import (
    create_directories,
    create_topology,
    send_evicts,
)


#
# Declare caches used by the protocol
#
class L1Cache(RubyCache):
    pass


class L2Cache(RubyCache):
    pass


def define_options(parser):
    return


def getNeighbourCache(curr, size):
    north = None
    south = None
    east = None
    west = None

    row = curr // size
    col = curr % size

    if row > 0:
        north = (row - 1) * size + col

    if row < size - 1:
        south = (row + 1) * size + col

    if col > 0:
        west = row * size + (col - 1)

    if col < size - 1:
        east = row * size + (col + 1)
    
    return north,south,east,west



def create_system_old(
    options, full_system, system, dma_ports, bootmem, ruby_system, cpus
):
    if buildEnv["PROTOCOL"] != "MESI_Two_Level":
        fatal("This script requires the MESI_Two_Level protocol to be built.")

    cpu_sequencers = []

    #
    # The ruby network creation expects the list of nodes in the system to be
    # consistent with the NetDest list.  Therefore the l1 controller nodes must be
    # listed before the directory nodes and directory nodes before dma nodes, etc.
    #
    l1_cntrl_nodes = []
    l2_cntrl_nodes = []
    dma_cntrl_nodes = []

    #
    # Must create the individual controllers before the network to ensure the
    # controller constructors are called before the network constructor
    #
    l2_bits = int(math.log(options.num_l2caches, 2))    #0
    block_size_bits = int(math.log(options.cacheline_size, 2))  #6

    for i in range(options.num_cpus):
        #
        # First create the Ruby objects associated with this cpu
        #
        l1i_cache = L1Cache(
            size=options.l1i_size,
            assoc=options.l1i_assoc,
            start_index_bit=block_size_bits,
            is_icache=True,
            ruby_system = ruby_system,
            # replacement_policy= BRRIPRP(),
        )
        l1d_cache = L1Cache(
            size=options.l1d_size,
            assoc=options.l1d_assoc,
            start_index_bit=block_size_bits,
            is_icache=False,
            ruby_system = ruby_system,
            # replacement_policy= BRRIPRP(),
        )

        prefetcher = RubyPrefetcher()

        clk_domain = cpus[i].clk_domain
        print("hello")

        l1_cntrl = L1Cache_Controller(
            version=i,
            L1Icache=l1i_cache,
            L1Dcache=l1d_cache,
            l2_select_num_bits=l2_bits,
            send_evictions=send_evicts(options),
            prefetcher=prefetcher,
            ruby_system=ruby_system,
            clk_domain=clk_domain,
            transitions_per_cycle=options.ports,
            enable_prefetch=True,
        )

        cpu_seq = RubySequencer(
            version=i,
            dcache=l1d_cache,
            clk_domain=clk_domain,
            ruby_system=ruby_system,
        )

        l1_cntrl.sequencer = cpu_seq
        exec("ruby_system.l1_cntrl%d = l1_cntrl" % i)

        # Add controllers and sequencers to the appropriate lists
        cpu_sequencers.append(cpu_seq)
        l1_cntrl_nodes.append(l1_cntrl)

        # Connect the L1 controllers and the network
        l1_cntrl.mandatoryQueue = MessageBuffer()
        l1_cntrl.requestFromL1Cache = MessageBuffer()
        l1_cntrl.requestFromL1Cache.out_port = ruby_system.network.in_port
        l1_cntrl.responseFromL1Cache = MessageBuffer()
        l1_cntrl.responseFromL1Cache.out_port = ruby_system.network.in_port
        l1_cntrl.unblockFromL1Cache = MessageBuffer()
        l1_cntrl.unblockFromL1Cache.out_port = ruby_system.network.in_port

        #New Buffers
        l1_cntrl.n_RequestFromL1Cache = MessageBuffer()
        l1_cntrl.n_RequestFromL1Cache.out_port = ruby_system.network.in_port
        l1_cntrl.n_ResponseFromL1Cache = MessageBuffer()
        l1_cntrl.n_ResponseFromL1Cache.out_port = ruby_system.network.in_port
        

        l1_cntrl.optionalQueue = MessageBuffer()

        l1_cntrl.requestToL1Cache = MessageBuffer()
        l1_cntrl.requestToL1Cache.in_port = ruby_system.network.out_port
        l1_cntrl.responseToL1Cache = MessageBuffer()
        l1_cntrl.responseToL1Cache.in_port = ruby_system.network.out_port

        #New
        l1_cntrl.n_RequestToL1Cache = MessageBuffer()
        l1_cntrl.n_RequestToL1Cache.in_port = ruby_system.network.out_port
        l1_cntrl.n_ResponseToL1Cache = MessageBuffer()
        l1_cntrl.n_ResponseToL1Cache.in_port = ruby_system.network.out_port

    l2_index_start = block_size_bits + l2_bits  #6

    for i in range(options.num_l2caches):
        #
        # First create the Ruby objects associated with this cpu
        #
        l2_cache = L2Cache(
            size=options.l2_size,
            assoc=options.l2_assoc,
            start_index_bit=l2_index_start,
            ruby_system = ruby_system,
            # replacement_policy= BRRIPRP(),
        )

        l2_cntrl = L2Cache_Controller(
            version=i,
            L2cache=l2_cache,
            transitions_per_cycle=options.ports,
            ruby_system=ruby_system,
        )

        exec("ruby_system.l2_cntrl%d = l2_cntrl" % i)
        l2_cntrl_nodes.append(l2_cntrl)

        # Connect the L2 controllers and the network
        l2_cntrl.DirRequestFromL2Cache = MessageBuffer()
        l2_cntrl.DirRequestFromL2Cache.out_port = ruby_system.network.in_port
        l2_cntrl.L1RequestFromL2Cache = MessageBuffer()
        l2_cntrl.L1RequestFromL2Cache.out_port = ruby_system.network.in_port
        l2_cntrl.responseFromL2Cache = MessageBuffer()
        l2_cntrl.responseFromL2Cache.out_port = ruby_system.network.in_port

        l2_cntrl.unblockToL2Cache = MessageBuffer()
        l2_cntrl.unblockToL2Cache.in_port = ruby_system.network.out_port
        l2_cntrl.L1RequestToL2Cache = MessageBuffer()
        l2_cntrl.L1RequestToL2Cache.in_port = ruby_system.network.out_port
        l2_cntrl.responseToL2Cache = MessageBuffer()
        l2_cntrl.responseToL2Cache.in_port = ruby_system.network.out_port

    # Run each of the ruby memory controllers at a ratio of the frequency of
    # the ruby system
    # clk_divider value is a fix to pass regression.
    ruby_system.memctrl_clk_domain = DerivedClockDomain(
        clk_domain=ruby_system.clk_domain, clk_divider=3
    )

    mem_dir_cntrl_nodes, rom_dir_cntrl_node = create_directories(
        options, bootmem, ruby_system, system
    )
    dir_cntrl_nodes = mem_dir_cntrl_nodes[:]
    if rom_dir_cntrl_node is not None:
        dir_cntrl_nodes.append(rom_dir_cntrl_node)
    for dir_cntrl in dir_cntrl_nodes:
        # Connect the directory controllers and the network
        dir_cntrl.requestToDir = MessageBuffer()
        dir_cntrl.requestToDir.in_port = ruby_system.network.out_port
        dir_cntrl.responseToDir = MessageBuffer()
        dir_cntrl.responseToDir.in_port = ruby_system.network.out_port
        dir_cntrl.responseFromDir = MessageBuffer()
        dir_cntrl.responseFromDir.out_port = ruby_system.network.in_port
        dir_cntrl.requestToMemory = MessageBuffer()
        dir_cntrl.responseFromMemory = MessageBuffer()

    for i, dma_port in enumerate(dma_ports):
        # Create the Ruby objects associated with the dma controller
        dma_seq = DMASequencer(
            version=i, ruby_system=ruby_system, in_ports=dma_port
        )

        dma_cntrl = DMA_Controller(
            version=i,
            dma_sequencer=dma_seq,
            transitions_per_cycle=options.ports,
            ruby_system=ruby_system,
        )

        exec("ruby_system.dma_cntrl%d = dma_cntrl" % i)
        dma_cntrl_nodes.append(dma_cntrl)

        # Connect the dma controller to the network
        dma_cntrl.mandatoryQueue = MessageBuffer()
        dma_cntrl.responseFromDir = MessageBuffer(ordered=True)
        dma_cntrl.responseFromDir.in_port = ruby_system.network.out_port
        dma_cntrl.requestToDir = MessageBuffer()
        dma_cntrl.requestToDir.out_port = ruby_system.network.in_port

    all_cntrls = (
        l1_cntrl_nodes + l2_cntrl_nodes + dir_cntrl_nodes + dma_cntrl_nodes
    )

    # Create the io controller and the sequencer
    if full_system:
        io_seq = DMASequencer(version=len(dma_ports), ruby_system=ruby_system)
        ruby_system._io_port = io_seq
        io_controller = DMA_Controller(
            version=len(dma_ports),
            dma_sequencer=io_seq,
            ruby_system=ruby_system,
        )
        ruby_system.io_controller = io_controller

        # Connect the dma controller to the network
        io_controller.mandatoryQueue = MessageBuffer()
        io_controller.responseFromDir = MessageBuffer(ordered=True)
        io_controller.responseFromDir.in_port = ruby_system.network.out_port
        io_controller.requestToDir = MessageBuffer()
        io_controller.requestToDir.out_port = ruby_system.network.in_port

        all_cntrls = all_cntrls + [io_controller]

    ruby_system.network.number_of_virtual_networks = 5
    topology = create_topology(all_cntrls, options)
    return (cpu_sequencers, mem_dir_cntrl_nodes, topology)


def create_system(options, full_system, system, dma_ports, bootmem, ruby_system, cpus):
    if buildEnv["PROTOCOL"] != "MESI_Two_Level":
        fatal("This script requires the MESI_Two_Level protocol to be built.")
    

    cpu_sequencers = []
    l1_cntrl_nodes = []
    l2_cntrl_nodes = []
    dma_cntrl_nodes = []

    l2_bits = int(math.log(options.num_l2caches, 2))
    block_size_bits = int(math.log(options.cacheline_size, 2))


    # Store caches for each core (index 0 to 15)
    l1_caches = []

    # First loop: Create all L1 caches (both L1I and L1D) for each core
    for i in range(options.num_cpus):
        l1i_cache = L1Cache(
            size=options.l1i_size,
            assoc=options.l1i_assoc,
            start_index_bit=block_size_bits,
            is_icache=True,
            ruby_system=ruby_system,
        )
        l1d_cache = L1Cache(
            size=options.l1d_size,
            assoc=options.l1d_assoc,
            start_index_bit=block_size_bits,
            is_icache=False,
            ruby_system=ruby_system,
        )

        # Store the caches as a tuple (L1Icache, L1Dcache)
        l1_caches.append((l1i_cache, l1d_cache))

    # Second loop: Assign caches and neighbors to L1Cache_Controller
    for i in range(options.num_cpus):
        

        l1i_cache, l1d_cache = l1_caches[i]  # Get current core's caches
        prefetcher = RubyPrefetcher()
        clk_domain = cpus[i].clk_domain

        north, south, east, west = getNeighbourCache(i, 4)
        if north == None:
            north = 0
        if south == None:
            south = 0
        if east == None:
            east = 0
        if west == None:
            west = 0
        

        # Get the neighbor L1D and L1I caches, or None if no neighbor exists
        north_L1D = l1_caches[int(north)][1] if north is not None else None
        north_L1I = l1_caches[int(north)][0] if north is not None else None
        south_L1D = l1_caches[int(south)][1] if south is not None else None
        south_L1I = l1_caches[int(south)][0] if south is not None else None
        east_L1D = l1_caches[int(east)][1] if east is not None else None
        east_L1I = l1_caches[int(east)][0] if east is not None else None
        west_L1D = l1_caches[int(west)][1] if west is not None else None
        west_L1I = l1_caches[int(west)][0] if west is not None else None

        # Create the L1Cache_Controller with the neighboring caches
        l1_cntrl = L1Cache_Controller(
            version=i,
            L1Icache=l1i_cache,
            L1Dcache=l1d_cache,
            l2_select_num_bits=l2_bits,
            send_evictions=send_evicts(options),
            prefetcher=prefetcher,
            ruby_system=ruby_system,
            clk_domain=clk_domain,
            transitions_per_cycle=options.ports,
            enable_prefetch=True,
            # Neighbors
            north_L1D=north_L1D,
            north_L1I=north_L1I,
            south_L1D=south_L1D,
            south_L1I=south_L1I,
            east_L1D=east_L1D,
            east_L1I=east_L1I,
            west_L1D=west_L1D,
            west_L1I=west_L1I,
        )

        # Create the RubySequencer for the current core
        cpu_seq = RubySequencer(
            version=i,
            dcache=l1d_cache,
            clk_domain=clk_domain,
            ruby_system=ruby_system,
        )

        # Link the sequencer to the controller
        l1_cntrl.sequencer = cpu_seq
        exec(f"ruby_system.l1_cntrl{i} = l1_cntrl")

        # Add controllers and sequencers to the appropriate lists
        cpu_sequencers.append(cpu_seq)
        l1_cntrl_nodes.append(l1_cntrl)

        # Connect the L1 controllers and the network
        l1_cntrl.mandatoryQueue = MessageBuffer()
        l1_cntrl.requestFromL1Cache = MessageBuffer()
        l1_cntrl.requestFromL1Cache.out_port = ruby_system.network.in_port
        l1_cntrl.responseFromL1Cache = MessageBuffer()
        l1_cntrl.responseFromL1Cache.out_port = ruby_system.network.in_port
        l1_cntrl.unblockFromL1Cache = MessageBuffer()
        l1_cntrl.unblockFromL1Cache.out_port = ruby_system.network.in_port

        # New Buffers for additional network requests/responses
        l1_cntrl.n_RequestFromL1Cache = MessageBuffer()
        l1_cntrl.n_RequestFromL1Cache.out_port = ruby_system.network.in_port
        l1_cntrl.n_ResponseFromL1Cache = MessageBuffer()
        l1_cntrl.n_ResponseFromL1Cache.out_port = ruby_system.network.in_port

        l1_cntrl.optionalQueue = MessageBuffer()

        l1_cntrl.requestToL1Cache = MessageBuffer()
        l1_cntrl.requestToL1Cache.in_port = ruby_system.network.out_port
        l1_cntrl.responseToL1Cache = MessageBuffer()
        l1_cntrl.responseToL1Cache.in_port = ruby_system.network.out_port

        # New Buffers for additional requests/responses
        l1_cntrl.n_RequestToL1Cache = MessageBuffer()
        l1_cntrl.n_RequestToL1Cache.in_port = ruby_system.network.out_port
        l1_cntrl.n_ResponseToL1Cache = MessageBuffer()
        l1_cntrl.n_ResponseToL1Cache.in_port = ruby_system.network.out_port

    # Continue with the rest of the system setup (L2 cache, directories, etc.)
    l2_index_start = block_size_bits + l2_bits  #6

    for i in range(options.num_l2caches):
        #
        # First create the Ruby objects associated with this cpu
        #
        l2_cache = L2Cache(
            size=options.l2_size,
            assoc=options.l2_assoc,
            start_index_bit=l2_index_start,
            ruby_system = ruby_system,
            # replacement_policy= BRRIPRP(),
        )

        l2_cntrl = L2Cache_Controller(
            version=i,
            L2cache=l2_cache,
            transitions_per_cycle=options.ports,
            ruby_system=ruby_system,
        )

        exec("ruby_system.l2_cntrl%d = l2_cntrl" % i)
        l2_cntrl_nodes.append(l2_cntrl)

        # Connect the L2 controllers and the network
        l2_cntrl.DirRequestFromL2Cache = MessageBuffer()
        l2_cntrl.DirRequestFromL2Cache.out_port = ruby_system.network.in_port
        l2_cntrl.L1RequestFromL2Cache = MessageBuffer()
        l2_cntrl.L1RequestFromL2Cache.out_port = ruby_system.network.in_port
        l2_cntrl.responseFromL2Cache = MessageBuffer()
        l2_cntrl.responseFromL2Cache.out_port = ruby_system.network.in_port

        l2_cntrl.unblockToL2Cache = MessageBuffer()
        l2_cntrl.unblockToL2Cache.in_port = ruby_system.network.out_port
        l2_cntrl.L1RequestToL2Cache = MessageBuffer()
        l2_cntrl.L1RequestToL2Cache.in_port = ruby_system.network.out_port
        l2_cntrl.responseToL2Cache = MessageBuffer()
        l2_cntrl.responseToL2Cache.in_port = ruby_system.network.out_port

    # Run each of the ruby memory controllers at a ratio of the frequency of
    # the ruby system
    # clk_divider value is a fix to pass regression.
    ruby_system.memctrl_clk_domain = DerivedClockDomain(
        clk_domain=ruby_system.clk_domain, clk_divider=3
    )

    mem_dir_cntrl_nodes, rom_dir_cntrl_node = create_directories(
        options, bootmem, ruby_system, system
    )
    dir_cntrl_nodes = mem_dir_cntrl_nodes[:]
    if rom_dir_cntrl_node is not None:
        dir_cntrl_nodes.append(rom_dir_cntrl_node)
    for dir_cntrl in dir_cntrl_nodes:
        # Connect the directory controllers and the network
        dir_cntrl.requestToDir = MessageBuffer()
        dir_cntrl.requestToDir.in_port = ruby_system.network.out_port
        dir_cntrl.responseToDir = MessageBuffer()
        dir_cntrl.responseToDir.in_port = ruby_system.network.out_port
        dir_cntrl.responseFromDir = MessageBuffer()
        dir_cntrl.responseFromDir.out_port = ruby_system.network.in_port
        dir_cntrl.requestToMemory = MessageBuffer()
        dir_cntrl.responseFromMemory = MessageBuffer()

    for i, dma_port in enumerate(dma_ports):
        # Create the Ruby objects associated with the dma controller
        dma_seq = DMASequencer(
            version=i, ruby_system=ruby_system, in_ports=dma_port
        )

        dma_cntrl = DMA_Controller(
            version=i,
            dma_sequencer=dma_seq,
            transitions_per_cycle=options.ports,
            ruby_system=ruby_system,
        )

        exec("ruby_system.dma_cntrl%d = dma_cntrl" % i)
        dma_cntrl_nodes.append(dma_cntrl)

        # Connect the dma controller to the network
        dma_cntrl.mandatoryQueue = MessageBuffer()
        dma_cntrl.responseFromDir = MessageBuffer(ordered=True)
        dma_cntrl.responseFromDir.in_port = ruby_system.network.out_port
        dma_cntrl.requestToDir = MessageBuffer()
        dma_cntrl.requestToDir.out_port = ruby_system.network.in_port

    all_cntrls = (
        l1_cntrl_nodes + l2_cntrl_nodes + dir_cntrl_nodes + dma_cntrl_nodes
    )

    # Create the io controller and the sequencer
    if full_system:
        io_seq = DMASequencer(version=len(dma_ports), ruby_system=ruby_system)
        ruby_system._io_port = io_seq
        io_controller = DMA_Controller(
            version=len(dma_ports),
            dma_sequencer=io_seq,
            ruby_system=ruby_system,
        )
        ruby_system.io_controller = io_controller

        # Connect the dma controller to the network
        io_controller.mandatoryQueue = MessageBuffer()
        io_controller.responseFromDir = MessageBuffer(ordered=True)
        io_controller.responseFromDir.in_port = ruby_system.network.out_port
        io_controller.requestToDir = MessageBuffer()
        io_controller.requestToDir.out_port = ruby_system.network.in_port

        all_cntrls = all_cntrls + [io_controller]

    ruby_system.network.number_of_virtual_networks = 5
    topology = create_topology(all_cntrls, options)
    return (cpu_sequencers, mem_dir_cntrl_nodes, topology)
