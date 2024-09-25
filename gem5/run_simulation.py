import os
import subprocess


"""
Address bits : 64

L1Dcache: Assoc : 2, Block size : 64, Number of sets: 256, set bits: 8.
L1Icache: Assoc : 2, Block size : 64, Number of sets: 256, set bits: 8.
L2cache: Assoc : 8, Block size : 64, Number of sets: 8192, set bits: 13.

L1 address : Tag (50) : Set Index (8) : Block Offset (6)
L2 address : Tag (45) : Set Index (13) : Block Offset (6)


CacheRubyInfoFlag
"""



# Define the command components
gem5_binary = "build/X86_MESI_Two_Level/gem5.opt"
debug_flags = "--debug-flags=CacheRubyInfoFlag"
debug_file = "--debug-file=output.out.gz"
config_script = "configs/deprecated/example/se.py"

# Define the list of workloads
workloads = [
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",

    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",

    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
    "tests/test-progs/hello/bin/x86/linux/hello",
 
]

# Construct the command with the workloads
cmd = ";".join(workloads)  # Join the workloads with semicolons

# Additional command options
options = "--options=\"\""
num_cpus = "-n 16"
num_dirs = "--num-dirs=16"
ruby = "--ruby"
network = "--network=garnet"
topology = "--topology=Mesh_XY"
mesh_rows = "--mesh-rows=4"
link_latency = "--link-latency=1"
router_latency = "--router-latency=1"
cpu_type = "--cpu-type=O3CPU"
l1d_size = "--l1d_size=32kB"
l1i_size = "--l1i_size=32kB"
l2_size = "--l2_size=4MB"
mem_size = "--mem-size=512MB"
l1d_assoc = "--l1d_assoc=2"
l1i_assoc = "--l1i_assoc=2"
l2_assoc = "--l2_assoc=8"
caches = "--caches"

# Construct the full command
command = [
    gem5_binary,
    debug_flags,
    debug_file,
    config_script,
    f'--cmd="{cmd}"',
    options,
    num_cpus,
    # caches,
    num_dirs,
    ruby,
    network,
    topology,
    mesh_rows,
    link_latency,
    router_latency,
    cpu_type,
    l1d_size,
    l1i_size,
    l2_size,
    mem_size,
    l1d_assoc,
    l1i_assoc,
    l2_assoc,
]

# Run the command
subprocess.run(" ".join(command), shell=True)
