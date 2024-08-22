import os
import subprocess

# Define the command components
gem5_binary = "gem5/build/X86/gem5.opt"
debug_flags = "--debug-flags=RubyCache"
debug_file = "--debug-file=output.out.gz"
config_script = "gem5/configs/deprecated/example/se.py"

# Define the list of workloads
workloads = [
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/a.out",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/fft/fft",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/heat/heat",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/dct",

    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/a.out",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/fft/fft",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/heat/heat",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/dct",

    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/a.out",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/fft/fft",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/heat/heat",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/dct",

    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/a.out",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/fft/fft",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/heat/heat",
    "/home/jagadeesh/Cache_simulation/gem5/SPEC.Small/dct/dct",
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
l1d_size = "--l1d_size=64kB"
l1i_size = "--l1i_size=64kB"
l2_size = "--l2_size=512kB"
mem_size = "--mem-size=512MB"
l1d_assoc = "--l1d_assoc=2"
l1i_assoc = "--l1i_assoc=2"
l2_assoc = "--l2_assoc=8"

# Construct the full command
command = [
    gem5_binary,
    debug_flags,
    debug_file,
    config_script,
    f'--cmd="{cmd}"',
    options,
    num_cpus,
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
