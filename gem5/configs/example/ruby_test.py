import argparse
import os
import sys

import m5
from m5.defines import buildEnv
from m5.objects import *
from m5.util import addToPath

addToPath("../")

from common import Options
from ruby import Ruby

"""Set Up Argument Parsing
"""

parser = argparse.ArgumentParser()
Options.addNoISAOptions(parser)

parser.add_argument(
    "--maxloads", metavar="N", default=100, help="Stop after N loads"
)
parser.add_argument(
    "-f",
    "--wakeup_freq",
    metavar="N",
    default=10,
    help="Wakeup every N cycles",
)

# Add Ruby and protocol specific options
Ruby.define_options(parser)

args = parser.parse_args()

"""Configure Cache Sizes and Associativities
"""

args.l1d_size = "32kB"
args.l1i_size = "32kB"
args.l2_size = "256kB"
args.l1d_assoc = 4
args.l1i_assoc = 4
args.l2_assoc = 8


"""Create the Ruby Tester
"""

check_flush = False
if buildEnv["PROTOCOL"] == "MOESI_hammer":
    check_flush = True

tester = RubyTester(
    check_flush=check_flush,
    checks_to_complete=args.maxloads,
    wakeup_frequency=args.wakeup_freq,
)


"""Create the M5 System
"""

system = System(cpu=tester, mem_ranges=[AddrRange(args.mem_size)])

system.voltage_domain = VoltageDomain(voltage=args.sys_voltage)
system.clk_domain = SrcClockDomain(
    clock=args.sys_clock, voltage_domain=system.voltage_domain
)


"""Create CPU List and Ruby System
"""

cpu_list = [system.cpu] * args.num_cpus
Ruby.create_system(args, False, system, cpus=cpu_list)

system.ruby.clk_domain = SrcClockDomain(
    clock=args.ruby_clock, voltage_domain=system.voltage_domain
)

assert args.num_cpus == len(system.ruby._cpu_ports)
tester.num_cpus = len(system.ruby._cpu_ports)


"""Configure Randomization and Tie Ruby Ports
"""

system.ruby.randomization = True

for ruby_port in system.ruby._cpu_ports:
    if ruby_port.support_data_reqs and ruby_port.support_inst_reqs:
        tester.cpuInstDataPort = ruby_port.in_ports
    elif ruby_port.support_data_reqs:
        tester.cpuDataPort = ruby_port.in_ports
    elif ruby_port.support_inst_reqs:
        tester.cpuInstPort = ruby_port.in_ports

    ruby_port.no_retry_on_stall = True
    ruby_port.using_ruby_tester = True


"""Set Up the Simulation
"""

root = Root(full_system=False, system=system)
root.system.mem_mode = "timing"

m5.ticks.setGlobalFrequency("1ns")
m5.instantiate()

exit_event = m5.simulate(args.abs_max_tick)
print("Exiting @ tick", m5.curTick(), "because", exit_event.getCause())

