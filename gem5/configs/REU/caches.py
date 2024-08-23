from m5.objects import Cache


class L1Cache(Cache):
    assoc = 2
    tag_latency = 2
    data_latency = 2
    response_latency = 2
    mshrs = 4
    tgts_per_mshr = 20

    def connectCPU(self, cpu):
        # need to define this in a base class!
        raise NotImplementedError

    def connectBus(self, bus):
        self.mem_side = bus.cpu_side_ports


class L1ICache(L1Cache):
    size = "32kB"

    def connectCPU(self, cpu):
        self.cpu_side = cpu.icache_port


class L1DCache(L1Cache):
    size = "32kB"

    def connectCPU(self, cpu):
        self.cpu_side = cpu.dcache_port

    def prefetcher_initialize(
        self, cpu, cache_line_size
    ):  # Added cache_line_size argument
        self.nextline_prefetcher = NextLinePrefetcher(
            self, cpu, cache_line_size
        )

    def prefetcher_access(self, addr):
        self.nextline_prefetcher.prefetch(addr)


class L2Cache(Cache):
    size = "512kB"
    assoc = 8
    tag_latency = 20
    data_latency = 20
    response_latency = 20
    mshrs = 20
    tgts_per_mshr = 12

    def connectCPUSideBus(self, bus):
        self.cpu_side = bus.mem_side_ports

    def connectMemSideBus(self, bus):
        self.mem_side = bus.cpu_side_ports


class NextLinePrefetcher:
    def __init__(self, cache, cpu, cache_line_size):
        self.cache = cache
        self.cpu = cpu
        self.cache_line_size = cache_line_size

    def prefetch(self, addr):
        next_line_addr = (
            addr + self.cache_line_size
        )  # Calculate the next cache line address
        self.cpu.initiateMemRead(next_line_addr)  # Prefetch the next line