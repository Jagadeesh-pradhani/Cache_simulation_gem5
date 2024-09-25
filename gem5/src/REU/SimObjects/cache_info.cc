#include "REU/SimObjects/cache_info.hh"
#include "mem/ruby/system/RubySystem.hh"
#include "mem/ruby/structures/CacheMemory.hh"
#include "mem/ruby/network/garnet/Router.hh"


#include "base/trace.hh"
#include "debug/CacheInfoFlag.hh"


namespace gem5
{

CacheInfo::CacheInfo(const Params &params):
    SimObject(params)
{
    DPRINTF(CacheInfoFlag, "%s: Hello World! From a "
                    "SimObject (constructor).\n", __func__);
}

void 
CacheInfo::getNeighbourID() {
    // int id = get_id();

}



} // namespace gem5