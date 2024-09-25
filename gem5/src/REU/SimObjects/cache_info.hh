
#ifndef __CACHE_INFO__
#define __CACHE_INFO__

#include "mem/ruby/system/RubySystem.hh"
#include "params/CacheInfo.hh"
#include "sim/sim_object.hh"

namespace gem5
{

class CacheInfo : public SimObject
{
  private:
  void getNeighbourID();
    



  public:
    PARAMS(CacheInfo);
    CacheInfo(const Params& params);
};

} // namespace gem5

#endif // __CACHE_INFO__