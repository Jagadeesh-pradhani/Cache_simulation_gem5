#ifndef __TUTORIAL_MY_SIMPLE_OBJECT_HH__
#define __TUTORIAL_MY_SIMPLE_OBJECT_HH__

#include "params/MySimpleObject.hh"
#include "sim/sim_object.hh"

namespace gem5
{

class MySimpleObject : public SimObject
{
  public:
    PARAMS(MySimpleObject);
    MySimpleObject(const Params &p);
};

} // namespace gem5

#endif // __TUTORIAL_MY_SIMPLE_OBJECT_HH__