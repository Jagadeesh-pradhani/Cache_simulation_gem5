#include "mem/ruby/structures/my_simple_object.hh"

#include <iostream>

namespace gem5
{

MySimpleObject::MySimpleObject(const Params &params) :
    SimObject(params)
{
    std::cout << "Hello World! From a SimObject!" << std::endl;
}

} // namespace gem5