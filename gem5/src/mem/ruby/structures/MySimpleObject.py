from m5.params import *
from m5.SimObject import SimObject

class MySimpleObject(SimObject):
    type = "MySimpleObject"
    cxx_header = "mem/ruby/structures/my_simple_object.hh"
    cxx_class = "gem5::MySimpleObject"