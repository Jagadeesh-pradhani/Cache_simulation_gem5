from gem5.resources.resource import obtain_resource

resource = obtain_resource("x86-ubuntu-18.04-img")

print(f"The resource is available at {resource.get_local_path()}")