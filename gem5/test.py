def getNorthCache(curr, size):
    north = None
    south = None
    east = None
    west = None

    row = curr / size
    col = curr % size

    if row > 0:
        north = (row - 1) * size + col

    if row < size - 1:
        south = (row + 1) * size + col

    if col > 0:
        west = row * size + (col - 1)

    if col < size - 1:
        east = row * size + (col + 1)
    
    return north,south,east,west

north, south, east, west = getNorthCache(4,4)

print(f"North : {north}, East : {east}, West : {west}, South : {south}")