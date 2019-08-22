# from http://www.roguebasin.com/index.php?title=Bresenham%27s_Line_Algorithm
def get_line(start, end):
    # Setup initial conditions
    x1, y1, z1 = start
    x2, y2, z2 = end
    dx = x2 - x1
    dy = y2 - y1
 
    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)
 
    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2
 
    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        z1, z2 = z2, z1
        swapped = True
 
    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
 
    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1
 
    zstep = dz / max(dx, 1)
    
    # Iterate over bounding box generating points between start and end
    y = y1
    z = z1
    for x in range(x1, x2 + 1):
        coord = (y, x, z) if is_steep else (x, y, z)
        yield coord
        error -= abs(dy)
        z += zstep
        if error < 0:
            y += ystep
            error += dx
 