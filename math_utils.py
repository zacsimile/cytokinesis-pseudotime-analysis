import numpy as np

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(hxy, z)
    az = np.arctan2(y, x)
    return az, el, r

def cart2cyl(x,y,z):
    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(y,x)

    return r, theta, z
