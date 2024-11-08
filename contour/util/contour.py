import numpy as np
from ..contour import Contour

def find_connections(lines, eps=1e-12):
    components = np.ones(lines.shape[0],dtype=int)*-1
    component = 0

    for _i in np.arange(lines.shape[0]):
        if components[_i] == -1:
            components[_i] = component
            component += 1
        for _j in np.arange(lines.shape[0]):
            if _i == _j:
                continue
            is_neighbor = np.prod(np.abs(lines[_i,0,:]-lines[_j,1,:])<eps)|np.prod(np.abs(lines[_i,1,:]-lines[_j,0,:])<eps)|np.prod(np.abs(lines[_i,0,:]-lines[_j,0,:])<eps)|np.prod(np.abs(lines[_i,1,:]-lines[_j,1,:])<eps)
            if is_neighbor:
                if components[_j] != -1:
                    min_component = np.minimum(components[_i], components[_j])
                    components[components==components[_i]] = min_component
                    components[components==components[_j]] = min_component
                else:
                    components[_j] = components[_i]
                    
    return components

def find_contours(lines, eps=1e-12):
    components = find_connections(lines, eps)
    contours = []
    for component in np.unique(components):
        contours.append(Contour(lines[components==component]))
    return contours