import numpy as np
# Quadtree based on David Baddeley's _octree code (https://github.com/python-microscopy/python-microscopy)

QUAD_SHIFT = np.zeros((4,2))

for n in range(4):
    QUAD_SHIFT[n,0] = 2*(n&1)-1
    QUAD_SHIFT[n,1] = (n&2)-1

INITIAL_NODES = 1000
NODE_DTYPE = [('depth', 'i4'), ('children', '4i4'), ('parent', 'i4'), ('center', '2f4'), ('density', 'f4')]

class ImageQuadtree():
    def __init__(self, threshold, density_function=np.std, **kwargs):
        """
        Quadtree based on pixel values.
        
        Parameters
        ----------
            threshold : float
                Isocontour threshold. Currently std of pixel values, but can change.
            density_function : function 
                Function applied to subdivided region of the image. Maps MxM subregion
                to a single value, which we call density. These values are compared
                to the threshold.
        """
        self._threshold = threshold
        self._density_function = density_function
        self._image = None  # Assumes 8-bit
        self._bounds = None  # Quadtree bounds
        self._xwidth = None
        self._ywidth = None
        self.max_depth = 5
        self._nodes = np.zeros(INITIAL_NODES, NODE_DTYPE)
        self._next_node = 1  # Index of next node to create
        self._resize_limit = INITIAL_NODES  # Maximum number of nodes
        
        # Set non-default values
        for k, v in kwargs.items():
            setattr(self, k, v)
            
    def _resize(self):
        old_nodes = self._nodes
        new_size = int(self._nodes.shape[0]*1.5 + 0.5)
        self._nodes = np.zeros(new_size, NODE_DTYPE)
        self._nodes[:self._next_node] = old_nodes
        self._resize_limit = new_size
        
    def _add_node(self, depth, parent, center):
        if self._next_node == self._resize_limit:
            self._resize()
        self._nodes['depth'][self._next_node] = depth
        self._nodes['parent'][self._next_node] = parent
        self._nodes['center'][self._next_node] = center
        self._nodes['density'][self._next_node] = self._density(self._next_node)
        self._next_node += 1
    
    def add_image(self, image):
        self._image = image
        if not self._bounds:
            self.set_bounds(0, self._image.shape[1], 0, self._image.shape[0])
        self._nodes['center'][0] = [(self._bounds[0]+self._bounds[1])/2.0,
                                    (self._bounds[2]+self._bounds[3])/2.0]
        self._nodes['density'] = np.mean(self._image)
        self._divide()
        
    def set_bounds(self, xl, xu, yl, yu):
        """
        Sets bounds on the quadtree. l is lower bound, u is upper bound, 
        x or y is axis.
        """
        self._bounds = [xl, xu, yl, yu]
        self._xwidth = xu-xl
        self._ywidth = yu-yl
        
    def _get_bounds(self, node_idx):
        node = self._nodes[node_idx]
        xw, yw = self._box_size(node['depth'])
        xs = xw/2.0
        ys = yw/2.0
        return int(node['center'][0]-xs), int(node['center'][0]+xs), int(node['center'][1]-ys), int(node['center'][1]+ys)

    def _box_size(self, depth):
        scale = 2**depth    
        return self._xwidth/scale, self._ywidth/scale
    
    def _density(self, node_idx):
        """
        Density of the quadtree node, used to compare to threshold.
        """
        xl, xu, yl, yu = self._get_bounds(node_idx)
        return self._density_function(self._image[yl:yu,xl:xu])
        
    def _divide(self):
        """
        The engine room. Actually subdivides the image.
        """
        if self._image is None:
            raise KeyError('No image to divide.')
        
        node_idx = 0
        while node_idx < self._next_node:
            node = self._nodes[node_idx]
            
            if (node_idx > 0) and (node['depth'] == 0):
                # We've somehow hit the empty node zone (we shouldn't be able to do this)
                print('Made it to the other world.')
                break
            
            node_idx += 1
            
            if node['density'] <= self._threshold:
                # We don't need to subdivide
                continue
            
            if node['depth'] >= self.max_depth:
                # We do not need to subdivide this node
                continue
            
            for _i in np.arange(4):
                # Subdivide
                new_center = node['center'] + 0.5*QUAD_SHIFT[_i]*self._box_size(node['depth']+1)
                node['children'][_i] = self._next_node
                self._add_node(node['depth']+1, node_idx-1, new_center)