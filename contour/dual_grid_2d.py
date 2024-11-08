import numpy as np

class DualGrid2D():
    def __init__(self, qt, **kwargs):
        """
        Construct a dual grid from a quad tree.
        """
        self._qt = qt  # quadtree
        
        self._vertices = []  # Dual grid vertices
        self._values = []  # Dual grid vertex values
        
        # Set non-default values
        for k, v in kwargs.items():
            setattr(self, k, v)
        
        # Construct the dual grid
        self._node_proc(0)
        
        # Cast to numpy array for ease of use
        self._vertices = np.array(self._vertices).squeeze()
        self._values = np.array(self._values).squeeze()
        
    @property
    def vertices(self):
        return self._vertices
    
    @property
    def values(self):
        return self._values
            
    def _subdivided(self, nodes):
        children = self._qt._nodes['children'][nodes]
        if len(children.shape) > 1:
            return (children.sum(1) != 0)
        return (children.sum() != 0)
            
    def _node_proc(self, node0):
        if self._subdivided(node0):
            children = self._qt._nodes['children'][node0]
            for child in children:
                self._node_proc(child)
        
            self._edge_proc_x(children[0], children[1])
            self._edge_proc_x(children[2], children[3])

            self._edge_proc_y(children[0], children[2])
            self._edge_proc_y(children[1], children[3])

            self._vert_proc(*children)
    
    def _edge_proc_x(self, node0, node1):
        children = [node0, node1, node0, node1]
        s0 = self._subdivided(node0)
        s1 = self._subdivided(node1)
        if s0:
            c0 = self._qt._nodes['children'][node0]
            children[0] = c0[1]
            children[2] = c0[3]
        if s1:
            c1 = self._qt._nodes['children'][node1]
            children[1] = c1[0]
            children[3] = c1[2]
        
        if s0 or s1:
            self._edge_proc_x(children[0], children[1])
            self._edge_proc_x(children[2], children[3])

            self._vert_proc(*children)

    def _edge_proc_y(self, node0, node1):
        children = [node0, node0, node1, node1]
        s0 = self._subdivided(node0)
        s1 = self._subdivided(node1)
        if s0:
            c0 = self._qt._nodes['children'][node0]
            children[0] = c0[2]
            children[1] = c0[3]
        if s1:
            c1 = self._qt._nodes['children'][node1]
            children[2] = c1[0]
            children[3] = c1[1]
        
        if s0 or s1:
            self._edge_proc_y(children[0], children[2])
            self._edge_proc_y(children[1], children[3])

            self._vert_proc(*children)
        
    def _vert_proc(self, node0, node1, node2, node3):
        children = [node0, node1, node2, node3]
        leaf_nodes = ~self._subdivided(children)
    
        if np.all(leaf_nodes):
            # All nodes have no children
            self._vertices.append(np.vstack([self._qt._nodes['center'][node0],
                                             self._qt._nodes['center'][node1],
                                             self._qt._nodes['center'][node2],
                                             self._qt._nodes['center'][node3]]))
            self._values.append(np.vstack([self._qt._nodes['density'][node0],
                                           self._qt._nodes['density'][node1],
                                           self._qt._nodes['density'][node2],
                                           self._qt._nodes['density'][node3]]))
            return
        
        # We need to recurse further
        if not leaf_nodes[0]:
            children[0] = self._qt._nodes['children'][node0][3]
        
        if not leaf_nodes[1]:
            children[1] = self._qt._nodes['children'][node1][2]
            
        if not leaf_nodes[2]:
            children[2] = self._qt._nodes['children'][node2][1]
            
        if not leaf_nodes[3]:
            children[3] = self._qt._nodes['children'][node3][0]
        
        self._vert_proc(*children)