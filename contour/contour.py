import numpy as np

VERTEX_DTYPE = [('position', '2f4'), ('segment', 'i4'), ('curvature', 'f4')]
SEGMENT_DTYPE = [('v0', 'i4'), ('v1', 'i4'), ('next', 'i4'), ('prev', 'i4')]

class Contour():
    def __init__(self, lines, **kwargs):
        self._vertices = None
        self._segments = None
        self._centroid = None
        self._bbox = None
        self._lines = None
        self._diameter = None
        self._curvature = False
        self._ordered_vertices = []
                
        self._initialize_contour(lines)
        
        # Set non-default values
        for k, v in kwargs.items():
            setattr(self, k, v)
            
    def _initialize_contour(self, lines):
        # Get the unique vertices
        lines_sorted = lines.reshape(2*lines.shape[0],2)
        v, idx, inv, cnt = np.unique(lines_sorted,axis=0, 
                                    return_index = True,
                                    return_inverse=True, 
                                    return_counts=True)
        self._vertices = np.zeros(v.shape[0],dtype=VERTEX_DTYPE)
        self._vertices['position'] = v
        self._vertices['segment'] = np.repeat(np.arange(lines.shape[0]),2)[idx]
        
        # Set up the line segments
        self._segments=np.zeros(lines.shape[0],dtype=SEGMENT_DTYPE)
        self._segments['v0'] = inv[::2]
        self._segments['v1'] = inv[1::2]
        self._segments['next'] = -1
        self._segments['prev'] = -1

        # Get line segment ordering
        curr = self._vertices['segment'][np.argmin(cnt)]  # if the loop isn't closed, pick the an endpoint
        visited = np.zeros(lines.shape[0], dtype=bool)  # mark vertices we've already been to
        safety_counter = 0  # we shouldn't need this
        while True:
            # prevent an infinite loop
            safety_counter += 1
            if safety_counter > v.shape[0]:
                break

            # mark whichever vertex we are at as accounted for
            visited[curr] = True

            # if we've been to every segement on the contour, we're done
            if np.all(visited):
                break

            # Loop over all other segments and find the one
            # adjacent to this one
            curr_seg = self._segments[curr]
            for i, seg in enumerate(self._segments):
                if visited[i]:
                    # Don't allow us to move backwards
                    continue
                # print(seg)
                if (curr_seg['v0'] == seg['v0']) \
                    | (curr_seg['v0'] == seg['v1']) \
                    | (curr_seg['v1'] == seg['v0']) \
                    | (curr_seg['v1'] == seg['v1']):

                    # now check the winding
                    # we want curr_seg['v1'] = seg['v0']
                    # first order curr_seg
                    if (curr_seg['v0'] == seg['v0']) or (curr_seg['v0'] == seg['v1']):
                        v0 = curr_seg['v0']
                        curr_seg['v0'] = curr_seg['v1']
                        curr_seg['v1'] = v0
                    
                    # then order seg
                    if curr_seg['v1'] == seg['v1']:
                        v0 = seg['v0']
                        seg['v0'] = seg['v1']
                        seg['v1'] = v0

                    self._ordered_vertices.extend([curr_seg['v0'], curr_seg['v1']])

                    # step forward
                    curr_seg['next'] = i
                    seg['prev'] = curr
                    curr = i
                    
                    break

    @property
    def bbox(self):
        if self._bbox is None:
            self._bbox = [np.min(self._vertices['position'][:,0]), 
                          np.max(self._vertices['position'][:,0]), 
                          np.min(self._vertices['position'][:,1]), 
                          np.max(self._vertices['position'][:,1])]
        return self._bbox
    
    @property
    def diameter(self):
        if self._diameter is None:
            self._diameter = np.sqrt((self.bbox[1]-self.bbox[0])**2+(self.bbox[3]-self.bbox[2])**2)
        return self._diameter
    
    @property
    def centroid(self):
        if self._centroid is None:
            self._centroid = [np.mean(self.vertices[:,0]), 
                              np.mean(self.vertices[:,1])]
        return self._centroid
    
    @property
    def lines(self):
        if self._lines is None:
            self._lines = self._vertices['position'][
                np.array([self._segments['v0'],
                          self._segments['v1']])
                ].swapaxes(0,1)
        return self._lines
    
    @property
    def vertices(self):
        return self._vertices[self._ordered_vertices]['position']
    
    @property
    def perimeter(self):
        return np.sqrt((self.lines[:,0]-self.lines[:,1])**2).sum()
    
    @property
    def curvature(self):
        if not self._curvature:
            for i, v in enumerate(self._vertices):
                # get the two edges
                seg0 = self._segments[v['segment']]
                
                if seg0['v0'] == i:
                    if seg0['prev'] == -1:
                        v['curvature'] = np.nan
                        continue
                    seg0 = self._segments[seg0['prev']]
                    seg1 = self._segments[v['segment']]
                elif seg0['v1'] == i:
                    if seg0['next'] == -1:
                        v['curvature'] = np.nan
                        continue
                    seg1 = self._segments[seg0['next']]
                else:
                    raise UserWarning(f"Segment {v['segment']} does not contain vertex {i}.")

                # get their normals
                seg0v0 = self._vertices[seg0['v0']]['position']
                seg0v1 = self._vertices[seg0['v1']]['position']
                seg1v0 = self._vertices[seg1['v0']]['position']
                seg1v1 = self._vertices[seg1['v1']]['position']
                n0 = (seg0v1 - seg0v0)[::-1]*np.array([-1,1])
                n1 = (seg1v1 - seg1v0)[::-1]*np.array([-1,1])
                nmag = np.sqrt(((n1-n0)**2).sum())

                # get their centroids
                c0 = (seg0v0+seg0v1)/2
                c1 = (seg1v0+seg1v1)/2
                cmag = np.sqrt(((c1-c0)**2).sum())

                # get their curvature
                v['curvature'] = nmag/cmag

            self._curvature = True
                
        return self._vertices['curvature']