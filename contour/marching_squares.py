import numpy as np

EDGES = np.array([[0,1],[0,2],[1,3],[2,3]])
# Cases wind segments in down, right
CASES = np.array([[-1,-1,-1,-1],
                  [0,1,-1,-1],
                  [0,2,-1,-1],
                  [1,2,-1,-1],
                  [1,3,-1,-1],
                  [0,3,-1,-1],
                  [0,1,2,3],
                  [2,3,-1,-1],
                  [2,3,-1,-1],
                  [0,2,1,3],
                  [0,3,-1,-1],
                  [1,3,-1,-1],
                  [1,2,-1,-1],
                  [0,2,-1,-1],
                  [0,1,-1,-1],
                  [-1,-1,-1,-1]])
CASE_MASK = np.array([1,2,4,8])

class MarchingSquares():
    def __init__(self, isolevel=1, **kwargs):
        """
        Implementation of marching squares algorithm.
        
        Parameters
        ----------
            isolevel : float
                Threshold for contouring
            vertices : np.array
                N x 4 x 2 vertices defining a series of squares
            values : np.array
                N x 4 vertex values for comparison to isolevel
        """
        self._vertices = None
        self._values = None
        self._isolevel = isolevel
        self._lines = []
        
        # Set non-default values
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    def add_points(self, vertices, values):
        self._vertices = vertices
        self._values = values
    
    def _case(self, values):
        if len(values.shape) > 1:
            return ((values < self._isolevel)*CASE_MASK).sum(1)
        return ((values < self._isolevel)*CASE_MASK).sum()
    
    def march(self, eps=1e-12):
        self._lines = []  # empty
        for _sq_idx in np.arange(self._vertices.shape[0]):
            case = self._case(self._values[_sq_idx])
            if (case == 0) or (case == 15):
                continue
            cases = CASES[case]
            cases_mask = (cases != -1)
            edges_idx = EDGES[cases[cases_mask]]
            
            # Check for saddle ambiguities
            if (case == 6) or (case == 9):
                vals = self._values[_sq_idx,edges_idx]
                center_val = vals.mean()
                if (center_val < self._isolevel):
                    if (case == 6):
                        case = 9
                    if (case == 9):
                        case = 6
                    cases = CASES[case]
                    cases_mask = (cases != -1)
                    edges_idx = EDGES[cases[cases_mask]]
            
            vert = self._vertices[_sq_idx,edges_idx]  # M (2 x # edges) x 2 (v1, v2) x 2 (x, y)
            vals = self._values[_sq_idx,edges_idx]  # M (2 x # edges) x 2 (v1, v2)
            # # Sort radially from 0
            # r = np.argsort((vert*vert).sum(2),axis=1)
            # _, row = np.meshgrid(np.arange(r.shape[1]),np.arange(r.shape[0]))
            # try:
            #     vert = vert[row,r]
            # except(IndexError):
            #     print(vert.shape, r.shape, row.shape)
            #     print(vert)
            #     print(r)
            #     print(row)
            #     vals = vals[row,r]
            denom = vals[:,1]-vals[:,0]  # M
            vert_interp = vert[:,0,:] + ((self._isolevel-vals[:,0])/denom)[:,None]*(vert[:,1,:]-vert[:,0,:])  # M (2 x # edges) x 2 (x, y)
            self._lines.append(vert_interp)
            idxs = (np.abs(denom)==0)
            vert_interp[idxs] = vert[idxs,0,:]  # div by zero case
            idxs = (np.abs(vals[:,0]-vals[:,1])<eps)
            vert_interp[idxs] = vert[idxs,0,:]  # close to v0
            idxs = (np.abs(vals[:,0]-self._isolevel)<eps)
            vert_interp[idxs] = vert[idxs,0,:]  # close to v0
            idxs = (np.abs(vals[:,1]-self._isolevel)<eps)
            vert_interp[idxs] = vert[idxs,1,:]  # close to v1
        self._lines = np.vstack(self._lines)
        self._lines = self._lines.reshape(-1,2,2,order='C')  # P (# lines) x 2 (# vertices/line) x 2 (x, y)