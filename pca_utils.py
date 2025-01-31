import numpy as np
from sklearn.decomposition import PCA

def dpoly(p):
    """
    Take the derivative of a polynomial of the form returned by np.polyfit,

    Parameters
    ----------
    p : list
        Polynomial coefficients in the order p(x) = p[0] * x**deg + ... + p[deg]
        where deg = len(p)-1

    Returns
    -------
    list
        Polynomial coefficients for derivative of p. Degree is deg-1.
    
    """
    N = len(p)
    dp = []
    for i, coef in enumerate(p[:-1]):
        dp.append((N-i-1)*coef)
    return dp

def point_poly_dist(x0, y0, p, return_coords=True):
    """ 

    Calculate the distance from a point to a polynomial p. Approximates this using 
    the point-line distance where the line is calculated from the derivative of the
    polynomial at x0.
    
    Parameters
    ----------
        x0, y0 : float
            (x0i, y0i) is a point, where i is an index into x0 and y0
        p : list
            Poly as defined in np.polyfit 
        return_coords : bool
            Return the coordinates of the point projection onto p.
    
    Returns
    -------
        dist : float
            Distance from the points (x0, y0) to the polynomial p.
        xp, yp : float, optional
            The projection of the points (x0, y0) onto the polynomial.
    """

    # find the derivative line at x
    # y = m*x+b
    z = np.poly1d(p)
    dfit = np.poly1d(dpoly(p))
    b = z(x0)-dfit(x0)*x0
    m = dfit(x0)

    dist = np.abs(b+m*x0-y0)/(1+m*m)

    if return_coords:
        xp = (x0 + m*y0 - m*b)/(m*m+1)
        yp = m*xp+b

        return (dist, xp, yp)

    return dist

def pca(features, n_components=2, transform=False, fit=False):
    """
    Normalize and compute pca on features

    Parameters
    ----------
    features : np.typing.ArrayLike
        (n_samples, n_features) matrix.
    
    Returns
    -------
    pca : sklearn.decomposition.PCA
        PCA object
    """
    # z-norm features
    features = (features - features.mean(0))/features.std(0)

    # Compute PCA
    pca = PCA(n_components=n_components)
    pca.fit(features)

    if transform:
        xx, yy = pca.transform(features).T[:2,:]
        if fit:
            fit = np.polyfit(xx, yy, 3)
            return pca, (xx, yy), fit
        return pca, (xx, yy)

    return pca


def sort_by_point_plane_dist(xx, yy, fit, nbins=None):
    """

    Sort features by point-plane distance from (xx,yy) to
    the line described by fit.

    Parameters
    ----------
    xx, yy : np.typing.ArrayLike
        2D point locations to map onto line
    fit : np.typing.ArrayLike
        Poly as defined in np.polyfit 

    Returns
    -------
    permutation : np.typing.ArrayLike
        Vector of indices sorting points
    """

    # Sort along line after collapsing via point-plane distance
    _, xp, _ = point_poly_dist(xx, yy, fit)
    permutation = np.argsort(xp)
    
    # Rebin as needed
    if nbins is not None:
        permutation = np.array_split(permutation, nbins)
    
    return permutation
