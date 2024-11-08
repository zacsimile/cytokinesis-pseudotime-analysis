import numpy as np

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
