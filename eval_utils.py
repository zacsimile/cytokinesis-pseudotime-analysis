import numpy as np

def min_chi2_dist(x01,x02,x11,x12):
    """
    Measure how close two point pairs are to one another. Compute the chi2 distance
    between sets of points. We find the pairing between (x01, x02) and (x11, x12)
    that minimizes the chi2 distance and report this. Useful for measuring how 
    accurately we predict two points.

    Paramaters
    ----------
    x01 : np.typing.ArrayLike
        First point in first set
    x02 : np.typing.ArrayLike
        Second point in first set
    x11 : np.typing.ArrayLike
        First point in second set
    x12 : np.typing.ArrayLike
        Second point in second set
    
    Returns
    -------
    chi2_dist : npt.typing.ArrayLike
        The sum of the squares of the minimum distances between the points in the first
        set and the points in the second set.
    """
    x1x1 = (x01-x11)**2
    x1x2 = (x01-x12)**2
    x2x1 = (x02-x11)**2
    x2x2 = (x02-x12)**2

    # compute the sum of the squares of the minimum distance pairings
    chi2_dist = -1*np.ones(x1x1.shape[0])

    # x1 is closer to x12, so we combine this with x2 to x22
    chi2_dist[x1x1 < x1x2] = x1x1[x1x1 < x1x2] + x2x2[x1x1 < x1x2]
    chi2_dist[x1x1 > x1x2] = x1x2[x1x1 > x1x2] + x2x1[x1x1 > x1x2]

    return chi2_dist