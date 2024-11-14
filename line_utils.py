import numpy as np
import pandas as pd
import skimage.measure as measure 
import scipy.special as spec

# ------------------------------ Line profile loading ------------------------------ #


# Subsequent sheets come in two flavors. Let's stick with the standard (not (2))
def load_line_profiles(workbook_path, 
                       sheet_name, 
                       expected_colnames=['Distance_(microns)', 'tub', 'RacGAP1', 'septin', 'DAPI'], 
                       res=dict()):
    """ 
    Load line profiles from Nadja's Excel sheets into a dictionary
    of pandas DataFrames. Each data frame will have the same columns
    and will be named for the Excel column preceeding the data.

    Parameters
    ----------
    workbooth_path : str
        Location of the Excel workbook.
    sheet_name : str
        Sheet name within the Excel workbook.
    expected_colnames : list, optional
        List of strings containing the names of the data column
        headers.
    res : dict, optional

    Returns
    -------
    res : dict
        A dictionary of pandas DataFrames. A pre-populated dictionary
        can be optionally passed to append the results of loading the
        line profiles from one sheet to another.
    """

    n_col = len(expected_colnames)
    
    # Open the sheet
    sheet = pd.read_excel(workbook_path, sheet_name=sheet_name, header=0)

    # Find the row containing "MAXIMUM"
    stop_row = sheet[sheet.iloc[:,0] == "MAXIMUM"].index[0]

    # Loop over columns until we do not find a block of expected columns anymore
    curr_col = 0
    while True:
        slice_col = slice(curr_col+1,curr_col+n_col+1)
        curr_colnames = sheet.columns[slice_col]
        
        # Need this brutal check because pandas does not allow duplicate column names
        if len(curr_colnames) != n_col or any([not isinstance(curr_colnames[i], str) for i in range(n_col)]) or not all([curr_colnames[i].startswith(expected_colnames[i]) for i in range(n_col)]):
            # print([not isinstance(curr_colnames[i], str) for i in range(n_col)])
            break

        # Add this set to the list, drop the empty rows, and ensure consistent column names
        res[sheet.columns[curr_col]] = sheet.iloc[:stop_row,slice_col].dropna(
            axis=0).rename(columns={k:v for k,v in zip(curr_colnames, expected_colnames)})

        # Push out to the next block
        curr_col += n_col + 1
    
    return res


def merge_df_information(source, target, id_key="Label", mapped_keys=None):
    """ 
    Merge information from source DataFrame into target DataFrame.

    Parameters
    ----------
    source : pd.DataFrame
        Source data frame
    target : pd.DataFrame
        Target data frame
    id_key : str
        Unique key to register information from source and target. Must be shared.
    mapped_keys : list
        List of keys to map from source to target.
    """
    if mapped_keys is None:
        return target

    for key in mapped_keys:
        target[key] = np.nan
        
    for i, ml in target.iterrows():
        for j, tl in source.iterrows():
            try:
                if tl[id_key].startswith(ml[id_key]):
                    for key in mapped_keys:
                        target.loc[i,key] = tl[key]
                    break
            except (AttributeError, TypeError):
                pass
    
    target = target.dropna()

    return target


# ------------------------------ Normalization ------------------------------ #

def normalize_line_profile(line_df):
    # // 2 instead of /2 deviates from Nadja's original analysis
    line_df_norm = line_df.join(line_df['Distance_(microns)']-line_df['Distance_(microns)'].max()//2,rsuffix='_norm')
    line_df_norm = line_df_norm.join(line_df[line_df.columns.drop('Distance_(microns)')].div(
        line_df[line_df.columns.drop('Distance_(microns)')].max(axis=0),axis=1),rsuffix='_norm')
    return line_df_norm


def normalize_line_profiles(line_dfs, inplace=True):
    if inplace:
        for k in line_dfs.keys():
            line_dfs[k] = normalize_line_profile(line_dfs[k])
    else:
        line_dfs_norm = {}
        for k in line_dfs.keys():
            line_dfs_norm[k] = normalize_line_profile(line_dfs[k])
        return line_dfs_norm
    

# ------------------------------ Display ------------------------------ #

def get_line_profile_endpoints(x, y, angle, length):
    # get line profile start and end coordinates for x, y, angle, length
    angle_rad = angle*np.pi/180
    length2 = length/2
    dx = length2*np.cos(angle_rad)
    dy = length2*np.sin(angle_rad)
    xl, xu = x - dx, x + dx
    yl, yu = y - dy, y + dy

    # Note that the ImageJ convention seems to have y flipped and 
    # corresponds to (xl, yu) to (xu, yl)
    return (xl, xu, yl, yu)

def imshow_with_profile(ax, im, xl, xu, yl, yu):
    # if im.shape[0] > 1:
    #     im = im.max(0)
    im_min = im.min(-1).min(-1)
    im = (im.T-im_min)/(im.max(-1).max(-1)-im_min)
    ax.imshow(im)#, cmap="gray")
    ax.set_axis_off()
    # ax.plot([xl, xu],[yu,yl],c='r')
    ax.plot([yu,yl],[xl, xu],c='r')

def image_with_profiles(fig, im, xl, xu, yl, yu, channel_names=["MTs", "septin2", "DAPI"], linewidth=25,
                        precalc_channel_profiles=None, precalc_peaks=None, precalc_fits=None,
                        precalc_means=None, precalc_midpoints=None, title=None):
    """
    Plot an image, a line profile, and the resulting plots from its line profile.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or matplotlib.figure.SubFigure
        Figure/subfigure within which to generate this figure
    im : np.typing.ArrayLike
        (C, N, M) image
    xl, xu, yl, yu : float
        Line profile start and stop coordinates in x and y
    channel_names : list
        List of strings containing human-readable names for profiles
    linewidth : int, optional
        Width of the scan, perpendicular to the line
    precalc_channel_profiles : np.typing.ArrayLike, optional
        (K, C) array of C channel profiles, which will overlay the calculated profiles.
    precalc_peaks : list
        List of lists of peaks calculated for the line profiles
    precalc_fist : list
        List of tuples of the form (function, [coefficients]) containing fits for
        line profiles.


    """
    subfigs = fig.subfigures(1,2)
    axs_left = subfigs[0].subplots(1,1)

    imshow_with_profile(axs_left, im, xl, xu, yl, yu)

    chs = measure.profile_line(im.T, 
                              [xl, yu], 
                              [xu, yl], 
                              linewidth=linewidth)
    
    nch = min(len(channel_names), chs.shape[1])
    axs_right = subfigs[1].subplots(nch,1,sharex=True)

    # colors = ['m', 'g', 'c']
    colors = ['r', 'g', 'b']
    for i in range(nch):
        ch = chs[:,i]
        axs_right[i].plot(ch, c=colors[i])

        if channel_names is not None:
            axs_right[i].annotate(channel_names[i], (0,np.mean(ch)), c=colors[i])
        
        if precalc_channel_profiles is not None:
            axs_right[i].plot(precalc_channel_profiles[:,i], c='k', linestyle='--')

        if precalc_peaks is not None and precalc_peaks[i] is not None:
            colors2 = ["orange", "cyan"]
            marker = ["x", "+"]
            if not isinstance(precalc_peaks[i], list):
                precalc_peaks[i] = [precalc_peaks[i]]
            for j, peaks in enumerate(precalc_peaks[i]):
                axs_right[i].plot(peaks, ch[peaks], marker[j], c=colors2[j], mew=2, ms=4)

        if precalc_means is not None and precalc_means[i] is not None:
            # axs_right[i].axhline(np.mean(ch), linestyle='--', color='gray')
            axs_right[i].axhline(precalc_means[i], linestyle='--', color='gray')
        
        if precalc_midpoints is not None and precalc_midpoints[i] is not None:
            # axs_right[i].axvline(len(ch)/2, linestyle='--', color='gray')
            axs_right[i].axvline(precalc_midpoints[i], linestyle='--', color='gray')

        if precalc_fits is not None and precalc_fits[i] is not None:
            func, coefs = precalc_fits[i]
            xx = np.arange(len(ch))
            axs_right[i].plot(xx, func(xx, *coefs), c="orange")
            axs_right[i].axhline(np.mean(ch), linestyle='--', color='gray')

    axs_right[i].set_xlabel("Distance (pixels)")

    if title is not None:
        subfigs[0].suptitle(title, c="white")

# ------------------------------ Profile fitting ------------------------------ #

def gauss(x, a=1, mu=0, sig=1, b=0):
    return a*np.exp(-(x-mu)**2/(sig**2)) + b

def top_hat(x, a, mu, sig, b=0):
    return b + a*(np.abs(x-mu)<sig).astype(int)

def shape_res(p0, fun, x, y):
    return fun(x,*p0) - y

def _gauss_convolved_semicircle_approx(r, t, sig):
    # Borrowed from https://github.com/bewersdorflab/nep-fitting/blob/master/nep_fitting/core/models.py
    
    return np.real((r*np.exp((-5*np.pi*(2j*r*t + 5*np.pi*sig**2))/(2.*r**2))*(12*spec.jv(1,5*np.pi)*(spec.erf((r*(r + t) - 5j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r**2 - r*t + 5j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + (spec.erf((r**2 - r*t - 5j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r*(r + t) + 5j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((10j*np.pi*t)/r)) + 60*spec.jv(1,np.pi)*((spec.erf((r**2 - r*t - 1j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r*(r + t) + 1j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((6*np.pi*(1j*r*t + 2*np.pi*sig**2))/r**2) + (spec.erf((r*(r + t) - 1j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r**2 - r*t + 1j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((4*np.pi*(1j*r*t + 3*np.pi*sig**2))/r**2)) + 20*spec.jv(1,3*np.pi)*((spec.erf((r**2 - r*t - 3j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r*(r + t) + 3j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((8*np.pi*(1j*r*t + np.pi*sig**2))/r**2) + (spec.erf((r*(r + t) - 3j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r**2 - r*t + 3j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((2*np.pi*(1j*r*t + 4*np.pi*sig**2))/r**2)) + 30*np.pi*spec.erf((r - t)/(np.sqrt(2)*sig))*np.exp((5*np.pi*(2j*r*t + 5*np.pi*sig**2))/(2.*r**2)) + 30*np.pi*spec.erf((r + t)/(np.sqrt(2)*sig))*np.exp((5*np.pi*(2j*r*t + 5*np.pi*sig**2))/(2.*r**2)) + 30*spec.jv(1,2*np.pi)*((spec.erf((r**2 - r*t - 2j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r*(r + t) + 2j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((7*np.pi*(2j*r*t + 3*np.pi*sig**2))/(2.*r**2)) + (spec.erf((r*(r + t) - 2j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r**2 - r*t + 2j*np.pi*sig**2)/(np.sqrt(2)*r*sig)))*np.exp((3*np.pi*(2j*r*t + 7*np.pi*sig**2))/(2.*r**2))) + 15*spec.jv(1,4*np.pi)*((spec.erf((r**2 - r*t - 4j*np.pi*sig**2)/(np.sqrt(2)*r*sig)) + spec.erf((r + t + (4j*np.pi*sig**2)/r)/(np.sqrt(2)*sig)))*np.exp((9*np.pi*(2j*r*t + np.pi*sig**2))/(2.*r**2)) + (-spec.erf((-r + t - (4j*np.pi*sig**2)/r)/(np.sqrt(2)*sig)) + spec.erf((r + t - (4j*np.pi*sig**2)/r)/(np.sqrt(2)*sig)))*np.exp((np.pi*(2j*r*t + 9*np.pi*sig**2))/(2.*r**2)))))/240.)

def gauss_convolved_annulus_approx(x, a, mu, sig, b, r_outer, r_inner):
    # Borrowed from https://github.com/bewersdorflab/nep-fitting/blob/master/nep_fitting/core/models.py

    t = x - mu
    # sig = psf_fwhm / 2.3548200450309493  # (2*np.sqrt(2*np.log(2)))
    
    outer = _gauss_convolved_semicircle_approx(r_outer, t, sig)
    inner = _gauss_convolved_semicircle_approx(r_inner, t, sig)

    return a*(outer - inner) + b