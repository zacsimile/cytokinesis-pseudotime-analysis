import os

import numpy as np
import scipy.ndimage as ndi
import tifffile as tf

import math_utils

def pad_rot_and_trans_im(im, angle, x, y, N=2048):
    """
    Pad, rotate and translate an image based on passed sizes, angles and positions.

    Parameters
    ----------
    im : np.typing.ArrayLike
        Image (CYX).
    angle : float
        Angle to rotate the image in degrees.
    x, y : float
        (x,y) is the coordinate on which to center the image.
    N : int
        Target size of the image. Expands im to be padded with zeros out to target
        size. Original image is centered in the padded image.
    """
    if np.any([N < x for x in im.shape]):
        raise Exception("Padding is smaller than image size. Increase N.")
    if len(im.shape) != 3 and len(im.shape) != 4:
        raise Exception("This only works for images of dimension 3 or 4. TODO: Extend.")

    # Pad the image
    if len(im.shape) == 3:
        im_pad = np.zeros((im.shape[0],N,N), dtype=im.dtype)
    elif len(im.shape) == 4:
        im_pad = np.zeros((im.shape[0],im.shape[1],N,N), dtype=im.dtype)
    x_off, y_off = (N-im.shape[-2])//2, (N-im.shape[-1])//2
    im_pad[...,x_off:(x_off+im.shape[-2]),y_off:(y_off+im.shape[-1])] = im

    # Now translate
    x_shift = int(im.shape[-2]//2-y)
    y_shift = int(im.shape[-1]//2-x)
    im_trans = np.roll(np.roll(im_pad, x_shift, axis=-2), y_shift, axis=-1)

    # And rotate
    im_rot = ndi.rotate(im_trans, -angle, axes=(-1,-2), reshape=False)

    return im_rot

def normalize_image(im):
    """
    Normalize image to range [0,1] for plotting with matplotlib.

    Expects im axis order CYX.
    """
    im_min = im.min(-1).min(-1)
    im = (im.T-im_min)/(im.max(-1).max(-1)-im_min)

    return im

def image_cylindrical_coordinates(im, xc, yc, zc, ang_xy, dx=1, dy=1, dz=1):
    """
    Compute the cylindirical coordinates of a 3D image along a line
    with the centroid (xc, yc, zc) at an angle ang_xy in the xy plane.

    Parameters
    ----------
    im : np.typing.ArrayLike
        CZYX or ZYX
    xc, yc, zc : float
        Line centroid in Euclidean space
    ang_xy : float
        Angle of the line in xy space
    dx, dy, dz : float
        pixel size along x, y and z
    """
    if len(im.shape) > 3:
        shape = im.shape[1:]
    else:
        shape = im.shape

    Z, Y, X = np.meshgrid(*[range(n) for n in shape], indexing='ij')
    Z, Y, X = dz*(Z - zc), dy*(Y - yc), dx*(X - xc)  # center

    ang = ang_xy*np.pi/180  # convert to radians

    # establish cylindrical coordinate space along line 
    r, theta, z = math_utils.cart2cyl(Z, 
                                      X*np.sin(ang)+Y*np.cos(ang), 
                                      X*np.cos(ang)-Y*np.sin(ang))

    return r, theta, z


def extract_channel_targets_from_filename(fn, wvls=[488, 568, 647], return_binned=True):
    """ 
    For each wavelength, return the target name. 

    Parameters
    ----------
    fn : str
        Filename
    wvls : list
        List of wavelengths to check for. Can be multiple wavelengths per list 
        element (binning).

    Returns
    -------
    chs_dict : dict
        Dictionary of targets keyed on wavelength
    """
    binned_wvl = []
    chs_dict = {}
    for wvl in wvls:
        if isinstance(wvl, list):
            for w in wvl:
                if str(w) in fn:
                    # choose wvl[0] as the bin
                    if wvl[0] not in binned_wvl:
                        binned_wvl.append(wvl[0])
                    target = fn[:(fn.index(str(w))-1)].split('_')[-1]
                    chs_dict[str(wvl[0])] = target
                    # break
        else:
            if str(wvl) in fn:
                if wvl not in binned_wvl:
                    binned_wvl.append(wvl)
                target = fn[:(fn.index(str(wvl))-1)].split('_')[-1]
                chs_dict[str(wvl)] = target

    if return_binned:
        return chs_dict, binned_wvl

    return chs_dict

class Image:
    def __init__(self, image_path, dc=1, dz=1, dy=1, dx=1):
        self._image_path = image_path

        self._shape_c = 1
        self._shape_z = 1
        self._shape_y = 1
        self._shape_x = 1

        self._dc = dc
        self._dz = dz
        self._dy = dy
        self._dx = dx

        self._metadata = {}

    @property
    def shape(self):
        """ CZYX """
        return (self._shape_c, self._shape_z, self._shape_y, self._shape_x)
    
    @property
    def voxel_size(self):
        """ Return voxel size in um """
        return (self._dc, self._dz, self._dy, self._dx)
    
    def __getitem__(self, keys):
        raise NotImplementedError("Implemented in a derived class.")

class NDImage(Image):
    def __init__(self, image_path, dc=1, dz=1, dy=1, dx=1, load_sorted=True):
        """Read ND files with TIFs.

        Example format:
            "NDInfoFile", Version 2.0
            "Description", Multi Dimensions Experiment
            "StartTime1", 20240719 14:19:15
            "DoTimelapse", FALSE
            "NTimePoints", 1
            "DoStage", FALSE
            "DoWave", TRUE
            "NWavelengths", 3
            "WaveName1", "CSU561"
            "WaveDoZ1", TRUE
            "WaveName2", "CSU491"
            "WaveDoZ2", TRUE
            "WaveName3", "CSU405 QUAD"
            "WaveDoZ3", TRUE
            "DoZSeries", TRUE
            "NZSteps", 66
            "ZStepSize", 1
            "WaveInFileName", TRUE
            "NEvents", 0
            "EndFile"

        """

        super().__init__(image_path, dc, dz, dy, dx)

        with open(self._image_path, "r") as fp:
            data = fp.read().splitlines()
            # Get rid of EndFile and everything after
            data = data[:data.index('"EndFile"')]
            data_list = [line.split(',') for line in data]
            self._metadata = {el[0].strip().replace('"',""): el[1].strip().replace('"',"") for el in data_list}

        # self.shape = im.shape[1:]
        self._shape_c = int(self._metadata["NWavelengths"])
        self._shape_z = int(self._metadata["NZSteps"])
        self._images = {}
        if self._shape_c > 0:
            self.channel_names = [self._metadata[f"WaveName{n+1}"] for n in range(self._shape_c)]

            if load_sorted:
                """ Load image so channels are sorted high to low by wavelength. """
                sorted_ch = np.argsort([int(x.split('CSU')[1].split(' ')[0]) for x in self.channel_names])[::-1]
            else:
                sorted_ch = range(len(self.channel_names))

            print("  loading ", self.channel_names, sorted_ch)

            base_path = os.path.splitext(image_path)[0]
            for i, ch in enumerate(sorted_ch):
                image_path_ch = f"{base_path}_w{ch+1}{self.channel_names[ch]}.TIF"
                try:
                    self._images[i] = tf.imread(image_path_ch)
                except FileNotFoundError:
                    # See if it saved the space in channel name with an underscore
                    image_path_ch = f"{base_path}_w{ch+1}{self.channel_names[ch].replace(' ', '_')}.TIF"
                    self._images[i] = tf.imread(image_path_ch)

            z, y, x = self._images[0].shape
            assert z == self._shape_z
            self._shape_y = y
            self._shape_x = x

    def __getitem__(self, keys):
        if not isinstance(keys, tuple):
            chs = range(self.shape[0])[keys]
        else:
            chs = np.arange(self.shape[0])[keys[0]]

        new_stack = []        
        for ch in chs:
            if isinstance(keys, tuple):
                new_stack.append(self._images[ch][keys[1:]])
            else:
                new_stack.append(self._images[ch])
    
        return np.stack(new_stack, axis=0)
    
