import os

import numpy as np
import scipy.ndimage as ndi
import tifffile as tf

import math_utils

DEG_TO_RAD = np.pi / 180.0

# Store already-calculated coordinates
CACHED_IMAGE_CYLINDRICAL_COORDS = {}

def pad_rot_and_trans_im(im, angle, x, y, N=2048, crop_length=None):
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

    # If crop_length is provided, crop to this
    if crop_length is not None:
        crop_length2 = crop_length // 2
        xc, yc = im_trans.shape[-1]//2, im_trans.shape[-2]//2
        im_trans = im_trans[...,(yc-crop_length2):(yc+crop_length2),(xc-crop_length2):(xc+crop_length2)]

    # And rotate
    im_rot = ndi.rotate(im_trans, -angle, axes=(-1,-2), reshape=False, order=0)

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

    # Unique key for coordinates that match this case
    cache_key = f"{'_'.join([str(x) for x in shape])}_{xc}_{yc}_{zc}_{ang_xy}_{dx}_{dy}_{dz}"

    try:
        r, theta, z = CACHED_IMAGE_CYLINDRICAL_COORDS[cache_key]
    except KeyError:
        Z, Y, X = np.meshgrid(*[range(n) for n in shape], indexing='ij')
        Z, Y, X = dz*(Z - zc), dy*(Y - yc), dx*(X - xc)  # center

        ang = ang_xy*DEG_TO_RAD  # convert to radians
        sin_ang, cos_ang = np.sin(ang), np.cos(ang)

        # establish cylindrical coordinate space along line 
        r, theta, z = math_utils.cart2cyl(Z, 
                                        X*sin_ang+Y*cos_ang, 
                                        X*cos_ang-Y*sin_ang)
        
        # Store for later
        CACHED_IMAGE_CYLINDRICAL_COORDS[cache_key] = (r, theta, z)

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

def target_names(targets, key):
    """ 
    Construct a list of the key name plus any aliases. Important for searching
    through file names.
    """
    try:
        target = targets[key]
    except KeyError:
        return []
    names = [key]
    try:
        names.extend(target["alias"])
    except KeyError:
        pass
    return names

def get_channel_orders(fn, wvls, n_ch, targets, desired_channel_order):
    """
    Get the sorting from the image to the 
    """

    # find wavelengths in file name and sort from high to low
    wvls_dict, binned_wvls = extract_channel_targets_from_filename(fn, wvls=wvls)

    # Establish target names in this data set and sort from high to low to match image load
    channel_targets = [wvls_dict[str(wvl)] for wvl in sorted(binned_wvls)[::-1]]

    # the last channel is always DAPI, if unknown
    if len(channel_targets) < n_ch:
        channel_targets.append("DAPI") 

    # Now find the resorting of the channels according to their target position
    channel_order = []
    group_channel_order = []
    mt_channel = 0
    for j, ch in enumerate(desired_channel_order):
        for opt in target_names(targets, ch):
            try:
                channel_order.append(channel_targets.index(opt))
                group_channel_order.append(j)
                if ch == "MTs":
                    mt_channel = channel_targets.index(opt)
            except ValueError:
                pass
    assert len(channel_order) == n_ch #len(desired_channel_order)
    print(f"  channel_targets: {channel_targets} channel_order: {channel_order} group_channel_order: {group_channel_order}")

    return channel_order, group_channel_order, mt_channel

def radial_projection(im, nbins, bin_size, x, y, angle, z_coord, dx=0.09, dy=0.09, dz=1):
    r, _, _ = image_cylindrical_coordinates(im, x, y, z_coord, angle, dx=dx, dy=dy, dz=dz)
    # bin_size = min(min(dx, dy), dz) #r.max()/nbins
    rbins = np.arange(nbins+1) * bin_size

    proj = np.zeros((im.shape[0], nbins, im.shape[-1]), dtype=im.dtype)

    for i in range(nbins):
        proj[:,i,:] = (im*((r>=rbins[i])&(r<rbins[i+1]))[None,...]).sum(-3).sum(-2)

    return proj

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
    
