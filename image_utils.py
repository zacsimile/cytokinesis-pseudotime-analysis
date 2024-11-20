import os

import numpy as np

import scipy.ndimage as ndi

from aicsimageio import AICSImage
import tifffile as tf

def pad_rot_and_trans_im(im, angle, x, y, N=2048):
    """
    Pad, rotate and translate an image based on passed sizes, angles and positions.

    Parameters
    ----------
    im : np.typing.ArrayLike
        Image (three-color).
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
    if len(im.shape) != 3:
        raise Exception("This only works for images of size 3. TODO: Extend.")

    # Pad the image
    im_pad = np.zeros((3,N,N), dtype=im.dtype)
    x_off, y_off = (N-im.shape[1])//2, (N-im.shape[2])//2
    im_pad[:,x_off:(x_off+im.shape[1]),y_off:(y_off+im.shape[2])] = im

    # Now translate
    x_shift = int(im.shape[1]//2-y)
    y_shift = int(im.shape[2]//2-x)
    im_trans = np.roll(np.roll(im_pad, x_shift, axis=1), y_shift, axis=2)

    # And rotate
    im_rot = ndi.rotate(im_trans, -angle, axes=(2,1), reshape=False)

    return im_rot

def normalize_image(im):
    """
    Normalize image to range [0,1] for plotting with matplotlib.

    Expects im axis order CYX.
    """
    im_min = im.min(-1).min(-1)
    im = (im.T-im_min)/(im.max(-1).max(-1)-im_min)

    return im

class NDImage:
    def __init__(self, image_path):
        im = AICSImage(image_path)

        self.shape = im.shape[1:]
        self._images = {}
        self.channel_names = []

        base_path = os.path.splitext(image_path)[0]
        for i, ch in enumerate(im.channel_names):
            image_path_ch = f"{base_path}_w{i+1}{ch}.TIF"
            self._images[i] = tf.imread(image_path_ch)
            self.channel_names.append(ch)

    def __getitem__(self, keys):
        if not isinstance(keys, tuple):
            chs = range(self.shape[0])[keys]
        else:
            chs = range(self.shape[0])[keys[0]]
        sliced_ds = np.empty((len(chs), *self.shape[1:]))
        for ch in chs:
            if isinstance(keys, tuple):
                sliced_ds[ch,...] = self._images[ch][keys[1:]]
            else:
                sliced_ds[ch,...] = self._images[ch]
    
        return sliced_ds
    
