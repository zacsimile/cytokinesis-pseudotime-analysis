import numpy as np
import scipy.ndimage as ndi

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