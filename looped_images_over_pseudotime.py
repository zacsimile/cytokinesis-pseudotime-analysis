# Imports
import os
import time
import logging

import tifffile as tf
import numpy as np

import scipy.ndimage as ndi

import line_utils
import image_utils
import file_utils

logger = logging.getLogger('pseudotime')
logging.basicConfig(
    filename='pseudotime_run.log',
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')
logger.addHandler(logging.StreamHandler())

targets = {}  # Start with an empty dictionary DO NOT DELETE

# Describe all of our target proteins here
# Any protein that does not have a specific workbook and image directory associated
# will be considered a general marker, available across all workbooks
targets = file_utils.load_targets('targets.yaml')

# Order of stages
time_key = "Stage"
# time_order = ["CF", "RC", "CS", "RS", "SM", "BA", "A"]
time_order = ["RC", "CS", "RS", "SM", "BA", "A"]


# time_order = [ "RC", "CS", "RS", "SM0", "SM1", "SM2", "SM3", "BA0", "BA1", "BA2", "BA3", "BA4", "A"]

# Don't fit the septin ring locations at these time points
# Added CF here because I expect no septin peaks at the furrow
time_do_not_fit = ["CF", "BA", "BA0", "BA1", "BA2", "BA3", "BA4", "A"]

# Channels per image (TODO: Auto detect)
n_ch = 4

# wavelengths to be found in the file names
# Sublists are grouped. First element of the sublist is a group name.
# NOTE: First element must be a number!
wvls = [488,[568, "orange"],[646,647,657]]

# desired channel order, specified by keys in targets
# must include MTs, septin and DAPI
desired_channel_order = ["MTs", "septin", "DAPI", "CellMask", "MKLP1", "RacGAP1", "PRC1", "Cit-K", "anillin", "myoIIA", "myoIIB", "actin", "Septin7", "Septin11", "Septin9", "BORG4", "Tsg101", "Tsg101-ab83m", "ALIXrb", "ALIXm", "IST1", "CHMP4B"]
# desired_channel_order = ["MTs", "septin", "DAPI", "CellMask", "MKLP1", "RacGAP1", "PRC1", "Cit-K", "anillin", "myoIIA", "myoIIB", "actin", "Septin7", "Septin11", "Septin9", "BORG4", "Tsg101", "ALIXrb", "ALIXm", "IST1", "CHMP4B"]
# desired_channel_order = ["MTs", "septin", "DAPI", "IST1", "CHMP4B"]#, "Cit-K", "PRC1"]

# desired_channel_order = ["MTs", "septin", "DAPI", "MKLP1", "RacGAP1", "anillin", "myoIIA", "myoIIB", "Cit-K", "CellMask", "PRC1", "actin"]
# desired_channel_order = ["MTs", "septin", "DAPI", "MKLP1", "RacGAP1", "anillin", "myoIIA", "Cit-K", "CellMask", "PRC1", "actin"]
# Length of cropped pseudotime region (should be roughly the line length)
length = 500
# length = 1000

# Choose the mode, one of "z-stack", "mean-proj", "mean-proj-individual" or "radial-proj"
# They are mostly self-explanatory, but mean-proj-individual produces one slice per mean projection
# Note: we can only do one at a time at the moment
# mode = "mean-proj"
# mode = "mean-proj-individual"
mode = "radial-proj"
# mode = "z-stack"

# In the worst case (one tubule is sandwiched at the top of the stack and the other
# at the bottom), this must be 2*<max stack length>-1 
num_planes = 100

# pixel sizes (we assume they are constant)
dx, dy, dz = 0.09, 0.09, 1

# Before we do anything, let's make sure all of our targets exist
for key in desired_channel_order:
    try:
        targets[key]
    except KeyError:
        raise KeyError(f"Element {key} does not exist in targets dictionary!")
    
# ...and let's make sure our mode is supported
assert mode in ["z-stack", "mean-proj", "mean-proj-individual", "radial-proj"]

# Load data from the workbooks
metrics = file_utils.load_workbooks(targets, desired_channel_order)
metrics = metrics[metrics[time_key].isin(time_order)]

# Set this to true if we want to return a z-stack
z_stack = mode == "z-stack"
radial_proj = mode == "radial-proj"
mean_proj_individual = mode == "mean-proj-individual"

groups = metrics.groupby(time_key)

plot_stack = None
n_groups = len(groups)
l2 = length // 2
rot_width = int((np.sqrt(2) * length) + 1)

if mean_proj_individual:
    num_planes = metrics.value_counts("Stage").max()
elif not z_stack:
    num_planes = 1

if radial_proj:
    group_img = np.zeros((n_groups, len(desired_channel_order), num_planes, length//8, length)).squeeze()
elif num_planes == 1:
    group_img = np.zeros((n_groups, len(desired_channel_order), length, length))
else:
    group_img = np.zeros((n_groups, len(desired_channel_order), num_planes, length, length))


# Establish columns for septin peaks (X12, X22) and distance between them (dX2)
metrics['dX2'], metrics['X12'], metrics['X22'] = np.nan, np.nan, np.nan

for group, tup in enumerate(groups):
    name, entries = tup
    n_group = len(entries)
    logger.info(f"{name}: {n_group} averaged")
    im_proj = {}

    # In a first pass, fit the septin ring distances for registration
    for i, ml in entries.iterrows():
        if ml[time_key] not in time_do_not_fit:
            # If we are in a class where it makes sense...
            logger.debug(f"  Septin ring fit for {os.path.basename(ml['filename'])}")

            # Get the image associated with this row and load it with the channels sorted from high to low
            im = image_utils.NDImage(ml["filename"], load_sorted=True)

            # get x, y, angle for this row
            x, y, angle = ml[["X", "Y", "Angle"]]

            # find wavelengths in file name and sort from high to low
            wvls_dict, binned_wvls = image_utils.extract_channel_targets_from_filename(ml["filename"], wvls=wvls)

            # Establish target names in this data set and sort from high to low to match image load
            channel_targets = [wvls_dict[str(wvl)] for wvl in sorted(binned_wvls)[::-1]]

            # the last channel is always DAPI, if unknown
            if len(channel_targets) < n_ch:
                channel_targets.append("DAPI") 

            # ... get the septin peaks
            mt_ch = [i for i, t in enumerate(channel_targets) if any([t == n for n in image_utils.target_names(targets, "MTs")])][0]
            septin_ch = [i for i, t in enumerate(channel_targets) if any([t == n for n in image_utils.target_names(targets, "septin")])][0]
            p0, p1, dX2 = line_utils.find_septin_peaks(im[:].mean(1).squeeze(), x, y, angle, length,
                                                        mt_ch=mt_ch, 
                                                        septin_ch=septin_ch)

            metrics.loc[i,['X12','X22','dX2']] = [p0, p1, dX2]


    # Now compute the average distance
    mean_dX2 = entries['dX2'].mean()

    # In our second pass, average these images
    for t, tup2 in enumerate(entries.groupby("target")):
        name2, entries2 = tup2
        n_target = len(entries2)
        logger.info(f"  {name2}: {n_target} averaged")
        if mean_proj_individual:
            k = 0
        for i, ml in entries2.iterrows():
            logger.info(f"  Analyzing {os.path.basename(ml['filename'])}")

            start = time.time()

            # Get the image associated with this row
            im = image_utils.NDImage(ml["filename"], load_sorted=True)

            stop = time.time()
            duration = stop-start
            logger.debug(f"  time to load image: {duration:.2f} s")
            start = time.time()

            channel_order, group_channel_order, mt_channel = image_utils.get_channel_orders(ml["filename"], 
                                                                                            wvls,
                                                                                            n_ch, 
                                                                                            targets, 
                                                                                            desired_channel_order)

            # get x, y, angle for this row
            x, y, angle = ml[["X", "Y", "Angle"]]

            # Rotate the image  # CYX
            im_rot = image_utils.pad_rot_and_trans_im(im[:], angle, x, y, crop_length=rot_width)

            stop = time.time()
            duration = stop-start
            logger.debug(f"  time to rotate image: {duration:.2f} s")
            start = time.time()

            if z_stack or radial_proj:
                im_min = im_rot.min(-1).min(-1).min(-1)
                im_rot = (im_rot - im_min[:, None, None, None]) / (im_rot.max(-1).max(-1).max(-1) - im_min)[:, None, None, None]
            else:
                im_rot = im_rot.mean(1).squeeze()

                im_min = im_rot.min(-1).min(-1)
                im_rot = (im_rot - im_min[:, None, None]) / (im_rot.max(-1).max(-1) - im_min)[:, None, None]

            stop = time.time()
            duration = stop-start
            logger.debug(f"  time to normalize image: {duration:.2f} s")
            start = time.time()

            # Grab coordinates
            xc, yc = im_rot.shape[-1]//2, im_rot.shape[-2]//2

            if z_stack or radial_proj:
                # if ml[time_key] not in time_do_not_fit:
                #     # Grab the z-coordinate of the central bit of the tubule, cast to integer
                #     # z_coord = int(round(line_utils.find_central_pos(im[:].max(2).squeeze(), ml["X"], ch=mt_channel)))
                #     z_coord = int(round(line_utils.find_central_pos(im[:].max(2).squeeze(), xc, ch=mt_channel)))
                # else:
                #     z_coord = im[:].shape[-3] // 2
                z_coord = int(round(line_utils.find_central_pos(im_rot[...,(yc-12):(yc+13),:].sum(2).squeeze(), xc, ch=mt_channel)))
                metrics.loc[i, "z_coord"] = z_coord
                logger.info(f"  im.shape: {im_rot.shape} projection shape: {im_rot[:].max(2).squeeze().shape} z_coord: {z_coord}")

            stop = time.time()
            duration = stop-start
            logger.debug(f"  time to get centroid: {duration:.2f} s")
            start = time.time()

            logger.info(f"  xc: {xc}  yc: {yc}  length: {length}")

            # Crop the image
            im_crop = im_rot[...,(yc-l2):(yc+l2),(xc-l2):(xc+l2)]

            logger.debug(f"  im_crop shape: {im_crop.shape}")

            stop = time.time()
            duration = stop-start
            logger.debug(f"  time to crop image: {duration:.2f} s")
            start = time.time()

            # rescale the image
            # if np.isnan(ml["dX (pxl)"]):
            if np.isnan(ml["dX2"]):
                mag = 1
                # im_zoom = im_crop
            else:
                # mag = ml["dX (pxl)"]/mean_dX
                mag = ml["dX2"]/mean_dX2

            if z_stack or radial_proj:
                zmag = dz // dx if radial_proj else 1
                im_zoom = ndi.zoom(im_crop, (1,zmag,1,mag), order=0)
            else:
                im_zoom = ndi.zoom(im_crop, (1,1,mag), order=0)

            stop = time.time()
            duration = stop-start
            logger.debug(f"  time to zoom image: {duration:.2f} s")
            start = time.time()

            # Crop the image again
            xc, yc = im_zoom.shape[-1]//2, im_zoom.shape[-2]//2
            im_crop2 = im_zoom[...,(yc-l2):(yc+l2),(xc-l2):(xc+l2)]
            logger.debug(f"  im_crop2 shape: {im_crop2.shape}")

            if radial_proj:
                num_planes2_crop = im_crop.shape[-3] // 2
                num_planes2 = im_crop2.shape[-3] // 2
                new_z_coord = num_planes2 - (num_planes2_crop - z_coord)*zmag
                im_crop2 = image_utils.radial_projection(im_crop2, length//8, 1,
                                                         l2, l2, 0, 
                                                         new_z_coord, dx=1, dy=1, dz=1)
                
                stop = time.time()
                duration = stop-start
                logger.debug(f"  time to get radial proj image: {duration:.2f} s")

            # Add the image with a weighting 1/length of the group 
            if z_stack:
                z_length = im_crop2.shape[-3]
                num_planes2 = num_planes // 2
                zl, zu = num_planes2 - z_coord - 1, num_planes2 + z_length - z_coord - 1

                group_img[group,group_channel_order,zl:zu,...] += (im_crop2[channel_order]/np.array([n_group, n_group, n_group, n_target])[:,None,None,None])
            elif mean_proj_individual:
                logger.debug(f"  group_img shape: {group_img[group,group_channel_order,k,...].shape}")
                group_img[group,group_channel_order,k,...] += (im_crop2[channel_order]/np.array([n_group, n_group, n_group, n_target])[:,None,None])
                k += 1
            else:
                group_img[group,group_channel_order,...] += (im_crop2[channel_order]/np.array([n_group, n_group, n_group, n_target])[:,None,None])
                
group_order = list(groups[time_key].unique().keys())
group_img_sorted = [group_order.index(g) for g in time_order if g in group_order]
logger.debug(group_img_sorted)

stack_fn = f'pseudotime_images_{mode}_{"_".join(["".join([x[0:2],x[-1]]) for x in desired_channel_order])}.ome.tif'
if z_stack or mean_proj_individual:
    tf.imwrite(stack_fn, group_img[group_img_sorted,...], metadata={'axes': 'TCZYX'}, dtype=group_img.dtype)
else:
    tf.imwrite(stack_fn, group_img[group_img_sorted,...], metadata={'axes': 'TCYX'}, dtype=group_img.dtype)
