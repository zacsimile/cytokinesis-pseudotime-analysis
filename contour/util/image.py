import numpy as np

def crop_image_from_contour(image, contour):
    bbox = np.round(contour.bbox).astype(np.int)
    return image[bbox[2]:bbox[3],bbox[0]:bbox[1]]
