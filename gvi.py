"""
* Green Visibility Index Script
"""
import time
from math import exp, hypot
from os import makedirs
from os.path import exists
from uuid import uuid1

import numpy as np
from numba import njit, prange
from numba.typed import List
from numpy import zeros, column_stack
from rasterio import open as rio_open
from rasterio.transform import rowcol
from skimage.draw import disk


def coords_2_array(a, x, y):
	"""
	* convert between coords and array position
	*  returns row,col (y,x) as expected by rasterio
	"""
	r, c = rowcol(a, x, y)
	return int(r), int(c)


@njit(fastmath=True)
def line_of_sight(pixels: np.ndarray, visible_height_arr: np.ndarray, output: np.ndarray):
    """
    This method calculates the line of sight for each pixel in a given array.

    Parameters:
    - pixels (np.ndarray): An array containing coordinates of pixels. Each pixel is represented by a tuple (r, c), where r is the row index and c is the column index.
    - visible_height_arr (np.ndarray): An array representing the height of each pixel.
    - output (np.ndarray): An array to store the line of sight information. For each pixel, if it has line of sight, the corresponding element in output array will be set to 1.

    """
    max_dydx = -500
    for r1, c1 in pixels:
        cur_dydx = visible_height_arr[r1, c1]
        if cur_dydx > max_dydx:
            max_dydx = cur_dydx
            output[r1, c1] = 1


@njit
def numba_line(r0, c0, r1, c1):
    """
    Compute the coordinates of a line between two given points using the Bresenham's line algorithm.

    Parameters:
    - r0: int, the row coordinate of the starting point.
    - c0: int, the column coordinate of the starting point.
    - r1: int, the row coordinate of the ending point.
    - c1: int, the column coordinate of the ending point.

    Returns:
    - rr: numpy.ndarray, the row coordinates of the line.
    - cc: numpy.ndarray, the column coordinates of the line.
    """
    # Compute differences and absolute differences
    dr = abs(r1 - r0)
    dc = abs(c1 - c0)
    sr = -1 if r0 > r1 else 1
    sc = -1 if c0 > c1 else 1
    err = dr - dc

    # Set up output arrays
    num_pts = max(dr, dc) + 1
    rr = np.zeros(num_pts, dtype=np.int64)
    cc = np.zeros(num_pts, dtype=np.int64)
    idx = 0

    # Main loop
    while r0 != r1 or c0 != c1:
        # Add point to output arrays
        rr[idx] = r0
        cc[idx] = c0
        idx += 1

        # Update error term and move to next pixel
        e2 = 2 * err
        if e2 > -dc:
            err -= dc
            r0 += sr
        if e2 < dr:
            err += dr
            c0 += sc

    # Add final point to output arrays
    rr[idx] = r0
    cc[idx] = c0

    return rr[:idx + 1], cc[:idx + 1]


@njit
def numba_circle_perimeter(r, c, radius, shape=None):
    """
    Generate indices of pixels of the circle perimeter.

    Parameters
    ----------
    r, c : int
        Row and column coordinates of the center of the circle.
    radius : int
        Radius of the circle.
    shape : tuple, optional
        Image shape which is used to determine the maximum extent of output pixel coordinates. If None, the shape is not
        used to determine the maximum extent of output pixel coordinates. Default is None.

    Returns
    -------
    rr, cc : ndarray of int
        Indices of pixels of the circle perimeter.

    Notes
    -----
    The perimeter coordinates are in row-column (i.e. y-x) coordinates.

     """
    if radius == 0:
        return np.array([r]), np.array([c])

    rr = []
    cc = []
    r = int(r)
    c = int(c)
    for i in range(radius + 1):
        j = int(np.sqrt(radius ** 2 - i ** 2) + 0.5)
        rr.append(i)
        cc.append(j)
        rr.append(-i)
        cc.append(j)
        rr.append(i)
        cc.append(-j)
        rr.append(-i)
        cc.append(-j)
        rr.append(j)
        cc.append(i)
        rr.append(-j)
        cc.append(i)
        rr.append(j)
        cc.append(-i)
        rr.append(-j)
        cc.append(-i)

    rr = np.asarray(rr, dtype=np.int64)
    cc = np.asarray(cc, dtype=np.int64)
    angles = np.arctan2(cc, rr)

    # Sort the points by angle
    idx = np.argsort(angles)
    rr = np.asarray(rr)[idx]
    cc = np.asarray(cc)[idx]
    if shape is not None:
        rr += r
        cc += c
        valid_idx = (rr >= 0) & (rr < shape[0]) & (cc >= 0) & (cc < shape[1])
        rr = rr[valid_idx]
        cc = cc[valid_idx]

    return rr + r, cc + c


@njit(fastmath=True)
def viewshed(pixel_line_list, visible_height_arr, output):
    """
    * Use Bresenham's Circle / Midpoint algorithm to determine endpoints for viewshed
    """
    for pixels in pixel_line_list:
        line_of_sight(pixels, visible_height_arr, output)


def distance_matrix(size, r, c, resolution):
    # Create an empty distance matrix of the same size as the input matrix
    dist_matrix = np.zeros((size, size), dtype=np.float32)

    # Iterate through each cell of the matrix and calculate the distance to the center cell
    for i in range(size):
        for j in range(size):
            dist_matrix[i, j] = hypot(i - r, j - c) * resolution
    return dist_matrix


def process_part(mask):
    """
    * main function for running with parallel.py
    """
    t = time.time()
    # create an output array at the same dimensions as data for output
    gvi = zeros(mask["dsm"].shape, dtype=np.float32)

    # radius in pixels
    radius_px = int(mask["options"]["radius"] // mask['meta']['transform'][0])

    # build weighting mask

    # get pixel references for aoi extents
    min_r, min_c = coords_2_array(mask["meta"]["transform"], mask["aoi"].bounds[0], mask["aoi"].bounds[3])
    max_r, max_c = coords_2_array(mask["meta"]["transform"], mask["aoi"].bounds[2], mask["aoi"].bounds[1])

    pixel_line_list = List(mask["pixel_line_list"])

    # Set up the kernel configuration.
    dtm_o_height = mask["dtm"] + mask["options"]["o_height"]
    r_slices = List([slice(r - radius_px, r + radius_px + 1) for r in range(radius_px, max_c - min_c + radius_px + 1)])
    c_slices = List([slice(c - radius_px, c + radius_px + 1) for c in range(radius_px, max_c - min_c + radius_px + 1)])
    calculate_green_visibility_index_for_part(gvi,
                                              dtm_o_height,  # observer height
                                              mask["dsm"],  # dsm dataset
                                              mask["green"],
                                              max_c - min_c + radius_px,
                                              max_r - min_r + radius_px,
                                              radius_px, mask["weights"], pixel_line_list, mask['distance_arr'],
                                              r_slices, c_slices)

    # clip gvi to aoi bounds
    gvi = gvi[radius_px:gvi.shape[0] - radius_px, radius_px:gvi.shape[1] - radius_px]

    # check that tmp folder exists
    if not exists('./tmp/'):
        makedirs('tmp')

    # make unique filename
    filepath = f'./tmp/{str(uuid1())}.tif'

    # output file with updated dimensions and transform
    with rio_open(filepath, 'w',
                  driver=mask["meta"]["driver"],
                  height=gvi.shape[0],
                  width=gvi.shape[1],
                  count=mask["meta"]["count"],
                  dtype='float32',
                  crs=mask["meta"]["crs"],
                  transform=mask['meta']['transform'],
                  tiled=True,
                  compress='deflate'
                  ) as dst:
        dst.write(gvi, 1)
    # print("part done", time.time() - t)
    # return the filepath to the result
    return filepath


def create_los_lines(radius_px):
    los_path_rr, los_path_cc = numba_circle_perimeter(radius_px, radius_px, radius_px)
    pixel_line_list = []
    previous_pixels = None
    i = 0
    for r, c in column_stack((los_path_rr, los_path_cc)):
        # calculate line of sight to each pixel
        pixels = column_stack(numba_line(radius_px, radius_px, r, c))[1:]
        pixels = pixels.astype(np.int32)
        if previous_pixels is not None:
            new_part = pixels.copy()
            for pp, pn in zip(previous_pixels, pixels):
                if np.alltrue(np.equal(pp, pn)):
                    partial_line_mask = ~((new_part[:, 0] == pp[0]) & (new_part[:, 1] == pp[1]))
                    new_part = new_part[partial_line_mask]
                else:
                    break
            if len(new_part) != len(pixels):
                new_part = np.vstack((pixels[len(pixels) - len(new_part) - 1], new_part))

            pixel_line_list.append(new_part)
        else:
            pixel_line_list.append(pixels)
        # previous_pixels = pixels  #deactivate double line skipping, because the code is memory bandwith bound, so calculating extra pixels is faster then keeping track of already done pixels
    print(sum([len(a) for a in pixel_line_list]))
    return pixel_line_list


def create_weighting_mask(resolution, radius_px):
    weighting_mask = zeros((radius_px * 2 + 1, radius_px * 2 + 1), dtype=np.float32)
    for r, c in column_stack(disk((radius_px, radius_px), radius_px, shape=weighting_mask.shape)):
        weighting_mask[(r, c)] = exp(-0.0003 * (hypot(radius_px - c, radius_px - r) * resolution))
    return weighting_mask


@njit(fastmath=True, parallel=True)
def calculate_green_visibility_index_for_part(gvi: np.ndarray, dtm_o_height: np.ndarray, dsm: np.ndarray, green: np.ndarray, max_c: int, max_r: int,
                                              radius_px: int, weighting_mask: np.ndarray, pixel_line_list: np.ndarray, distance_arr: np.ndarray,
                                              r_slices, c_slices):
    # Preallocate memory for output
    all_outputs = [np.zeros((2 * radius_px + 1, 2 * radius_px + 1), dtype=np.byte)
                   for i in range(radius_px, max_r + 1)]

    for r in prange(radius_px, max_r + 1):
        # Directly pick the buffer based on the loop index
        output = all_outputs[r - radius_px]

        for c in range(radius_px, max_c + 1):
            # Calculate the slice indices only once
            r_slice = r_slices[r - radius_px]
            c_slice = c_slices[c - radius_px]

            # Extract the necessary slices
            dsm_data = dsm[r_slice, c_slice] - dtm_o_height[r, c]
            visible_height_arr = dsm_data / distance_arr
            output.fill(0)
            # set the start location as visible automatically
            output[radius_px, radius_px] = 1
            # Calculate viewshed and update output and reached_height_arr
            viewshed(pixel_line_list, visible_height_arr, output)

            # Extract the viewshed data from the output surface and apply weighting mask
            visible = output * weighting_mask

            # Multiply extract of (weighted) viewshed with extract of (weighted) green dataset
            visible_green = visible * green[r_slice, c_slice] * weighting_mask

            # Calculate the ratio for greenness in the view
            gvi[r, c] = np.sum(visible_green) / np.sum(visible)
    return gvi


"""
* Do not call this script directly
"""
if __name__ == '__main__':
    print("please call this script using parallel.py")
