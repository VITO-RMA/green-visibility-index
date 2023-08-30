"""
* Divide and Conquour script to divide, process in parallel and merge again without having
*	to worry about edge effects.
* This version is intended to be used with call gvi.py
"""
import argparse
import os
from datetime import datetime
import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import numpy as np
import rasterio
from affine import Affine
from rasterio import open as rio_open
from rasterio.merge import merge as rio_merge
from rasterio.windows import from_bounds
from shapely.geometry import Polygon
from tqdm import tqdm

from gvi import process_part, create_weighting_mask, create_los_lines, distance_matrix


def extractPolygon(src, poly):
    """
    * extract a subset from a raster according to the specified bounds
    * returns the dataset (numpy array) and metadata object
    """

    # Get the extent of the geometry
    geometry_bounds = poly[0].bounds
    resolution = src.res[0]
    geometry_bounds_rounded = (np.floor(geometry_bounds[0] / resolution) * resolution,
                               np.floor(geometry_bounds[1] / resolution) * resolution,
                               np.ceil(geometry_bounds[2] / resolution) * resolution,
                               np.ceil(geometry_bounds[3] / resolution) * resolution)

    transform = src.transform  # get the coordinate transformation

    geometry_window = from_bounds(*geometry_bounds_rounded, transform=transform)
    # data, transform = rio_mask(src, poly, crop=True, all_touched=True)
    data = src.read(1, window=geometry_window)
    clipped_data = np.full((max(int(geometry_window.height), data.shape[0]), max(data.shape[1], int(geometry_window.width))), fill_value=0, dtype=np.float32)
    transform = rasterio.windows.transform(geometry_window, transform)
    if data.shape == clipped_data.shape:
        clipped_data = data
    else:

        # Clip the raster using the geometry
        # Update the portion of clipped_data that overlaps with the geometry
        window = rasterio.windows.from_bounds(*geometry_bounds_rounded, transform=transform)
        row_offset = abs(int(window.row_off)) if window.row_off < 0 else 0
        col_offset = abs(int(window.col_off)) if window.col_off < 0 else 0
        clipped_data[row_offset:row_offset + int(data.shape[0]), col_offset:col_offset + int(data.shape[1])] = data

    # create the metadata for the dataset
    out_meta = src.meta
    out_meta.update({
        "height": clipped_data.shape[0],
        "width": clipped_data.shape[1],
        "transform": transform
    })

    # return the mask and the metadata
    return clipped_data.astype(np.float32), out_meta


def run(res, padding, landbouw=True, blauw=True, grid='grid_vl', total_parts=1, part_nr=0, dsm_path=None, green_path=None, grid_path=None,
        output_name_postfix="",year=2015, overwrite=False):
    global time, meta, bounds, transform, path
    landbouw_str = "_landbouw"
    if not landbouw:
        landbouw_str = ""
    blauw_str = "_blauw"
    if not blauw:
        blauw_str = ""
    if dsm_path is None:
        dsm_path = rf"input_data/{year}/DSM_{res}m.tif"
    dtm_path = rf"input_data/2015/DTM_{res}m.tif"
    if green_path is None:
        green_path = rf"input_data/{year}/green_01{landbouw_str}{blauw_str}_{res}m.tif"
    if grid_path is None:
        grid_path = rf"input_data/{grid}.gpkg"
    out_path = rf"output/{year}/green_vis_vl_{res}m_{padding}m{landbouw_str}{blauw_str}_{part_nr}{output_name_postfix}.tif"
    print(out_path)
    if os.path.exists(out_path):
        if overwrite:
            os.remove(out_path)
        else:
            print(f"output already exists for this part_nr {part_nr}. skipping!")
            exit()
    # log start time and log start
    time = datetime.now()
    print(f"processes, started at {time}.")
    # initialise masks array for the results
    masks = []
    grid_df = gpd.read_file(grid_path)
    # get aoi bounds from shapefile
    # read in raster data
    grid_split_df = np.array_split(grid_df, total_parts)
    grid_part_df = grid_split_df[part_nr]
    with rio_open(dtm_path) as dtm_input:
        with rio_open(dsm_path) as dsm_input:
            with rio_open(green_path) as green_input:
                meta = dtm_input.meta
                if (dsm_input.width == dtm_input.width == green_input.width) == False or \
                        (dsm_input.height == dtm_input.height == green_input.height) == False:

                    print("rasters do not match!")
                    print("width: \t\t", dsm_input.width == dtm_input.width == green_input.width)
                    print("height: \t", dsm_input.height == dtm_input.height == green_input.height)
                    if (dsm_input.res[0] == dtm_input.res[0] == green_input.res[0]) == False:
                        print("resolution: \t", dsm_input.res[0] == dtm_input.res[0] == green_input.res[0])
                        exit()
                # build weighting mask
                resolution = meta['transform'][0]
                radius_px = int(padding // resolution)
                weighting_mask = create_weighting_mask(resolution, radius_px)
                pixel_line_list = create_los_lines(radius_px)
                distance_arr = distance_matrix((radius_px * 2) + 1, radius_px, radius_px, resolution)

                for i, row in tqdm(grid_part_df.iterrows(), total=len(grid_part_df)):
                    # verify raster dimensions and resolution
                    bounds = row.geometry.bounds
                    xmin, ymin, xmax, ymax = bounds
                    transform = dtm_input.transform  # get the coordinate transformation
                    raster_extent = rasterio.windows.bounds(rasterio.windows.Window(col_off=0, row_off=0, width=dtm_input.width, height=dtm_input.height), transform=transform)
                    pixel_size = dtm_input.transform[0]
                    x1 = np.floor((xmin - raster_extent[0]) / pixel_size) * pixel_size + raster_extent[0]
                    y1 = raster_extent[3] - (np.ceil((raster_extent[3] - ymin) / pixel_size) * pixel_size)
                    x2 = np.ceil((xmax - raster_extent[0]) / pixel_size) * pixel_size + raster_extent[0]
                    y2 = raster_extent[3] - (np.floor((raster_extent[3] - ymax) / pixel_size) * pixel_size)
                    bounds = (x1, y1, x2, y2)
                    # construct a Shapely polygon for use in processing
                    polygon = Polygon([
                        (bounds[0] - (padding), bounds[1] - (padding)),  # bl
                        (bounds[0] - (padding), bounds[3] + (padding + resolution)),  # tl
                        (bounds[2] + (padding + resolution), bounds[3] + (padding + resolution)),  # tr
                        (bounds[2] + (padding + resolution), bounds[1] - (padding))  # br
                    ])

                    # extract the polygon and append to masks list
                    dtm, meta = extractPolygon(dtm_input, [polygon])
                    dsm, _ = extractPolygon(dsm_input, [polygon])
                    green, _ = extractPolygon(green_input, [polygon])

                    # fix dsm for ocean
                    dsm = np.where(dsm < dtm, dtm, dsm)

                    xmin, ymin, xmax, ymax = bounds
                    x_res, y_res = green_input.res
                    new_transform = Affine(x_res, 0, xmin, 0, -y_res, ymax)

                    # Create a new window object based on the new bounds
                    new_window = from_bounds(xmin, ymin, xmax, ymax, new_transform)

                    # Update the metadata for the new raster
                    new_meta = green_input.meta.copy()
                    new_meta.update({
                        'width': new_window.width,
                        'height': new_window.height,
                        'transform': new_transform})
                    # make result object and append to masks list
                    masks.append({
                        "dtm": dtm,
                        "dsm": dsm,
                        "green": green,
                        "meta": new_meta,
                        "aoi": row.geometry,
                        "weights": weighting_mask,
                        "distance_arr": distance_arr,
                        "pixel_line_list": pixel_line_list,
                        "options": {
                            "radius": padding,  # viewshed radius
                            "o_height": 1.7,  # observer height
                            "t_height": 0  # target height
                        }
                    })

                    # print(masks[0])
                    # exit()

            # output files for debugging purposes
            # outputFiles(dtm_input.crs, masks)
            # exit()
    # make as many processes as are required and launch them
    results = []
    for m in tqdm(masks, total=len(masks)):
        results.append(process_part(m))
    # with Pool(processes) as p:
    #     results = p.map(process_part, masks)
    # open all result files in read mode
    files = []
    for filepath in results:
        files.append(rio_open(filepath, 'r'))
    # merge result files
    merged, out_transform = rio_merge(files)
    # update the metadata
    out_meta = files[0].meta.copy()
    out_meta.update({
        "height": merged.shape[1],
        "width": merged.shape[2],
        "transform": out_transform,
        "dtype": 'float32'
    })
    # create a raster file
    with rio_open(out_path, 'w', **out_meta) as dst:
        dst.write(merged[0], 1)
    # use GDAL binary to calculate histogram and statistics
    # call(["gdalinfo", "-stats", out_path], stdout=DEVNULL)
    # close all of the files
    for file in files:
        file.close()
    for path in results:
        os.remove(path)
    # print how long it took
    print(datetime.now() - time)
    print("done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run green visibility index calculation.')
    parser.add_argument('total_parts', metavar='-tp', type=int, help='number of parts to split')
    parser.add_argument('part_nr', type=int, metavar='-p', help='part to process')
    parser.add_argument('-l', '--landbouw', action='store_true', help='enable landbouw')
    parser.add_argument('-b', '--blauw', action='store_true', help='enable blauw')

    args = parser.parse_args()
    total_parts = args.total_parts
    part_nr = args.part_nr
    landbouw = args.landbouw
    blauw = args.blauw
    year = 2022
    res = 1
    view_distances = [800]
    overwrite=True
    # total_parts = 10
    # part_nr = 0
    print(landbouw, blauw)
    grid = 'grid_vl'
    for view_distance in view_distances:
        run(res, view_distance, landbouw=landbouw, blauw=blauw, grid=grid, total_parts=total_parts, part_nr=part_nr, output_name_postfix=grid,year=year, overwrite=overwrite)
