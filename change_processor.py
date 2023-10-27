import os

import numpy as np

import geopandas as gpd
import rasterio
from rasterio.features import rasterize

os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
from parallel import extractPolygon, run
from rasterio import open as rio_open
import rasterio
import numpy as np
from rasterio.features import rasterize


def create_sphere_array(radius, height):
    size = radius * 2 + 1
    sphere = np.zeros((size, size))
    center = radius

    for x in range(size):
        for y in range(size):
            distance = np.sqrt((x - center) ** 2 + (y - center) ** 2)
            if distance <= radius:
                sphere[x, y] = height * np.sqrt(1 - (distance / radius) ** 2)

    return sphere


def create_geometry_raster(geom_df, height, width, transform):
    # Obtain the transform and dimensions of the raster
    # Step 2: Convert the GeoPandas DataFrame to GeoJSON-like format
    geoms = geom_df.geometry
    shapes = ((geom, 1) for geom in geoms)  # Create iterable of (geometry, value) pairs

    # Step 3: Rasterize the GeoJSON-like data
    rasterized = rasterize(shapes=shapes, out_shape=(height, width), transform=transform)

    # Step 4: Optionally, write the rasterized image to a new GeoTIFF file
    return rasterized


# Further processing or analysis with the rasterized image
# ...

def place_small_array(dsm_arr, dtm_arr, small_array, center_x, center_y):
    small_height, small_width = small_array.shape
    big_height, big_width = dsm_arr.shape

    start_x = center_x - small_width // 2
    start_y = center_y - small_height // 2 - 1

    end_x = start_x + small_width
    end_y = start_y + small_height

    # Ensure the start and end indices are within bounds
    start_x = max(0, start_x)
    start_y = max(0, start_y)
    end_x = min(end_x, big_width)
    end_y = min(end_y, big_height)

    # Calculate the indices for placing the small array
    small_start_x = max(0, -start_x)
    small_start_y = max(0, -start_y)
    small_end_x = small_start_x + (end_x - start_x)
    small_end_y = small_start_y + (end_y - start_y)

    # Place the small array into the big array
    small_arr_extract = small_array[small_start_y:small_end_y, small_start_x:small_end_x]
    # Place the small array into the big array
    dsm_arr[start_y:end_y, start_x:end_x] = np.where(small_arr_extract != 0,
                                                     small_arr_extract + dtm_arr[start_y:end_y, start_x:end_x],
                                                     dsm_arr[start_y:end_y, start_x:end_x])

    return dsm_arr


def add_tree(x, y, height, size, res, view_distance, landbouw=False, blauw=False, year=2015):
    landbouw_str = "_landbouw"
    if not landbouw:
        landbouw_str = ""
    blauw_str = "_blauw"
    if not blauw:
        blauw_str = ""
    dsm_path = rf"input_data/{year}/DSM_{res}m.tif"
    dtm_path = rf"input_data/{year}/DTM_{res}m.tif"
    green_path = rf"input_data/{year}/green_01{landbouw_str}{blauw_str}_{res}m.tif"
    with rio_open(dtm_path) as src:
        row, col = src.index(x, y)
        cell_centroid = rasterio.transform.xy(src.transform, row, col, offset='center')
        snapped_x, snapped_y = cell_centroid
    df = gpd.GeoDataFrame(geometry=gpd.points_from_xy([snapped_x], [snapped_y], crs='epsg:31370'))
    grid = gpd.points_from_xy([snapped_x], [snapped_y], crs='epsg:31370').buffer(view_distance * 2.05 + res / 2).envelope
    df['geometry'] = df.buffer(view_distance + res / 2).envelope
    sphere_arr = create_sphere_array(int(size / res), height)
    with rio_open(dtm_path) as dtm_input:
        with rio_open(dsm_path) as dsm_input:
            with rio_open(green_path) as green_input:
                dsm_arr, t_poly = extractPolygon(dsm_input, grid)
                dtm_arr, _ = extractPolygon(dtm_input, grid)
                green_arr, _ = extractPolygon(green_input, grid)
                row, col = rasterio.transform.rowcol(t_poly['transform'], snapped_x, snapped_y)
    new_dsm_arr = place_small_array(dsm_arr, dtm_arr, sphere_arr, int(row), int(col))
    new_green_arr = place_small_array(green_arr, np.zeros_like(green_arr), np.where(sphere_arr > 0, 1, 0), int(row), int(col))
    # Define the output filename
    output_grid_filename = 'input_data/tree_box.gpkg'
    output_dsm_filename = 'input_data/dsm_tree.tif'
    output_green_filename = 'input_data/green_tree.tif'
    df.to_file(output_grid_filename)

    # Open the output GeoTIFF file in write mode
    with rasterio.open(output_dsm_filename, 'w', driver='GTiff', height=new_dsm_arr.shape[0],
                       width=new_dsm_arr.shape[1], count=1, dtype=new_dsm_arr.dtype, crs='EPSG:31370', transform=t_poly['transform']) as dst:
        # Write the array data to the file
        dst.write(new_dsm_arr, 1)
    with rasterio.open(output_green_filename, 'w', driver='GTiff', height=new_green_arr.shape[0],
                       width=new_green_arr.shape[1], count=1, dtype=new_green_arr.dtype, crs='EPSG:31370', transform=t_poly['transform']) as dst:
        # Write the array data to the file
        dst.write(new_green_arr, 1)

    return output_grid_filename, output_dsm_filename, output_green_filename


def change_building(geometry_df, height, add, res, view_distance, year):
    landbouw_str = "_landbouw"
    if not landbouw:
        landbouw_str = ""
    blauw_str = "_blauw"
    if not blauw:
        blauw_str = ""
    dsm_path = rf"input_data/{year}/DSM_{res}m.tif"
    dtm_path = rf"input_data/{year}/DTM_{res}m.tif"
    green_path = rf"input_data/{year}/green_01{landbouw_str}{blauw_str}_{res}m.tif"
    grid = geometry_df.geometry.buffer(view_distance * 2.05 + res / 2).envelope

    # Convert the buffered polygon's exterior coordinates to pixel coordinates

    with rio_open(dtm_path) as dtm_input:
        with rio_open(dsm_path) as dsm_input:
            with rio_open(green_path) as green_input:
                dsm_arr, t_poly = extractPolygon(dsm_input, grid)
                dtm_arr, _ = extractPolygon(dtm_input, grid)
                green_arr, _ = extractPolygon(green_input, grid)
                building_arr = create_geometry_raster(geometry_df, dtm_arr.shape[0], dtm_arr.shape[1], t_poly['transform'])
    geometry_df['geometry'] = geometry_df.buffer(view_distance + res / 2).envelope
    new_dsm_arr = np.where(building_arr == 1, dtm_arr + height if add else dtm_arr, dsm_arr)
    new_green_arr = np.where(building_arr == 1, 0, green_arr)
    output_grid_filename = 'input_data/building_box.gpkg'
    output_dsm_filename = 'input_data/dsm_building.tif'
    output_green_filename = 'input_data/green_building.tif'
    geometry_df.to_file(output_grid_filename)
    with rasterio.open(output_dsm_filename, 'w', driver='GTiff', height=new_dsm_arr.shape[0],
                       width=new_dsm_arr.shape[1], count=1, dtype=new_dsm_arr.dtype, crs='EPSG:31370', transform=t_poly['transform']) as dst:
        # Write the array data to the file
        dst.write(new_dsm_arr, 1)
    with rasterio.open(output_green_filename, 'w', driver='GTiff', height=new_green_arr.shape[0],
                       width=new_green_arr.shape[1], count=1, dtype=new_green_arr.dtype, crs='EPSG:31370', transform=t_poly['transform']) as dst:
        # Write the array data to the file
        dst.write(new_green_arr, 1)
    return output_grid_filename, output_dsm_filename, output_green_filename


if __name__ == '__main__':
    x = 182482.4
    y = 207534.2
    size = 10
    height = 20
    res = 5
    view_distance = 800
    landbouw = True
    blauw = True
    year = 2015
    building_df = gpd.read_file('input_data/test_building.gpkg')
    building_df.to_crs(crs='epsg:31370', inplace=True)
    add = False
    # output_grid_filename, output_dsm_filename, output_green_filename = add_tree(x, y, height, size, res, view_distance, landbouw=landbouw, blauw=blauw, year=year)
    output_grid_filename, output_dsm_filename, output_green_filename = change_building(building_df, height, add, res, view_distance, year)
    run(res, view_distance, landbouw=landbouw, blauw=blauw, grid='grid_vl', total_parts=1, part_nr=0,
        dsm_path=output_dsm_filename, green_path=output_green_filename, grid_path=output_grid_filename, output_name_postfix="_building")
