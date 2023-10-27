import rasterio
from pysdx import raster_utils as ru
import numpy as np
import os

os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
from osgeo import gdal
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
ru.show_progressbar(True)

def funct(in_arr, limit):
    return np.where((in_arr <= limit) & (in_arr != 0), 1, 0)


def cut_grid_with_vl():
    global vl_df
    grid_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\grid_all.gpkg")
    vl_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\vlaanderen.gpkg")
    vl_df.crs = grid_df.crs
    grid = gpd.sjoin(grid_df, vl_df, how="inner", op="intersects")
    grid.to_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\grid_vl.gpkg")


def filter_green(input_arr, serres_arr):
    input_arr = np.where(serres_arr == 1, 4, input_arr)
    return input_arr


def filter_green_blauw(input_arr, serres_arr, water_arr, zee_arr):
    input_arr = np.where(serres_arr == 1, 4, input_arr)
    input_arr = np.where(water_arr == 1, 2, input_arr)
    input_arr = np.where(zee_arr == 1, 2, input_arr)
    return input_arr

import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
import fiona



def make_green_binary_maps():
    res = 1
    year= 2022

    ru.rasterize_shape_match_raster(input_file=rf'C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\serres.shp',
                                    output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\serres.tif",
                                    match_raster=fr"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif")

    ru.rasterize_shape_match_raster(rf'C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\Wtz.shp',
                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\water.tif",
                                    fr"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif")

    ru.rasterize_shape_match_raster(input_file=rf'C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\bcs.shp',
                                    output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\zee.tif",
                                    match_raster=fr"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif")
    ru.process_raster_list_with_function_in_chunks(input_file_list=[fr"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                                                    fr"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\serres.tif"],
                                                   output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_filter_{res}m.tif",
                                                   like_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                                   function_to_apply=filter_green)
    ru.process_raster_list_with_function_in_chunks(input_file_list=[rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\serres.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\water.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\zee.tif", ],
                                                   output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_blauw_filter_{res}m.tif",
                                                   like_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                                   function_to_apply=filter_green_blauw)

    ru.process_raster_with_function_in_chunks(input_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_filter_{res}m.tif",
                                              output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_01_{res}m.tif",
                                              like_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                              function_to_apply=funct,
                                              function_arguments={'limit': 2},
                                              dtype='uint8')

    ru.process_raster_with_function_in_chunks(input_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_filter_{res}m.tif",
                                              output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_01_landbouw_{res}m.tif",
                                              like_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                              function_to_apply=funct,
                                              function_arguments={'limit': 3},
                                              dtype='uint8')

    ru.process_raster_with_function_in_chunks(input_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_blauw_filter_{res}m.tif",
                                              output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\green_01_landbouw_blauw_{res}m.tif",
                                              like_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\{year}\total_{res}m.tif",
                                              function_to_apply=funct,
                                              function_arguments={'limit': 3},
                                              dtype='uint8')


def cleanup_green():
    make_green_binary_maps()

def fill_water(dsm_arr,dtm_arr,water_arr,zee_arr, strand_arr):
    return np.where((water_arr == 1) | (zee_arr == 1) | (strand_arr == 1), dtm_arr, dsm_arr)
def fix_dsm():
    # ru.convert_raster_to_likeraster(rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\zee.tif",
    #                                 rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2015\DTM_1m.tif",
    #                                 rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\zee_fix.tif")
    # ru.rasterize_shape_match_raster(rf'R:\Landgebruikskaart_2022\postgis_output\methodeD\BWK_natuur.shp',
    #                                 rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\strand_fix.tif",
    #                                 fr"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\zee_fix.tif",
    #                                 where='''"natcat" = 'strand' ''')

    ru.process_raster_list_with_function_in_chunks(input_file_list=[rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_1m.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2015\DSM_1m.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\water_fix.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\zee_fix.tif",
                                                                    rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\strand_fix.tif",
                                                                    ],
                                              output_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_1m_s.tif",
                                              like_file=rf"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_1m.tif",
                                              function_to_apply=fill_water, chunks = 40)

def to_5m():
    translateoptions = gdal.TranslateOptions(gdal.ParseCommandLine("-tr 5 5 -r average -of Gtiff -co COMPRESS=DEFLATE -co TILED=YES -co BIGTIFF=YES -co NUM_THREADS=ALL_CPUS"))
    gdal.Translate(fr'C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_5m_s.tif',
                   r'C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_1m_s.tif', options=translateoptions)
if __name__ == '__main__':
    cleanup_green()
    fix_dsm()
    to_5m()
    # translateoptions = gdal.TranslateOptions(gdal.ParseCommandLine("-of COG -co COMPRESS=LERC_DEFLATE -co MAX_Z_ERROR=0.01 -co BIGTIFF=YES -a_nodata -9999 -co NUM_THREADS=ALL_CPUS"))
    # gdal.Translate(rf'C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_1m_lerc.tif', r'C:\RMAbuild\Projects\green-visibility-index\input_data\2022\DSM_1m.tif', options=translateoptions)