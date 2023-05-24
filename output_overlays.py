import glob
import os
import os

os.environ['USE_PYGEOS'] = '0'
import pandas as pd
from pysdx import raster_utils as ru
import numpy as np
import geopandas as gpd
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


# for f in glob.glob(r"Y:\Unit_RMA\GIS\Flanders\DSM\DSM_5m\e\*\*\*.tif"):
#     print(f)
#     shutil.move(f,os.path.join(r"Y:\Unit_RMA\GIS\Flanders\DSM\DSM_5m\GeoTiff",os.path.basename(f)))

def funct(in_arr):
    return np.where(in_arr < 3, 1, 0)


def maak_wegen_punten():
    import geopandas as gpd
    from shapely.geometry import LineString, Point
    from tqdm import tqdm
    line_data = gpd.read_file('C:\RMAbuild\Projects\green-visibility-index\input_data\Wegsegment.shp')
    print(line_data.columns)
    point_data = []

    interval_distance = 5
    for index, line in tqdm(line_data.iterrows(), total=len(line_data)):
        # Extract additional column values from the line input
        column1_val = line['LSTRNMID']
        column2_val = line['LSTRNM']
        column3_val = line['WEGCAT']
        column4_val = line['MORF']

        # Calculate the length of the line
        line_length = line.geometry.length

        # Generate points at regular intervals along the line
        for distance in range(0, int(line_length), interval_distance):
            # Get the point geometry at the specified distance
            point = line.geometry.interpolate(distance)

            # Create a GeoSeries from the point geometry and append it to the point_data GeoDataFrame
            point_data.append({'geometry': point, 'LSTRNMID': column1_val, 'LSTRNM': column2_val, 'WEGCAT': column3_val, 'MORF': column4_val})
    point_data = gpd.GeoDataFrame(point_data, crs=line_data.crs)
    # Save the resulting points as a new shapefile
    point_data.to_file('C:\RMAbuild\Projects\green-visibility-index\input_data\weg_points.gpkg')
    import rasterio
    from rasterio.features import geometry_mask

    keep_cols = list(set(point_data.columns.values))
    point_data = ru.get_zonal_statistics_for_dataframe(input_df=point_data,
                                                       input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw.tif",
                                                       prefix="green_vis_5m_800m_landbouw_", ignore_value=0)
    keep_cols.append("green_vis_5m_800m_landbouw_mean")
    point_data = point_data[keep_cols]
    point_data = ru.get_zonal_statistics_for_dataframe(input_df=point_data,
                                                       input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m.tif",
                                                       prefix="green_vis_5m_800m_", ignore_value=0)
    keep_cols.append("green_vis_5m_800m_mean")
    point_data = point_data[keep_cols]

    point_data.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points.gpkg")


def add_buurtgroen():
    global buildings_edge_df
    buildings_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0.gpkg", engine='pyogrio')
    buildings_df['unique_id'] = buildings_df.index
    buildings_edge_df = buildings_df.copy()
    buildings_edge_df.geometry = buildings_edge_df.buffer(400)
    buildings_edge_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_df,
                                                              input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_landbouw_mean.tif",
                                                              prefix="buurtgroen_5m_800m_landbouw_", ignore_value=0)
    buildings_edge_df = pd.merge(buildings_df, buildings_edge_df[['unique_id', 'buurtgroen_5m_800m_landbouw_mean', 'buurtgroen_5m_800m_landbouw_sum']], on='unique_id')
    buildings_edge_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_df,
                                                              input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_mean.tif",
                                                              prefix="buurtgroen_5m_800m_", ignore_value=0)
    r = pd.merge(buildings_df, buildings_edge_df[['unique_id', 'buurtgroen_5m_800m_mean', 'buurtgroen_5m_800m_sum']], on='unique_id')
    r.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0_buurtgroen.gpkg", engine='pyogrio')


def add_fields():
    global buildings_df2
    mask_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\grid_vl.gpkg")
    buildings_df = gpd.read_file(r"Y:\Unit_RMA\GIS\Flanders\GRB\2017\shapefile_20170109\Shapefile\Gbg.shp", mask=mask_df)
    buildings_df['idx'] = buildings_df.index
    buildings_df2 = buildings_df.copy()
    buildings_df2.geometry = buildings_df2.geometry.centroid
    other_df = gpd.read_file(r"Y:\Unit_RMA\GIS\Flanders\GRB\2017\shapefile_20170109\Shapefile\Adp.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["adp_CAPAKEY"] = joined_gdf['CAPAKEY']
    other_df = gpd.read_file(r"Y:\Unit_RMA\GIS\Flanders\GRB\2017\shapefile_20170109\Shapefile\Adp_1.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["adp_CAPAKEY2"] = joined_gdf['CAPAKEY']
    buildings_df["adp_CAPAKEY"] = np.where(pd.isna(buildings_df['adp_CAPAKEY']), buildings_df['adp_CAPAKEY'], buildings_df['adp_CAPAKEY2'])
    buildings_df.drop(columns=['adp_CAPAKEY2'], inplace=True)
    print('adp_CAPAKEY')
    other_df = gpd.read_file(r"R:\Landgebruikskaart_2016\Indicatoren\KLV\percelen_classified_with_onbebouwd.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["bebouwingstype"] = joined_gdf['TYPOLOGIE']
    print('bebouwingstype')
    other_df = gpd.read_file(r"Y:\Unit_RMA\RDM\_Projecten\1910056 - Reftaak VPO\3 Werkdocumenten\Berekeningen\stedelijke_gebieden\Versie_oktober2020_RURA2.0\Versted_rand_land_vlaa_2016_versie2.shp",
                             mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["verstedelijkt-randstelijk-landelijk"] = joined_gdf['type']
    buildings_df["NISCODE"] = joined_gdf['NISCODE']
    buildings_df["CODSEC"] = joined_gdf['CODSEC']
    print('CODSEC')
    other_df = gpd.read_file(r"R:\Landgebruikskaart_2016\Werkplan_2018\Verwerking\Verstedelijking\residentiële_percelen\residentiële_percelen_v2.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["residentieel"] = joined_gdf['resPerc']
    print('residentieel')
    other_df = gpd.read_file(r"R:\Landgebruikskaart_2016\postgis_output\beb_percelen_voorzieningen_lu35.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["lu35_sector"] = joined_gdf['sector']
    print('lu35_sector')
    other_df = gpd.read_file(r"R:\Landgebruikskaart_2016\postgis_output\methodeD\vkbo_end_result_max_sector.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["vitosector"] = joined_gdf['vitosector']
    buildings_df["tot_tewerk"] = joined_gdf['tot_tewerk']
    print('vitosector')
    other_df = gpd.read_file(r"R:\Landgebruikskaart_2016\postgis_output\beb_percelen_voorzieningen_lu33.shp", mask=mask_df)
    joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
    joined_gdf = joined_gdf.groupby('idx').first()
    buildings_df["lu33_sector"] = joined_gdf['sector']
    print('lu33_sector')
    buildings_df['area'] = buildings_df.area
    keep_cols = list(set(buildings_df.columns.values))
    buildings_df['unique_id'] = buildings_df.index
    buildings_df.geometry = buildings_df.buffer(0)  # makevalid
    buildings_df.to_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_joined.gpkg')
    buildings_df = gpd.read_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_joined.gpkg')
    keep_cols = list(set(buildings_df.columns.values))
    buildings_edge_df = buildings_df.copy()
    buildings_edge_df.geometry = buildings_edge_df.buffer(1).boundary
    buildings_edge_no_overlap_df = buildings_edge_df.overlay(buildings_df, how='difference')
    # add edges that have no overlap at all back ot dataframe
    # buildings_edge_no_overlap_df = buildings_edge_no_overlap_df.append(buildings_edge_df.loc[~buildings_edge_df['unique_id'].isin(buildings_edge_no_overlap_df['unique_id'].unique()), :])
    buildings_edge_no_overlap_df = pd.concat([buildings_edge_no_overlap_df, buildings_edge_df.loc[~buildings_edge_df['unique_id'].isin(buildings_edge_no_overlap_df['unique_id'].unique()), :]])
    buildings_edge_df = None
    buildings_edge_no_overlap_df['len_tot'] = buildings_edge_no_overlap_df.geometry.length
    buildings_edge_no_overlap_df.reset_index(inplace=True, drop=True)
    # buildings_edge_no_overlap_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_no_overlap_df,
    #                                                                      input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_1m_800m_landbouw.tif",
    #                                                                      prefix="green_vis_1m_800m_landbouw_", ignore_value=0)
    # keep_cols.append("green_vis_1m_800m_landbouw_mean")
    # buildings_edge_no_overlap_df = buildings_edge_no_overlap_df[keep_cols]
    buildings_edge_no_overlap_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_no_overlap_df,
                                                                         input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw.tif",
                                                                         prefix="green_vis_5m_800m_landbouw_", ignore_value=0)
    keep_cols.append("green_vis_5m_800m_landbouw_mean")
    buildings_edge_no_overlap_df = buildings_edge_no_overlap_df[keep_cols]
    buildings_edge_no_overlap_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_no_overlap_df,
                                                                         input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m.tif",
                                                                         prefix="green_vis_5m_800m_", ignore_value=0)
    keep_cols.append("green_vis_5m_800m_mean")
    buildings_edge_no_overlap_df = buildings_edge_no_overlap_df[keep_cols]
    buildings_edge_no_overlap_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_no_overlap_df,
                                                                         input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_5m_200m_landbouw.tif",
                                                                         prefix="green_vis_5m_200m_landbouw_", ignore_value=0)
    keep_cols.append("green_vis_5m_200m_landbouw_mean")
    buildings_edge_no_overlap_df = buildings_edge_no_overlap_df[keep_cols]
    buildings_edge_no_overlap_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_no_overlap_df,
                                                                         input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_1m_200m_landbouw.tif",
                                                                         prefix="green_vis_1m_200m_landbouw_", ignore_value=0)
    keep_cols.append("green_vis_1m_200m_landbouw_mean")
    buildings_edge_no_overlap_df = buildings_edge_no_overlap_df[keep_cols]
    print(buildings_edge_no_overlap_df.columns)
    buildings_edge_no_overlap_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0.gpkg")
    buildings_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\gbg_green_vis_zonder_0.gpkg")
    buildings_df2 = buildings_df.copy()
    buildings_df2.geometry = buildings_df2.geometry.centroid
    adp_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\adp_green_vis_zonder_0.gpkg")
    joined_gdf = gpd.sjoin(buildings_df2, adp_df, how="inner", op='intersects')
    buildings_df["adp_CAPAKEY"] = joined_gdf['CAPAKEY']
    buildings_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\gbg_green_vis_zonder_0_adp.gpkg")


def make_adp():
    adp_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\adp.gpkg")
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_1m_200m_landbouw.tif",
                                                   prefix="green_vis_1m_200m_landbouw_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_1m_800m_landbouw.tif",
                                                   prefix="green_vis_1m_800m_landbouw_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_5m_200m_landbouw.tif",
                                                   prefix="green_vis_5m_200m_landbouw_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_5m_800m.tif",
                                                   prefix="green_vis_5m_800m_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_5m_800m_landbouw.tif",
                                                   prefix="green_vis_5m_800m_landbouw_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\green_vis_5m_800m_landbouw_10xweight.tif",
                                                   prefix="green_vis_5m_800m_landbouw_10xweight_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)
    adp_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\adp_green_vis_zonder_0.gpkg")


def cut_grid_with_vl():
    global vl_df
    grid_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\grid_all.gpkg")
    vl_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\vlaanderen.gpkg")
    vl_df.crs = grid_df.crs
    grid = gpd.sjoin(grid_df, vl_df, how="inner", op="intersects")
    grid.to_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\grid_vl.gpkg")


def make_green():
    ru.process_raster_with_function_in_chunks(input_file=r"C:\RMAbuild\Projects\green-visibility-index\input_data\total.tif",
                                              output_file=r"C:\RMAbuild\Projects\green-visibility-index\input_data\green_01_1m.tif",
                                              like_file=r"C:\RMAbuild\Projects\green-visibility-index\input_data\total.tif",
                                              function_to_apply=funct,
                                              dtype='uint8')


if __name__ == '__main__':
    # make_green()
    # cut_grid_with_vl()
    # make_adp()
    # add_fields()
    # maak_wegen_punten()
    add_buurtgroen()
