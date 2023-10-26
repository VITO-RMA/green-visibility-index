from time import process_time_ns

from osgeo import gdal
import os
import uuid
import glob

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
import numpy as np
import rasterio
import geopandas as gpd
def get_point(geom):
    centroid = geom.centroid
    if geom.contains(centroid):
        return centroid
    else:
        return geom.representative_point()


def sample_raster(raster_ds, pts, band=1, value_col='raster_value'):
    """
    Providing a rasterio raster and numpy 2-D array of coordinates, returns
    a GeoDataFrame with the values at each point as an extra column. Assumes CRS
    of raster and coordinates are the same.
    Args:
        raster_ds: Rasterio raster
        pts: Numpy array of dimensions (n,2) or GeoPandas DataFrame
        band: The raster band number to be evaluated, defaults to 1
        value_col: Name of the column to store raster values in the GeoDataFrame
    Returns:
        GeoDataFrame: A new GeoDataFrame with the added column for raster values
    Example:
        r = rasterio.open("raster.tif")
        pts = np.genfromtxt("points.csv", delimiter=",", skip_header=1, usecols=(1,2))
        gdf = sampleRaster(r, pts)
    """
    ras = raster_ds.read(band)

    if isinstance(pts, gpd.GeoDataFrame):
        # Extract coordinates from the GeoDataFrame
        filtered_gdf = pts.cx[raster_ds.bounds[0]:raster_ds.bounds[2], raster_ds.bounds[1]:raster_ds.bounds[3]]
        filtered_g = filtered_gdf.geometry.apply(lambda geom: (geom.x, geom.y)).values
    inPts = np.array([np.array(p) for p in filtered_g])
    inPts_filtered = inPts[
        np.where(
            (inPts[:, 0] >= raster_ds.bounds[0])
            & (inPts[:, 0] < raster_ds.bounds[2])
            & (inPts[:, 1] >= raster_ds.bounds[1])
            & (inPts[:, 1] < raster_ds.bounds[3])
        )
    ]
    originX = raster_ds.bounds[0]
    originY = raster_ds.bounds[3]
    cellSize = raster_ds.transform[0]
    col = ((inPts_filtered[:, 0] - originX) / cellSize).astype(int)
    row = ((originY - inPts_filtered[:, 1]) / cellSize).astype(int)
    res = ras[row, col]

    # Create a new GeoDataFrame with the additional column
    new_gdf = gpd.GeoDataFrame(filtered_gdf)
    new_gdf[value_col] = res

    return new_gdf

def maak_wegen_punten():
    from shapely.geometry import LineString, Point
    from tqdm import tqdm
    line_data = gpd.read_file(r'\\vito.local\VITO\Unit_RMA\GIS\Flanders\Wegenregister\2022\Wegenregister_SHAPE_20220616_correctie\Shapefile\Wegsegment.shp', engine='pyogrio')
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
    # # Save the resulting points as a new shapefile
    point_data.to_file(r'C:\RMAbuild\Projects\green-visibility-index\input_data\2022\weg_points.gpkg', engine='pyogrio')
    point_data = gpd.read_file(r'C:\RMAbuild\Projects\green-visibility-index\input_data\2022\weg_points.gpkg', engine='pyogrio')
    # point_data = ru.get_zonal_statistics_for_dataframe(input_df=point_data,
    #                                                    input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw_lerc.tif",
    #                                                    prefix="green_vis_5m_800m_landbouw_", ignore_value=0)
    # point_data = sample_raster(rasterio.open(r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw_lerc.tif"), point_data, value_col='green_vis_5m_800m_landbouw_mean')
    point_data = sample_raster(rasterio.open(r"C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_2022_5m_800m_landbouw_blauw.tif"), point_data, value_col='green_vis_vl_2022_5m_800m_landbouw_blauw_mean')
    # point_data = sample_raster(rasterio.open(r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_lerc.tif"), point_data, value_col='green_vis_5m_800m_mean')
    # point_data = ru.get_zonal_statistics_for_dataframe(input_df=point_data,
    #                                                    input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw_blauw_lerc.tif",
    #                                                    prefix="green_vis_5m_800m_landbouw_", ignore_value=0)
    # keep_cols.append("green_vis_5m_800m_landbouw_mean")
    # point_data = point_data[keep_cols]
    # point_data = ru.get_zonal_statistics_for_dataframe(input_df=point_data,
    #                                                    input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_lerc.tif",
    #                                                    prefix="green_vis_5m_800m_", ignore_value=0)
    # keep_cols.append("green_vis_5m_800m_mean")
    # point_data = point_data[keep_cols]

    point_data.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points.gpkg", engine='pyogrio')
    ru.rasterize_shape_match_raster(r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points.gpkg",
                                    output_file=r'C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_landbouw_blauw_mean.tif',
                                    match_raster=r"C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_2022_5m_800m_landbouw_blauw.tif",
                                    attribute='green_vis_vl_2022_5m_800m_landbouw_blauw_mean')
    # ru.rasterize_shape_match_raster(r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points.gpkg",
    #                                 output_file=r'C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_landbouw_mean.tif',
    #                                 match_raster=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw_blauw_lerc.tif",
    #                                 attribute='green_vis_5m_800m_landbouw_mean')
    # ru.rasterize_shape_match_raster(r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points.gpkg",
    #                                 output_file=r'C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_mean.tif',
    #                                 match_raster=r"C:\RMAbuild\Projects\green-visibility-index\output\green_vis_vl_5m_800m_landbouw_blauw_lerc.tif",
    #                                 attribute='green_vis_5m_800m_mean')


def add_buurtgroen():
    buildings_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0_adp.gpkg", engine='pyogrio')
    buildings_df['unique_id'] = buildings_df.index
    buildings_edge_df = buildings_df.copy()
    buildings_edge_df.geometry = buildings_edge_df.geometry.buffer(0.1).buffer(400)
    buildings_edge_df['geometry'] = buildings_edge_df['geometry'].simplify(0.1)
    buildings_edge_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0_buurtgroen_buffers.gpkg", engine='pyogrio')
    buildings_edge_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_df,
                                                              input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_2022_5m_800m_landbouw_blauw.tif",
                                                              prefix="buurtgroen_5m_800m_landbouw_blauw_", ignore_value=0)
    buildings_df = pd.merge(buildings_df, buildings_edge_df[['unique_id', 'buurtgroen_5m_800m_landbouw_blauw_median']], on='unique_id')
    # buildings_edge_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_df,
    #                                                           input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_landbouw_mean.tif",
    #                                                           prefix="buurtgroen_5m_800m_landbouw_", ignore_value=0)
    # buildings_df = pd.merge(buildings_df, buildings_edge_df[['unique_id', 'buurtgroen_5m_800m_landbouw_mean', 'buurtgroen_5m_800m_landbouw_sum']], on='unique_id')
    # buildings_edge_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_df,
    #                                                           input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\wegen_points_mean.tif",
    #                                                           prefix="buurtgroen_5m_800m_", ignore_value=0)
    # r = pd.merge(buildings_df, buildings_edge_df[['unique_id', 'buurtgroen_5m_800m_mean', 'buurtgroen_5m_800m_sum']], on='unique_id')
    buildings_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0_buurtgroen.gpkg", engine='pyogrio')


def add_fields(clean=False, overlap=False):
    if clean:
        mask_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\grid_vl.gpkg", engine='pyogrio')
        buildings_df = gpd.read_file(r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Gbg.shp", bbox=mask_df, engine='pyogrio')
        buildings_df['idx'] = buildings_df.index
        buildings_df2 = buildings_df.copy()
        buildings_df2.geometry = buildings_df2.geometry.apply(get_point)
        other_df = gpd.read_file(r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["adp_CAPAKEY"] = joined_gdf['CAPAKEY']
        other_df = gpd.read_file(r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp_1.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["adp_CAPAKEY2"] = joined_gdf['CAPAKEY']
        buildings_df["adp_CAPAKEY"] = np.where(pd.isna(buildings_df['adp_CAPAKEY']), buildings_df['adp_CAPAKEY2'], buildings_df['adp_CAPAKEY'])
        buildings_df.drop(columns=['adp_CAPAKEY2'], inplace=True)
        print('adp_CAPAKEY')
        other_df = gpd.read_file(r"R:\Landgebruikskaart_2022\Indicatoren\KLV\2022_v3\percelen_classified_with_onbebouwd.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["bebouwingstype"] = joined_gdf['typologie']
        print('bebouwingstype')
        other_df = gpd.read_file(r"Y:\Unit_RMA\RDM\_Projecten\1910056 - Reftaak VPO\3 Werkdocumenten\Berekeningen\stedelijke_gebieden\Versie_oktober2020_RURA2.0\Versted_rand_land_vlaa_2019_versie2.shp", engine='pyogrio',
                                 bbox=mask_df)
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["verstedelijkt-randstelijk-landelijk"] = joined_gdf['type']
        buildings_df["NISCODE"] = joined_gdf['NISCODE']
        buildings_df["CODSEC"] = joined_gdf['CODSEC']
        print('CODSEC')
        other_df = gpd.read_file(r"R:\Landgebruikskaart_2022\postgis_output\methodeD\resident_percelen.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["residentieel"] = joined_gdf['resperceel']
        print('residentieel')
        other_df = gpd.read_file(r"R:\Landgebruikskaart_2022\postgis_output\methodeD\beb_percelen_voorzieningen_lu35.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["lu35_sector"] = joined_gdf['sector']
        print('lu35_sector')
        other_df = gpd.read_file(r"R:\Landgebruikskaart_2022\postgis_output\methodeD\vkbo_end_result_max_sector.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["vitosector"] = joined_gdf['vitosector']
        buildings_df["tot_tewerk"] = joined_gdf['tot_tewerk']
        print('vitosector')
        other_df = gpd.read_file(r"R:\Landgebruikskaart_2022\postgis_output\methodeD\beb_percelen_voorzieningen_lu33.shp", bbox=mask_df, engine='pyogrio')
        joined_gdf = gpd.sjoin(buildings_df2, other_df, how="left", op='intersects')
        joined_gdf = joined_gdf.groupby('idx').first()
        buildings_df["lu33_sector"] = joined_gdf['sector']
        print('lu33_sector')
        buildings_df['area'] = buildings_df.area

        buildings_df['unique_id'] = buildings_df.index
        buildings_df.geometry = buildings_df.buffer(0)  # makevalid
        buildings_df['gid'] = buildings_df.index
        buildings_df.to_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_joined.gpkg', engine='pyogrio')
    if overlap:
        buildings_df = gpd.read_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_joined.gpkg', engine='pyogrio')
        keep_cols = list(set(buildings_df.columns.values))
        buildings_edge_df = buildings_df.copy()
        buildings_edge_df.geometry = buildings_edge_df.buffer(1).boundary.simplify(0.1)
        buildings_edge_no_overlap_df = buildings_edge_df.overlay(buildings_df, how='difference')
        # add edges that have no overlap at all back ot dataframe
        # buildings_edge_no_overlap_df = buildings_edge_no_overlap_df.append(buildings_edge_df.loc[~buildings_edge_df['unique_id'].isin(buildings_edge_no_overlap_df['unique_id'].unique()), :])
        buildings_edge_no_overlap_df = pd.concat([buildings_edge_no_overlap_df, buildings_edge_df.loc[~buildings_edge_df['unique_id'].isin(buildings_edge_no_overlap_df['unique_id'].unique()), :]])
        buildings_edge_df = None
        buildings_edge_no_overlap_df['geometry'] = buildings_edge_no_overlap_df['geometry'].simplify(0.1)
        buildings_edge_no_overlap_df['len_tot'] = buildings_edge_no_overlap_df.geometry.length
        buildings_edge_no_overlap_df.reset_index(inplace=True, drop=True)
        buildings_edge_no_overlap_df.to_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_edge_no_overlap_df.gpkg', engine='pyogrio')
    buildings_edge_no_overlap_df = gpd.read_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_edge_no_overlap_df.gpkg', engine='pyogrio')
    keep_cols = list(set(buildings_edge_no_overlap_df.columns.values))

    buildings_edge_no_overlap_df = ru.get_zonal_statistics_for_dataframe(input_df=buildings_edge_no_overlap_df,
                                                                         input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_2022_5m_800m_landbouw_blauw.tif",
                                                                         prefix="green_vis_vl_2022_5m_800m_landbouw_blauw_", ignore_value=0)
    keep_cols.append("green_vis_vl_2022_5m_800m_landbouw_blauw_median")
    buildings_edge_no_overlap_df = buildings_edge_no_overlap_df[keep_cols]
    buildings_edge_no_overlap_df['geometry'] = buildings_edge_no_overlap_df['geometry'].simplify(0.1)
    buildings_edge_no_overlap_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0.gpkg", engine='pyogrio')
    buildings_df = gpd.read_file(r'C:\RMAbuild\Projects\green-visibility-index\output\buildings_joined.gpkg', engine='pyogrio')
    buildings_df['gid'] = buildings_df.index
    buildings_df['geometry'] = buildings_df['geometry'].simplify(0.1)
    buildings_df2 = buildings_df.copy()
    buildings_df2.geometry = buildings_df2.geometry.apply(get_point)
    dfs = []

    # Read each shapefile and append to the list
    for shapefile in [r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp.shp", r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp_1.shp"]:
        df = gpd.read_file(shapefile, engine='pyogrio')
        dfs.append(df)

    # Concatenate the individual GeoDataFrames into one
    adp_df = gpd.GeoDataFrame(pd.concat(dfs, ignore_index=True))
    buildings_df.drop(columns=['adp_CAPAKEY'],inplace=True)
    joined_gdf = gpd.sjoin(buildings_df2, adp_df, how="inner", op='intersects')
    buildings_df = pd.merge(buildings_edge_no_overlap_df, joined_gdf[['gid','CAPAKEY']], on='gid')
    buildings_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0_adp.gpkg", engine='pyogrio')


def make_adp():
    dfs = []
    for shapefile in [r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp.shp", r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp_1.shp"]:
        df = gpd.read_file(shapefile, engine='pyogrio')
        dfs.append(df)
    adp_df = gpd.GeoDataFrame(pd.concat(dfs, ignore_index=True))
    adp_df = ru.get_zonal_statistics_for_dataframe(input_df=adp_df,
                                                   input_raster_path=r"C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_2022_5m_800m_landbouw_blauw.tif",
                                                   prefix="green_vis_vl_2022_5m_800m_landbouw_blauw_", ignore_value=0)
    adp_df.drop(columns=['index'], inplace=True)

    adp_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\testgebied\adp_green_vis_zonder_0.gpkg", engine='pyogrio')

def add_inwoners():
    # inw_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\input_data\2022\inwaant_clean.gpkg", engine='pyogrio')
    # dfs=[]
    # print('.')
    # for shapefile in [r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp.shp", r"\\vito.local\VITO\Unit_RMA\GIS\Flanders\GRB\2022\GRB_20220919_Shapefile\Shapefile\Adp_1.shp"]:
    #     df = gpd.read_file(shapefile, engine='pyogrio')
    #     dfs.append(df)
    # adp_df = gpd.GeoDataFrame(pd.concat(dfs, ignore_index=True))
    # print('.')
    # adp_join_df = gpd.sjoin(inw_df, adp_df[['geometry','CAPAKEY']], how='left', op='within')
    # adp_join_df.drop(columns=['index_left', 'index_right'], errors='ignore', inplace=True)
    # buildings_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\buildings_joined.gpkg", engine='pyogrio')
    #
    # bld_join_df = gpd.sjoin(adp_join_df, buildings_df[['geometry','UIDN']], how='left', op='within')
    # bld_join_df.drop(columns=['index_left', 'index_right'], errors='ignore', inplace=True)
    #
    # bld_join_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\inwaant_clean_capakey_bld_id.gpkg", engine='pyogrio')
    adp_join_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\inwaant_clean_capakey_bld_id.gpkg", engine='pyogrio')
    print('..')
    buildings_df = gpd.read_file(r"C:\RMAbuild\Projects\green-visibility-index\output\gbg_green_vis_zonder_0_buurtgroen.gpkg", engine='pyogrio')
    buildings_df.drop(columns=['geometry'],inplace=True)
    print('.')
    # Sort the DataFrame
    df_sorted = buildings_df.sort_values(['CAPAKEY', 'TYPE', 'LENGTE'], ascending=[True, True, False])

    # Group by 'CAPAKEY' and aggregate
    result_df = df_sorted.groupby('CAPAKEY').agg({
        'LENGTE': 'first',
        'TYPE': 'first',
        'green_vis_vl_2022_5m_800m_landbouw_blauw_median': 'mean',
        'buurtgroen_5m_800m_landbouw_blauw_median': 'mean'
    }).reset_index()

    print('..')
    # Finally, merge with adp_join_df
    inw_merge_df = adp_join_df.merge(result_df, on='CAPAKEY', how='left')

    print('..')
    # print(inw_merge_df.columns)
    # Finally, merge with adp_join_df
    # inw_merge_df = inw_merge_df.merge(result_df, on='UIDN', how='left')
    print('.')
    # Convert back to GeoDataFrame if needed and set geometry column

    # Identify rows where no match was found (i.e., some columns are NaN)
    no_match_df = inw_merge_df[inw_merge_df['buurtgroen_5m_800m_landbouw_blauw_median'].isna()]
    #
    df_sorted = buildings_df.sort_values(['UIDN', 'TYPE', 'LENGTE'], ascending=[True, True, False])

    # Group by 'CAPAKEY' and aggregate
    result_df = df_sorted.groupby('UIDN').agg({
        'LENGTE': 'first',
        'TYPE': 'first',
        'green_vis_vl_2022_5m_800m_landbouw_blauw_median': 'mean',
        'buurtgroen_5m_800m_landbouw_blauw_median': 'mean'
    }).reset_index()
    no_match_df.drop(columns=['LENGTE','TYPE','green_vis_vl_2022_5m_800m_landbouw_blauw_median','buurtgroen_5m_800m_landbouw_blauw_median'], inplace=True)
    # # Perform another merge on those rows with a different field
    no_match_merged_df = no_match_df.merge(result_df, on='UIDN', how='left')
    #
    # # Drop the rows where the first merge didn't find a match
    inw_merge_df.dropna(subset=['buurtgroen_5m_800m_landbouw_blauw_median'], inplace=True)

    # Concatenate the two results
    final_df = pd.concat([inw_merge_df, no_match_merged_df])
    inw_merge_df = gpd.GeoDataFrame(final_df, geometry='geometry')
    inw_merge_df.to_file(r"C:\RMAbuild\Projects\green-visibility-index\output\inwaant_clean_capakey_buurtgroen_mean_hfdg.gpkg", engine='pyogrio')


def combine_outputs():
    blauw = ""
    landbouw = ""
    # blauw = "blauw_"
    # landbouw = "landbouw_"
    tif_list = []
    for f in glob.glob(rf'C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_5m_800m_landbouw_blauw_*grid_vl.tif'):
        tif_list.append(f)
    filename_vrt = str(uuid.uuid4())
    vrt_options = gdal.BuildVRTOptions(srcNodata=0, VRTNodata=0)

    vrt = gdal.BuildVRT(rf"c:/temp/{filename_vrt}.vrt", tif_list, options=vrt_options)
    translateoptions = gdal.TranslateOptions(gdal.ParseCommandLine("-of COG -co COMPRESS=LERC_DEFLATE -co MAX_Z_ERROR=0.001 -co BIGTIFF=YES -a_nodata 0"))
    gdal.Translate(rf'C:\RMAbuild\Projects\green-visibility-index\output\2022\green_vis_vl_5m_800m_landbouw_blauw_2022.tif', vrt, options=translateoptions)
    vrt = None
    os.remove(f"c:/temp/{filename_vrt}.vrt")
if __name__ == '__main__':
    # make_adp()
    # add_fields(False,False)
    # maak_wegen_punten()
    # add_buurtgroen()
    # combine_outputs()
    add_inwoners()
