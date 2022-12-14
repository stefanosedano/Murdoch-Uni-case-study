import urllib
from netCDF4 import Dataset, num2date
from bs4 import BeautifulSoup
import urllib.request
import os
import numpy as np
import pandas as pd

from osgeo import gdal

import sys
sys.path.append('..')
from General.geo import get_GADM, writeGeoTiff_v3, clip_to_extent, get_country_mask, submatnanmean, get_list_of_countries,get_country_mask_projected,clip_to_extent_resolution,clip_to_extent_resolution_inproj_outproj

from osgeo import osr, ogr, gdal
from operator import itemgetter, attrgetter


GHS_POP = gdal.Open("C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_POP/{}/index.vrt".format(2000))

res = GHS_POP.GetGeoTransform()[1]
wkt = GHS_POP.GetProjection()
epsg = "ESRI:54009"

minx, miny, maxx, maxy, provincemask, new_geotrans, names_code= get_country_mask_projected("PAK",res,epsg,wkt)
extent = minx, miny, maxx, maxy
provincemask = np.array(provincemask.ReadAsArray(), dtype=np.int32)

SMOD_codes = [
[30,"Urban Centre grid cell"],
[23,"Dense Urban Cluster grid cell"],
[22,"Semi-dense Urban Cluster grid cell"],
[21,"Suburban or per-urban grid cell"],
[13,"Rural cluster grid cell"],
[12,"Low Density Rural grid cell"],
[11,"Very low density rural grid cell"],
[10,"Water grid cell"],
]
SMOD_codes = pd.DataFrame(SMOD_codes,columns=["CODE","NAME"])


GHS_YEARS = [2000,2005,2010,2015,2020,2025,2030]
lis_of_values_per_smod_class_per_countries = []

smod_classes = [10, 11, 12, 13, 21, 22, 23, 30]
country_ids = np.unique(provincemask).tolist()
country_ids.remove(0)

if not os.path.exists("pakistan_stats_pop_smod.csv"):
    for year in GHS_YEARS:
        print("processing {}".format(year))

        ## surface stats
        GHS_SMOD = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_SMOD/{}/index.vrt".format(year))

        GHS_SMOD = clip_to_extent_resolution(GHS_SMOD, extent,res)
        GHS_SMOD[provincemask==0] = 0

        class_base_layer = (GHS_SMOD + provincemask)
        class_base_layer[provincemask == 0] = 0

        for country in country_ids:
            lis_of_values_per_smod_class = []
            lis_of_values_per_smod_class.append(names_code.loc[names_code["CODE"]==country,"NAME_1"].values[0])
            lis_of_values_per_smod_class.append(year)
            lis_of_values_per_smod_class.append("num_pixel_in_class")
            for smod_class in smod_classes:
                id_search = country+smod_class
                num_pixel = np.sum(class_base_layer==id_search)
                if num_pixel < 0:
                    num_pixel = 0
                lis_of_values_per_smod_class.append(num_pixel)
            lis_of_values_per_smod_class_per_countries.append(lis_of_values_per_smod_class)

        ## pop stats
        GHS_POP = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_POP/{}/index.vrt".format(year))
        GHS_POP = clip_to_extent_resolution(GHS_POP, extent, res)
        GHS_POP[provincemask == 0] = 0

        for country in country_ids:
            lis_of_values_per_smod_class = []
            lis_of_values_per_smod_class.append(names_code.loc[names_code["CODE"]==country,"NAME_1"].values[0])
            lis_of_values_per_smod_class.append(year)
            lis_of_values_per_smod_class.append("pop_in_class")
            for smod_class in smod_classes:
                id_search = country+smod_class
                pop_pixel = np.sum(GHS_POP[class_base_layer==id_search])
                if pop_pixel < 0:
                    pop_pixel = 0
                lis_of_values_per_smod_class.append(pop_pixel)
            lis_of_values_per_smod_class_per_countries.append(lis_of_values_per_smod_class)

    columns=[
        "country",
        "year",
        "variable",
        SMOD_codes.loc[SMOD_codes["CODE"]==10,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==11,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==12,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==13,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==21,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==22,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==23,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==30,"NAME"].values[0]
    ]
    out=pd.DataFrame(lis_of_values_per_smod_class_per_countries,columns=columns)
    out.to_csv("pakistan_stats_pop_smod.csv")


if not os.path.exists("pakistan_stats_flood.csv"):
    lis_of_values_per_smod_class_per_countries = []
    for year in range(2003,2019):
        print("processing flood {}".format(year))
        arr = np.asarray(GHS_YEARS)
        i = (np.abs(arr - year)).argmin()
        closest_year = arr[i]

        ## flood stats
        FLOOD = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GLOBAL_FLOOD_DB/{}/index.vrt".format(year))
        FLOOD_ = clip_to_extent_resolution_inproj_outproj(FLOOD, extent, res,"EPSG:4326",epsg)
        FLOOD = FLOOD_
        FLOOD = FLOOD>0
        FLOOD[provincemask == 0] = 0

        ## surface
        GHS_SMOD = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_SMOD/{}/index.vrt".format(closest_year))

        GHS_SMOD = clip_to_extent_resolution(GHS_SMOD, extent,res)
        GHS_SMOD[provincemask==0] = 0

        class_base_layer = (GHS_SMOD + provincemask)
        class_base_layer[provincemask == 0] = 0
        class_base_layer[FLOOD == 0] = 0

        ## pop
        GHS_POP = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_POP/{}/index.vrt".format(closest_year))
        GHS_POP = clip_to_extent_resolution(GHS_POP, extent, res)
        GHS_POP[provincemask == 0] = 0
        #GHS_POP[FLOOD == 0] = 0
        GHS_POP = GHS_POP * FLOOD_


        for country in country_ids:
            lis_of_values_per_smod_class = []
            lis_of_values_per_smod_class.append(names_code.loc[names_code["CODE"]==country,"NAME_1"].values[0])
            lis_of_values_per_smod_class.append(year)
            lis_of_values_per_smod_class.append("num_pixel_flood_in_class")
            for smod_class in smod_classes:
                id_search = country+smod_class
                #num_pixel = np.sum(class_base_layer==id_search)

                num_pixel = np.sum(FLOOD_ * (class_base_layer==id_search))
                if num_pixel < 0:
                    num_pixel = 0
                lis_of_values_per_smod_class.append(num_pixel)
            lis_of_values_per_smod_class_per_countries.append(lis_of_values_per_smod_class)


        for country in country_ids:
            lis_of_values_per_smod_class = []
            lis_of_values_per_smod_class.append(names_code.loc[names_code["CODE"]==country,"NAME_1"].values[0])
            lis_of_values_per_smod_class.append(year)
            lis_of_values_per_smod_class.append("pop_flood_in_class")
            for smod_class in smod_classes:
                id_search = country+smod_class
                pop_pixel = np.sum(GHS_POP[class_base_layer==id_search])
                if pop_pixel < 0:
                    pop_pixel = 0
                lis_of_values_per_smod_class.append(pop_pixel)
            lis_of_values_per_smod_class_per_countries.append(lis_of_values_per_smod_class)


    columns=[
        "country",
        "year",
        "variable",
        SMOD_codes.loc[SMOD_codes["CODE"]==10,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==11,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==12,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==13,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==21,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==22,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==23,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==30,"NAME"].values[0]
    ]
    out=pd.DataFrame(lis_of_values_per_smod_class_per_countries,columns=columns)
    out.to_csv("pakistan_stats_flood.csv")




if not os.path.exists("nightlight_CCNL.csv"):
    lis_of_values_per_smod_class_per_countries = []
    for year in range(2000,2014):
        print("processing night {}".format(year))
        arr = np.asarray(GHS_YEARS)
        i = (np.abs(arr - year)).argmin()
        closest_year = arr[i]


        NIGHT = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/nightlight_CCNL/{}/index.vrt".format(year))
        NIGHT = clip_to_extent_resolution_inproj_outproj(NIGHT, extent, res,"EPSG:4326",epsg)
        NIGHT[provincemask == 0] = 0
        NIGHT[NIGHT< 0] = 0

        ## surface
        GHS_SMOD = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_SMOD/{}/index.vrt".format(closest_year))

        GHS_SMOD = clip_to_extent_resolution(GHS_SMOD, extent,res)
        GHS_SMOD[provincemask==0] = 0

        class_base_layer = (GHS_SMOD + provincemask)
        class_base_layer[provincemask == 0] = 0
        class_base_layer[NIGHT == 0] = 0

        for country in country_ids:
            lis_of_values_per_smod_class = []
            lis_of_values_per_smod_class.append(names_code.loc[names_code["CODE"]==country,"NAME_1"].values[0])
            lis_of_values_per_smod_class.append(year)
            lis_of_values_per_smod_class.append("nightlight_CCNL_in_class")
            for smod_class in smod_classes:
                id_search = country+smod_class
                pop_pixel = np.mean(NIGHT[class_base_layer==id_search])
                if pop_pixel < 0:
                    pop_pixel = 0
                lis_of_values_per_smod_class.append(pop_pixel)
            lis_of_values_per_smod_class_per_countries.append(lis_of_values_per_smod_class)


    columns=[
        "country",
        "year",
        "variable",
        SMOD_codes.loc[SMOD_codes["CODE"]==10,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==11,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==12,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==13,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==21,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==22,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==23,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==30,"NAME"].values[0]
    ]
    out=pd.DataFrame(lis_of_values_per_smod_class_per_countries,columns=columns)
    out.to_csv("nightlight_CCNL.csv")


if not os.path.exists("nightlight_VCMCFG.csv"):
    lis_of_values_per_smod_class_per_countries = []
    for year in range(2012,2023):
        print("processing night {}".format(year))
        arr = np.asarray(GHS_YEARS)
        i = (np.abs(arr - year)).argmin()
        closest_year = arr[i]


        NIGHT = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/nightlight_VCMCFG/{}/index.vrt".format(year))
        NIGHT = clip_to_extent_resolution_inproj_outproj(NIGHT, extent, res,"EPSG:4326",epsg)
        NIGHT[provincemask == 0] = 0
        NIGHT[NIGHT< 0] = 0

        ## surface
        GHS_SMOD = gdal.Open(
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GHS_SMOD/{}/index.vrt".format(closest_year))

        GHS_SMOD = clip_to_extent_resolution(GHS_SMOD, extent,res)
        GHS_SMOD[provincemask==0] = 0

        class_base_layer = (GHS_SMOD + provincemask)
        class_base_layer[provincemask == 0] = 0
        class_base_layer[NIGHT == 0] = 0

        for country in country_ids:
            lis_of_values_per_smod_class = []
            lis_of_values_per_smod_class.append(names_code.loc[names_code["CODE"]==country,"NAME_1"].values[0])
            lis_of_values_per_smod_class.append(year)
            lis_of_values_per_smod_class.append("nightlight_VCMCFG_in_class")
            for smod_class in smod_classes:
                id_search = country+smod_class
                pop_pixel = np.mean(NIGHT[class_base_layer==id_search])
                if pop_pixel < 0:
                    pop_pixel = 0
                lis_of_values_per_smod_class.append(pop_pixel)
            lis_of_values_per_smod_class_per_countries.append(lis_of_values_per_smod_class)


    columns=[
        "country",
        "year",
        "variable",
        SMOD_codes.loc[SMOD_codes["CODE"]==10,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==11,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==12,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==13,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==21,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==22,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==23,"NAME"].values[0],
        SMOD_codes.loc[SMOD_codes["CODE"]==30,"NAME"].values[0]
    ]
    out=pd.DataFrame(lis_of_values_per_smod_class_per_countries,columns=columns)
    out.to_csv("nightlight_VCMCFG.csv")




