# import requests module
import ee
import subprocess
import io
import math
import os
import requests
import numpy as np
from geojson import Polygon
from shapely.geometry import shape
import os, glob
#ee.Authenticate()
credentials = ee.ServiceAccountCredentials(
            "test1-landsat8@geelandsat8.iam.gserviceaccount.com",
            "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/google_keys/geelandsat8-6af86334d7ec.json",
        )
ee.Initialize(credentials)
from osgeo import gdal, osr

from General.gee_data import gee

myImg = gee()
myImg.CRS = "EPSG:4326"
myImg.REGION = ee.Geometry.BBox(60.8742484882, 23.6919650335, 77.8374507995, 37.1330309108)
myImg.BANDS = ['flooded']
myImg.SATELLITE = "GLOBAL_FLOOD_DB/MODIS_EVENTS/V1"
myImg.SCALE = 250

subset_list = [
        ee.Geometry.Polygon([[[63.8934151, 31.0970116], [66.8934151, 31.0970116], [66.8934151, 28.0970116],
                              [63.8934151, 28.0970116], [63.8934151, 31.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[63.8934151, 28.0970116], [66.8934151, 28.0970116], [66.8934151, 25.0970116],
                              [63.8934151, 25.0970116], [63.8934151, 28.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[72.8934151, 34.0970116], [75.8934151, 34.0970116], [75.8934151, 31.0970116],
                              [72.8934151, 31.0970116], [72.8934151, 34.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[72.8934151, 31.0970116], [75.8934151, 31.0970116], [75.8934151, 28.0970116],
                              [72.8934151, 28.0970116], [72.8934151, 31.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[69.8934151, 25.0970116], [72.8934151, 25.0970116], [72.8934151, 22.0970116],
                              [69.8934151, 22.0970116], [69.8934151, 25.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[72.8934151, 37.0970116], [75.8934151, 37.0970116], [75.8934151, 34.0970116],
                              [72.8934151, 34.0970116], [72.8934151, 37.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[69.8934151, 31.0970116], [72.8934151, 31.0970116], [72.8934151, 28.0970116],
                              [69.8934151, 28.0970116], [69.8934151, 31.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[69.8934151, 28.0970116], [72.8934151, 28.0970116], [72.8934151, 25.0970116],
                              [69.8934151, 25.0970116], [69.8934151, 28.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[69.8934151, 37.0970116], [72.8934151, 37.0970116], [72.8934151, 34.0970116],
                              [69.8934151, 34.0970116], [69.8934151, 37.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[69.8934151, 34.0970116], [72.8934151, 34.0970116], [72.8934151, 31.0970116],
                              [69.8934151, 31.0970116], [69.8934151, 34.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[75.8934151, 37.0970116], [78.8934151, 37.0970116], [78.8934151, 34.0970116],
                              [75.8934151, 34.0970116], [75.8934151, 37.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[63.8934151, 34.0970116], [66.8934151, 34.0970116], [66.8934151, 31.0970116],
                              [63.8934151, 31.0970116], [63.8934151, 34.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[60.8934151, 28.0970116], [63.8934151, 28.0970116], [63.8934151, 25.0970116],
                              [60.8934151, 25.0970116], [60.8934151, 28.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[60.8934151, 31.0970116], [63.8934151, 31.0970116], [63.8934151, 28.0970116],
                              [60.8934151, 28.0970116], [60.8934151, 31.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[66.8934151, 28.0970116], [69.8934151, 28.0970116], [69.8934151, 25.0970116],
                              [66.8934151, 25.0970116], [66.8934151, 28.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[66.8934151, 25.0970116], [69.8934151, 25.0970116], [69.8934151, 22.0970116],
                              [66.8934151, 22.0970116], [66.8934151, 25.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[66.8934151, 34.0970116], [69.8934151, 34.0970116], [69.8934151, 31.0970116],
                              [66.8934151, 31.0970116], [66.8934151, 34.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[66.8934151, 31.0970116], [69.8934151, 31.0970116], [69.8934151, 28.0970116],
                              [66.8934151, 28.0970116], [66.8934151, 31.0970116]]]).bounds(),
        ee.Geometry.Polygon([[[63.8934151, 25.0970116], [66.8934151, 25.0970116], [66.8934151, 22.0970116],
                              [63.8934151, 22.0970116], [63.8934151, 25.0970116]]]).bounds()
    ]
#subset_list = [
#    ee.Geometry.BBox(60.8742484882, 23.6919650335, 68, 37.1330309108),
#    ee.Geometry.BBox(68, 23.6919650335, 77.8374507995, 37.1330309108)
#]
for year in range(2000, 2019):

    myImg.START_DATE = "{}-01-01".format(year)
    myImg.END_DATE = "{}-12-31".format(year)

    myImg.image = myImg.get_image_by_country_name_max()



    counter = 0
    for subset in subset_list:
        try:
            myImg.SUB_REGION = subset
            myImg.OUTPATH = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GLOBAL_FLOOD_DB/"
            myImg.OUTPATH = "{}/{}".format(myImg.OUTPATH, year)
            if not os.path.exists(myImg.OUTPATH):
                os.makedirs(myImg.OUTPATH)

            myImg.OUTPATH = "{}/subpart_{}.tif".format(myImg.OUTPATH, counter)
            if not os.path.exists(myImg.OUTPATH):
                myImg.get_image_to_file()
        except:
            pass

        counter = counter + 1
        if counter == len(subset_list):
            #build vrt
            outfile = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GLOBAL_FLOOD_DB/{}/index.vrt".format(year)
            datapath = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/GLOBAL_FLOOD_DB/{}/*.tif".format(
                year)
            out = subprocess.run(
                ["C:/Program Files/QGIS 3.16/bin/gdalbuildvrt.exe", outfile, datapath], shell=True)

            out = subprocess.run(
                ["C:/Program Files/QGIS 3.16/bin/gdalwarp.exe", "-t_srs" ,"ESRI:54009", "-co","COMPRESS=LZW", outfile, outfile.replace(".vrt",".tif")], shell=True)

