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
myImg.BANDS = ['avg_rad']
myImg.SATELLITE = "NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG"
myImg.SCALE = 500

for year in range(2012, 2023):

    myImg.START_DATE = "{}-01-01".format(year)
    myImg.END_DATE = "{}-12-31".format(year)

    myImg.image = myImg.get_image_by_country_name()

    subset_list = [
        ee.Geometry.BBox(60.8742484882, 23.6919650335, 68, 37.1330309108),
        ee.Geometry.BBox(68, 23.6919650335, 77.8374507995, 37.1330309108)
    ]

    counter = 0
    for subset in subset_list:
        myImg.REGION = subset
        myImg.OUTPATH = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/nightlight_VCMCFG/"
        myImg.OUTPATH = "{}/{}".format(myImg.OUTPATH, year)
        if not os.path.exists(myImg.OUTPATH):
            os.makedirs(myImg.OUTPATH)

        myImg.OUTPATH = "{}/subpart_{}.tif".format(myImg.OUTPATH, counter)
        if not os.path.exists(myImg.OUTPATH):
            myImg.get_image_to_file()

        counter = counter + 1
        if counter == 2:
            #build vrt
            outfile = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/nightlight_VCMCFG/{}/index.vrt".format(year)
            datapath = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/nightlight_VCMCFG/{}/*.tif".format(
                year)
            out = subprocess.run(
                ["C:/Program Files/QGIS 3.16/bin/gdalbuildvrt.exe", outfile, datapath], shell=True)
