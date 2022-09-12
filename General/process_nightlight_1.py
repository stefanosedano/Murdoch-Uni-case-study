# import requests module
import ee
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
            "google_keys/geelandsat8-6af86334d7ec.json",
        )
ee.Initialize(credentials)
from osgeo import gdal, osr

from General.gee_data import gee

myImg = gee()
myImg.CRS = "EPSG:4326"
myImg.REGION = ee.Geometry.BBox(60.8742484882, 23.6919650335, 77.8374507995, 37.1330309108)
myImg.BANDS = ['avg_rad']
myImg.SATELLITE = "NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG"
myImg.SCALE = 500

for year in range(2014, 2023):

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
        myImg.OUTPATH = "CA_Pakistan/nightlight_VCMSLCFG"
        myImg.OUTPATH = "{}/{}".format(myImg.OUTPATH, year)
        if not os.path.exists(myImg.OUTPATH):
            os.makedirs(myImg.OUTPATH)

        myImg.OUTPATH = "{}/subpart_{}.tif".format(myImg.OUTPATH, counter)
        myImg.get_image_to_file()
        counter = counter + 1
