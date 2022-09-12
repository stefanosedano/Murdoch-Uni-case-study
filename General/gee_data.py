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


def writeGeoTiff2(npArray, geoTrans, outpath, dataType: object, wkt, subset, projection):
    if len(npArray.shape) == 3:
        driver = gdal.GetDriverByName('GTiff')
        dst = driver.Create(outpath, subset, subset, npArray.shape[0], gdal.GDT_UInt16,
                            ['COMPRESS=DEFLATE', 'TILED=YES'])
        dst.SetGeoTransform(geoTrans)
        dst.SetProjection(projection)
        for i in range(1, npArray.shape[0] + 1):
            dst.GetRasterBand(i).SetNoDataValue(0)
            dst.GetRasterBand(i).WriteArray(npArray[i - 1, 0:subset, 0:subset])
        dst.FlushCache()
        dst = None

class gee:
    def __init__(self):
        # start the session
        self.COORDS = None
        self.SATELLITE = None
        self.SCALE = None
        self.CRS = None
        self.BANDS = None
        self.START_DATE = None
        self.END_DATE = None
        self.CLOUD_FILTER = 60
        self.CLD_PRB_THRESH = 50
        self.NIR_DRK_THRESH = 0.15
        self.CLD_PRJ_DIST = 1
        self.ADM0_NAME = None
        self.REGION=None

        self.OUTPATH = None
        self.BUFFER = None
        self.max_NDVI = None
        self.image=None


    def get_image_by_country_name(self):






        image = ee.ImageCollection(self.SATELLITE) \
            .filterDate(self.START_DATE, self.END_DATE) \
            .filterBounds(self.REGION)

        return image.median().clip(self.REGION)



    def get_image_to_file(self):
        response = self.get_image_url()
        with open(self.OUTPATH, 'wb') as fd:
            fd.write(response.content)


    def get_image_url(self):
        import requests
        image = self.image
        path = image.getDownloadUrl({
            'bands': self.BANDS,
            'scale': self.SCALE,
            'crs': self.CRS,
            'region': self.REGION,
            'format': 'GEO_TIFF'
        })
        response = requests.get(path)
        return response




if __name__ == "__main__":

    myImg = gee()
    myImg.CRS="EPSG:4326"
    myImg.REGION  = ee.Geometry.BBox(60.8742484882, 23.6919650335, 77.8374507995, 37.1330309108)
    myImg.BANDS = ['avg_rad']
    myImg.SATELLITE = "NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG"
    myImg.SCALE = 500

    for year in range(2014,2023):

        myImg.START_DATE ="{}-01-01".format(year)
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
            myImg.OUTPATH = "{}/{}".format(myImg.OUTPATH,year)
            if not os.path.exists(myImg.OUTPATH):
                os.makedirs(myImg.OUTPATH)

            myImg.OUTPATH = "{}/subpart_{}.tif".format(myImg.OUTPATH, counter)
            myImg.get_image_to_file()
            counter = counter + 1






