import os
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from osgeo import osr, ogr, gdal



def clip_to_extent_and_resoultion(in_src, extent, resolution,resampling_method,epsg):

    # Resample with GDAL warp
    outDs = gdal.Warp('',
                      in_src,
                      outputBounds=extent,
                      dstSRS=epsg,
                      outputType=gdal.GDT_Float32,
                      xRes=resolution, yRes=resolution,
                      resampleAlg=resampling_method,
                      format="MEM")



    array_out = outDs.GetRasterBand(1).ReadAsArray()

    return array_out


def clip_to_extent_resolution_inproj_outproj(in_src, extent,resolution,inproj,outproj):

    # Resample with GDAL warp
    outDs = gdal.Warp('',
                      in_src,
                      srcSRS=inproj,
                      dstSRS=outproj,
                      outputBounds=extent,
                      xRes=resolution, yRes=resolution,
                      format="MEM")

    array_out = outDs.GetRasterBand(1).ReadAsArray()

    return array_out

def clip_to_extent_resolution(in_src, extent,resolution):

    # Resample with GDAL warp
    outDs = gdal.Warp('',
                      in_src,
                      outputBounds=extent,
                      xRes=resolution, yRes=resolution,
                      format="MEM")

    array_out = outDs.GetRasterBand(1).ReadAsArray()

    return array_out

def clip_to_extent(in_src, extent):

    # Resample with GDAL warp
    outDs = gdal.Warp('',
                      in_src,
                      outputBounds=extent,
                      format="MEM")

    array_out = outDs.GetRasterBand(1).ReadAsArray()

    return array_out


# a function to find the index of the point closest pt
# (in squared distance) to give lat/lon value.
def getclosest_ij(lats,lons,latpt,lonpt):
  # find squared distance of every point on grid
  latsq = (lats-latpt)**2
  lonsq = (lons-lonpt)**2

  return latsq.argmin(), lonsq.argmin()


def writeGeoTiff_v3(npArray,geoTrans,outpath,dataType,wkt,MEM):

    if (MEM == 1):
        driver = gdal.GetDriverByName('MEM')
    else:
        driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(
         outpath,
         npArray.shape[1],
         npArray.shape[0],
         1,
         dataType)
         #['COMPRESS=LZW']


    if (geoTrans != ''):
        dataset.SetGeoTransform(geoTrans)
    dataset.GetRasterBand(1).WriteArray(npArray)
    if (wkt != ''):
        dataset.SetProjection(wkt)

    dataset.FlushCache()

    return dataset


def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and
        new axes must divide old ones.

    Example
    -------
    #>>> m = np.arange(0,100,1).reshape((10,10))
    #>>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    #>>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray


def submatsum(data,n,m):
    # return a matrix of shape (n,m)
    bs = data.shape[0]//n,data.shape[1]//m  # blocksize averaged over
    return np.reshape(np.array([np.sum(data[k1*bs[0]:(k1+1)*bs[0],k2*bs[1]:(k2+1)*bs[1]]) for k1 in range(n) for k2 in range(m)]),(n,m))

def submatnansum(data,n,m):
    # return a matrix of shape (n,m)
    bs = data.shape[0]//n,data.shape[1]//m  # blocksize averaged over
    return np.reshape(np.array([np.nansum(data[k1*bs[0]:(k1+1)*bs[0],k2*bs[1]:(k2+1)*bs[1]]) for k1 in range(n) for k2 in range(m)]),(n,m))

def submatnanmean(data,n,m):
    # return a matrix of shape (n,m)
    bs = data.shape[0]//n,data.shape[1]//m  # blocksize averaged over
    return np.reshape(np.array([np.nanmean(data[k1*bs[0]:(k1+1)*bs[0],k2*bs[1]:(k2+1)*bs[1]]) for k1 in range(n) for k2 in range(m)]),(n,m))



def get_GADM():
    """
    Returns a geopandas dataframe of GADM
    """
    import geopandas as gpd
    import shapely
    shapely.speedups.disable()


    dataset = gpd.read_file("C:/Users/email/Documents/German FFO/CODE/reference_datasets/GADM_UPDATE_SPARSE_20210922.gpkg", layer="GADM_ADM1_SPARSE")
    return dataset


def get_GAUL():
    """
    Returns a geopandas dataframe of GAUL
    """

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pathGeo = os.path.join(base_dir, "reference_datasets/Gaul/gaul1_asap.shp")
    dataset = gpd.read_file(pathGeo)
    return dataset


def gaul1_to_gadm1(GADM, GAUL, asap_id):

    return GADM[
        GADM.intersects(
            GAUL[GAUL["asap1_id"] == asap_id].representative_point().values[0]
        )
    ]

def build_gid_0_lookup_table(GID_0_list, res):

    #build region mask

    # Burkina Faso, Chad, Mali, Mauritania, Niger, Nigeria and Cameroon
    gadm = get_GADM()
    if (len(GID_0_list)>0):
        gadm = gadm.loc[
            gadm["GID_0"].isin(GID_0_list)]

    else:
        gadm = gadm.loc[
            gadm["NAME_0"].isin(["Burkina Faso", "Chad", "Mali", "Mauritania", "Niger", "Nigeria", "Cameroon"])]


    gadm = gadm.dissolve(by="GID_0")
    gadm.insert(0, 'GID_0_id', range(0, 0 + len(gadm)))

    gadm.to_file("../reference_datasets/gadm_selected.shp")

    #input_shp = ogr.Open("../reference_datasets/gadm_selected.shp")
    #shp_layer = input_shp.GetLayer()

    minx, miny, maxx, maxy = gadm.geometry.total_bounds
    minx, miny, maxx, maxy = round(minx - res, 1), round(miny - res, 1), round(maxx + res, 1), round(maxy + res, 1)

    gadm_tif = gdal.Rasterize("../reference_datasets/gadm_selected.tif", "../reference_datasets/gadm_selected.shp", xRes=res, yRes=res,
                        attribute="GID_0_id", outputBounds=[minx, miny, maxx, maxy],
                        outputType=gdal.GDT_Int32)

    path_XYZ = "output/region_lookup_mask.csv"
    format = "XYZ"
    driver = gdal.GetDriverByName(format)

    # Output to new format
    XYZ_region_mask = driver.CreateCopy(path_XYZ, gadm_tif, 0)

    # BUILD XYZ dataset
    a = gadm.reset_index()
    a = a[["GID_0","GID_0_id"]]
    xyz = pd.read_csv(path_XYZ, sep=" ", names=["lon", "lat", "value"], header=None)

    xyz = xyz.merge(a, left_on="value",right_on="GID_0_id", how="left")

    xyz["value"] = xyz["GID_0"]

    xyz = xyz[["lon", "lat", "value"]]



    xyz.loc[xyz['value'] == "-nan(ind)", "value"] = ""



    xyz.loc[xyz['value'] == 0, "value"] = ""

    xyz["lon"] = xyz["lon"].round(2)
    xyz["lat"] = xyz["lat"].round(2)
#    xyz["value"] = pd.to_numeric(xyz["value"]).round(4)

    xyz = xyz.dropna()

    xyz = xyz.loc[ xyz["lon"] >= -180]
    xyz = xyz.loc[xyz["lat"] >= -90]
    xyz = xyz.loc[xyz["lon"] <= 180]
    xyz = xyz.loc[xyz["lat"] <= 90]

    xyz.to_csv(path_XYZ, index=False)


    return minx, miny, maxx, maxy, gadm_tif


def build_gid_1_lookup_table(GID_0_list, res):

    #build region mask

    # Burkina Faso, Chad, Mali, Mauritania, Niger, Nigeria and Cameroon
    gadm = get_GADM()
    if (len(GID_0_list)>0):
        gadm = gadm.loc[
            gadm["GID_0"].isin(GID_0_list)]

    else:
        gadm = gadm.loc[
            gadm["NAME_0"].isin(["Burkina Faso", "Chad", "Mali", "Mauritania", "Niger", "Nigeria", "Cameroon"])]


    gadm = gadm.dissolve(by="GID_1")
    gadm.insert(0, 'GID_1_id', range(0, 0 + len(gadm)))

    gadm.to_file("../reference_datasets/gadm_selected.shp")

    #input_shp = ogr.Open("../reference_datasets/gadm_selected.shp")
    #shp_layer = input_shp.GetLayer()

    minx, miny, maxx, maxy = gadm.geometry.total_bounds
    minx, miny, maxx, maxy = round(minx - res, 1), round(miny - res, 1), round(maxx + res, 1), round(maxy + res, 1)

    gadm_tif = gdal.Rasterize("../reference_datasets/gadm_selected_gid_1.tif", "../reference_datasets/gadm_selected.shp", xRes=res, yRes=res,
                        attribute="GID_1_id", outputBounds=[minx, miny, maxx, maxy],
                        outputType=gdal.GDT_Int32)

    path_XYZ = "output/region_lookup_mask_gid_1.csv"
    format = "XYZ"
    driver = gdal.GetDriverByName(format)

    # Output to new format
    XYZ_region_mask = driver.CreateCopy(path_XYZ, gadm_tif, 0)

    # BUILD XYZ dataset
    a = gadm.reset_index()
    a = a[["GID_1","GID_1_id"]]
    xyz = pd.read_csv(path_XYZ, sep=" ", names=["lon", "lat", "value"], header=None)

    xyz = xyz.merge(a, left_on="value",right_on="GID_1_id", how="left")

    xyz["value"] = xyz["GID_1"]

    xyz = xyz[["lon", "lat", "value"]]



    xyz.loc[xyz['value'] == "-nan(ind)", "value"] = ""



    xyz.loc[xyz['value'] == 0, "value"] = ""

    xyz["lon"] = xyz["lon"].round(2)
    xyz["lat"] = xyz["lat"].round(2)
#    xyz["value"] = pd.to_numeric(xyz["value"]).round(4)

    xyz = xyz.dropna()

    xyz = xyz.loc[ xyz["lon"] >= -180]
    xyz = xyz.loc[xyz["lat"] >= -90]
    xyz = xyz.loc[xyz["lon"] <= 180]
    xyz = xyz.loc[xyz["lat"] <= 90]

    xyz.to_csv(path_XYZ, index=False)


    return minx, miny, maxx, maxy, gadm_tif

def get_list_of_countries():
    gadm = get_GADM()
    return np.unique(gadm["GID_0"])


def get_list_of_provinces():
    gadm = get_GADM()
    return np.unique(gadm["GID_1"])


def get_country_mask_projected(gid_0,res,crs,wkt):
    gadm = get_GADM()
    gadm = gadm.loc[
        gadm["GID_0"].isin([gid_0])]

    gadm["dis"] = 1
    gadm=gadm.to_crs(crs)
    gadm = gadm.dissolve(by="dis")
    geom = gadm.geometry.to_wkt().values[0]

    minx, miny, maxx, maxy = gadm.geometry.total_bounds
    minx, miny, maxx, maxy = round(minx - res, 1), round(miny - res, 1), round(maxx + res, 1), round(maxy + res, 1)

    gadm = get_GADM()
    gadm = gadm.loc[
        gadm["GID_0"].isin([gid_0])]

    gadm['CODE'] = np.arange(len(gadm))
    gadm['CODE'] = (gadm['CODE'] + 1)*100000
    name_code = gadm[["CODE","NAME_1"]]
    gadm = gadm.to_crs(crs)

    gadm.to_file("gadm_selected.shp")
    driver = ogr.GetDriverByName('ESRI Shapefile')
    gadm = driver.Open("gadm_selected.shp", 0)  # 0 means read-only. 1 means writeable.
    lyr = gadm.GetLayer()


    target_ds = gdal.GetDriverByName('MEM').Create('', round((maxx - minx) / res), round((maxy - miny) / res), 1,
                                                   gdal.GDT_Int32)
    new_geotrans = (minx, res, 0, maxy, 0, -res)
    target_ds.SetGeoTransform(new_geotrans)
    target_ds.SetProjection(wkt)

    gdal.RasterizeLayer(target_ds, [1], lyr, options=['ATTRIBUTE=CODE'])


    return minx, miny, maxx, maxy, target_ds,new_geotrans,name_code


def clip_to_extent_resolution_projection(in_src, extent,resolution,epsg):
    # Resample with GDAL warp
    outDs = gdal.Warp('',
                      in_src,
                      dstSRS=epsg,
                      outputBounds=extent,
                      xRes=resolution, yRes=resolution,
                      format="MEM")

    array_out = outDs.GetRasterBand(1).ReadAsArray()

    return array_out


def get_country_mask(GID_0_list, res):


    #build region mask


    # Burkina Faso, Chad, Mali, Mauritania, Niger, Nigeria and Cameroon

    gadm = 0
    if ("all" in GID_0_list ):
        gadm = get_GADM()
        gadm = gadm


    elif (len(GID_0_list)>0):
        gadm = get_GADM()
        gadm = gadm.loc[
            gadm["GID_0"].isin(GID_0_list)]


    if (gadm.shape[0] ==0):
        gadm = get_GADM()
        gadm = gadm.loc[
            gadm["GID_1"].isin(GID_0_list)]

    #else:
    #    gadm = get_GADM()
    #    gadm = gadm.loc[
    #        gadm["NAME_0"].isin(["Burkina Faso", "Chad", "Mali", "Mauritania", "Niger", "Nigeria", "Cameroon"])]

    gadm["dis"] = 1
    gadm = gadm.dissolve(by="dis")
    geom = gadm.geometry.to_wkt().values[0]
    #gadm.to_file("../reference_datasets/sahel.shp")

    #input_shp = ogr.Open("../reference_datasets/sahel.shp")
    #shp_layer = input_shp.GetLayer()

    minx, miny, maxx, maxy = gadm.geometry.total_bounds
    minx, miny, maxx, maxy = round(minx - res, 1), round(miny - res, 1), round(maxx + res, 1), round(maxy + res, 1)

    # Create a memory layer to rasterize from.
    rast_ogr_ds = ogr.GetDriverByName('Memory').CreateDataSource('wrk')
    sr_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
    sr = osr.SpatialReference(sr_wkt)
    rast_mem_lyr = rast_ogr_ds.CreateLayer('poly', srs=sr)

    # Add a polygon.
    wkt_geom = geom

    feat = ogr.Feature(rast_mem_lyr.GetLayerDefn())
    feat.SetGeometryDirectly(ogr.Geometry(wkt=wkt_geom))

    rast_mem_lyr.CreateFeature(feat)


    target_ds = gdal.GetDriverByName('MEM').Create('', round((maxx - minx)/res), round((maxy - miny)/res), 1,gdal.GDT_Byte)
    target_ds.SetGeoTransform((minx, res, 0, maxy, 0, -res))
    target_ds.SetProjection(sr_wkt)



    gdal.RasterizeLayer(target_ds, [1], rast_mem_lyr, burn_values=[255])

    target_ds.FlushCache()


    return minx, miny, maxx, maxy, target_ds




if __name__ == "__main__":

    countries_of_interest = ["MRT", "BFA", "MLI", "NER", "TCD"]
    dates = ["2000-01-01", "2021-06-30"]
