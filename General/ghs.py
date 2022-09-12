# import requests module
import os
import requests
import ssl
import zipfile
import subprocess
ssl._create_default_https_context = ssl._create_unverified_context


import subprocess
def run_win_cmd(cmd):
    result = []
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    for line in process.stdout:
        result.append(line)
    errcode = process.returncode
    for line in result:
        print(line)
    if errcode is not None:
        raise Exception('cmd %s failed, see above for details', cmd)

class ghs:
    def __init__(self):
        # start the session
        self.COORDS = None
        self.SATELLITE = None
        self.SCALE = None
        self.CRS = None

    def build_vrt(self):
        datapath = "{}{}/{}/*.tif".format(self.workingdir,self.product,self.year)
        outfile = "{}{}/{}/index.vrt".format(self.workingdir,self.product,self.year)
        out = subprocess.run(
            ["C:/Program Files/QGIS 3.16/bin/gdalbuildvrt.exe", outfile, datapath], shell=True)


    def get_GHS(self):
        self.out_dir = "{}{}/{}".format(self.workingdir,self.product,self.year)
        self.out_zip_file ="{}_{}{}_GLOBE_{}_54009_{}_V1_0_{}.zip".format(self.product, self.type, self.year, self.version, self.res,self.row_col)
        self.full_out_zip_file_path = "{}/{}".format(self.out_dir,self.out_zip_file)

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        if not os.path.exists(self.full_out_zip_file_path):
            self.downloadIt()

        if not os.path.exists(self.full_out_zip_file_path.replace(".zip",".tif")):
            with zipfile.ZipFile(self.full_out_zip_file_path, "r") as zip_ref:
                zip_ref.extractall(self.out_dir)

        #test if 100m file exists
        out_dir = "{}{}/{}".format(self.workingdir, self.product, self.year)
        out_zip_file = "{}_{}{}_GLOBE_{}_54009_{}_V1_0_{}.zip".format(self.product, self.type, self.year,
                                                                      self.version, 100, self.row_col)
        outfile = "{}/{}".format(self.out_dir, self.out_zip_file).replace(".zip", ".tif")
        infile = self.full_out_zip_file_path.replace(".zip", ".tif")
        #if not os.path.exists(outfile):

        #    if self.product == "GHS_POP":
        #        resample = "average"
        #    else:
        #        resample = "near"

        #        out = subprocess.run(["C:/Program Files/QGIS 3.16/bin/gdalwarp.exe", "-t_srs", "EPSG:4326","-r",resample, "-co","compress=LZW" , infile,  outfile ],shell=True)
        #        #print(out)

    def downloadIt(self):
        import urllib.request
        # Download the file from `url` and save it locally under `file_name`:
        urllib.request.urlretrieve(self.url, self.full_out_zip_file_path)




if __name__ == "__main__":

    myGHS = ghs()
    myGHS.rows_cols = ["R5_C24","R5_C25","R6_C24","R6_C25","R7_C25"]
    myGHS.products = ["GHS_POP","GHS_SMOD"]
    myGHS.years=["2000","2005","2010","2015","2020","2025","2030"]
    myGHS.version = "R2022A"
    myGHS.workingdir = "C:/Users/email/Documents/German FFO/Murdoch-Uni-case-study/DATA/"


    for product in myGHS.products:
        for year in myGHS.years:
            for row_col in myGHS.rows_cols:
                if int(year) > 2022:
                    type = "P"
                else:
                    type = "E"
                if product == "GHS_POP":
                    res = 100
                if product == "GHS_SMOD":
                    res = 1000

                myGHS.res = res
                myGHS.type = type
                myGHS.product=product



                url = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/{}_GLOBE_{}/{}_{}{}_GLOBE_{}_54009_{}/V1-0/tiles/{}_{}{}_GLOBE_{}_54009_{}_V1_0_{}.zip".format(
                    product, myGHS.version, product, type, year, myGHS.version, res, product, type, year, myGHS.version, res,row_col)
                myGHS.product = product
                myGHS.year = year
                myGHS.row_col = row_col
                myGHS.url=url
                myGHS.get_GHS()

    #build VRT

    for product in myGHS.products:
        for year in myGHS.years:
            myGHS.product = product
            myGHS.year = year
            myGHS.build_vrt()



