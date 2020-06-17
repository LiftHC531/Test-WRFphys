# -*- coding: utf-8 -*-
#from __future__ import print_function
import wrf_dim_info #local function
from time import time as cpu_time
import copy
import numpy as np
import glob
import Ngl; import Nio
import netCDF4 as nc
import xarray as xr
from wrf import getvar,ALL_TIMES,get_pyngl

def get_rain(ff,time):
    rain_exp = getvar(ff, "RAINNC", timeidx=time) #proudced by microphysics Par.
    rain_con = getvar(ff, "RAINC", timeidx=time)  #produced by cumulus Par.
    print("Variable: {}, Units: {}\nVariable: {}, Units: {} "\
         .format(rain_exp.description,rain_exp.units,rain_con.description,rain_con.units))
    rain_con = rain_con + rain_exp
    res = get_pyngl(rain_exp) #Set some map options based on information in WRF output file
    res.tfDoNDCOverlay = True # required for native projection
    return rain_con, res
    
def wks_setting(time,res):
    time = time[0:19].replace("-",":").replace("T",":")
    time = time.split(":")
    time = time[0]+time[1]+time[2]+"_"+time[3]+time[4]
    image_name = "QPF_"+time
    print("---Plot:",image_name,"---")
    wks_type = "png"
    rlist = Ngl.Resources()
    rlist.wkWidth  = 3000 #page resolutio
    rlist.wkHeight = 3000
    wks = Ngl.open_wks(wks_type,image_name,rlist)
    Ngl.define_colormap(wks, wrf_dim_info.cwbrain_colorbar())
    res.nglMaximize = True
    res.nglDraw  = False; res.nglFrame = False
    #res.vpWidthF = 0.74; res.vpHeightF = 0.54
    LeftString = "Valid Time: "+time+" UTC"
    return wks, res, LeftString

def add_shapefile_polylines(wks,plot,color="black",thick=10.0):
    """ Attach shapefile polylines to map """
    path = "/home/WRF/shapefile/Shimen/"
    ff = "COUNTY_MOI_1080617.shp" #"TWN_adm2.shp"
    f_shap = Nio.open_file(path+ff, "r")
    lon = f_shap.variables["x"][:] #np.ravel()
    lat = f_shap.variables["y"][:]
    lnres = Ngl.Resources()
    lnres.gsLineColor = color 
    lnres.gsLineThicknessF = thick
    lnres.gsSegments = f_shap.variables["segments"][:,0]
    return Ngl.add_polyline(wks, plot, lon, lat, lnres)

#--------------------------------------------------------------------
start_time = cpu_time()
undef = 9.96920996839e+36
#--------------------------------------------------------------------
Path = "./"
FileID = glob.glob(Path+"*wrfout*")
for i, Filename in enumerate(FileID): print("({}){}".format(i,Filename))
print("Input File ID:")
File = FileID[np.int(input())]; print("The File u pick:"+File)
ff = nc.Dataset(File); #ff = Nio.open_file(File+".nc")
(ntimes, nlev, nlat, nlon, times) = wrf_dim_info.info(ff)
print("Domain:{}".format(ntimes)+"|{}".format(nlev)+\
    "|{}".format(nlat)+"|{}".format(nlon)," .......Get Dimensions!!")

t = 0 
dtime = str(times[t].values)[0:19]; print("\t"+dtime)
(rain,res) = get_rain(ff,t)
""" Plot """
(wks,res,leftname) = wks_setting(dtime,res)
    #resl = copy.deepcopy(res) #for plot line
res.cnLineLabelFormat = "@*+^sg"
#res.cnFillMode        = "RasterFill"        # These two resources
#res.trGridType        = "TriangularMesh"    # can speed up plotting.
res.cnLevelSelectionMode = "ExplicitLevels" #"ManualLevels"
res.cnLevels = np.array([1,2,6,10,15,20,30,40,50,70,90,110,130,150,200,300], 'f')
res.cnFillOn       = True
res.cnLineLabelsOn = False # label the contour line
res.cnLinesOn      = False # draw the contour line

res.mpFillOn        = False
res.mpOutlineOn     = False
res.mpGridAndLimbOn = False
res.tmXTOn = False; res.tmYROn = False
res.tmXBLabelFont = 26; res.tmYLLabelFont = res.tmXBLabelFont
res.tmXBLabelFontHeightF = 0.015; res.tmYLLabelFontHeightF = res.tmXBLabelFontHeightF
res.tmXBMajorThicknessF  = 25.0 ; res.tmYLMajorThicknessF = res.tmXBMajorThicknessF
res.tmBorderThicknessF   = 25.0 

res.lbOrientation = "horizontal"
res.lbBoxSeparatorLinesOn = False
res.lbTitleOn          = True
res.lbTitleString      = "~F25~mm"
res.lbTitlePosition    = "Top"
res.lbTitleFontHeightF = 0.020    # make title smaller
res.lbTitleDirection   = "Across" # title direction
res.lbLabelFontHeightF = 0.013; res.lbLabelFont = 26
#res_pn.lbLabelAlignment = "BoxCenters" #label orientation
res.pmLabelBarHeightF = 0.07; res.pmLabelBarWidthF  = 0.65
res.pmLabelBarOrthogonalPosF = 0.02 #-0.02 #+D -U

res.tiMainString = ""

plot = []
plot.append(Ngl.contour_map(wks,rain[:,:],res))


shp_lines = add_shapefile_polylines(wks,plot[0])
wrf_dim_info.ngl_Strings(wks, plot[0], left=leftname, center="", right="")
Ngl.draw(plot[0])
Ngl.frame(wks)
Ngl.destroy(wks); del res, plot



Ngl.end()
end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("plt_rain.py has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")
