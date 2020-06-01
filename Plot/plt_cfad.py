# -*- coding: utf-8 -*-
import wrf_dim_info #local function
from time import time as cpu_time
import copy
import numpy as np
import Ngl; import Nio
import netCDF4 as nc
import xarray as xr
from wrf import getvar,ALL_TIMES,interplevel, interp1d
import multiprocessing #, threading
import concurrent.futures
from numba import jit, vectorize, int32, float32, types
#from joblib import Parallel, delayed, dump, load
#import dill as pickle
#from multiprocessing.dummy import Pool as ThreadPool
#from pathos.multiprocessing import ProcessingPool as ThreadPool

#print(Nio.__version__)
##._FillValue = 9.96920996839e+36

#@jit(float32[:,:]( int32, int32, float32[:,:,:]), nopython=True, nogil=True)
def interp_var(lev):
    global dbz, z
    result = interplevel(dbz, z, lev)
    return result

def delta_z():
    kk = 100.0 #m bottom
    nk = 150   #top 15.0 km
    lev = np.zeros(nk,dtype=np.float32)
    for k in range(nk):
        lev[k] = kk; #print(lev[k])
        kk += 100.0
    return lev, nk

def interp(ff,time):
    global nlat, nlon, lev, nk
    var_m = np.zeros(shape=(nk,nlat,nlon),dtype=np.float32)
    """
    for k,ll in enumerate(lev):
        var_m[k,:,:] = interp_var(lev[k]) 
    """
    #with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
    with concurrent.futures.ProcessPoolExecutor\
         (max_workers=multiprocessing.cpu_count()) as executor:
         for k,output in enumerate(executor.map(interp_var, lev)):
             var_m[k,:,:] = output; #print(k)
   
    return var_m 

@jit(types.Tuple((float32[:], float32[:,:]))\
(int32, int32, int32, float32[:,:,:]), nopython=True, nogil=True)
def cfad(nz,jc,ic,mdb):
    nni = 65+1 #bin: 0~65 dBZ
    delta = np.arange(nni,dtype=np.float32); #print(delta) #dBZ range !! not range(ni)
    countd = np.zeros(shape=(nni,nz),dtype=np.float32)
    total = np.zeros(shape=(nz),dtype=np.float32)
    for k in range(nz):
        for j in range(jc-10, jc+11):
            for i in range(ic-10, ic+11):
                total[k] += 1.0
                for ni,db in enumerate(delta):
                    if delta[ni] == 0.0 and mdb[k,j,i] <= delta[ni] and \
                       mdb[k,j,i] > -1.0 :
                       countd[ni,k] += 1.0
                    elif delta[ni] == 65.0 and mdb[k,j,i] <= delta[ni] and \
                         mdb[k,j,i] > 64.0 :
                       countd[ni,k] += 1.0
                    else:
                       if mdb[k,j,i] <= delta[ni] and \
                          mdb[k,j,i] > delta[ni-1] :
                          countd[ni,k] += 1.0

    for k in range(nz):
        for ni in range(nni):
            countd[ni,k] /= total[k]
            countd[ni,k] *= 100.0 

    return delta, countd

def stat_kernel(dbz,medid):
    global nk, undef 
    med_list = []  
    for k in range(nk): 
        med_list.append(dbz[k,medid[k]])

    med = np.array(med_list, dtype=np.float32)
    med = np.where(np.isnan(med), undef, med) 
    return med    

def wks_setting(time):
    time = time[0:19].replace("-",":").replace("T",":")
    time = time.split(":") 
    time = time[0]+time[1]+time[2]+"_"+time[3]+time[4]
    image_name = "CFAD_"+time
    print("---Plot:",image_name,"---")
    wks = Ngl.open_wks(wks_type,image_name,rlist)
    Ngl.define_colormap(wks,"matlab_jet")
    res = Ngl.Resources()
    res.nglMaximize = True
    res.nglDraw  = False; res.nglFrame = False
    res.vpWidthF = 0.74; res.vpHeightF = 0.54
    LeftString = "Valid Time: "+time+" UTC"
    return wks, res, LeftString   

def cmap():
    colors = np.array([[0.10, 0.10, 0.10], [0.15, 0.15, 0.15], \
                    [0.20, 0.20, 0.20], [0.25, 0.25, 0.25], \
                    [0.30, 0.30, 0.30], [0.35, 0.35, 0.35], \
                    [0.40, 0.40, 0.40], [0.45, 0.45, 0.45], \
                    [0.50, 0.50, 0.50], [0.55, 0.55, 0.55], \
                    [0.60, 0.60, 0.60], [0.65, 0.65, 0.65], \
                    [0.70, 0.70, 0.70], [0.75, 0.75, 0.75], \
                    [0.80, 0.80, 0.80], [0.85, 0.85, 0.85]],'f')
    return colors

#--------------------------------------------------------------------
start_time = cpu_time()
undef = 9.96920996839e+36
ic = 68-1 #analysis center
jc = 34-1
#for j in range(jc-10, jc+11): print(j)
(lev,nk) = delta_z()
#--------------------------------------------------------------------
Path = "/home/mlhchen/2020_rainmaking/20191205/Sensitivity/"
File = Path + "WDM6log_wrfout_d04_2019-12-05_00:00:00"
ff = nc.Dataset(File) 
(ntimes, nlev, nlat, nlon, times) = wrf_dim_info.info(ff)
print("Domain:{}".format(ntimes)+"|{}".format(nlev)+\
    "|{}".format(nlat)+"|{}".format(nlon)," .......Get Dimensions!!")
#File_out = "dbz_result.bin"
#newFile = open(File_out,"wb")
 
wks_type = "png"
rlist = Ngl.Resources()
rlist.wkWidth  = 3000 #page resolutio
rlist.wkHeight = 3000
for t in range(ntimes):
    dtime = str(times[t].values)[0:19]; print("\t"+dtime)
    dbz = getvar(ff, "REFL_10CM", timeidx=t)
    z = getvar(ff, "z", timeidx=t)
    dbz_m = interp(ff,t); #print(dbz_m)
    #newFile.write(bytearray(dbz_m)); del dbz, z 
    (delta,countd) = cfad(nk,jc,ic,dbz_m) 
    countd = np.transpose(countd) #due to plot, (bin,lev)->(lev,bin)
    #med = np.median(countd,axis=1); print(med.shape) 
    #med = np.argsort(countd[:,:],axis=1)[:, len(countd[0,:])//2]; #print(med.shape)
    dbz_stat = np.reshape(dbz_m[:,jc-10:jc+11,ic-10:ic+11],(nk,21*21))
    medid = np.argsort(dbz_stat,axis=1)[:, len(dbz_stat[0,:])//2]; #print(med.shape)
    IQ1id = np.argsort(dbz_stat,axis=1)[:, len(dbz_stat[0,:])//4]; #print(med.shape)
    IQ3id = np.argsort(dbz_stat,axis=1)[:, len(dbz_stat[0,:])*3//4]; #print(med.shape)
    #print(medid.shape)
    med = stat_kernel(dbz_stat,medid)
    IQ1 = stat_kernel(dbz_stat,IQ1id)
    IQ3 = stat_kernel(dbz_stat,IQ3id)
    """ Plot """
    (wks,res,leftname) = wks_setting(dtime)
    ##._FillValue = 9.96920996839e+36
    resl = copy.deepcopy(res) #for plot line
    res.cnLineLabelFormat="@*+^sg"
    res.cnLevelSelectionMode = "ManualLevels"
    res.cnMinLevelValF       =  0.5 
    res.cnMaxLevelValF       = 10.0 #5.25
    res.cnLevelSpacingF      =  0.5 #0.25 
    res.cnFillColors = np.concatenate(((0,), np.arange(2,65,3)), axis=0)
    res.cnFillOn             = True
    res.cnLineLabelsOn       = False # label the contour line
    res.cnLinesOn            = False # draw the contour line

    res.trXMinF = 0; res.trXMaxF = 65.0
    res.trYMinF = -1; res.trYMaxF = nk-1
    res.tmXBMode = "Explicit"
    res.tmXBValues = np.arange(0,65+1,5) 
    res.tmXBLabels = res.tmXBValues
    res.tmXBMinorValues = np.arange(0,65+1,1) 
    res.tiXAxisString = "~F25~dBZ"
    #res.sfYArray = np.arange(15) # X and Y
    res.tmYLMode = "Explicit"
    res.tmYLValues = np.arange(9,nk,10); #print(np.arange(9,nk,10))
    res.tmYLLabels = lev[res.tmYLValues]/1000.0
    res.tmYLMinorValues = np.arange(4,nk,10) 
    res.tiYAxisString = "~F25~Height (km)"

    res.lbTitleFontHeightF= 0.015
    res.lbTitleString = "~F26~(%)"; res.lbTitlePosition  = "Top"
    res.lbTitleOffsetF = -0.03
    res.lbLabelFont        = 26; res.lbLabelFontHeightF = 0.013
    res.pmLabelBarHeightF  = 0.65; res.pmLabelBarWidthF   = 0.07
    res.pmLabelBarOrthogonalPosF = 0.03 #L- & R+ 
    res.pmLabelBarParallelPosF = 0.55   #B- & T+  

    res.tmXBLabelFont = 26; res.tmYLLabelFont = res.tmXBLabelFont 
    res.tmXBMajorThicknessF = 10; res.tmYLMajorThicknessF = res.tmXBMajorThicknessF
    res.tmXBMinorThicknessF = 10; res.tmYLMinorThicknessF = res.tmXBMinorThicknessF
    res.tmXBLabelFontHeightF = 0.015; res.tmYLLabelFontHeightF = res.tmXBLabelFontHeightF 
    res.tmBorderThicknessF = 15.0

    resl.caXMissingV = undef;resl.caYMissingV = resl.caXMissingV
    resl.xyLineThicknesses = 20 
    resl.xyLineColors = "gray76"

    #res.cnFillPalette = cmap()
    plot = []
    plot.append(Ngl.contour(wks,countd[:,:],res))  
    plot.append(Ngl.xy(wks,med,np.arange(nk),resl))    
    resl.xyDashPatterns = 1
    plot.append(Ngl.xy(wks,IQ1,np.arange(nk),resl))    
    resl.xyDashPatterns = 2
    plot.append(Ngl.xy(wks,IQ3,np.arange(nk),resl))    
    Ngl.overlay(plot[0],plot[1])
    Ngl.overlay(plot[0],plot[2])
    Ngl.overlay(plot[0],plot[3])

     
    wrf_dim_info.ngl_Strings(wks, plot[0], left=leftname, center="", right="")
    Ngl.draw(plot[0])
    Ngl.frame(wks)
    Ngl.destroy(wks); del res, resl, plot

Ngl.end()
end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("plt_cfad.py has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")

