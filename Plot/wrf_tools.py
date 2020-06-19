import numpy as np
from wrf import getvar,ALL_TIMES,get_pyngl,latlon_coords,to_np
import Ngl, Nio
import glob

#Li-Hsin Chen 2020

def info(ff):
    var = getvar(ff, "QCLOUD", timeidx=0)
    nz = np.int32(len(var[:,0,0]))
    ny = np.int32(len(var[0,:,0]))
    nx = np.int32(len(var[0,0,:])) 
    (lat,lon) = latlon_coords(var) # Get the latitude and longitude points
    lat = to_np(lat); lon = to_np(lon)
    #Set some map options based on information in WRF output file
    res = get_pyngl(var); del var 
    res.tfDoNDCOverlay = True # required for native projection
    t = getvar(ff, "times", timeidx=ALL_TIMES) #-1 isn't work
    nt = np.int32(len(t))
    print("Domain:{}".format(nt)+"|{}".format(nz)+\
    "|{}".format(ny)+"|{} ".format(nx)+"."*6+"Get Dimensions!! \
    \nGet WRF map-related resources and coordinates(lat. & long.)")
    return nt, nz, ny, nx, t, lat, lon, res

def search_wrf_file(Path="./"):
    FileID = glob.glob(Path+"*wrf*:*") #just like "ls *wrf*:*"
    for i, Filename in enumerate(FileID): print("\033[92m({}){}\033[0m".format(i,Filename))
    ID = input("\033[92mInput File ID:\t\033[0m")
    if ID.replace(".","",1).isdigit() and float(ID) >= 0.0: #integer or float True
       File = FileID[int(round(float(ID),0))]
    else:
       File = FileID[0] #Default: 0
       print("\033[91m\tCan't find ur input ID (Default: 0)\033[0m")
    print("\033[103mThe choosen file: "+File+"\033[0m")
    return File

def add_shapefile_polylines(ff,wks,plot,color="black",thick=10.0):
    """ Attach shapefile polylines to map """
    f_shap = Nio.open_file(ff, "r")
    lon = f_shap.variables["x"][:] #np.ravel()
    lat = f_shap.variables["y"][:]
    lnres = Ngl.Resources()
    lnres.gsLineColor = color
    lnres.gsLineThicknessF = thick
    lnres.gsSegments = f_shap.variables["segments"][:,0]
    return Ngl.add_polyline(wks, plot, lon, lat, lnres)

def target_area(res,lat,lon,jc=33,ic=67,jrg=15,irg=15):
    plt_target_area = False
    check_plt_area = input("\033[92mPlot target area ? (Y/N) \033[0m")
    check_plt_area = np.where(check_plt_area == "y","Y",check_plt_area)
    if check_plt_area == "Y": plt_target_area = True
    if plt_target_area:
       pt = [jc,ic] 
       jj = np.array([pt[0]-jrg, pt[0]+jrg]).astype(int) 
       ii = np.array([pt[1]-irg, pt[1]+irg]).astype(int)
    else:
       jj = np.array([0, len(lat[:,0])-1]).astype(int)
       ii = np.array([0, len(lat[0,:])-1]).astype(int)
    # Zoom in on map area of interest
    res.mpLimitMode = "Corners" 
    res.mpLeftCornerLatF  = lat[jj[0],ii[0]] #bottom y left
    res.mpLeftCornerLonF  = lon[jj[0],ii[0]] #bottom x left
    res.mpRightCornerLatF = lat[jj[1],ii[1]] #Top y left
    res.mpRightCornerLonF = lon[jj[1],ii[1]] #Top x left
    print("start point: ({:3d}){:.3f}, ({:3d}){:.3f}\
          \n  end point: ({:3d}){:.3f}, ({:3d}){:.3f}"\
     .format(jj[0],res.mpLeftCornerLatF,ii[0],res.mpLeftCornerLonF, \
             jj[1],res.mpRightCornerLatF,ii[1],res.mpRightCornerLonF))
        
    return res, ii, jj


def cwbrain_colorbar():
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [1.000,1.000,1.000], [0.608,1.000,1.000], \
     [0.000,0.812,1.000], [0.039,0.596,1.000], [0.039,0.396,1.000], \
     [0.188,0.600,0.039], [0.196,1.000,0.000], \
     [0.973,1.000,0.000], [1.000,0.796,0.000], [1.000,0.603,0.000], \
     [0.980,0.012,0.000], [0.800,0.000,0.012], [0.627,0.000,0.000], \
     [0.596,0.000,0.604], [0.765,0.016,0.800], \
     [0.973,0.020,0.953], [0.996,0.796,1.000], \
    ], 'f')
    return colors

#-- Routine ngl_Strings: draw left, right and/or center string
def ngl_Strings(wks, plot, left='', center='', right=''):
    """
       ngl_Strings(wks, plot, left='', center='', right='')

       Add annotations
	- left, right or center string above plot
			
       Correspond to NCL's 
	  gsnLeftString, gsnCenterString, gsnRightString'
    """
    assert str(getattr(wks,'__class__')  == "<class 'int'>"), 'ERROR - 1st parameter is not a Ngl wks'
    assert str(getattr(plot,'__class__') == "<class 'ngl.PlotIds'>"), 'ERROR - 2nd parameter is not a Ngl plot'
	
    vpx = Ngl.get_float(plot,"vpXF")         #-- retrieve value of res.vpXF from plot
    vpy = Ngl.get_float(plot,"vpYF")         #-- retrieve value of res.vpYF from plot
    vpw = Ngl.get_float(plot,"vpWidthF")     #-- retrieve value of res.vpWidthF from plot
	
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.018              #-- font size for left, center and right string
    txres.txFont = 26
	
    y = vpy + 0.025                          #-- y-position
	
    if(left != ''):
       txres.txJust = "CenterLeft"           #-- text justification
       x = vpx                               #-- x-position
       Ngl.text_ndc(wks, left, x, y, txres)  #-- add text to wks
	   
    if(center != ''):
      txres.txJust = "CenterCenter"          #-- text justification	   x = vpx + vpw/2
      x = vpx + vpw/2
      Ngl.text_ndc(wks, center, x, y, txres) #-- add text to wks

    if(right != ''):
       txres.txJust = "CenterRight"          #-- text justification
       x = vpx+vpw                           #-- x-position
       Ngl.text_ndc(wks, right, x, y, txres) #-- add text to wks

def Date_string(yy,mm,dd,hh,mi,TW_LST=False):
    def det_mon_yr_add(y,m,d,h):
        #string to integer
        yy = int(y); mm = int(m)
        dd = int(d); hh = int(h)
        months_day = [31,28,31,30,31,30,31,31,30,31,30,31]
        # leap year or ordinary year
        if yy % 4 == 0 and yy % 100 != 0:
           months_day[2-1] = 29
        elif yy % 400 == 0:
           months_day[2-1] = 29
        #update hour & day & month
        if hh >= 24:
           hh  = hh - 24
           dd = dd + 1
           if dd > months_day[mm-1]:
              dd = dd - months_day[mm-1]
              mm = mm + 1
              if mm > 12:
                 mm = mm - 12
                 yy = yy + 1
        return str(yy), format(mm,'02d'), \
               format(dd,'02d'), format(hh,'02d')
    #--------------------------------------
    time_cord = "UTC"
    if TW_LST:
       time_cord = "LST"
       hh = str(int(hh) + 8) # UTC to LST, GMT+8
       (yy,mm,dd,hh) = det_mon_yr_add(yy,mm,dd,hh)
       yy = str(int(yy) - 1911) # the year of the Republic Era
       #print("{} LST {} {}, {}".format(hh,dd,mm,yy))

    year = int(yy) #string to integer
    mon  = int(mm)
    day  = int(dd)
    hr   = int(hh)
    months = ["Jan.","Feb.","Mar.","Apr.","May","June", \
              "July","Aug.","Sept.","Oct.","Nov.","Dec."]
    for i, month in enumerate(months):
        mm = np.where(mon == i + 1, month, mm)

    days_str = ["st","nd","rd"]
    if day <= 3:
       for i, day_str in enumerate(days_str):
           if day == i+1: dd = dd+"~S~"+day_str+"~N~"
    else:
       dd = dd+"~S~th~N~"

    date = "{}{} {} {} {}, {}".format(hh,mi,time_cord,dd,mm,yy)

    print(date)
    return date

if __name__ == '__main__':
   info()
   search_wrf_file()
   add_shapefile_polylines()
   cwbrain_colorbar()
   ngl_Strings()
   Date_string()
