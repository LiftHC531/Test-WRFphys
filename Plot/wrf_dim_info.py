import numpy as np
from wrf import getvar,ALL_TIMES
import Ngl

def info(ff):
    wa = getvar(ff, "wa", timeidx=0)
    nk = np.int32(len(wa[:,0,0]))
    ny = np.int32(len(wa[0,:,0]))
    nx = np.int32(len(wa[0,0,:])); del wa
    t = getvar(ff, "times", timeidx=ALL_TIMES) #-1 isn't work
    nt = np.int32(len(t))
    return nt, nk, ny, nx, t

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
      txres.txJust = "CenterCenter"         #-- text justification	   x = vpx + vpw/2
      x = vpx + vpw/2
      Ngl.text_ndc(wks, center, x, y, txres)#-- add text to wks

    if(right != ''):
       txres.txJust = "CenterRight"          #-- text justification
       x = vpx+vpw                           #-- x-position
       Ngl.text_ndc(wks, right, x, y, txres) #-- add text to wks

if __name__ == '__main__':
   info()
   ngl_Strings()
