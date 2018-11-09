## solar_wind_to_dst.py

## This is a code for showing an example how to convert
## solar wind observations (here OMNI2 data) to the geomagnetic Dst index
## https://github.com/cmoestl/solar_wind_to_dst
## Author: C. Moestl, IWF Graz, Austria

## needs file dst_module.py

## If results made with this code are used for scientific publications, 
## please contact me before submission: christian.moestl@oeaw.ac.at 
## or twitter @chrisoutofspace

## latest update: November 2018
## tested in python 3.5.2 with sunpy and seaborn installed

## available methods: 
## Burton et al. 1975 doi:10.1029/JA080i031p04204  
## OBrien & McPherron 2000 doi: 10.1029/1998JA000437
## Temerin and Li 2002 doi: 10.1029/2001JA007532

## MIT LICENSE
## Copyright 2018, Christian Moestl 
## Permission is hereby granted, free of charge, to any person obtaining a copy of this 
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, 
## merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
## permit persons to whom the Software is furnished to do so, subject to the following 
## conditions:
## The above copyright notice and this permission notice shall be included in all copies 
## or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
## PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
## CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
## OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import sys
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
import seaborn as sns
import os
import copy
import pdb
import urllib

import dst_module
from dst_module import get_omni_data
from dst_module import plot_omni
from dst_module import convert_omni_time
from dst_module import make_dst_from_wind
from dst_module import interslice_omni



################################### main program ########################################

#closes all plots
plt.close('all')

########### load data (only needed 1 time)
#if omni2 hourly data is not here, download:
if not os.path.exists('omni2_all_years.dat'):
  #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
  print('download OMNI2 data from')
  omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
  print(omni2_url)
  try: urllib.request.urlretrieve(omni2_url, 'omni2_all_years.dat')
  except urllib.error.URLError as e:
      print(' ', http_sta_pla_file_str[p],' ',e.reason)

#if omni2 hourly data is not yet converted and saved as pickle
if not os.path.exists('omni2_all_years_pickle.p'):
 #load OMNI2 dataset from .dat file
 o=get_omni_data()
 #o.time,day,hour,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp
 
 
 #save for faster loading:
 pickle.dump(o, open('omni2_all_years_pickle.p', 'wb') )

else:  o=pickle.load(open('omni2_all_years_pickle.p', 'rb') )
print('loaded OMNI2 data')
#######################


#plot example of a variable in the OMNI2 data
#plt.plot_date(otime,o.btot,'-')
#plot omni data for a given time range
#plot_omni(o,'2009-Jan-1','2017-Dec-31',365.25)
#for manipulating plot
#plt.show(block='true')


slice_start='2017-Jan-1'
slice_end='2017-Dec-31'



#slice the data including linear interpolation over NaNs 
#(this is needed for make_dst_from_wind
so=interslice_omni(o,slice_start,slice_end)


plot_omni(so,slice_start,slice_end,20)


print('Calculate Dst')



# Use the function make_dst_from_wind, defined in dst_module.py, to
# calculate Dst with 3 methods
[dst_burton, dst_obrien, dst_temerin_li]=make_dst_from_wind(so.btot,so.bx, \
 so.bygsm, so.bzgsm, so.speed, so.speedx, so.den,so.time)

   
#def make_dst_from_wind(btot_in,bx_in, by_in,bz_in,v_in,vx_in,density_in,time_in):
 #this makes from synthetic or observed solar wind the Dst index	
 #all nans in the input data must be removed prior to function call
 #3 models are calculated: Burton et al., OBrien/McPherron, and Temerin/Li
 #btot_in IMF total field, in nT, GSE or GSM (they are the same)
 #bx_in - the IMF Bx field in nT, GSE or GSM (they are the same)
 #by_in - the IMF By field in nT, GSM
 #bz_in - the IMF Bz field in nT, GSM
 #v_in - the speed in km/s
 #vx_in - the solar wind speed x component (GSE is similar to GSM) in km/s
 #time_in - the time in matplotlib date format 
 
 
 
############################# plot solar wind to Dst conversion results

sns.set_context("talk")     
sns.set_style("darkgrid")  
fig=plt.figure(2,figsize=(10,6))
wide=1
fsize=10

plt.suptitle('Dst prediction from solar wind speed and magnetic field', fontsize=15)

ax1 = fig.add_subplot(411)
plt.plot_date(so.time,so.btot, 'k', linewidth=wide, label='B')
plt.ylabel('B [nT]',fontsize=fsize)
plt.yticks(fontsize=15) 
plt.tick_params(labelbottom=False)

#ax1 = fig.add_subplot(412)
#plt.plot_date(timesi,thetai_deg, 'k', linewidth=wide, label='theta')
#plt.ylabel('IMF clock angle [deg]')

ax1 = fig.add_subplot(412)
plt.plot_date(so.time,so.bzgsm, 'b', linewidth=wide, label='Bz GSM')
plt.ylabel('Bz component [nT]',fontsize=fsize)
plt.tick_params(labelbottom=False)
plt.yticks(fontsize=fsize) 

ax2 = fig.add_subplot(413)
plt.plot_date(so.time,so.speed, 'r', linewidth=wide, label='V')
plt.ylabel('V [km/s]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.tick_params(labelbottom=False)

ax3 = fig.add_subplot(414)
plt.plot_date(so.time,so.dst, 'ok', markersize=wide+1,linewidth=1) #, label='Observed hourly Dst')
#plt.plot_date(so.time,dst_burton, 'b-', linewidth=wide, label='Burton et al. 1975')
#plt.plot_date(so.time,dst_obrien, 'r-', linewidth=wide, label='OBrien & McPherron 2000')
plt.plot_date(so.time,dst_temerin_li, 'g-', linewidth=1, label='Temerin & Li 2002')
plt.legend(loc=3,fontsize=fsize-2)
#ax3.set_ylim([-130,40])
plt.ylabel('Dst index [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.tight_layout()

#save plot
filename='test_coupling.eps'
plt.savefig(filename, dpi=300)
filename='wind_to_dst_coupling_test.png'
plt.savefig(filename)



#################################### Metrics:

print('Metrics for interval ', slice_start,'  to ',slice_end)
print()
print('Absolute difference model to observed Dst for this time interval:')
abserr=(sum((abs(abs(so.dst)-abs(dst_burton))))/len(so.dst))
print('Burton:', round(abserr,1), ' nT')
abserr=(sum((abs(abs(so.dst)-abs(dst_obrien))))/len(so.dst))
print('OBrien:', round(abserr,1), ' nT')
abserr=(sum((abs(abs(so.dst)-abs(dst_temerin_li))))/len(so.dst))
print('TemerinLi:', round(abserr,1), ' nT')

print()
print('RMS error TL model to observed Dst for this time interval: ')
rms=(sum((so.dst-dst_burton)**2)/len(so.dst))**0.5
print('Burton:', round(rms,1), ' nT')
rms=(sum((so.dst-dst_obrien)**2)/len(so.dst))**0.5
print('OBrien:', round(rms,1), ' nT')
rms=(sum((so.dst-dst_temerin_li)**2)/len(so.dst))**0.5
print('TemerinLi:', round(rms,1), ' nT')

print()
print('prediction efficiency:')
residual=so.dst-dst_burton
effic=(1-np.std(residual)**2/np.std(so.dst)**2)*100
print('Burton: ', int(effic), ' %')
residual=so.dst-dst_obrien
effic=(1-np.std(residual)**2/np.std(so.dst)**2)*100
print('OBrien: ', int(effic), ' %')
residual=so.dst-dst_temerin_li
effic=(1-np.std(residual)**2/np.std(so.dst)**2)*100
print('TemerinLi: ', int(effic), ' %')

