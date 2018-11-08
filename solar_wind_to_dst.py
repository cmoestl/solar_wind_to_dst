## solar_wind_to_dst.py
##
##
## This is a code for showing an example how to convert
## solar wind observations (here OMNI2 data) to the geomagnetic Dst index
## https://github.com/cmoestl/solar_wind_to_dst
##
## Author: C. Moestl, IWF Graz, Austria
##
## If results made with this code are used for scientific publications, 
## please contact me before submission: christian.moestl@oeaw.ac.at or twitter @chrisoutofspace
##
## latest update: November 2018
## tested in python 3.5.2 with sunpy and seaborn installed
##
## available methods: 
## Burton et al. 1975 doi:10.1029/JA080i031p04204  
## OBrien & McPherron 2000 doi: 10.1029/1998JA000437
## Temerin and Li 2002 doi: 10.1029/2001JA007532
##
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

import predstorm_module

from predstorm_module import get_omni_data
from predstorm_module import convert_omni_time
from predstorm_module import make_dst_from_wind

#from predstorm_module import time_to_num_cat

#from predstorm_module import convert_GSE_to_GSM


############################################################
def plotwind(start, end, tickdays):
 
 plotstartdate=matplotlib.dates.date2num(sunpy.time.parse_time(start))
 plotenddate=matplotlib.dates.date2num(sunpy.time.parse_time(end))
  
     
 #eigenes tick spacing for major ticks
 #jeder 5. Tag 1 tick
 #meinex_majorticks=np.arange(plotstartdate,plotenddate,366) 
 meinex_majorticks=np.arange(plotstartdate,plotenddate,tickdays) 


 fig=plt.figure(1)

 sns.set_context("talk")     
 sns.set_style("darkgrid")  
 
 wide=1.1

 #############
 ax1 = fig.add_subplot(411)
 
 plt.title('Near Earth solar wind and Dst (OMNI2) at 1 hour resolution')
 
 plt.plot_date(times1,np.zeros(dataset),'k--',linewidth=wide)
 plt.plot_date(times1,btot,'k',linewidth=wide, label='B total')
 #plt.plot_date(times1,bx,'r',linewidth=wide,label='Bx')
 #plt.plot_date(times1,by,'g',linewidth=wide, label='By')
 
 plt.yticks(fontsize=10) 
 plt.ylabel('IMF in GSE [nT]')
 plt.ylim((0,np.nanmax(btot)+10)) 
 plt.xlim((plotstartdate, plotenddate))
 plt.legend(loc=1,fontsize=6)
 plt.xticks(meinex_majorticks,fontsize=10,rotation=17) 
 plt.tick_params(labelbottom=False)
 plt.grid()
 
 #################
 ax2 = fig.add_subplot(412)
 plt.plot_date(times1,np.zeros(dataset),'k--',linewidth=wide)
 plt.plot_date(times1,bz,'b',linewidth=wide)
 plt.ylabel('Bz [nT]')
 plt.ylim((np.nanmin(bz)-10,np.nanmax(bz)+10))
 plt.xlim((plotstartdate, plotenddate))
 plt.xticks(meinex_majorticks,fontsize=10,rotation=17) 
 plt.tick_params(labelbottom=False)
 plt.yticks(fontsize=10) 
 plt.grid()
 
 #################
 ax3 = fig.add_subplot(413)
 plt.plot_date(times1,np.zeros(dataset),'k--',linewidth=wide)
 plt.plot_date(times1,speed,'b',linewidth=wide)
 plt.ylabel('V [$\mathrm{km} \mathrm{\; s^{-1}}}$]')
 plt.ylim((np.nanmin(speed)-20,np.nanmax(speed)+100))
 plt.xlim((plotstartdate, plotenddate))
 plt.xticks(meinex_majorticks,fontsize=10,rotation=17) 
 plt.tick_params(labelbottom=False)
 plt.yticks(fontsize=10) 
 plt.grid()
   
 ################## 
 ax4 = fig.add_subplot(414)
 plt.plot_date(times1,np.zeros(dataset),'k--',linewidth=wide)
 plt.plot_date(times1,dst,'k',linewidth=wide,label='Dst')
 
 #set major ticks from own array
 plt.xticks(meinex_majorticks,fontsize=10,rotation=17) 
  
 #dateformat
 #myformat = mdates.DateFormatter('%Y %m %d %H:%M')
 myformat = mdates.DateFormatter('%Y %m %d')
 ax4.xaxis.set_major_formatter(myformat)
 
   
 plt.ylabel('Dst [nT]')
 plt.ylim((np.nanmin(dst)-50,np.nanmax(dst)+20))
 plt.xlim((plotstartdate, plotenddate))
 plt.yticks(fontsize=10) 
 plt.grid()
 
 ####################################################
 
 print('Bmax for selected interval', start, ' to ', end)
 #get total Bmax for selected interval
 
 start_ind=np.where(plotstartdate==times1)[0][0]
 end_ind=np.where(plotenddate==times1)[0][0]
 time_interval=times1[start_ind:end_ind]
 btot_interval=btot[start_ind:end_ind]
 print('Maximum total B',np.nanmax(btot_interval), 'nT')
 bmax_ind=np.nanargmax(btot_interval) #index of max
 timebtotmax=matplotlib.dates.num2date(time_interval[bmax_ind]) #time of minimum
 print('happened on', sunpy.time.break_time(timebtotmax))
 print(' ')



################################### main program ########################################

#closes all plots
plt.close('all')

#if omni2 hourly data is not here, download:
if not os.path.exists('omni2_all_years.dat'):
  #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
  omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
  try: urllib.request.urlretrieve(omni2_url, 'omni2_all_years.dat')
  except urllib.error.URLError as e:
      print(' ', http_sta_pla_file_str[p],' ',e.reason)



#load OMNI2 dataset from .dat file
[year,day,hour,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp]=get_omni_data()
#convert year, day, hour to matplotlib format
otime=convert_omni_time(year,day,hour)

#plot example:
plt.plot_date(otime,btot,'-')
#for manipulating plot
#plt.show(block='true')




#plotwind('1995-Jan-1','2016-Jul-10',365)


########## slice data for comparison of solar wind to Dst conversion

#test time range
start='2015-Jan-20'
ndays=500

#smaller ndays array with hourly values
btoti=np.zeros(24*ndays)
bxi=np.zeros(24*ndays)
bygsmi=np.zeros(24*ndays)
bzgsmi=np.zeros(24*ndays)
speedi=np.zeros(24*ndays)
speedix=np.zeros(24*ndays)
timesi=np.zeros(24*ndays)
dsti=np.zeros(24*ndays)
kpi=np.zeros(24*ndays)

deni=np.zeros(24*ndays)
pdyni=np.zeros(24*ndays)


s=mdates.date2num(sunpy.time.parse_time(start))
#"i" stands for interval
ind=np.where(s==otime)[0][0]
btoti=btot[ind:ind+24*ndays]
bxi=bx[ind:ind+24*ndays]
bygsei=by[ind:ind+24*ndays]
bzgsei=bz[ind:ind+24*ndays]

bygsmi=bygsm[ind:ind+24*ndays]
bzgsmi=bzgsm[ind:ind+24*ndays]
speedi=speed[ind:ind+24*ndays]
speedix=speedx[ind:ind+24*ndays]
deni=den[ind:ind+24*ndays]
timesi=otime[ind:ind+24*ndays]
dsti=dst[ind:ind+24*ndays]
kpi=kp[ind:ind+24*ndays]

#IMF clock angle
thetai=np.arctan2(bygsmi,bzgsmi)
thetai_deg=thetai*180/np.pi

################  linearly interpolate over nan

time_new=np.arange(timesi[0],timesi[-1]+1/24,1/24)

good = np.where(np.isfinite(btoti))
btoti=np.interp(time_new,timesi[good],btoti[good])

good = np.where(np.isfinite(bxi))
bxi=np.interp(time_new,timesi[good],bxi[good])

good = np.where(np.isfinite(bygsmi))
bygsmi=np.interp(time_new,timesi[good],bygsmi[good])

good = np.where(np.isfinite(bzgsmi))
bzgsmi=np.interp(time_new,timesi[good],bzgsmi[good])

good = np.where(np.isfinite(speedi))
speedi=np.interp(time_new,timesi[good],speedi[good])

good = np.where(np.isfinite(speedix))
speedix=np.interp(time_new,timesi[good],speedix[good])

good = np.where(np.isfinite(deni))
deni=np.interp(time_new,timesi[good],deni[good])



######### Use function for calculation
[dst_burton, dst_obrien, dst_temerin_li]=make_dst_from_wind(btoti,bxi,  bygsmi, bzgsmi, speedi, speedix, deni, timesi)



############################# plot solar wind to Dst conversion results

sns.set_context("talk")     
sns.set_style("darkgrid")  
fig=plt.figure(2,figsize=(10,6))
wide=1
fsize=10

plt.suptitle('Dst prediction from solar wind speed and magnetic field', fontsize=15)

ax1 = fig.add_subplot(411)
plt.plot_date(timesi,btoti, 'k', linewidth=wide, label='B')
plt.ylabel('B [nT]',fontsize=fsize)
plt.yticks(fontsize=15) 
plt.tick_params(labelbottom=False)

#ax1 = fig.add_subplot(412)
#plt.plot_date(timesi,thetai_deg, 'k', linewidth=wide, label='theta')
#plt.ylabel('IMF clock angle [deg]')

ax1 = fig.add_subplot(412)
plt.plot_date(timesi,bzgsmi, 'b', linewidth=wide, label='Bz GSM')
plt.ylabel('Bz component [nT]',fontsize=fsize)
plt.tick_params(labelbottom=False)
plt.yticks(fontsize=fsize) 

ax2 = fig.add_subplot(413)
plt.plot_date(timesi,speedi, 'r', linewidth=wide, label='V')
plt.ylabel('V [km/s]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.tick_params(labelbottom=False)

ax3 = fig.add_subplot(414)
plt.plot_date(timesi,dsti, 'ok', markersize=wide+1,linewidth=1) #, label='Observed hourly Dst')
#plt.plot_date(timesi,dst_burton, 'b-', linewidth=wide, label='Burton et al. 1975')
#plt.plot_date(timesi,dst_obrien, 'r-', linewidth=wide, label='OBrien & McPherron 2000')
plt.plot_date(timesi,dst_temerin_li, 'g-', linewidth=1, label='Temerin & Li 2002')
plt.legend(loc=3,fontsize=fsize-2)
#ax3.set_ylim([-130,40])
plt.ylabel('Dst index [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 

plt.tight_layout()


#save plot

#filename='test_coupling.eps'
#plt.savefig(filename, dpi=300)
filename='wind_to_dst_coupling_test.png'
plt.savefig(filename)



#################################### Metrics:

print('Metrics for ',ndays,' days interval starting with ',start)
print()
print('Absolute difference model to observed Dst for this time interval:')
abserr=(sum((abs(abs(dsti)-abs(dst_burton))))/len(dsti))
print('Burton:', round(abserr,1), ' nT')
abserr=(sum((abs(abs(dsti)-abs(dst_obrien))))/len(dsti))
print('OBrien:', round(abserr,1), ' nT')
abserr=(sum((abs(abs(dsti)-abs(dst_temerin_li))))/len(dsti))
print('TemerinLi:', round(abserr,1), ' nT')

print()
print('RMS error TL model to observed Dst for this time interval: ')
rms=(sum((dsti-dst_burton)**2)/len(dsti))**0.5
print('Burton:', round(rms,1), ' nT')
rms=(sum((dsti-dst_obrien)**2)/len(dsti))**0.5
print('OBrien:', round(rms,1), ' nT')
rms=(sum((dsti-dst_temerin_li)**2)/len(dsti))**0.5
print('TemerinLi:', round(rms,1), ' nT')

print()
print('prediction efficiency:')
residual=dsti-dst_burton
effic=(1-np.std(residual)**2/np.std(dsti)**2)*100
print('Burton: ', int(effic), ' %')
residual=dsti-dst_obrien
effic=(1-np.std(residual)**2/np.std(dsti)**2)*100
print('OBrien: ', int(effic), ' %')
residual=dsti-dst_temerin_li
effic=(1-np.std(residual)**2/np.std(dsti)**2)*100
print('TemerinLi: ', int(effic), ' %')

