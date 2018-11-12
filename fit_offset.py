## fit_offset.py

## determines the offset terms in Temerin and Li 2002 for a given timerange in OMNI2 data
## Author: C. Moestl, IWF Graz, Austria
##
## latest update: July 10 2018
## tested in ipython 3.5 with sunpy and seaborn installed


## If results made with this code are used for scientific publications, 
## please contact me before submission: christian.moestl@oeaw.ac.at or twitter @chrisoutofspace


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
import scipy
import math



import dst_module
from dst_module import get_omni_data
from dst_module import plot_omni
from dst_module import make_dst_from_wind
from dst_module import interslice_omni



########################## initialize

#get current directory
os.system('pwd')
#closes all plots
plt.close('all')

 
### FUNCTION FOR SOLAR WIND TO DST CONVERSION optimization

def optimize_offset(sarr):

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

 s1=sarr[0]
 s2=sarr[1]
 s3=sarr[2]


 
 ######## model 3: Xinlin Li LASP Colorado and Mike Temerin
 
  
 #2002 version 
 
 #define all terms
 dst1=np.zeros(len(bz_in))
 dst2=np.zeros(len(bz_in))
 dst3=np.zeros(len(bz_in))
 pressureterm=np.zeros(len(bz_in))
 directterm=np.zeros(len(bz_in))
 offset=np.zeros(len(bz_in))
 dst_temerin_li_out=np.zeros(len(bz_in))
 bp=np.zeros(len(bz_in))
 bt=np.zeros(len(bz_in))
 
 
 #define inital values (needed for convergence, see Temerin and Li 2002 note)
 dst1[0:10]=-15
 dst2[0:10]=-13
 dst3[0:10]=-2
 
 

 #define all constants
 p1=0.9
 p2=2.18e-4
 p3=14.7
 

 
 a1=6.51e-2
 a2=1.37
 a3=8.4e-3  
 a4=6.053e-3
 a5=1.12e-3
 a6=1.55e-3
 
 tau1=0.14 #days
 tau2=0.18 #days
 tau3=9e-2 #days
 
 b1=0.792
 b2=1.326
 b3=1.29e-2
 
 c1=-24.3
 c2=5.2e-2

 #Note: vx has to be used with a positive sign throughout the calculation
 
 
 
 #----------------------------------------- loop over each timestep
 for i in np.arange(1,len(bz_in)-1):

      
  #t time in days since beginning of 1995   #1 Jan 1995 in Julian days
  t1=sunpy.time.julian_day(mdates.num2date(time_in[i]))-sunpy.time.julian_day('1995-1-1 00:00')
  
 
  yearli=365.24 
  tt=2*np.pi*t1/yearli
  ttt=2*np.pi*t1
  alpha=0.078
  beta=1.22
  cosphi=np.sin(tt+alpha)*np.sin(ttt-tt-beta)*(9.58589*1e-2)+np.cos(tt+alpha)*(0.39+0.104528*np.cos(ttt-tt-beta))
 
  #equation 1 use phi from equation 2
  sinphi=(1-cosphi**2)**0.5
  
  pressureterm[i]=(p1*(btot_in[i]**2)+density_in[i]*((p2*((v_in[i])**2)/(sinphi**2.52))+p3))**0.5
  
  #2 directbzterm 
  directterm[i]=0.478*bz_in[i]*(sinphi**11.0)

  #3 offset term - the last two terms were cut because don't make sense as t1 rises extremely for later years
  offset[i]=s1+s2*np.sin(2*np.pi*t1/yearli+s3)
  #or just set it constant
  #offset[i]=-5
  bt[i]=(by_in[i]**2+bz_in[i]**2)**0.5  
  #mistake in 2002 paper - bt is similarly defined as bp (with by bz); but in Temerin and Li's code (dst.pro) bp depends on by and bx
  bp[i]=(by_in[i]**2+bx_in[i]**2)**0.5  
  #contains t1, but in cos and sin 
  dh=bp[i]*np.cos(np.arctan2(bx_in[i],by_in[i])+6.10) * ((3.59e-2)*np.cos(2*np.pi*t1/yearli+0.04)-2.18e-2*np.sin(2*np.pi*t1-1.60))
  theta_li=-(np.arccos(-bz_in[i]/bt[i])-np.pi)/2
  exx=1e-3*abs(vx_in[i])*bt[i]*np.sin(theta_li)**6.1
  #t1 and dt are in days
  dttl=sunpy.time.julian_day(mdates.num2date(time_in[i+1]))-sunpy.time.julian_day(mdates.num2date(time_in[i]))

 
  #4 dst1 
  #find value of dst1(t-tau1) 
  #time_in is in matplotlib format in days: 
  #im time_in den index suchen wo time_in-tau1 am nächsten ist
  #und dann bei dst1 den wert mit dem index nehmen der am nächsten ist, das ist dann dst(t-tau1)
  #wenn index nicht existiert (am anfang) einfach index 0 nehmen
  #check for index where timesi is greater than t minus tau
  
  indtau1=np.where(time_in > (time_in[i]-tau1))
  dst1tau1=dst1[indtau1[0][0]]
  #similar search for others  
  dst2tau1=dst2[indtau1[0][0]]
  th1=0.725*(sinphi**-1.46)
  th2=1.83*(sinphi**-1.46)
  fe1=(-4.96e-3)*  (1+0.28*dh)*  (2*exx+abs(exx-th1)+abs(exx-th2)-th1-th2)*  (abs(vx_in[i])**1.11)*((density_in[i])**0.49)*(sinphi**6.0)
  dst1[i+1]=dst1[i]+  (a1*(-dst1[i])**a2   +fe1*   (1+(a3*dst1tau1+a4*dst2tau1)/(1-a5*dst1tau1-a6*dst2tau1)))  *dttl
  
  #5 dst2    
  indtau2=np.where(time_in > (time_in[i]-tau2))
  dst1tau2=dst1[indtau2[0][0]]
  df2=(-3.85e-8)*(abs(vx_in[i])**1.97)*(btot_in[i]**1.16)*(np.sin(theta_li)**5.7)*((density_in[i])**0.41)*(1+dh)
  fe2=(2.02*1e3)*(sinphi**3.13)*df2/(1-df2)
  dst2[i+1]=dst2[i]+(b1*(-dst2[i])**b2+fe2*(1+(b3*dst1tau2)/(1-b3*dst1tau2)))*dttl
  
  #6 dst3  
  indtau3=np.where(time_in > (time_in[i]-tau3))
  dst3tau3=dst3[indtau3[0][0]]
  df3=-4.75e-6*(abs(vx_in[i])**1.22)*(bt[i]**1.11)*np.sin(theta_li)**5.5*((density_in[i])**0.24)*(1+dh)
  fe3=3.45e3*(sinphi**0.9)*df3/(1-df3)
  dst3[i+1]=dst3[i]+  (c1*dst3[i]   + fe3*(1+(c2*dst3tau3)/(1-c2*dst3tau3)))*dttl
  
   
  #print(dst1[i], dst2[i], dst3[i], pressureterm[i], directterm[i], offset[i])
  #debugging
  #if i == 30: pdb.set_trace()
  #for debugging
  #print()
  #print(dst1[i])
  #print(dst2[i])
  #print(dst3[i])
  #print(pressureterm[i])
  #print(directterm[i])


  #add time delays: ** to do
  
  #The dst1, dst2, dst3, (pressure term), (direct IMF bz term), and (offset terms) are added (after interpolations) with time delays of 7.1, 21.0, 43.4, 2.0, 23.1 and 7.1 min, respectively, for comparison with the ‘‘Kyoto Dst.’’ 

  #dst1
  
  dst_temerin_li_out[i]=dst1[i]+dst2[i]+dst3[i]+pressureterm[i]+directterm[i]+offset[i]
  
  
  
 rms_TL=(sum((dsti-dst_temerin_li_out)**2)/len(dsti))**0.5
 
 print(round(s1,2),round(s2,2),round(s3,2), '      rms: ',round(rms_TL,2))

 return rms_TL









####################################################### MAIN PROGRAM
#closes all plots
plt.close('all')


################################## (1) LOAD OMNI2 DATA

# if not here download OMNI2 data (only needed first time running the program, currently 155 MB)
if not os.path.exists('omni2_all_years.dat'):

  #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
  print('download OMNI2 data from')
  omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
  print(omni2_url)
  try: urllib.request.urlretrieve(omni2_url, 'omni2_all_years.dat')
  except urllib.error.URLError as e:
      print(' ', omni2_url,' ',e.reason)

#if omni2 hourly data is not yet converted and saved as pickle, do it:
if not os.path.exists('omni2_all_years_pickle.p'):
 #load OMNI2 dataset from .dat file with a function from dst_module.py
 o=get_omni_data()
 #contains: o. time,day,hour,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp
 #save for faster loading later
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


################################# (2) Optimize offset terms for a given interval

#test different intervals

#slice_start='2017-Jan-1'
#slice_end='2017-Jan-31'

slice_start='2017-Jan-1'
slice_end='2017-Dec-31'

#slice_start='2017-Jul-1' 
#slice_end='2017-Jul-31'

#slice_start='2017-Dec-1' 
#slice_end='2017-Dec-31'


#slice the data, and linear interpolate over NaNs (this is needed for the function make_dst_from_wind later)
so=interslice_omni(o,slice_start,slice_end)

plot_omni(so,slice_start,slice_end,20)


#these are global variables to be used in the function optimize_offset, 
#containing the Temerin/Li 2002 model

btot_in=so.btot
bx_in=so.bx
by_in=so.bygsm
bz_in=so.bzgsm
v_in=so.speed
vx_in=so.speedx
density_in=so.den
time_in=so.time
dsti=so.dst


#s1 s2 s3
x0=[4.2,5.95,-3.93]

#smin=scipy.optimize.fmin(optimize_offset, x0)
#smin=scipy.optimize.minimize(optimize_offset, x0, method='Powell')



#2017 data: 2.11 -1.22 -3.93       rms:  7.84
#2017 data: 2.11 -1.55 -4.62       rms:  7.79


 
rms=optimize_offset(x0)
print(rms)


######### Use function for calculation
#[dst_burton, dst_obrien, dst_temerin_li]=make_dst_from_wind(btoti,bxi,  bygsmi, bzgsmi, speedi, speedix, deni, timesi)

#optimize s1-4 of offset function

 # these need to be found with a fit for 1-2 years before calculation
 # taken from the TL code:    offset_term_s1 = 6.70       ;formerly named dsto
 #   offset_term_s2 = 0.158       ;formerly hard-coded     2.27 for 1995-1999
 #   offset_term_s3 = -0.94       ;formerly named phasea  -1.11 for 1995-1999
 #   offset_term_s4 = -0.00954    ;formerly hard-coded
 #   offset_term_s5 = 8.159e-6    ;formerly hard-coded
 
 #x0=[6.7, 0.15, -0.94, -5]

 #s1=6.7
 #s2=0.158
 #s3=-0.94
 #set by myself as a constant in the offset term
 #s4=-3
 

 #s4 and s5 as in the TL 2002 paper are not used due to problems with the time
 #s4=-1.054*1e-2
 #s5=8.6e-6

#x0=[6.7, 0.15, -0.9, -3]

#x0=[-2.79, 1.44, -0.92, -7]
