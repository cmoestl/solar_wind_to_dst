## solar_wind_to_dst.py

## available methods: Burton et al. 1975, OBrien & McPherron 2000, and Temerin and Li 2002 
## Author: C. Moestl, IWF Graz, Austria
##
## latest update: July 10 2018
## tested in ipython 3.5 with sunpy and seaborn installed

## Copyright <2018> <Christian Moestl> MIT LICENSE
## Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## If results made with this code are used for scientific publications, 
## please contact me before submission: christian.moestl@oeaw.ac.at or twitter @chrisoutofspace


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


########################## initialize

#get current directory
os.system('pwd')
#closes all plots
plt.close('all')

#define global variables from OMNI2 dataset
#see http://omniweb.gsfc.nasa.gov/html/ow_data.html

dataset=473376;
#global Variables
spot=np.zeros(dataset) 
btot=np.zeros(dataset) #floating points
bx=np.zeros(dataset) #floating points
by=np.zeros(dataset) #floating points
bz=np.zeros(dataset) #floating points
bzgsm=np.zeros(dataset) #floating points
bygsm=np.zeros(dataset) #floating points

speed=np.zeros(dataset) #floating points
speedx=np.zeros(dataset) #floating points
speed_phi=np.zeros(dataset) #floating points
speed_theta=np.zeros(dataset) #floating points

dst=np.zeros(dataset) #float
kp=np.zeros(dataset) #float

den=np.zeros(dataset) #float
pdyn=np.zeros(dataset) #float
year=np.zeros(dataset)
day=np.zeros(dataset)
hour=np.zeros(dataset)
t=np.zeros(dataset) #index time
times1=np.zeros(dataset) #datetime time



#############################################################   
#reads OMNI2 data from attached savefile:
def getdata():

 #FORMAT(2I4,I3,I5,2I3,2I4,14F6.1,F9.0,F6.1,F6.0,2F6.1,F6.3,F6.2, F9.0,F6.1,F6.0,2F6.1,F6.3,2F7.2,F6.1,I3,I4,I6,I5,F10.2,5F9.2,I3,I4,2F6.1,2I6,F5.1)
 #1963   1  0 1771 99 99 999 999 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 9999999. 999.9 9999. 999.9 999.9 9.999 99.99 9999999. 999.9 9999. 999.9 999.9 9.999 999.99 999.99 999.9  7  23    -6  119 999999.99 99999.99 99999.99 99999.99 99999.99 99999.99  0   3 999.9 999.9 99999 99999 99.9
 
 j=0
 print('start reading variables from file')
 with open('omni2_all_years.dat') as f:
  for line in f:
   line = line.split() # to deal with blank 
   #print line #41 is Dst index, in nT
   dst[j]=line[40]
   kp[j]=line[38]
   
   if dst[j] == 99999: dst[j]=np.NaN
   #40 is sunspot number
   spot[j]=line[39]
   #if spot[j] == 999: spot[j]=NaN

   #25 is bulkspeed F6.0, in km/s
   speed[j]=line[24]
   if speed[j] == 9999: speed[j]=np.NaN
 
   #get speed angles F6.1
   speed_phi[j]=line[25]
   if speed_phi[j] == 999.9: speed_phi[j]=np.NaN

   speed_theta[j]=line[26]
   if speed_theta[j] == 999.9: speed_theta[j]=np.NaN
   #convert speed to GSE x see OMNI website footnote
   speedx[j] = - speed[j] * np.cos(np.radians(speed_theta[j])) * np.cos(np.radians(speed_phi[j]))

   #9 is total B  F6.1 also fill ist 999.9, in nT
   btot[j]=line[9]
   if btot[j] == 999.9: btot[j]=np.NaN

   #GSE components from 13 to 15, so 12 to 14 index, in nT
   bx[j]=line[12]
   if bx[j] == 999.9: bx[j]=np.NaN
   by[j]=line[13]
   if by[j] == 999.9: by[j]=np.NaN
   bz[j]=line[14]
   if bz[j] == 999.9: bz[j]=np.NaN
 
   #GSM
   bygsm[j]=line[15]
   if bygsm[j] == 999.9: bygsm[j]=np.NaN
 
   bzgsm[j]=line[16]
   if bzgsm[j] == 999.9: bzgsm[j]=np.NaN 	
 
   #24 in file, index 23 proton density /ccm
   den[j]=line[23]
   if den[j] == 999.9: den[j]=np.NaN
 
   #29 in file, index 28 Pdyn, F6.2, fill values sind 99.99, in nPa
   pdyn[j]=line[28]
   if pdyn[j] == 99.99: pdyn[j]=np.NaN 		
 
   year[j]=line[0]
   day[j]=line[1]
   hour[j]=line[2]
   j=j+1     


 print('done reading variables from file')
 print(j, ' datapoints')   #for reading data from OMNI file
 
 
 
############################################################# 
def converttime(): #to matplotlib format

 #http://docs.sunpy.org/en/latest/guide/time.html
 #http://matplotlib.org/examples/pylab_examples/date_demo2.html

 print('convert time start')
 for index in range(0,dataset):
      #first to datetimeobject 
      timedum=datetime.datetime(int(year[index]), 1, 1) + datetime.timedelta(day[index] - 1) +datetime.timedelta(hours=hour[index])
      #then to matlibplot dateformat:
      times1[index] = matplotlib.dates.date2num(timedum)
      #print time
      #print year[index], day[index], hour[index]
 print('convert time done')   #for time conversion




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



 
######################################################################################## 
###  THIS IS THE MAIN FUNCTION FOR SOLAR WIND TO DST CONVERSION  
########################################################################################
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
  
  

  


  #print(dst_temerin_li_out[i])
 
  #---------------- loop over
 #plt.plot(dst1)
 #plt.plot(dst2)
 #plt.plot(dst3)
 #plt.plot(pressureterm)
 #plt.plot(directterm)
 #plt.plot(offset) 
  
 rms_TL=(sum((dsti-dst_temerin_li_out)**2)/len(dsti))**0.5
 
 print(round(s1,2),round(s2,2),round(s3,2), '      rms: ',round(rms_TL,2))

 return rms_TL
   
   












####################################################### MAIN PROGRAM


########################################## get data
#read in data from omni file -> 1 , from save_file -> 0
data_from_omni_file = 0 #

if data_from_omni_file == 1:
 getdata()
 converttime()
 pickle.dump([spot,btot,bx,by,bz,bygsm,bzgsm,speed,speedx, dst,kp, den,pdyn,year,day,hour,times1], open( "omni2save.p", "wb" ) ) 
else: [spot,btot,bx,by,bz,bygsm, bzgsm,speed,speedx, dst,kp,den,pdyn,year,day,hour,times1]= pickle.load( open( "omni2save.p", "rb" ) )


###################################### for plotting a specific interval of the data in an extra plot

#plotwind('1995-Jan-1','2016-Jul-10',365)


################################## slice data for comparison of solar wind to Dst conversion

#test time range
start='2015-Jan-20'
ndays=300

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
ind=np.where(s==times1)[0][0]
btoti=btot[ind:ind+24*ndays]
bxi=bx[ind:ind+24*ndays]
bygsei=by[ind:ind+24*ndays]
bzgsei=bz[ind:ind+24*ndays]

bygsmi=bygsm[ind:ind+24*ndays]
bzgsmi=bzgsm[ind:ind+24*ndays]
speedi=speed[ind:ind+24*ndays]
speedix=speedx[ind:ind+24*ndays]
deni=den[ind:ind+24*ndays]
timesi=times1[ind:ind+24*ndays]
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


btot_in=btoti
bx_in=bxi
by_in=bygsmi
bz_in=bzgsmi
v_in=speedi
vx_in=speedix
density_in=deni
time_in=timesi


#define variables
#dynamic pressure
pdyn1=np.zeros(len(bz_in))
protonmass=1.6726219*1e-27  #kg
#assume pdyn is only due to protons
pdyn1=density_in*1e6*protonmass*(v_in*1e3)**2*1e9  #in nanoPascal

 





#x0=[6.7, 0.15, -0.9, -3]

#x0=[-2.79, 1.44, -0.92, -7]

x0=[4.2,5.95,-3.93]

#smin=scipy.optimize.fmin(optimize_offset, x0)
smin=scipy.optimize.minimize(optimize_offset, x0, method='Powell')

sys.exit()

 
rms=optimize_offset(x0)
print(rms)

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

