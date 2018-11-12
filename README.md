# Conversion of solar wind to Dst index

In space weather prediction, it is difficult to infer geomagnetic indices such as Dst or Kp 
from knowledge of the solar wind at the Sun-Earth L1 point. This package provides 3 methods
for achieving this for Dst (Kp will be added in the future), namely: Temerin and Li (2002), 
O'Brien and McPherron (2000), and Burton et al. (1975).


Status: work in progress, November 2018 / issue with offsets needs to be resolved (fit_offset.py)


If you plan to use this code for generating results for 
peer-reviewed scientific publications, please contact me (see bio).


## Dependencies

* python 3 anaconda, non-standard packages: sunpy, seaborn
* I use python 3.5.2 / IPython 4.2.0 on MacOSX Mojave 

## Running the code

* After cloning the repository, run on the command line:

 $ python solar_wind_to_dst.py
 
* or in ipython
 
 $ run solar_wind_to_dst

* The code will download OMNI2 hourly data (158 MB currently) into the main folder, 
convert the time to the matplotlib time format, and save the data as a numpy recarray in a python pickle file, 
so when the structure 'omni' contains the data, the variables can be used as 
'omni.time','omni.btot', 'omni.speed' etc.

* The main function to convert a given solar wind to dst is 'make_dst_from_wind' in the dst_module.py file

* The main program solar_wind_to_dst.py creates a plot of the solar wind, the observed Dst and a Dst calculated from the solar wind 
for a time interval selected in solar_wind_to_dst.py
