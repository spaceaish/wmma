# This script reads the NDBC data using python panda and stores it in an array 
import sys
filename = sys.argv[1]
print(filename)
import os
import numpy as np
import pandas as pd
from io import StringIO
from numpy import genfromtxt

headers = ('station_id',
           'sensor_id',
           'latitude',
           'longitude',
           'date_time',
           'sea_surface_wave_significant_height',
           'sea_surface_wave_peak_period',
           'sea_surface_wave_mean_period',
           'sea_surface_swell_wave_significant_height',
           'sea_surface_swell_wave_period',
           'sea_surface_wind_wave_significant_height',
           'sea_surface_wind_wave_period',
           'sea_water_temperature',
           'sea_surface_wave_to_direction',
           'sea_surface_swell_wave_to_direction',
           'sea_surface_wind_wave_to_direction',
           'number_of_frequencies',
           'center_frequencies',
           'bandwidths',
           'spectral_energy',
           'mean_wave_direction',
           'principal_wave_direction',
           'polar_coordinate_r1',
           'polar_coordinate_r2',
           'calculation_method',
           'sampling_rate',)

# filename = 'waves.csv' 

df = pd.read_csv(filename, sep=',',skipinitialspace=True, index_col=False, names = headers,header = 1, dtype = 'S')

# station id
sid= df.station_id
count = 0 
sid3_list = []

for s in sid:
	sid2 = str(s)
	findwmo = sid2.find("wmo:")
	sidstart = findwmo + 4
	sid3 = (sid2[sidstart:sidstart+5])

	sid3_list.append(sid3) 
# print(sid3_list)

# sensor id
sen_id = df.sensor_id

# location
lat = pd.to_numeric(df.latitude, errors = 'coerce').fillna(0)
lon = pd.to_numeric(df.longitude, errors = 'coerce').fillna(0)

datetime = df.date_time
# convert datetime to a number
dt_list = []
for d in datetime:
	dt = str(d)
	dt_temp1 = dt[0:4]
	dt_temp2 = dt[5:7]
	dt_temp3 = dt[8:10]
	dt_temp4 = dt[11:13]
	dt_temp5 = dt[14:16]
	dt_temp6 = dt[17:19]
 	
	dt_list.append(dt_temp1 + dt_temp2 + dt_temp3 + dt_temp4 + dt_temp5 + dt_temp6)
dtval = pd.to_numeric(dt_list, errors = 'coerce')

# wave params
wsighgt = pd.to_numeric(df.sea_surface_wave_significant_height, errors = 'coerce').fillna(0)

wpeakperiod = pd.to_numeric(df.sea_surface_wave_peak_period, errors = 'coerce').fillna(0)

wmeanperiod = pd.to_numeric(df.sea_surface_wave_mean_period, errors = 'coerce').fillna(0)

swell_wsighgt = pd.to_numeric(df.sea_surface_swell_wave_significant_height, errors = 'coerce').fillna(0)	

swell_wperiod = pd.to_numeric(df.sea_surface_swell_wave_period, errors = 'coerce').fillna(0)


wnd_wsighgt = pd.to_numeric(df.sea_surface_wind_wave_significant_height, errors = 'coerce').fillna(0) 


wnd_wperiod = pd.to_numeric(df.sea_surface_wind_wave_period, errors = 'coerce').fillna(0)

sea_water_temp = pd.to_numeric(df.sea_water_temperature, errors = 'coerce').fillna(0)

sea_wave_dir = pd.to_numeric(df.sea_surface_wave_to_direction, errors = 'coerce').fillna(0)

swell_sea_wave_dir = pd.to_numeric(df.sea_surface_swell_wave_to_direction, errors = 'coerce').fillna(0)

sea_wnd_wave_dir = pd.to_numeric(df.sea_surface_wind_wave_to_direction, errors = 'coerce').fillna(0)

num_freq = pd.to_numeric(df.number_of_frequencies, errors = 'coerce').fillna(0)

center_freq = pd.to_numeric(df.center_frequencies, errors = 'coerce').fillna(0)

bandwidths = pd.to_numeric(df.bandwidths, errors = 'coerce').fillna(0)


spec_energy = pd.to_numeric(df.spectral_energy, errors = 'coerce').fillna(0)

mean_wave_dir = pd.to_numeric(df.mean_wave_direction, errors = 'coerce').fillna(0)

prin_wave_dir = pd.to_numeric(df.principal_wave_direction, errors = 'coerce').fillna(0)

polar_coordinate_r1 = pd.to_numeric(df.polar_coordinate_r1, errors = 'coerce').fillna(0)

polar_coordinate_r2 = pd.to_numeric(df.polar_coordinate_r2, errors = 'coerce').fillna(0)


cal_method = df.calculation_method

samp_rate = df.sampling_rate



# something to print
print('----------------------------')
print('Data testing') 
print('----------------------------')
print('Latitude')
print(lat[1:10])
print('Longitude')
print(lon[1:10])
print('station_id')
print(sid3_list[1:10])
print('datetime')
print(dtval[1:10])
print('wave significant height')
print(wsighgt[1:10])
print('NDBC buoy data read from % s; All done!', filename)
