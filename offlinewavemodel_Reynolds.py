import os
import sys

from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt
import numpy as np
from netCDF4 import Dataset 

nc = Dataset('/nobackup/1/adarmeno/wave_model/validation/run_c180_merra2/output/umwmout_2017-01-12_10:00:00.nc')

data1read = nc.variables['ust']
ust = np.squeeze(data1read)
data2read = nc.variables['swh']
swh = np.squeeze(data2read)
data3create = np.multiply(ust,swh)
vis  = 1.05e-6 

'''
viscocity at 20 degree C m^2 s-1

'''
Re = np.divide(data3create,vis)
wspd = np.squeeze(nc.variables['wspd'])
data5read = nc.variables['cd']
Cd  = np.squeeze(data5read)


'''
create U10 using Cd and U*

'''
U10 = np.divide(ust, (np.power(Cd,0.5)))
n, bins, patches = plt.hist(U10[~np.isnan(U10)],10)
'''
Gong 2003 dfdDp

print("Calculate white cap fractions for U10 bins"
WM_bins = np.multiply(3.84e-6, np.power(bins,3.41))
tune = 30
Dp = np.arange(0.01,10,0.2)
A = 4.7*(
'''

'''
calculate white cap fraction (Monahan and O'Muircheartaigh, 1980)
'''

WM = np.multiply(3.84e-6, np.power(U10,3.41))

readlat = nc.variables['lat'][:]
readlon = nc.variables['lon'][:]
lat = np.squeeze(readlat)
lon =np.squeeze(readlon)

''' 
plotting

'''

fig1a = plt.figure(figsize=(12, 6), edgecolor='w')
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmapboundary(fill_color='aqua')
map.drawmeridians(np.arange(-180.,181.,60.), labels = [True, False,False,True])
map.drawparallels(np.arange(-90.,91.,30.), labels = [False, True, True,False])
ny = Re.shape[0]
nx = Re.shape[1] 
lons,lats = map.makegrid(nx,ny)
x,y = map(lons,lats)
plt.pcolor(x,y,Re, cmap='RdBu',vmin = 0, vmax = 1e6)
map.colorbar(location='bottom',size="5%",pad="10%")
plt.title("Re" )
plt.savefig('Re.tif')
plt.clf()

fig1b = plt.figure(figsize=(12, 6), edgecolor='w')
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmapboundary(fill_color='aqua')
map.drawmeridians(np.arange(-180.,181.,60.), labels = [True, False,False,True])
map.drawparallels(np.arange(-90.,91.,30.), labels = [False, True, True,False])
plt.pcolor(x,y,np.power(ust,2.41), cmap='RdBu',vmin = 0, vmax = 1e6)
map.colorbar(location='bottom',size="5%",pad="10%")
plt.title("u$\mathregular{_{*}}$$\mathregular{^{2.41}}$(ms$\mathregular{^{-1}}$)" )
plt.clim(0, 0.005)
plt.savefig('Ust.tif')
plt.clf()






fig1c = plt.figure(figsize=(12, 6), edgecolor='w')
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmapboundary(fill_color='aqua')
map.drawmeridians(np.arange(-180.,181.,60.), labels = [True, False,False,True])
map.drawparallels(np.arange(-90.,91.,30.), labels = [False, True, True,False])
plt.pcolor(x,y,WM, cmap='RdBu',vmin = 0, vmax = 1e6)
map.colorbar(location='bottom',size="5%",pad="10%")
plt.title("Whitecap fraction" )
plt.clim(0, 0.05)
plt.savefig('WM.tif')
plt.clf()

fig1d= plt.figure(figsize=(12, 6), edgecolor='w')
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmapboundary(fill_color='aqua')
map.drawmeridians(np.arange(-180.,181.,60.), labels = [True, False,False,True])
map.drawparallels(np.arange(-90.,91.,30.), labels = [False, True, True,False])
plt.pcolor(x,y,swh, cmap='RdBu',vmin = 0, vmax = 1e6)
map.colorbar(location='bottom',size="5%",pad="10%")
plt.title("SWH (m)" )
plt.clim(0,10)
plt.savefig('SWH.tif')
plt.clf()

'''
Plot subset of data for scatterplot of u* Vs Re

'''

lat_bnds, lon_bnds = [15,45], [-80,-30]
lat_inds = np.array(np.logical_and((lat >= lat_bnds[0]) , (lat <= lat_bnds[1])))
lon_inds = np.array(np.logical_and((lon >= lon_bnds[0]) , (lon <= lon_bnds[1])))
latlon_inds = np.logical_and(lat_inds,lon_inds)
dataust_subset = ust[latlon_inds]
dataRe_subset = Re[latlon_inds]
datawspd_subset = wspd[latlon_inds]
dataWM_subset = WM[latlon_inds]


fig2 = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(np.power(dataust_subset,2.41),dataRe_subset,c="g",marker = 'o')
plt.tight_layout()
plt.xlabel(r'u$\mathregular{^{*}}$$\mathregular{^{2.41}}$(ms$\mathregular{^{-1}}$)', fontsize=18)
plt.ylabel(r'R$\mathregular{_{e}}$',fontsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax = plt.gca()

'''windspeed
ax.set_xlim(0,5000)
'''
ax.set_xlim(0,0.2)
ax.set_ylim(0,1e6)
plt.tight_layout()
plt.title('[15,45N];[-80,-30]')
plt.savefig('Ust_Vs_Re[15,45N;-80,-30].tif')
plt.clf()
'''
plot White cap fraction Vs Re
'''

fig3 = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(dataWM_subset,dataRe_subset,c="b",marker = 'o')
plt.tight_layout()
plt.xlabel('W', fontsize=18)
plt.ylabel(r'R$\mathregular{_{e}}$',fontsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax = plt.gca()
ax.set_xlim(0,0.015)
ax.set_ylim(0,1e6)
plt.tight_layout()
plt.title('[15,45N];[-80,-30]')
plt.savefig('W_Vs_Re[15,45N;-80,-30].tif')
plt.clf()
