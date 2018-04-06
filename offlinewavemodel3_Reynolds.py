import os
import sys
import math 
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt
import numpy as np
from netCDF4 import Dataset 
import matplotlib as mpl


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
dwp = np.squeeze(nc.variables['dwp'])
omegap = 2.0*180.0*dwp 



'''
create U10 using Cd and U*

'''
U10 = np.divide(ust, (np.power(Cd,0.5)))

N = 4
j=0
cind = np.array(np.random.rand(U10.shape[0],U10.shape[1]))
cind[np.where(np.logical_and(U10>=0, U10<3))] = 0
cind[np.where(np.logical_and(U10>=3, U10<6))] = 1
cind[np.where(np.logical_and(U10>=6, U10<10))] = 2
cind[np.where(U10>=10)] = 3
cind[np.where(np.isnan(U10))] = 3
cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0,N,N+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)




'''
calculate white cap fraction (Monahan and O'Muircheartaigh, 1980)
'''

WM = np.multiply(3.84e-6, np.power(U10,3.41))
readlat = nc.variables['lat'][:]
readlon = nc.variables['lon'][:]
lat = np.squeeze(readlat)
lon =np.squeeze(readlon)

'''
Aerosol flux parameterization
'''


'''
Limiting steep criterion
'''
# flux is a function of radius, significant wave height, and peak radial frequency
# white cap spray function 
r =  np.linspace(10e-3,12,50)
rad_func = 16.1 - (np.multiply(3.43,np.log10(r)) + np.multiply(2.49, np.power(np.log10(r),2)) + np.multiply(1.2, np.power(np.log10(r),3)))
fprod = np.exp(rad_func)


g = 9.81 

state_ls = np.divide((np.multiply(swh,np.power(omegap,2))),g) 
temp1_ls = -0.1933*-0.3048
state2_ls = np.divide(temp1_ls,np.power(state_ls,2))
state3_ls = np.exp(state2_ls)

countrx = 0
for rx in r:
	FE_LS = np.multiply(fprod[countrx],state3_ls)
	countrx+=1


  

'''
Threshold vertical acceleration
'''
phi = 2 ## change this ; probability integral 
state_tva = np.divide((np.multiply(swh,np.power(omegap,2))),g)
temp1_tva = 0.447
state2_tva = np.multiply(phi, (np.divide(temp1_tva,state_tva)))
state3_tva = 1 - state2_tva
country = 0
for ry in r:
	FE_TVA = np.multiply(fprod[country],state3_tva)
	country+=1


'''Callaghan and white'''

## Wind speed < 10.18 m/s

tau_DWM = 5.3 # timescale for discrete whitecap method
tuning_parm = 30 

U10_gt10 = np.where(U10 > 10.18)
U10_lt23 = np.where(np.logical_and((U10>=10.18),(U10<=23.09)))

## add these in here fopr wind regimes
for ru10_18 in r:
	temp_cw1 = np.power((U10 - 3.7),3) * ru10_18 *tau_DWM
	temp_cw2_1 = np.divide(93.55,temp_cw1)
	temp_cw3_2 =  (1 + (0.057 * np.power(ru10_18,3.45)))
	temp_cw4 = np.exp(-5.33 * np.power((0.433 - np.log10(ru10_18)),2))
	temp_cw5_3 = np.exp(3.68*temp_cw4)
	temp_cw6 = -0.17*ru10_18 - 1.44
	temp_cw7_4  = -4.7 * np.log(np.power((1 + tuning_parm *ru10_18),temp_cw6))
	FE_CW1 = (temp_cw2_1  * temp_cw3_2 * temp_cw5_3) + (temp_cw7_4)
	


## wind speed > 10.18 < 23.09
for ru10_18 in r:
    temp2_cw1 = np.power((U10 + 1.98),3) * ru10_18 *tau_DWM
    temp2_cw2_1 = np.divide(14.18,temp2_cw1)
    temp2_cw3_2 =  (1 + (0.057 * np.power(ru10_18,3.45)))
    temp2_cw4 = np.exp(-5.33 * np.power((0.433 - np.log10(ru10_18)),2))
    temp2_cw5_3 = np.exp(3.68*temp_cw4)
    temp2_cw6 = -0.17*ru10_18 - 1.44
    temp2_cw7_4  = -4.7 * np.log(np.power((1 + tuning_parm *ru10_18),temp_cw6))
    FE_CW2 = (temp2_cw2_1  * temp2_cw3_2 * temp2_cw5_3) + (temp2_cw7_4)

## Petelski 2003 based on dissipation energy 

A = 1.2e-6
B = -1.6e-7
Hs = np.power(swh,0.5) # m
U10 = np.power(U10,0.5)  # m/s
state_pt03 = np.power(np.multiply(Hs,U10),(2.0/3.0)) 

#Emission flux expressed in microgram/s/m2

FE_PT03 = np.multiply(A,state_pt03) + B 

## Ovadnevaite 202014

CMD = [0.018,0.041,0.09,0.23,0.83]
sig = [1.37,1.5,1.42,1.53,1.85]
Fi = [(104.5 * np.power((Re - 1e5),0.556)), (0.0442 * np.power((Re - 1e5),1.08)),(149.6 * np.power((Re - 1e5),0.545)),(2.96 * np.power((Re - 1e5),0.79)),(0.51 * np.power((Re - 2e5),0.87))]
mode = [0,1,2,3,4]

dFdlogD = np.zeros((361,576,50))
rad_count = 0
for rad in r:
   dFdlogD[:,:,rad_count] = 0.0
   for c in mode:
        temp_ovad_1 = np.divide(np.log((rad/CMD[c])),np.log(sig[c]))
        temp_ovad_2 = np.exp(-0.5*(np.power(temp_ovad_1,2)))
        temp_ovad_3 = np.power((2*3.14),0.5) * np.log(sig[c])
        temp_ovad_4 = Fi[c]
        dFdlogD[:,:,rad_count] = dFdlogD[:,:,rad_count] + np.divide((temp_ovad_4 * temp_ovad_2),temp_ovad_3)
   rad_count+=1



 







sys.exit()

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
'''plt.savefig('Re.tif')
'''
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
'''
plt.savefig('Ust.tif')
'''
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
'''
plt.savefig('WM.tif')
'''
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
'''
plt.savefig('SWH.tif')
'''
plt.clf()

'''
Plot subset of data for scatterplot of u* Vs Re

'''

lat_bnds1, lon_bnds1 = [-50, -40],[np.nanmin(lon),np.nanmax(lon)]
lat_bnds2, lon_bnds2 = [-23, 23],[np.nanmin(lon),np.nanmax(lon)]
lat_bnds3, lon_bnds3 = [40, 50],[np.nanmin(lon),np.nanmax(lon)]


'''
South 
'''
lat_inds1 = np.array(np.logical_and((lat >= lat_bnds1[0]) , (lat <= lat_bnds1[1])))
lon_inds1 = np.array(np.logical_and((lon >= lon_bnds1[0]) , (lon <= lon_bnds1[1])))
latlon_inds1 = np.logical_and(lat_inds1,lon_inds1)
dataust1_subset = ust[latlon_inds1]
dataRe1_subset = Re[latlon_inds1]
datawspd1_subset = wspd[latlon_inds1]
dataWM1_subset = WM[latlon_inds1]
colortags1_subset = cind[latlon_inds1]
'''
Equator
'''
lat_inds2 = np.array(np.logical_and((lat >= lat_bnds2[0]) , (lat <= lat_bnds2[1])))
lon_inds2 = np.array(np.logical_and((lon >= lon_bnds2[0]) , (lon <= lon_bnds2[1])))
latlon_inds2 = np.logical_and(lat_inds2,lon_inds2)
dataust2_subset = ust[latlon_inds2]
dataRe2_subset = Re[latlon_inds2]
datawspd2_subset = wspd[latlon_inds2]
dataWM2_subset = WM[latlon_inds2]
colortags2_subset = cind[latlon_inds2]
'''
North
'''
lat_inds3 = np.array(np.logical_and((lat >= lat_bnds3[0]) , (lat <= lat_bnds3[1])))
lon_inds3 = np.array(np.logical_and((lon >= lon_bnds3[0]) , (lon <= lon_bnds3[1])))
latlon_inds3 = np.logical_and(lat_inds3,lon_inds3)
dataust3_subset = ust[latlon_inds3]
dataRe3_subset = Re[latlon_inds3]
datawspd3_subset = wspd[latlon_inds3]
dataWM3_subset = WM[latlon_inds3]
colortags3_subset = cind[latlon_inds3]
'''
Plotting
'''
fig2a = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(np.power(dataust1_subset,2.41),dataRe1_subset,c=colortags1_subset,cmap=cmap,norm=norm)
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
plt.title('Zonal([50S ; 40S])')
plt.savefig('South_zonal_UstVsRe.tif')
plt.clf()


fig2b = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(np.power(dataust2_subset,2.41),dataRe2_subset,c=colortags2_subset,cmap=cmap,norm=norm)
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
plt.title('Zonal([23S; 23N])')
plt.savefig('Equator_zonal_UstVsRe.tif')
plt.clf()

fig2c = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(np.power(dataust3_subset,2.41),dataRe3_subset,c=colortags3_subset,cmap=cmap,norm=norm)
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
plt.title('Zonal([40N;50N])')
plt.savefig('North_zonal_UstVsRe.tif')
plt.clf()





'''
plot White cap fraction Vs Re
'''

fig3a = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(dataWM1_subset,dataRe1_subset,c=colortags1_subset,cmap=cmap,norm=norm)
plt.tight_layout()
plt.xlabel('W', fontsize=18)
plt.ylabel(r'R$\mathregular{_{e}}$',fontsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax = plt.gca()

ax.set_xlim(0,0.015)
ax.set_ylim(0,1e6)
plt.tight_layout()
plt.title('Zonal[(50S ; 40S)]')
plt.savefig('South_zonal_WMVsRe.tif')
plt.clf()

fig3b = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(dataWM2_subset,dataRe2_subset,c=colortags2_subset,cmap=cmap,norm=norm)
plt.tight_layout()
plt.xlabel('W', fontsize=18)
plt.ylabel(r'R$\mathregular{_{e}}$',fontsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax = plt.gca()
ax.set_xlim(0,0.015)
ax.set_ylim(0,1e6)
plt.tight_layout()
plt.title('Zonal[(23S ; 23N)]')
plt.savefig('Equator_zonal_WMVsVsRe.tif')
plt.clf()

fig3c = plt.figure(figsize=(8, 6), edgecolor = 'w')
plt.scatter(dataWM3_subset,dataRe3_subset,c=colortags3_subset,cmap=cmap,norm=norm)
plt.tight_layout()
plt.xlabel('W', fontsize=18)
plt.ylabel(r'R$\mathregular{_{e}}$',fontsize=18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
ax = plt.gca()
ax.set_xlim(0,0.015)
ax.set_ylim(0,1e6)
plt.tight_layout()
plt.title('Zonal[(40N ; 50N)]')
plt.savefig('North_zonal_WMVsRe.tif')
plt.clf()

