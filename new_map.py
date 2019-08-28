#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import plot_utils
import utils
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager

mpl.rc('font',serif='DejaVu Sans') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

def setupmap(central_lon, central_lat,handle):
    #handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    handle.set_extent(extent)

    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle


net ='N4'
comps = ['Z', 'R']
inv = plot_utils.get_dataless(net)
ax ={}
fig= plt.figure(figsize=(12,12))
for idx, comp in enumerate(comps):
    
    dic_all ={}
    files = glob.glob(net + '_results/Results_' + net + '_' + comp + '*_2019*.csv')

    for curfile in files:
        stuff = plot_utils.parse(curfile)
        for item in stuff:
            try:
                dic_all[item] += stuff[item]
            except:
                dic_all[item] = []
                dic_all[item] += stuff[item]


    print(dic_all['time'])

    key = 'corr'

    lats, lons, staval, parameter, stas = plot_utils.get_plot_params(dic_all, inv, key)
    
    key2 = 'snr'
    
    lats2, lons2, staval2, parameter2, stas2 = plot_utils.get_plot_params(dic_all, inv, key2)

    boxcoords=[min(lats) -1., min(lons)-1., max(lats) +1. , max(lons) + 1.]
    extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
    central_lon = np.mean(extent[:2])
    central_lat = np.mean(extent[2:])
    ax[str(idx)] = plt.subplot(2,1,idx+1, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    ax[str(idx)] = setupmap(central_lon, central_lat, ax[str(idx)])

    sc = ax[str(idx)].scatter(lons, lats,  c = staval,  transform=ccrs.PlateCarree())
    

    
    cbar = fig.colorbar(sc, orientation='vertical', shrink=0.7)
    cbar.ax.set_ylabel(parameter)
    if idx == 0:
        plt.text(min(lons)-10, max(lats)+1, '(a)', fontsize=26, transform=ccrs.PlateCarree())
    else:
        plt.text(min(lons)-10, max(lats)+1, '(b)', fontsize=26, transform=ccrs.PlateCarree())
    #cbar.ax.xaxis.set_label_position('top')
#for trip in zip(lons, lats, stas):
#    print(trip)
#    plt.text(trip[0]+0.1, trip[1]+0.1, trip[2], transform=ccrs.PlateCarree())
#cbar = plt.colorbar()
#cbar.ax.set_ylabel(parameter)


plt.savefig('map_' + key + '_' + net + '.png', format='PNG', dpi=400)
plt.show()
