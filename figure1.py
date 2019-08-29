#!/usr/bin/env python
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import utils
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

debug = True
client = Client()


def setupmap(central_lon, central_lat,handle):
    #handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    handle.set_extent(extent)

    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    #handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES)
    handle.add_feature(cfeature.RIVERS)
    #handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle

net = 'NE,IW'

stime = UTCDateTime('2009-204T00:00:00')
etime = UTCDateTime('2019-204T00:00:00')


fig= plt.figure(figsize=(12,6))

inv = client.get_stations(network=net, station="*",
                        channel = '*HZ', level="response",
                        starttime=stime, endtime = etime)

lats, lons, cols = [], [],[]
for cnet in inv:
    for sta in cnet:
        for chan in sta:
            lats.append(chan.latitude)
            lons.append(chan.longitude)
            if cnet.code =='NE':
                cols.append('C1')
            else:
                cols.append('C2')
             

mlon = np.mean(lons)
mlat = np.mean(lats)
boxcoords=[min(lats) -1., min(lons)-1., max(lats) +1. , max(lons) + 1.]
extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])

#ax1 = plt.subplot(1,2,1, projection=ccrs.Mercator(central_lon))
#ax1 = setupmap(central_lon, central_lat, ax1)
#sc = ax1.scatter(lons, lats, 200., color='C1', edgecolor='C1', marker="v", transform=ccrs.PlateCarree())
    

ax = plt.subplot(1,1,1, projection=ccrs.Mollweide(central_lon))
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.add_feature(cfeature.RIVERS)
ax.add_feature(cfeature.STATES, edgecolor='gray')
ax.set_global()












inv = client.get_stations(network='IW', station="*",
                        channel = '*HZ', level="response",
                        starttime=stime, endtime = etime)

lats, lons= [], []
for cnet in inv:
    for sta in cnet:
        for chan in sta:
            lats.append(chan.latitude)
            lons.append(chan.longitude)

sc = ax.scatter(lons, lats, 200., color='C1', edgecolor='k', marker="v", transform=ccrs.PlateCarree(), label='IW Stations')

mlon = np.mean(lons)
mlat = np.mean(lats)

paramdic = utils.get_parameters('P')
cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], 
                        latitude=mlat, longitude=mlon, maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])
evelats = []
evelons = []
for eve in cat:
    evelats.append(eve.origins[0].latitude)
    evelons.append(eve.origins[0].longitude)

print(evelons)

ax.scatter(evelons, evelats, 65, color="C2", edgecolor="C2", marker="o", zorder=3, alpha=0.5, transform=ccrs.PlateCarree(), label='IW Earthquakes')


inv = client.get_stations(network='NE', station="*",
                        channel = '*HZ', level="response",
                        starttime=stime, endtime = etime)

lats, lons= [], []
for cnet in inv:
    for sta in cnet:
        for chan in sta:
            lats.append(chan.latitude)
            lons.append(chan.longitude)

sc = ax.scatter(lons, lats, 200., color="C3", edgecolor='k', marker="v", transform=ccrs.PlateCarree(), label='NE Stations')

mlon = np.mean(lons)
mlat = np.mean(lats)

paramdic = utils.get_parameters('P')
cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], 
                        latitude=mlat, longitude=mlon, maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])
evelats = []
evelons = []
for eve in cat:
    evelats.append(eve.origins[0].latitude)
    evelons.append(eve.origins[0].longitude)

print(evelons)

ax.scatter(evelons, evelats, 65, color="C4", edgecolor="C4", alpha=0.5, marker="o", zorder=3, transform=ccrs.PlateCarree(), label='NE Earthquake')
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=5)

plt.savefig('figure1.png',format='PNG', dpi=400)
plt.show()

