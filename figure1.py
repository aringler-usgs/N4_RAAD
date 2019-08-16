#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import numpy as np
import utils



import matplotlib.pyplot as plt
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

debug = True
client = Client()

net = "IW"

stime = UTCDateTime('2018-001T00:00:00')
etime = UTCDateTime('2019-001T00:00:00')



def plot_stations(net):
    inv = client.get_stations(network=net, station="*",
                            channel = '*HZ', level="response",
                            starttime=stime, endtime = etime)

    lats, lons = [], []


    for net in inv:
        for sta in net:
            for chan in sta:
                coors = inv.get_coordinates(net.code + '.' + sta.code + '.' + chan.location_code + '.' + chan.code)
                lats.append(coors['latitude'])
                lons.append(coors['longitude'])
                break
            


    mlon = np.mean(lons)
    mlat = np.mean(lats)


    m = Basemap(projection='merc', lon_0=mlon, lat_0=mlat, llcrnrlon=mlon-7, 
                llcrnrlat=mlat-7, urcrnrlon=mlon+7, urcrnrlat=mlat+7, epsg=4269, resolution="i")
    m.drawcountries(linewidth=1., zorder=3.)
    m.fillcontinents(color='.9', lake_color='C0')
    m.drawmapboundary(fill_color='C0')
    m.drawstates(linewidth=1., zorder=3.)
    ax=plt.gca()
    x,y =m(lons, lats)
    m.scatter(x,y, 200, color="C1", edgecolor="C1", marker="v", zorder=3)
    return
    
def plot_events(net):
    inv = client.get_stations(network=net, station="*",
                            channel = '*HZ', level="response",
                            starttime=stime, endtime = etime)

    lats, lons = [], []


    for net in inv:
        for sta in net:
            for chan in sta:
                coors = inv.get_coordinates(net.code + '.' + sta.code + '.' + chan.location_code + '.' + chan.code)
                lats.append(coors['latitude'])
                lons.append(coors['longitude'])
                break
            


    mlon = np.mean(lons)
    mlat = np.mean(lats)
    m = Basemap(projection='ortho',lon_0=mlon,lat_0=mlat,resolution='l')
    m.drawcoastlines()
    m.fillcontinents(color='.9' ,lake_color='C0')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    m.drawmapboundary(fill_color='C0')
    x,y =m(lons, lats)
    m.scatter(x,y, 100, color="C1", edgecolor="C1", marker="v", zorder=3)

    paramdic = utils.get_parameters('P')

    cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], 
                        latitude=mlat, longitude=mlon, maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])
    evelats = []
    evelons = []
    for eve in cat:
        evelats.append(eve.origins[0].latitude)
        evelons.append(eve.origins[0].longitude)
    x, y = m(evelons, evelats)
    m.scatter(x, y, 50, color="C2", edgecolor="C2", marker="o", zorder=3)
    
    
    
    
    return    



fig = plt.figure(1, figsize=(12,12))

plt.subplot(2,2,1)
plot_stations("IW")
plt.title('Intermountain West (IW) Network')
plt.annotate('(a)', xy=(-0.1, 1.1), xycoords='axes fraction')
plt.subplot(2,2,2)
plot_events("IW")
plt.title('Events used for IW')
plt.annotate('(b)', xy=(-0.1, 1.1), xycoords='axes fraction')
plt.subplot(2,2,3)
plot_stations("NE")
plt.title('New England (NE) Network')
plt.annotate('(c)', xy=(-0.1, 1.1), xycoords='axes fraction')
plt.subplot(2,2,4)
plot_events("NE")
plt.title('Events used for NE')
plt.annotate('(d)', xy=(-0.1, 1.1), xycoords='axes fraction')


plt.savefig('figure1.png', format='PNG', dpi=400)
