#!/usr/bin/env python
from obspy.core import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
import matplotlib.pyplot as plt
import sys
import numpy as np
from obspy.signal.cross_correlation import xcorr
model = TauPyModel(model="iasp91")

import matplotlib as mpl
#mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)



debug = True
# We want to N4 check this station via the Rebecca, Andrew, Adam method
net, staGOOD = 'N4', 'Z38B'

client = Client("IRIS")

stime = UTCDateTime('2019-001T00:00:00')
etime = stime + 30*24*60*60


inv = client.get_stations(network=net, station=staGOOD,
                            channel = 'HHZ', level="response",
                            starttime=stime, endtime = etime)
if debug:
    print(inv)

# We know where our station is 
coors = inv.get_coordinates(net + '.' + staGOOD + '..HHZ')

stas = client.get_stations(stime, etime, network = net, channel="*HZ",
            longitude=coors['longitude'], latitude=coors['latitude'], maxradius = 3., level="channel")[0]

if debug:
    print(stas) 

# We will want to play with this to figure out the correct events for the phase of interest
cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=6.0, latitude=coors['latitude'], 
                        longitude=coors['longitude'], maxradius=85., minradius = 30.)
                        
for eve in cat:
    # We have our events we have our stations and our test station
    st = Stream()
    for sta in stas:
        for chan in sta:

            if chan.code == 'BHZ':
                if debug:
                    print(sta.code)
                coors = stas.get_coordinates(net + '.' + sta.code +'..HHZ')
                if debug:
                    print(coors)
    
                ## Time to ray trace
                (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], eve.origins[0].latitude,eve.origins[0].longitude)
                ## This might be a distance in m.
                dis = kilometer2degrees(dis/1000.)
                arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=dis, phase_list=['P'])
                for arrival in arrivals:

                    winlen = 20.
                
                    pstime = (eve.origins[0].time + arrival.time)-5.
                    petime = pstime + winlen
                    if debug:
                        print(net + ' ' + sta.code + ' ' + chan.code)
                        print(pstime)
                        print(petime) 
                    try:
                        st += client.get_waveforms(net, sta.code, '*', chan.code, pstime, petime, attach_response = True)
                    except:
                        continue
    
    
    # Now we have the data so lets cross-correlate and stack
    #st.plot()
    st.detrend('constant')
    st.taper(0.05)
    
    st.remove_response()
    st.filter('bandpass',freqmin=0.1, freqmax=5.)
    trRef = st.select(station=staGOOD)[0]
    
    
    fig = plt.figure(1)
    for tr in st:
        idx, val = xcorr(tr, trRef, 50)
        if debug:
            print(idx)
            print(val)
        t = np.arange(0, tr.stats.npts)/tr.stats.sampling_rate
        plt.plot(t - float(idx)/40., tr.data, label= tr.id + ' ' + str(val))
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (m/s)')
        tr.stats.starttime -= float(idx)/40.
        if 'stack' not in vars():
            stack = tr.data
        else:
            stack += tr.data
    stack /= float(len(st))
    plt.plot(t, stack, color = 'C1', linewidth=3, label ='Stack')
    plt.xlim((min(t), max(t)))
    plt.legend()
    plt.savefig('FULL_CODA.jpg', format='JPEG')
    plt.show()
    del stack
    sys.exit()
