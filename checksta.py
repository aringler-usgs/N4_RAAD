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

def get_data(nets, eve, phase, debug = True):
    st = Stream()
    for net in nets:
        for sta in net:
            for chan in sta:
                # Here we need to deal with different channel codes e.g. BH and HH should probably get compared
                if chan.code == 'LHZ':
                    if debug:
                        print(sta.code)
                    coors = nets.get_coordinates(net.code + '.' + sta.code + '.' + chan.location_code + '.' + chan.code )
                    if debug:
                        print(coors)
        
                    ## Time to ray trace
                    (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], eve.origins[0].latitude,eve.origins[0].longitude)
                    ## This might be a distance in m.
                    if phase == 'R':
                        # Using surface waves
                        dis /= 1000.
                        arrivtime = dis / 5.5
                        winlen = dis / 3.5
                        pstime = (eve.origins[0].time + arrivtime)
                        petime = (eve.origins[0].time + winlen)
                        
                    else:
                        dis = kilometer2degrees(dis/1000.)
                        arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=dis, phase_list=[phase])
                        if debug:
                            print(arrivals)
                        for arrival in arrivals:

                            winlen = 100.
                        
                            pstime = (eve.origins[0].time + arrival.time)-5.
                            petime = pstime + winlen
                    if debug:
                        print(net.code + '.' + sta.code + chan.location_code + '.' + chan.code)
                        print(pstime)
                        print(petime) 
                    try:
                    #if True:
                        st += client.get_waveforms(net.code, sta.code, chan.location_code, chan.code, pstime, petime, attach_response = True)
                    except:
                        print('No data for: ' + net.code + '.' + sta.code + chan.location_code + '.' + chan.code)
        return st





debug = True
# We want to N4 check this station via the Rebecca, Andrew, Adam method
net, staGOOD = 'IU', 'COR'

client = Client("IRIS")

stime = UTCDateTime('2018-001T00:00:00')
etime = UTCDateTime('2019-150T00:00:00')

# Need to correct channels for BH and HH at the stame time.  We also need to remove LH

inv = client.get_stations(network=net, station=staGOOD,
                            channel = 'LHZ', level="response",
                            starttime=stime, endtime = etime)
if debug:
    print(inv[0][0][0].code)

# We know where our station is 
coors = inv.get_coordinates(inv[0].code + '.' + inv[0][0].code + '.' + inv[0][0][0].location_code + '.' + inv[0][0][0].code)

nets = client.get_stations(stime, etime, channel="*HZ",
            longitude=coors['longitude'], latitude=coors['latitude'], maxradius = 2., level="channel")

if debug:
    print(nets) 

# We will want to play with this to figure out the correct events for the phase of interest

# We should pick our event parameters for the correct phase e.g. shallow for rayleigh

cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=6.5, latitude=coors['latitude'], 
                        longitude=coors['longitude'], maxradius=85., minradius = 30.)
  


times, amps, corrs = [], [], []
for eve in cat:
    print(eve)
    # We have our events we have our stations and our test station
    
    st = get_data(nets, eve, 'P')   
        
    if debug:
        print(st)
    sys.exit()
    # Now we have the data so lets cross-correlate and stack
    #st.plot()
    st.detrend('constant')
    st.taper(0.05)
    
    st.remove_response()
    st.filter('bandpass',freqmin=0.01, freqmax=0.5)
    
    print(st)
    trRef = st.select(station=staGOOD)
    print(trRef)
    trRef = trRef[0]
    
    # Need to add an off/on flag for plotting events
    #fig = plt.figure(1, figsize=(8,8))
    for tr in st:
        idx, val = xcorr(tr, trRef, 20)
        if debug:
            print(idx)
            print(val)
        t = np.arange(0, tr.stats.npts)/tr.stats.sampling_rate
        plt.plot(t - float(idx)/40., tr.data, label= tr.id + ' ' + str(val))
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (m/s)')
        tr.stats.starttime -= float(idx)/100.
        if 'stack' not in vars():
            stack = tr.data
        else:
            stack += tr.data
    stack /= float(len(st))
    
    idx, val = xcorr(trRef, stack, 20)
    amp = np.sqrt(np.sum(trRef.data**2)/np.sum(stack**2))
    amps.append(amp)
    times.append(float(idx)/100.)
    corrs.append(val)
    plt.plot(t, stack, color = 'C1', linewidth=3, label ='Stack')
    plt.xlim((min(t), max(t)))
    plt.legend()
    plt.savefig('FULL_CODA.jpg', format='JPEG')
    plt.show()
    del stack
    
print(amps)
print(times)
print(corrs)
