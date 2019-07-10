#!/usr/bin/env python
import os
from obspy.core import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
import matplotlib.pyplot as plt
import csv
import sys
import numpy as np
import matplotlib as mpl
from obspy.signal.cross_correlation import xcorr
model = TauPyModel(model="iasp91")


#mpl.rc('font',family='serif')
#mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
#mpl.rc('font',size=14)

#added site_dis as an argument since for R phase, winlens were varied, made dependent on only the reference site.
# make winlen and site_dis not required 
def get_data(nets, eve, phase, winlen, dis, debug = True):
    dis =  dis/1000
    st = Stream()
    for net in nets:
        for sta in net:
            for chan in sta:
                # Here we need to deal with different channel codes e.g. BH and HH should probably get compared
                if chan.code == 'HHZ':
                    if debug:
                        print(sta.code)
                    # grab station coordinates
                    coors = nets.get_coordinates(net.code + '.' + sta.code + '.' + chan.location_code + '.' + chan.code )
                    if debug:
                        print(coors)
                    ## Time to ray trace
                    ## (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], eve.origins[0].latitude,eve.origins[0].longitude)
                    ## This might be a distance in m.
                    if phase == 'R':
                        # Using surface waves
                        # is the assumption that the R speed falls between 3.5 and 5.5
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
                            # do we need to use P-S time to decide window length
                            # winlen = 50.
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



####################
# 1. What are the best frequencies bands for P, S, and Raleigh Waves?
# 2. What are the best parameters for earthquakes for P, S, and Raleigh (e.g depth restrictions, st-eq distance)
# 3. Test maxr (st-st radius) to see where program breaks down for P, S, Raleigh?
###################


# extra info --> command line argument
debug = True

# show plot
plot = True

# station --> command line argument
net, staGOOD = 'N4', 'O52A'

client = Client("IRIS")

# start and end times --> command line argument
stime = UTCDateTime('2019-001T00:00:00')
etime = UTCDateTime('2019-150T00:00:00')

# Need to correct channels for BH and HH at the stame time.  We also need to remove LH

inv = client.get_stations(network=net, station=staGOOD,
                            channel = 'HHZ', level="response",
                            starttime=stime, endtime = etime)
if debug:
    print(inv[0][0][0].code)

# grab test station coordinates
coors = inv.get_coordinates(inv[0].code + '.' + inv[0][0].code + '.' + inv[0][0][0].location_code + '.' + inv[0][0][0].code)

# max station-station degree distance
maxr = 1.5

# get all stations that fall within +/- maxr of staGOOD
try:
    nets = client.get_stations(starttime=stime, endtime=etime, channel="HHZ",
                longitude=coors['longitude'], latitude=coors['latitude'],
                maxradius = maxr, level="channel", network="N4")
except:
    sys.exit('No stations nearby {}_{}'.format(net, staGOOD))

if debug:
    print(nets) 

# We will want to play with this to figure out the correct events for the phase of interest

# We should pick our event parameters for the correct phase e.g. shallow for rayleigh
try:
    cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=6.0, maxmagnitude=7.5, latitude=coors['latitude'], 
                            longitude=coors['longitude'], maxradius=85., minradius = 30.)
except:
    sys.exit('No events available')

# describe frequency band --> command line argument
[fmin, fmax] = [0.01, 0.1]
#winlen = 1./fmin * 2
winlen = 50

# describe phase --> command line argument
phase = 'R'

# create file to store data if it does not exist
if not os.path.isfile('raad_earthquakes.csv'):
    with open('raad_earthquakes.csv', mode='w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(['Time', 'Longitude', 'Latitude', 'Depth', 'Magnitude', 'Station', 'Distance', 'Frequency', 'Phase', 'St-St Deg', 'Xcorr'])

# for each event run the analysis
times, amps, corrs = [], [], []
for i,eve in enumerate(cat):
    # test station distance for window length of phase R
    (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], eve.origins[0].latitude, eve.origins[0].longitude)
    # get data for each station for a given event (eve)
    st = get_data(nets, eve, phase, winlen, dis, True)
    if debug:
        print(st)
    # Now we have the data so lets cross-correlate and stack
    #st.plot()
    st.detrend('constant')
    st.taper(0.05)
    st.remove_response()
    st.filter('bandpass',freqmin=fmin, freqmax=fmax)
    # get trace for reference station
    trRef = st.select(station=staGOOD)
    print(trRef)
    try:
        trRef = trRef[0]
    except:
        continue
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
    # write results to csv for analysis
    with open('raad_earthquakes.csv', mode='a') as file:
        write = csv.writer(file, delimiter=',')
        write.writerow([eve['origins'][0]['time'],eve['origins'][0]['longitude'],
                eve['origins'][0]['latitude'], float(eve['origins'][0]['depth'])/1000,
                eve['magnitudes'][0]['mag'], '{}_{}'.format(net,staGOOD), round(site_dis/1000,2),
                '({},{})'.format(fmin,fmax), phase, maxr, round(val,5)])
    amp = np.sqrt(np.sum(trRef.data**2)/np.sum(stack**2))
    amps.append(amp)
    times.append(float(idx)/100.)
    corrs.append(val)
    plt.plot(t, stack, color = 'C1', linewidth=3, label ='Stack')
    plt.xlim((min(t), max(t)))
    plt.legend(loc='upper left')
    plt.savefig('RAAD_{}_{}_{}_{}_{}.png'.format(staGOOD, fmin, fmax, eve['origins'][0]['time'], eve['magnitudes'][0]['mag']), format='png')
    plt.close()
    if plot:
        plt.show()
    del stack
print(amps)
print(times)
print(corrs)
