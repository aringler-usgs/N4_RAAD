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


client = Client()


mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

#added site_dis as an argument since for R phase, winlens were varied, made dependent on only the reference site.
# make winlen and site_dis not required 
def get_data(nets, eve, phase, winlen, dis, debug = True):
    st = Stream()
    for net in nets:
        for sta in net:
            for chan in sta:
                # Here we need to deal with different channel codes e.g. BH and HH should probably get compared
                #if chan.code == 'HHZ':
                if True:
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
                        disdeg = kilometer2degrees(dis/1000.)
                        if debug:
                            print('Here is the distance')
                        arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=disdeg, phase_list=[phase])
                        if debug:
                            print(arrivals)
                        for arrival in arrivals:
                            # do we need to use P-S time to decide window length
                            # winlen = 50.
                            pstime = (eve.origins[0].time + arrival.time)-10.
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


def common_decimation(st):
    # Need to do the prime factors, but it was having issues with the set differences
    for tr in st:
        if tr.stats.sampling_rate == 100:
            tr.decimate(5)
        elif tr.stats.sampling_rate == 40:
            tr.decimate(2)
        elif tr.stats.sampling_rate == 200:
            tr.decimate(2)
            tr.decimate(5)
            
    return st


def remove_channels(st, debug = False):
    # This removed redundant channels for a given location code
    stGood = Stream()
    stas = list(set([(tr.id)[:-4] for tr in st]))
    if debug:
        print(stas)
    for sta in stas:
        # We now have the unique pieces
        net, sta, loc = sta.split('.')
        st2 = st.select(station=sta, location=loc, network=net)
        for tr in st2:
            if 'trGood' not in vars():
                trGood = tr
            else:
                if tr.stats.sampling_rate > trGood.stats.sampling_rate:
                    trGood = tr
        stGood += trGood
        del trGood
    return stGood

def choptocommon(st):
    """ A function to chop the data to a common time window. """
    stime = max([tr.stats.starttime for tr in st])
    etime = min([tr.stats.endtime for tr in st])
    st.trim(starttime=stime, endtime=etime)
    if debug:
        print 'starttime: '+str(stime)+' endtime: '+str(etime)
    return st


def get_parameters(phase):
    paramdic = {}
    if phase == 'P':
        
        paramdic['station_radius'] = 1.5
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 1./20.
        paramdic['fmax'] =1./5.
    elif phase == 'Rayleigh':
        paramdic['station_radius'] = 1.5
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 0.1
        paramdic['fmax'] =1.
    elif  phase == 'S':
        paramdic['station_radius'] = 1.5
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 0.1
        paramdic['fmax'] =1.
    else:
        pass
    
    
    return paramdic


def proc_sta(net, staGOOD, phase):

    paramdic = get_parameters(phase)

    stime = UTCDateTime('2017-001T00:00:00')
    #stime = UTCDateTime('2019-130T00:00:00')
    etime = UTCDateTime('2019-150T00:00:00')
    

    inv = client.get_stations(network=net, station=staGOOD,
                            channel = '*HZ', level="response",
                            starttime=stime, endtime = etime)

    
    if debug:
        print(inv[0][0][0].code)
    
    # grab test station coordinates
    coors = inv.get_coordinates(inv[0].code + '.' + inv[0][0].code + '.' + inv[0][0][0].location_code + '.' + inv[0][0][0].code)
    
    
    # get all stations that fall within paramdic['station_radius'] of staGOOD
    try:
        nets = client.get_stations(starttime=stime, endtime=etime, channel="*HZ",
                    longitude=coors['longitude'], latitude=coors['latitude'],
                    maxradius = paramdic['station_radius'], level="channel", network="IU,II,NE,IW,N4,GS,CU")
    except:
        sys.exit('No stations nearby {}_{}'.format(net, staGOOD))
    
    try:
        cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], latitude=coors['latitude'], 
                                longitude=coors['longitude'], maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])
    except:
        sys.exit('No events available')
    
    if debug:
        print(cat)
    
    
    # create file to store data if it does not exist
    if not os.path.isfile('raad_earthquakes.csv'):
        with open(staGOOD + '_raad_earthquakes.csv', mode='w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(['Time', 'Longitude', 'Latitude', 'Depth', 'Magnitude', 'Station', 'Distance', 'Frequency', 'Phase', 'St-St Deg', 'Xcorr', 'Amplitude', 'Lag', 'NumStas'])
    
    # for each event run the analysis
    times, amps, corrs = [], [], []
    for i,eve in enumerate(cat):
        if debug:
            print('In the main loop')
        
        # test station distance for window length of phase R
        (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], eve.origins[0].latitude, eve.origins[0].longitude)
        # get data for each station for a given event (eve)
        st = get_data(nets, eve, phase, paramdic['length'], dis, False)
        print('Here we are')
        if debug:
            print(st)
        st = remove_channels(st)
        try:
            st.remove_response()
        except:
            print('Bad metadata for:')
            print(st)
            
        st = common_decimation(st)
        print(st)
        try:
            st = choptocommon(st)
        except:
            continue
        
        # Now we have the data so lets cross-correlate and stack
        #st.plot()
        st.detrend('constant')
        st.taper(0.05)
        st.filter('bandpass',freqmin=paramdic['fmin'], freqmax=paramdic['fmax'])
        
        
        # get trace for reference station
        trRef = st.select(station=staGOOD)
        print(trRef)
        try:
            trRef = trRef[0]
        except:
            continue
        fig = plt.figure(1, figsize=(12,12))
        for tr in st:
            idx, val = xcorr(tr, trRef, 20)
            if debug:
                print(idx)
                print(val)
            if val < 0.8:
                st.remove(tr)
                continue
            
            t = np.arange(0, tr.stats.npts)/tr.stats.sampling_rate
            plt.plot(t - float(idx)/float(tr.stats.sampling_rate), tr.data*10**6, label= tr.id + ' ' + str(round(val,5)))
            plt.xlabel('Time (s)')
            plt.ylabel('Velocity ($\mu$m/s)')
            tr.stats.starttime -= float(idx)/tr.stats.sampling_rate
            if 'stack' not in vars():
                stack = tr.data
            else:
                stack += tr.data
        stack /= float(len(st))
        idx, val = xcorr(trRef, stack, 20)
        # write results to csv for analysis
        amp = np.sqrt(np.sum(trRef.data**2)/np.sum(stack**2))
        if len(st) <= 1:
            # Not enough data
            continue
            
        with open(staGOOD + '_raad_earthquakes.csv', mode='a') as file:
            write = csv.writer(file, delimiter=',')
            write.writerow([eve['origins'][0]['time'],eve['origins'][0]['longitude'],
                    eve['origins'][0]['latitude'], float(eve['origins'][0]['depth'])/1000,
                    eve['magnitudes'][0]['mag'], '{}_{}'.format(net,staGOOD), round(dis/1000,2),
                    '({},{})'.format(paramdic['fmin'], paramdic['fmax']), phase, paramdic['station_radius'], round(val,5), amp, float(idx)/float(tr.stats.sampling_rate), len(st)])
        
    
        plt.plot(t, stack*10**6, color = 'C1', linewidth=3, label ='Stack')
        plt.xlim((min(t), max(t)))
        plt.legend(loc=2)
        strtime = str(eve['origins'][0]['time'].year) + ' '
        strtime += str(eve['origins'][0]['time'].julday) + ' '
        strtime += str(eve['origins'][0]['time'].hour).zfill(2) + ':' + str(eve['origins'][0]['time'].minute).zfill(2)
        thand = strtime.replace(' ','_')
        thand = thand.replace(':','_')
        plt.title(phase + ' Comparison M' + str(eve['magnitudes'][0]['mag']) + ' ' + strtime)
        plt.savefig('RAAD_{}_{}_{}_{}_{}.png'.format(staGOOD, paramdic['fmin'], paramdic['fmax'], thand, eve['magnitudes'][0]['mag']), format='png')
        plt.clf()
        plt.close()
        del stack
    return





####################
# 1. What are the best frequencies bands for P, S, and Raleigh Waves?
# 2. What are the best parameters for earthquakes for P, S, and Raleigh (e.g depth restrictions, st-eq distance)
# 3. Test maxr (st-st radius) to see where program breaks down for P, S, Raleigh?
###################

# describe phase --> command line argument
phase = 'P'



# extra info --> command line argument
debug = True

# show plot
plot = True

# station --> command line argument


net, staGOOD = 'IW', 'FXWY'

stas = ['DLMT', 'FLWY', 'FXWY', 'IMW', 'LOHW', 'MFID', 'MOOW', 'PHWY', 'PLID', 'REDW', 'RRI2', 'RWWY', 'SMCO', 'SNOW', 'TPAW']

def proc_net(sta):
    proc_sta(net, sta, phase) 

net = sys.argv[1]
sta = sys.argv[2]

proc_net(sta)


#from multiprocessing import Pool
#for sta in stas:
#    proc_net(sta)
#pool = Pool(14)
#pool.map(proc_net, stas)
