#!/usr/bin/env python
import numpy as np
from obspy.core import Stream
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
from itertools import combinations

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)


def mean_coors(inv):
    lons, lats = [], []
    for net in inv:
        for sta in net:
            lons.append(sta.longitude)
            lats.append(sta.latitude)
    return np.mean(lons), np.mean(lats)
    
def get_parameters(phase):
    paramdic = {}
    if phase == 'P':
        
        paramdic['station_radius'] = 2.0
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 1./20.
        paramdic['fmax'] = 1./5.
        paramdic['phase'] = 'P'
        paramdic['winlength'] = 50.
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

def get_data(inv, eve, paramdic, model, client, debug=False):
    st = Stream()
    bad_stas = []
    for net in inv:
        for sta in net:
            for chan in sta:
                sncl = net.code + '.' + sta.code + '.' + chan.location_code + '.' + chan.code
                coors = inv.get_coordinates(sncl)
                (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], 
                                                   eve.origins[0].latitude, eve.origins[0].longitude)
                disdeg = kilometer2degrees(dis/1000.)
                arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000.,
                                                  distance_in_degree=disdeg, phase_list=paramdic['phase'])
                if len(arrivals) == 0:
                    break
                pstime = (eve.origins[0].time + arrivals[0].time)-10.
                petime = pstime + paramdic['winlength']
                try:
                    st += client.get_waveforms(net.code, sta.code, chan.location_code, 
                                                chan.code, pstime, petime, attach_response = False)
                    if debug:
                        print('Got data for: ' + sncl)
                except:
                    print('No data for: ' + sncl)
                    bad_stas.append(sncl)
    return st, bad_stas

def proc_data(st, inv, paramdic, eve):
    for tr in st:
        if tr.stats.sampling_rate <= 1.:
            st.remove(tr)
    st.detrend('linear')
    st.detrend('constant')
    st.remove_response(inventory=inv, output='VEL')
    st.filter('bandpass', freqmin=paramdic['fmin'], freqmax=paramdic['fmax'])
    for tr in st:
        coors = inv.get_coordinates(tr.id)
        (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], 
                                            eve.origins[0].latitude, eve.origins[0].longitude)
        tr.stats.back_azimuth = bazi
    st.rotate('->ZNE', inventory=inv)
    st.rotate('NE->RT', inventory=inv)
    return st



def comp_stack(st, comp):
    st2 = st.select(component=comp)
    comb = combinations(range(len(st2)),2)
    used, not_used = [], []
    results = {}
    for ele in list(comb):
        tr1 = st2[ele[0]].copy()
        tr2 = st2[ele[1]].copy()
        cc = correlate(tr1.data, tr2.data, 20)
        shift, value = xcorr_max(cc)
        if value <= 0.8:
            continue
        tr2.stats.starttime -= float(shift)/tr2.stats.sampling_rate
        if 'stack' not in vars():
            stack = tr2.data + tr1.data
            used.append(ele[1])
            used.append(ele[0])
        elif ele[1] not in used:
            stack += tr2.data
            used.append(ele[1])
    try:
        stack /= float(len(used))
    except:
        stack = []
    for idx in range(len(st)):
        if idx not in used:
            not_used.append(idx)
    return stack, not_used
        
    
def pretty_plot(st, stack, eve, not_used, comp, inv, paramdic):
    st2 = st.select(component=comp)
    diss = []
    # compute distances
    for tr in st2:
        coors = inv.get_coordinates(tr.id[:-1] + 'Z')
        (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], 
                                            eve.origins[0].latitude, eve.origins[0].longitude)
        disdeg = kilometer2degrees(dis/1000.)
        diss.append(disdeg)
    print(diss)
    mdiss = min(diss)
    Mdiss = max(diss)
    ptp = np.ptp(stack)
    ran = 0.3*(Mdiss-mdiss)*ptp
    fig = plt.figure(1,figsize=(12,12))
    tithand = st[0].stats.network + ' ' + paramdic['phase'] + '-Wave ' 
    if comp == 'R':
        tithand += ' Radial '
    elif comp == 'Z':
        tithand += ' Vertical '
    elif comp == 'T':
        tithand += ' Transverse '
    tithand += str(eve['origins'][0]['time'].year) + ' '
    tithand += str(eve['origins'][0]['time'].julday) + ' '
    tithand += str(eve['origins'][0]['time'].hour).zfill(2) + ':' + str(eve['origins'][0]['time'].minute).zfill(2)
    mag = eve.magnitudes[0].mag
    magstr = eve.magnitudes[0].magnitude_type
    if 'Lg' in magstr:
        magstr = 'mb_{Lg}'
    gmax, gmin = -100., 500.
    tithand += ' $' + magstr + '$=' + str(mag)
    plt.title(tithand)
    for pair in zip(diss, st2):
        t = pair[1].times()
        if max(pair[1].data/ran + pair[0]) > gmax:
            gmax = max(pair[1].data/ran + pair[0])
        if min(pair[1].data/ran + pair[0]) <  gmin:
            gmin = min(pair[1].data/ran + pair[0])
        p = plt.plot(t, pair[1].data/ran + pair[0])
        plt.text(min(t)+1., pair[0]-+.2, (pair[1].id)[:-4].replace('.',' '), color=p[0].get_color())
        plt.plot(t, stack/ran + pair[0], color='k', alpha=0.5, linewidth=3)
    plt.plot([10., 10.], [0., 2*Mdiss+ ran], color='k', linewidth=3)
    plt.ylim((gmin - 0.02*gmin, gmax + 0.02*gmax))
    plt.xlim((min(t),max(t)))
    plt.xlabel('Time (s)')
    plt.ylabel('Distance (deg)')
    plt.savefig(st[0].stats.network + '_' + comp + '_' + str(eve['origins'][0]['time'].year) +
                str(eve['origins'][0]['time'].julday) + '_' +  str(eve['origins'][0]['time'].hour).zfill(2) +
                str(eve['origins'][0]['time'].minute).zfill(2) + '.png', format='PNG', dpi=400)

    plt.clf()
    plt.close()
    return
    
    
def write_event_results(st, stack, eve, not_used, comp, inv, paramdic):
    st2 = st.select(component=comp)

    # we will make a csv file with the infor for each channel for the event
    filehand = 'Results_' + st2[0].stats.network + '_' + comp + '_' + paramdic['phase'] + \
                '_' + str(eve['origins'][0]['time'].year) + \
                str(eve['origins'][0]['time'].julday) + '_' + \
                str(eve['origins'][0]['time'].hour).zfill(2) + \
                str(eve['origins'][0]['time'].minute).zfill(2) + '.csv'
    
    f = open(filehand,'w')
    f.write('ID, dis, depth, mag, amp, shift, corr, used, ptp, snr \n')
    
    for idx, tr in enumerate(st2):
        f.write(tr.id + ', ')
        
        coors = inv.get_coordinates(tr.id[:-1] + 'Z')

        (dis,azi, bazi) = gps2dist_azimuth(coors['latitude'], coors['longitude'], 
                                            eve.origins[0].latitude, eve.origins[0].longitude)
        disdeg = kilometer2degrees(dis/1000.)
        f.write(str(disdeg) + ', ')
        f.write(str(float(eve['origins'][0]['depth'])/1000) + ', ')
        f.write(str(eve.magnitudes[0].mag) + ', ')
        amp = np.sqrt(np.sum(tr.data**2)/np.sum(stack**2))
        f.write(str(amp) + ', ')
        cc = correlate(tr.data, stack, 20)
        shift, value = xcorr_max(cc)
        f.write(str(shift/float(tr.stats.sampling_rate)) + ', ')
        f.write(str(round(value,5)) + ', ')
        if idx in not_used:
            f.write('Bad, ')
        else:
            f.write('Good, ' )
        tr2 = tr.copy()
        tr2.trim(tr2.stats.starttime, tr2.stats.starttime+5.)
        f.write(str(np.ptp(tr.data)))
        f.write(str(np.ptp(tr2.data)/np.ptp(tr.data)) + '\n')
    f.close()

    return
    
    
def proc_net(net):  
    model = TauPyModel(model="iasp91")
    client = Client()  
    stime = UTCDateTime('2017-001T00:00:00')
    etime = UTCDateTime('2019-150T00:00:00')
    #etime = stime + 5.*24*60.*60
    inv = client.get_stations(starttime=stime, endtime=etime, 
                          channel="*H*", network=net, level="response")
    
    mlon, mlat = mean_coors(inv)
    # We want to eventually scrub the inventory

    paramdic = get_parameters('P')

    cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], 
                        latitude=mlat, longitude=mlon, maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])

    for idx, eve in enumerate(cat):
        print('On event: ' + str(idx+1) + ' of  ' + str(len(cat)))
        # Make a flat file with the important info
        
        # Get the data 
        st, bad_stas = get_data(inv, eve, paramdic, model, client)
        st = proc_data(st, inv, paramdic, eve)
        for comp in ['Z', 'R']:
            stack, not_used = comp_stack(st, comp)
            if len(stack) == 0:
                print('Bad event')
                continue
                
            # We have a pretty plot showing all of the components
            pretty_plot(st, stack, eve, not_used, comp, inv, paramdic)
            
            write_event_results(st, stack, eve, not_used, comp, inv, paramdic)
    
    return
        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

