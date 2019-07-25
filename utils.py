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

UTCDateTime.DEFAULT_PRECISION = 1

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
        paramdic['min_radius'] = 2.
        paramdic['max_radius'] = 60.
        paramdic['min_mag'] = 5.0
        paramdic['max_mag'] = 8.0
        paramdic['length'] = 60
        paramdic['fmin'] = 1./15.
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
                pstime = (eve.origins[0].time + arrivals[0].time)-20.
                petime = pstime + paramdic['winlength']
                try:
                    st += client.get_waveforms(net.code, sta.code, chan.location_code, 
                                                chan.code, pstime, petime, attach_response = False)
                    if debug:
                        print('Got data for: ' + sncl)
                except:
                    print('No data for: ' + sncl)
                    bad_stas.append(sncl)
    #st = choptocommon(st)
    return st, bad_stas

def choptocommon(st):
    """ A function to chop the data to a common time window. """
    st2 = st.copy()
    stime = max([tr.stats.starttime for tr in st])
    etime = min([tr.stats.endtime for tr in st])
    st2.trim(starttime=stime, endtime=etime)
    if debug:
        print 'starttime: '+str(stime)+' endtime: '+str(etime)
    return st2


def proc_data(st, inv, paramdic, eve):
    for tr in st:
        if tr.stats.sampling_rate <= 1.:
            st.remove(tr)
        if tr.stats.sampling_rate == 200:
            tr.decimate(5)
            tr.decimate(2)
            tr.decimate(2)
        if tr.stats.sampling_rate == 40:
            tr.decimate(4)
        if tr.stats.sampling_rate == 100:
            tr.decimate(2)
            tr.decimate(5)    
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
    
    st.trim(st[0].stats.starttime+20., st[0].stats.endtime, pad=True, fill_value=0.)
    st.taper(0.05)
    newlen = min([tr.stats.npts for tr in st])
    for tr in st:
        tr.data[:newlen]

    return st



def comp_stack(st, comp):
    st2 = st.select(component=comp)
    print(st2)
    comb = combinations(range(len(st2)),2)
    used, not_used = [], []
    results = {}
    for ele in list(comb):
        tr1 = st2[ele[0]].copy()
        tr2 = st2[ele[1]].copy()
        cc = correlate(tr1.data, tr2.data, 20)
        shift, value = xcorr_max(cc)
        print(shift)
        print(value)
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
    st2 = st2.copy()
    
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
    
    # super nearby stations
    if Mdiss - mdiss < 2:
        smallplot = True
        diss = np.arange(float(len(st2)))
        for tr in st2:
            tr.data /= np.max(np.abs(stack))
        stack /= np.max(np.abs(stack))
        
    else:
        smallplot = False
                
    
    ptp = np.ptp(stack)
    
    print(diss)
    ran = 0.3*(Mdiss-mdiss)*ptp
    if smallplot:
        ran = 1.
    
    diss /= ran
    fig = plt.figure(1,figsize=(16,12))
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
    labs = []
    for pair in zip(diss, st2):
        labs.append((pair[1].id).replace('.',' '))
        t = pair[1].times()
        if max(pair[1].data/ran + pair[0]) > gmax:
            gmax = max(pair[1].data/ran + pair[0])
        if min(pair[1].data/ran + pair[0]) <  gmin:
            gmin = min(pair[1].data/ran + pair[0])
        if pair[1].max() > np.max(np.abs(stack))*3.:
            p = plt.plot(t, (pair[1].data)/(np.max(np.abs(stack))*3.) + pair[0])
            plt.text(min(t)+1., pair[0]+.2, (pair[1].id)[:-4].replace('.',' ') + ' gain', color=p[0].get_color())
        else:
            p = plt.plot(t, pair[1].data/ran + pair[0])
            plt.text(min(t)+1., pair[0]-+.2, (pair[1].id)[:-4].replace('.',' '), color=p[0].get_color())

        plt.plot(t, stack/ran + pair[0], color='k', alpha=0.5, linewidth=3)
    if smallplot:
        plt.yticks(diss, labs)
    plt.plot([10., 10.], [-1000., 1000.], color='k', linewidth=3)
    plt.ylim((gmin - 0.02*gmin, gmax + 0.02*gmax))
    if smallplot:
        plt.ylim((min(diss)-1,max(diss)+1))
    plt.xlim((min(t),max(t)))
    plt.xlabel('Time (s)')
    if smallplot:
        plt.ylabel('Station index')
    else:
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
        f.write(str(np.ptp(tr.data)) + ', ')
        f.write(str(np.ptp(tr.data)/np.ptp(tr2.data)) + '\n')
    f.close()

    return
    
    

        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

