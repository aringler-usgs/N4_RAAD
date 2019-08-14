#!/usr/bin/env python
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from obspy.clients.fdsn import Client




def parse(curfile, debug = True):
    f = open(curfile, 'r')
    for line in f:
        line = line.strip()
        if debug:
            print(line)
        if 'stuff' not in vars():
            stuff = {}
            params = line.split(',')
            for ele in params:
                stuff[ele] = []
        else:
            line = line.split(',')
            for ele in zip(params,line):
                try:
                    stuff[ele[0]].append(float(ele[1]))
                except:
                    stuff[ele[0]].append(ele[1])
    f.close()
    
    goodstuff = stuff
    del stuff
    
    return goodstuff
    
def plot_figure(net, sta):
    
    results = parse(net + '_' + sta + '_raad_earthquakes.csv')
    results['Time'] = [UTCDateTime(ctime) for ctime in results['Time']]
    times = [float(ctime.year) + float(ctime.julday)/365.25 for ctime in results['Time']]
    amps = results['Amplitude']
    lags = results['Lag']
    corrs = results['Xcorr']

    fig = plt.figure(1, figsize=(12,12))
    plt.subplot(3,1,1)
    plt.title(net + ' ' + sta)
    plt.plot(times, amps,'.', alpha=.5)
    plt.plot(times, np.convolve(amps, np.ones((5,))/5, mode='same'))
    plt.xlim((min(times), max(times)))
    plt.ylabel('Rel. Amplitude')
    plt.subplot(3,1,2)
    plt.plot(times, lags,'.', alpha=.5)
    plt.plot(times, np.convolve(lags, np.ones((5,))/5, mode='same'))
    plt.xlim((min(times), max(times)))
    plt.ylim((-2,2))
    plt.ylabel('Lag (s)')
    plt.subplot(3,1,3)
    plt.plot(times, corrs,'.', alpha=.5)
    plt.plot(times, np.convolve(corrs, np.ones((5,))/5, mode='same'))
    plt.xlim((min(times), max(times)))
    #plt.ylim((-2,2))
    plt.ylabel('Correlation')
    plt.savefig(net + '_' + sta + '_summary.png', format='PNG')
    plt.clf()
    plt.close()
    #plt.show()
    return
    
client = Client()
stime = UTCDateTime('2017-001T00:00:00')
etime = UTCDateTime('2019-196T00:01:00')

net = 'N4'
inv = client.get_stations(network=net, station="*",
                            channel = '*',location = '*', level="response",
                            starttime=stime, endtime = etime)
stas = []
for nets in inv:
    for sta in nets:
        try:
            plot_figure(net, sta.code)
        except:
            print('Bad station: ' + sta.code)
            plt.close()
            plt.clf()
