#!/usr/bin/env python
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from obspy.clients.fdsn import Client
import plot_utils
import glob
import sys


    
def plot_figure(net, sta, dic_list, sta_inv):
    comps = ['Z', 'R']
    fig = plt.figure(1, figsize=(12,12))
    for idxcomp, dic_all in enumerate(dic_list):
        amps, lags, corrs, times = [], [], [], []
        for idx, tsta in enumerate(dic_all['ID']):
            if sta in tsta:
                tamp = dic_all['amp'][idx]
                if tamp > 2.:
                    tamp = 2.
                amps.append(tamp)
                lags.append(dic_all['shift'][idx])
                corrs.append(dic_all['corr'][idx])
                times.append(float(dic_all['time'][idx].year) + float(dic_all['time'][idx].julday)/365.25)
                
        
        times = np.array(times)
        idx = np.argsort(times)
        lags = np.array(lags)[idx]
        corrs = np.array(corrs)[idx]
        amps = np.array(amps)[idx]
        times = times[idx]
        epoch_times =[]
        for chan in sta_inv:
            epoch_times.append(float(chan.start_date.year) + float(chan.start_date.julday)/365.25)
            
        
        plt.subplot(3,1,1)
        plt.title(net + ' ' + sta)
        plt.plot(times, amps, '.', alpha=.5, label=comps[idxcomp] + ' Comp.')
        plt.plot(times, np.convolve(amps, np.ones((5,))/5, mode='same'), label='5 pt ' + comps[idxcomp] + ' Comp.')
        for epoch in epoch_times:
            plt.plot([epoch, epoch], [-10, 10], label='Epoch Boundary')
        plt.xlim((min(times), max(times)))
        plt.ylabel('Rel. Amplitude')
        plt.ylim((0., 3.))
        plt.subplot(3,1,2)
        plt.plot(times, lags, '.', alpha=.5, label=comps[idxcomp])
        plt.plot(times, np.convolve(lags, np.ones((5,))/5, mode='same'), label='5 pt ' + comps[idxcomp] + ' Comp.')
        for epoch in epoch_times:
            plt.plot([epoch, epoch], [-10, 10])
        plt.xlim((min(times), max(times)))
        plt.ylim((-2,2))
        plt.ylabel('Lag (s)')
        plt.subplot(3,1,3)
        plt.plot(times, corrs,'.', alpha=.5, label=comps[idxcomp])
        plt.plot(times, np.convolve(corrs, np.ones((5,))/5, mode='same'), label='5 pt ' + comps[idxcomp])
        for epoch in epoch_times:
            plt.plot([epoch, epoch], [-10, 10])
        plt.xlim((min(times), max(times)))
        plt.ylim((-1,1))
        plt.ylabel('Correlation')
        plt.xlabel('Time (year)')
    
    ax = plt.subplot(3,1,3)
    handles, labels = ax.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc = 'lower center', ncol = 5, fontsize = 15)
    plt.savefig(net + '_' + sta + '_summary.png', format='PNG')
    #plt.show()
    plt.clf()
    plt.close()
    
    return
    
client = Client()
stime = UTCDateTime('2017-001T00:00:00')
etime = UTCDateTime('2019-196T00:01:00')

net = 'N4'
inv = plot_utils.get_dataless(net)
#comp ='Z'


dicZ = plot_utils.get_dic(net, 'Z')
dicR = plot_utils.get_dic(net, 'R')


stas = []
for nets in inv:
    if nets.code == net:
        for sta_inv in nets:
            try:
                
                plot_figure(net, sta_inv.code, [dicZ, dicR], sta_inv)
            
            except:
                print('Bad station: ' + sta_inv.code)
                plt.close()
                plt.clf()
