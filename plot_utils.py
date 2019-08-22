#!/usr/bin/env python
import numpy as np
import os
import pickle
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

def transform_sta(stas):
    goodstas = []
    for sta in stas:
        vals = sta.split('.')
        goodstas.append(vals[0] + '.' + vals[1] + '.00.LHZ')
    return goodstas


def parse(curfile, debug = False):
    f = open(curfile, 'r')
    for line in f:
        line = line.strip()
        if debug:
            print(line)
        if 'stuff' not in vars():
            stuff = {}
            line = line.replace(' ','')
            params = line.split(',')
            for ele in params:
                stuff[ele] = []
        else:
            line = line.split(',')
            for ele in zip(params,line):
                try:
                    stuff[ele[0]].append(float(ele[1]))
                except:
                    lineval = ele[1].replace(' ','')
                    stuff[ele[0]].append(lineval)
    f.close()
    goodstuff = stuff
    del stuff
    return goodstuff

def get_plot_params(dic_all, inv, key, debug = False):
    finalstas =[]
    allvals = np.asarray(dic_all[key])
    allstas = np.asarray(dic_all['ID'])
    cleanallstas = []
    for sta in allstas:
        cleanallstas.append(sta.split('.')[1])
    cleanallstas = np.asarray(cleanallstas)
    if debug:
        print(cleanallstas)
    lats ,lons, staval = [], [], []
    stas = list(set(cleanallstas))
    for sta in stas:
        vals = allvals[cleanallstas == sta]    
        if key == 'used':
            goods = len(vals[vals == 'Good'])
            bads = len(vals[vals == 'Bad'])
            if (goods == 0) and (bads == 0):
                continue
        for chan in allstas:
            chan = chan[:-1] + 'Z'
            if sta in chan:
                try:
                    coors = inv.get_coordinates(chan)
                    lats.append(coors['latitude'])
                    lons.append(coors['longitude'])
                    staval.append(100*goods/(goods+ bads))
                    finalstas.append(sta)
                    parameter = 'Percentage (\%)'
                    break
                except:
                    print('Can not find station: ' + sta)
                    pass
        
    staval = np.asarray(staval)
    lats = np.asarray(lats)
    lons = np.asarray(lons)
    finalstas = np.asarray(finalstas)
    return lats, lons, staval, parameter, finalstas

def get_dataless(net):
    if os.path.exists(net + '_metadata.pickle'):
        with open('N4_metadata.pickle', 'rb') as fhand:
            inv = pickle.load(fhand) 
    else:
        if net == 'N4':
            net2 = "N4,US,IU,II"
        else:
            net2 = net
        client = Client()  
        stime = UTCDateTime('2017-204T00:00:00')
        etime = UTCDateTime('2019-204T00:00:00')
        inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                                  channel="*H*", network=net2, level="response")
        inv = utils.scrub_inventory(inv)
        with open(net + '_metadata.pickle', 'wb') as fhand:
            pickle.dump(inv, fhand)
    return inv
