#!/usr/bin/env python
import os
import sys
from multiprocessing import Pool
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

net = 'N4'


def runcheck(sta):
    os.system('python checksta.py ' + net + ' ' + sta)
    
client = Client()
stime = UTCDateTime('2019-196T00:00:00')
etime = UTCDateTime('2019-196T00:01:00')

inv = client.get_stations(network="N4", station="*",
                            channel = '*',location = '*', level="response",
                            starttime=stime, endtime = etime)
stas = []
for nets in inv:
    for sta in nets:
        stas.append(str(sta.code))

print(stas)
stas = list(set(stas))
#stas = ['DLMT', 'FLWY', 'FXWY', 'IMW', 'LOHW', 'MFID', 'MOOW', 'PHWY', 'PLID', 'REDW', 'RRI2', 'RWWY', 'SMCO', 'SNOW', 'TPAW']
print(stas)
pool = Pool(10)
pool.map(runcheck, stas)
