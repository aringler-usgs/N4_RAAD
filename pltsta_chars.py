#!/usr/bin/env python
import glob
import sys
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt

import matplotlib as mpl
#Set font parameters using matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


debug = True

def parse(curfile):
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


net =''

dic_all ={}
files = glob.glob(net + '*.csv')
for curfile in files:
    stuff = parse(curfile)
    for item in stuff:
        try:
            dic_all[item] += stuff[item]
        except:
            dic_all[item] =[]
            dic_all[item] += stuff[item]

for item in dic_all:
    print(item)

item1 = 'Amplitude'
item2 = 'Xcorr'
fig = plt.figure(1,figsize=(9,9))
plt.plot(dic_all[item1], dic_all[item2], '.')
plt.xlabel(item1)
plt.ylabel(item2)
plt.savefig(net + '_' + item1 + '_' + item2 + '.png', format='PNG')
plt.show()
