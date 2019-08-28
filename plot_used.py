#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import plot_utils
#mpl.rc('font',family='serif')
#mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
#mpl.rc('font',size=16)

net ='N4'
comp = 'R'
dic_all ={}
files = glob.glob(net + '_results/Results_' + net + '_' + comp + '_*.csv')

for curfile in files:
    stuff = plot_utils.parse(curfile)
    for item in stuff:
        try:
            dic_all[item] += stuff[item]
        except:
            dic_all[item] = []
            dic_all[item] += stuff[item]


stas = list(set(dic_all['ID']))
allstas = np.asarray(dic_all['ID'])
used = np.asarray(dic_all['used'])
snr = np.asarray(dic_all['snr'])

goodsta = []
staratio = []
meansnr = []
for sta in stas:
    vals = used[allstas == sta]
    snrs = snr[allstas == sta]
    
    goods = len(vals[vals == 'Good'])
    bads = len(vals[vals == 'Bad'])
    goodsta.append(100*goods/(goods+ bads))
    meansnr.append(np.mean(snr))

goodsta = np.asarray(goodsta)
meansnr = np.asarray(meansnr)
stas = np.asarray(stas)
stas = stas[np.argsort(goodsta)]
meansnr = meansnr[np.argsort(goodsta)]
goodsta = goodsta[np.argsort(goodsta)]

fig = plt.figure(1, figsize=(12,12))
badstas = []
for idx, triple in enumerate(zip(stas, goodsta, meansnr)):
    if triple[1] < 50:
        plt.plot(triple[1], idx, 'k.')
        badstas.append(triple[0])
plt.yticks(range(len(badstas)), badstas, fontsize=8)
plt.show()

    
