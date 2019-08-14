#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

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
                    stuff[ele[0]].append(ele[1])
    f.close()
    
    goodstuff = stuff
    del stuff
    return goodstuff


def cull_for_plotting_2(dic_all, vari1, vari2, vari1mM, vari2mM, good_sta):
    good1, good2 = [], []
    for idx, sta in enumerate(dic_all['ID']):
        if good_sta in sta:
            if (dic_all[vari1][idx] > vari1mM[0]) and (dic_all[vari1][idx] < vari1mM[1]):
                if (dic_all[vari2][idx] > vari2mM[0]) and (dic_all[vari2][idx] < vari2mM[1]):
                    good1.append(dic_all[vari1][idx])
                    good2.append(dic_all[vari2][idx])
    return good1, good2


def cull_for_plotting_3(dic_all, vari1, vari2, vari3, vari1mM, vari2mM, vari3mM, good_sta):
    good1, good2, good3 = [], [], []
    for idx, sta in enumerate(dic_all['ID']):
        if good_sta in sta:
            print(dic_all['used'][idx])
            if dic_all['used'][idx] == ' Good':
                if (dic_all[vari1][idx] > vari1mM[0]) and (dic_all[vari1][idx] < vari1mM[1]):
                    if (dic_all[vari2][idx] > vari2mM[0]) and (dic_all[vari2][idx] < vari2mM[1]):
                        if (dic_all[vari3][idx] > vari3mM[0]) and (dic_all[vari3][idx] < vari3mM[1]):
                            good1.append(dic_all[vari1][idx])
                            good2.append(dic_all[vari2][idx])
                            good3.append(dic_all[vari3][idx])
    return good1, good2, good3


    

net ='*'
comp = '*'
dic_all ={}
files = glob.glob('Results_' + net + '_' + comp + '_*.csv')
for curfile in files:
    stuff = parse(curfile)
    for item in stuff:
        try:
            dic_all[item] += stuff[item]
        except:
            dic_all[item] = []
            dic_all[item] += stuff[item]

cols = list(set(dic_all['ID']))



goodsta = ''
vari1 = 'mag'
vari1mM = [3., 9.]
vari2 = 'amp'
vari2mM = [0., 2.]
vari3 = 'corr'
vari3mM = [-1., 1.]

good1, good2, good3 = cull_for_plotting_3(dic_all, vari1, vari2, vari3, vari1mM, vari2mM, vari3mM, goodsta)


print(good3)

fig = plt.figure(1)
plt.scatter(good1, good2, c=good3) 
#plt.legend(net + ' ' + comp + '-Component')
plt.xlabel(vari1)
plt.ylabel(vari2)
plt.xlim((vari1mM[0], vari1mM[1]))
plt.ylim((vari2mM[0], vari2mM[1]))
plt.colorbar()    
#plt.savefig(net + '_' + comp + '_summary.jpg', format='JPG',dpi=400)
plt.show()
    

