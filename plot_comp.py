#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import plot_utils
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)




def cull_for_plotting_2(dic_all, vari1, vari2, vari1mM, vari2mM, good_sta):
    good1, good2 = [], []
    for idx, sta in enumerate(dic_all['ID']):
        if good_sta in sta:
            if (dic_all[vari1][idx] > vari1mM[0]) and (dic_all[vari1][idx] < vari1mM[1]):
                if (dic_all[vari2][idx] > vari2mM[0]) and (dic_all[vari2][idx] < vari2mM[1]):
                    good1.append(dic_all[vari1][idx])
                    good2.append(dic_all[vari2][idx])
    good1 = np.asarray(good1)
    good2 = np.asarray(good2)
    return good1, good2


def cull_for_plotting_3(dic_all, vari1, vari2, vari3, vari1mM, vari2mM, vari3mM, good_sta):
    good1, good2, good3 = [], [], []
    for idx, sta in enumerate(dic_all['ID']):
        if good_sta in sta:
            if dic_all['used'][idx] == 'Good':
                if (dic_all[vari1][idx] > vari1mM[0]) and (dic_all[vari1][idx] < vari1mM[1]):
                    if (dic_all[vari2][idx] > vari2mM[0]) and (dic_all[vari2][idx] < vari2mM[1]):
                        if (dic_all[vari3][idx] > vari3mM[0]) and (dic_all[vari3][idx] < vari3mM[1]):
                            good1.append(dic_all[vari1][idx])
                            good2.append(dic_all[vari2][idx])
                            good3.append(dic_all[vari3][idx])
    good1 = np.asarray(good1)
    good2 = np.asarray(good2)
    good3 = np.asarray(good3)
    return good1, good2, good3


    

net ='N4'
comp = 'Z'
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


goodsta = ''
vari1 = 'mag'
vari1mM = [3., 9.]
vari2 = 'amp'
vari2mM = [0., 3.]
vari3 = 'corr'
vari3mM = [0.75, 1.]

good1, good2, good3 = cull_for_plotting_3(dic_all, vari1, vari2, vari3, vari1mM, vari2mM, vari3mM, goodsta)


distinct = list(set(good1))

mean =[]
std = []
for ele in distinct:
    print(ele)
    newamps = good2[good1 == ele]
    mean.append(np.mean(newamps))
    std.append(np.std(newamps))

fig = plt.figure(1)

plt.scatter(good1, good2, c=good3) 
plt.errorbar(distinct, mean,yerr=std, fmt='o')
#plt.legend(net + ' ' + comp + '-Component')
plt.xlabel(vari1)
plt.ylabel(vari2)
plt.xlim((min(good1)-0.1, max(good1)+0.1))
plt.ylim((min(good2)-0.1, max(good2)+0.1))
plt.colorbar()    
#plt.savefig(net + '_' + comp + '_summary.jpg', format='JPG',dpi=400)
plt.show()
    

