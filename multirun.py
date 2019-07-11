#!/usr/bin/env python
import os
import sys
from multiprocessing import Pool

net = 'IW'


def runcheck(sta):
    os.system('python checksta.py ' + net + ' ' + sta)
    

stas = ['DLMT', 'FLWY', 'FXWY', 'IMW', 'LOHW', 'MFID', 'MOOW', 'PHWY', 'PLID', 'REDW', 'RRI2', 'RWWY', 'SMCO', 'SNOW', 'TPAW']

pool = Pool(14)
pool.map(runcheck, stas)
