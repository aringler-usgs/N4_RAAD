#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel

import utils


    
debug = True
net ="GS"
model = TauPyModel(model="iasp91")
client = Client()  
stime = UTCDateTime('2019-188T00:00:00')
etime = UTCDateTime('2019-204T00:00:00')
#etime = stime + 5.*24*60.*60
inv = client.get_stations(starttime=stime, endtime=etime, station="CA*",
                          channel="HH*", network=net, level="response")
    
mlon, mlat = utils.mean_coors(inv)
# We want to eventually scrub the inventory

paramdic = utils.get_parameters('P')

cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], 
                        latitude=mlat, longitude=mlon, maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])


def proc_event(eve):
    # Get the data 
    try:
        st, bad_stas = utils.get_data(inv, eve, paramdic, model, client)
        print(st)
        st = utils.proc_data(st, inv, paramdic, eve)
    except:
        return
    for comp in ['Z', 'R']:
        stack, not_used = utils.comp_stack(st, comp)
        if len(stack) == 0:
            print('Bad event')
            continue
                
        # We have a pretty plot showing all of the components
        try:
        #if True:
            utils.pretty_plot(st, stack, eve, not_used, comp, inv, paramdic)
            utils.write_event_results(st, stack, eve, not_used, comp, inv, paramdic)
        except:
            print('Problem')
            print(eve)
            continue
    return



from multiprocessing import Pool
for eve in cat:
    proc_event(eve)
#pool = Pool(40)
#pool.map(proc_event, cat)
