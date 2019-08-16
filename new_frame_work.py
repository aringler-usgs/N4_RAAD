#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import numpy as np
import utils


    
debug = True
net ="N4"
model = TauPyModel(model="iasp91")
client = Client()  
stime = UTCDateTime('2019-154T00:00:00')
etime = UTCDateTime('2019-204T00:00:00')
small = False
inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                          channel="*H*", network=net, level="response")

inv = utils.scrub_inventory(inv)

if small:
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
                
                utils.pretty_plot(st, stack, eve, not_used, comp, inv, paramdic)
                utils.write_event_results(st, stack, eve, not_used, comp, inv, paramdic)
            except:
                print('Problem')
                print(eve)
                continue
        return



    from multiprocessing import Pool
    #for eve in cat:
        #proc_event(eve)
        #import sys
        #sys.exit()
    if 'pool' not in vars():
        pool = Pool(60)
    pool.map(proc_event, cat)
    
else:
    minlon, Mlon, minlat, Mlat = utils.min_max_coors(inv)
    # We want to eventually scrub the inventory
    
    paramdic = utils.get_parameters('P')
    mlon, mlat = utils.mean_coors(inv)
    cat = client.get_events(starttime=stime, endtime=etime, minmagnitude=paramdic['min_mag'], maxmagnitude=paramdic['max_mag'], 
                        latitude=mlat, longitude=mlon, maxradius=paramdic['max_radius'], minradius = paramdic['min_radius'])

    for lat in np.arange(minlat, Mlat, 1.):
        for lon in np.arange(minlon, Mlon, 1.):

            try:
                inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                          channel="*H*", network="IU,N4,US", level="response", latitude=lat,
                          longitude=lon, maxradius=paramdic['station_radius'])
            except:
                print('No stations for: ' + str(lat) + ' '  + str(lon))
                continue
            
            inv = utils.scrub_inventory(inv)

            if len(inv[0]) < 3:
                print('too few stations')
                continue

            def proc_event(eve):
                ## Get the data 
                try:
                    st, bad_stas = utils.get_data(inv, eve, paramdic, model, client)
                    st = utils.proc_data(st, inv, paramdic, eve)
                except:
                    return
                
                for comp in ['Z', 'R']:
                    if len(st.select(component=comp)) < 3:
                        continue
                    stack, not_used = utils.comp_stack(st, comp)
                    if len(stack) == 0:
                        print('Bad event')
                        continue
                    ## We have a pretty plot showing all of the components
                    try:
                    ##if True:
                        utils.pretty_plot(st, stack, eve, not_used, comp, inv, paramdic)
                        ## Add the lats and lons for book keeping.
                        utils.write_event_results(st, stack, eve, not_used, comp, inv, paramdic, lat, lon)
                    except:
                        print('Problem')
                        print(eve)
                        continue
                return

            from multiprocessing import Pool
            #for eve in cat:
                #proc_event(eve)
                #import sys
                #sys.exit()
            if 'pool' not in vars():
                pool = Pool(120)
            pool.map(proc_event, cat)

