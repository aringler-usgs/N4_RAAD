#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel


import utils





debug = True
net ="IW"

utils.proc_net(net)
