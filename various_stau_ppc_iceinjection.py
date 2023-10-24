#!/usr/bin/env python

####
# A script for the stau generation, the stau propagation through the ice, 
# and the photon propagation using ppc. 
####

import argparse
from os.path import expandvars
import os

usage = "usage: %prog [options] inputfile"
parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-o', '--out',default=".",
                  help='Write output to OUTFILE (.i3{.gz} format)')
parser.add_argument('-s', '--seed', default=-1, type=int,
                  help="Initial seed for the random number generator")
parser.add_argument('-r', '--run', default=0, type=int,
                  help="The run number for this simulation")
parser.add_argument('-d', '--debug', action='store_true', default=False, help="Debug mode")

### options for the ppc 
parser.add_argument('--domos', default=5, type=float)
parser.add_argument('--domeff', default=1., type=float)
parser.add_argument('--gcd', default=os.path.os.path.expandvars('$I3_DATA/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz'),
                    type=str,
                    help='gcd file')
parser.add_argument('--icemodel', choices=('SpiceMie', 'SpiceLea', 'Spice3.2', 'Spice3.2.1',
                                           'Spicebfr-v1', 'Spicebfr-v2', 'Spicebfr-v2_flat',
                                           'Spicebfr-v2_tilt'),
                    default='Spicebfr-v2', help='The ice model to use for photon propagation.')
parser.add_argument('--holeice',default='as.h2-50cm', choices=('as.flasher_p1_0.30_p2_0',
                                                               'as.flasher_p1_0.30_p2_-1',
                                                               'as.flasher_p1_0.30_p2_+1',
                                                               'as.flasher_p1_0.30_p2_-2',
                                                               'as.flasher_p1_0.30_p2_-3',
                                                               'as.h2-50cm',
                                                               'as.nominal',
                                                               'as.set0_p0=0.0_p1=0.0'), help="Holeice file")

parser.add_argument('--usecpu',action='store_true',default=False)
parser.add_argument('--infile',type=str,default="")

args = parser.parse_args()

from I3Tray import *
import sys

from icecube import icetray, dataclasses, dataio, phys_services, sim_services

from PropagateStaus import PropagateStaus

#
import math
import numpy
import random

from validation import *

icetray.set_log_level(icetray.I3LogLevel.LOG_TRACE)

import time

def stautrack_ppc(args, infilename):
    print("RNG being initiated...")
    # set up a random number generator
    if args.seed == -1:
        seed = int(random.uniform(0,10000))
    else:
        seed = 8888

    randomService = phys_services.I3SPRNGRandomService(
        seed = seed,
        nstreams = 100000000,
        streamnum = args.run)

    load('ppc')
    tray = I3Tray()

    tray.context['I3RandomService'] = randomService
    tray.Add('I3Reader', Filenamelist=[args.gcd]+[infilename])

    def frame_time_0(frame):
        frame['Time_0'] = dataclasses.I3String(f'{time.time()}')
    tray.AddModule(frame_time_0, 'frame_time_0',streams=[icetray.I3Frame.DAQ])

    if args.usecpu:
        import os
        os.putenv('OCPU','1')
    tray.Add('i3ppc','ppc',
            MCTree = "I3MCTree",
            infoName = "PPCInfoDict", # create a dict
            verbose=args.debug,
            photons=True,
            pseries=True)

    def frame_time_1(frame):
        frame['Time_1'] = dataclasses.I3String(f'{time.time()}')
    tray.AddModule(frame_time_1, 'frame_time_1',streams=[icetray.I3Frame.DAQ])

    tray.Add("Rename", "PPCRename", Keys=["MCPESeriesMap", "I3MCPESeriesMap"])

    outfilename = args.out + infilename.split('.i')[0].split('/')[-1] + '_ppc.i3'

    tray.Add('I3Writer', filename=outfilename, streams=[icetray.I3Frame.DAQ])
    tray.Execute()

if __name__=='__main__':
    os.makedirs(args.out,exist_ok=True)
    stautrack_ppc(args, args.infile)
