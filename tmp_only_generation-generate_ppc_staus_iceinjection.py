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

### options for the generation 
parser.add_argument('-o', '--out',default="test.i3",
                  help='Write output to OUTFILE (.i3{.gz} format)')
parser.add_argument('-s', '--seed', default=-1, type=int,
                  help="Initial seed for the random number generator")
parser.add_argument('-r', '--run', default=0, type=int,
                  help="The run number for this simulation")
parser.add_argument('-n', '--nevents', default=1, type=int,
                  help="The number of events per run")
parser.add_argument('-w', '--weightshow', action='store_true', default=False, 
                  help="weight showing to the out filename")
parser.add_argument('--minloge', default=5, type=float, 
                  help="Minimum LogE value")
parser.add_argument('--maxloge', default=12, type=float,
                  help="Maximum LogE value")
parser.add_argument('-d', '--debug', action='store_true', default=False, help="Debug mode")
parser.add_argument('--mass', default=150, type=int,
                  help="stau mass (only int supported now)")
parser.add_argument('--factor', default=1.5, type=float, help="Event generation weight: E^{factor}, default: 1.5")

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

parser.add_argument('--muon',action='store_true',default=False)

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

mass = args.mass 

pathdir = f'inputs/mass{mass}/fluxnpy'
npylogE = f'{pathdir}/logE.npy'
npyflux = f'{pathdir}/flux.npy'

if args.muon:
    pathdir='inputs/muon'
    npyflux = f'{pathdir}/flux_fullzen.npy'

factor = args.factor

if os.path.isfile(npylogE) * os.path.isfile(npyflux):
    if not os.path.isdir(pathdir):
        os.makedirs(pathdir)
    with open(npylogE, 'rb') as f:
        logEs = np.load(npylogE)
    with open(npyflux, 'rb') as f:
        fluxes = np.load(npyflux)
else:
    logEs, fluxes = makeicflux2d()
    with open(npylogE, 'wb') as f:
        np.save(f, logEs)
    with open(npyflux, 'wb') as f:
        np.save(f, fluxes)
weight = 1

class mySimpleStau(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.STauMinus)
        self.AddParameter("EnergyMin", "", 10000*I3Units.GeV)
        self.AddParameter("EnergyMax", "", 10000*I3Units.GeV)
        self.AddParameter("ZenithMin", "", 90*I3Units.degree)
        self.AddParameter("ZenithMax", "", 180*I3Units.degree)
        self.AddParameter("AzimuthMin", "", 0*I3Units.degree)
        self.AddParameter("AzimuthMax", "", 360*I3Units.degree)
        self.AddParameter("DiskRadius", "", 600 * I3Units.m)
        self.AddParameter("SphereRadius", "", 600 * I3Units.m)
        self.AddParameter("NEvents", "", 100)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.particleType = self.GetParameter("Type")
        self.energyMin = self.GetParameter("EnergyMin")
        self.energyMax = self.GetParameter("EnergyMax")
        self.zenithMin = self.GetParameter("ZenithMin")
        self.zenithMax = self.GetParameter("ZenithMax")
        self.azimuthMin = self.GetParameter("AzimuthMin")
        self.azimuthMax = self.GetParameter("AzimuthMax")
        self.diskRadius = self.GetParameter("DiskRadius")
        self.sphereRadius = self.GetParameter("SphereRadius")
        self.nEvents = self.GetParameter("NEvents")

    def DAQ(self, frame):
        cos_zen_low  = math.cos(self.zenithMin / I3Units.radian)
        cos_zen_high = math.cos(self.zenithMax / I3Units.radian )

        #energyGeV = self.rs.uniform(self.energyMin,self.energyMax)
        def inversePolyPDF(x,alpha,xmin,xmax):
            A = (-alpha+1)/(xmax**(-alpha+1)-xmin**(-alpha+1))
            C = -A/(-alpha+1)*xmin**(-alpha+1)
            return ((x-C)*(-alpha+1)/A)**(1/(-alpha+1))

        #log_energyGeV = self.rs.uniform(np.log10(self.energyMin),np.log10(self.energyMax))
        random_param = random.random()
        energyGeV = inversePolyPDF(random_param, factor, self.energyMin, self.energyMax)
        log_energyGeV = np.log10(energyGeV)

        azi = self.rs.uniform(self.azimuthMin,self.azimuthMax)
        zen = math.acos(self.rs.uniform(cos_zen_low,cos_zen_high))

        r = self.ConstructPerpVector(zen,azi) * math.sqrt(self.rs.uniform(0,self.diskRadius**2))

        diskCenter = self.sphereRadius * numpy.array([math.sin(zen) * math.cos(azi),\
                                                      math.sin(zen) * math.sin(azi),
                                                      math.cos(zen)])

        pos = diskCenter + r
        
        energy = energyGeV * I3Units.GeV

        daughter = dataclasses.I3Particle()
        daughter.type = self.particleType
        daughter.energy = energy
        daughter.pos = dataclasses.I3Position(pos[0], pos[1], pos[2])
        daughter.dir = dataclasses.I3Direction(zen,azi)
        daughter.time = 0.
        daughter.location_type = dataclasses.I3Particle.LocationType.InIce

        primary = dataclasses.I3Particle()
        primary.type = dataclasses.I3Particle.ParticleType.STauMinus
        primary.energy = energy
        primary.pos = dataclasses.I3Position(pos[0], pos[1], pos[2])
        primary.dir = dataclasses.I3Direction(zen,azi)
        primary.time = 0.
        primary.location_type = dataclasses.I3Particle.LocationType.Anywhere

        mctree = dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary,daughter)

        frame["I3MCTree_preStauProp"] = mctree

        weight = getTheWeight(log_energyGeV, int(180-zen/np.pi*180), logEs, fluxes)
        #print(weight[0]*10**(args.minloge-5))
        frame['EventWeight'] = dataclasses.I3String(f'{weight[0]}')

        self.PushFrame(frame)

    def ConstructPerpVector(self, zenith, azimuth):
        x = math.sin(zenith) * math.cos(azimuth)
        y = math.sin(zenith) * math.sin(azimuth)
        z = math.cos(zenith)
        
        v = numpy.array([x,y,z])
        
        # construct another vector in a random direction
        ru_azimuth = self.rs.uniform(0,2.* math.pi)
        ru_zenith = math.acos(self.rs.uniform(-1.0,1.0))
        
        xi = math.sin(ru_zenith) * math.cos(ru_azimuth)
        yi = math.sin(ru_zenith) * math.sin(ru_azimuth)
        zi = math.cos(ru_zenith)
        
        vi = numpy.array([xi,yi,zi])
        
        # calculate the displacement vector from the center of the disk
        # construct a vector in the disk plane
        v_cross_vi = numpy.cross(v,vi)
        # normalize the vector
        return v_cross_vi / math.sqrt(numpy.dot(v_cross_vi,v_cross_vi))


def stau_generation(args):
    tray = I3Tray()

    if args.seed == -1:
        seed = int(random.uniform(0,10000))
    else:
        seed = 8888
    
    # a random number generator
    randomService = phys_services.I3SPRNGRandomService(
        seed = seed,
        nstreams = 1000000,
        streamnum = args.run)
    
    tray.AddModule("I3InfiniteSource","streams",
                   Prefix=args.gcd,
                   Stream=icetray.I3Frame.DAQ)
    
    tray.AddModule("I3MCEventHeaderGenerator","gen_header",
                   Year=2009,
                   DAQTime=158100000000000000,
                   RunNumber=1,
                   EventID=1,
                   IncrementEventID=True)
    
    energyMin = 10**(args.minloge)
    energyMax = 10**(args.maxloge)

    particle = dataclasses.I3Particle.ParticleType.STauMinus
    if args.muon:
        particle = dataclasses.I3Particle.ParticleType.MuMinus
    
    tray.AddModule(mySimpleStau, "injectStau",
                   I3RandomService = randomService,
                   #Type = dataclasses.I3Particle.ParticleType.STauMinus,
                   Type = particle,
                   NEvents = args.nevents,
                   EnergyMin = energyMin*I3Units.GeV,
                   EnergyMax = energyMax*I3Units.GeV,
                   DiskRadius = 200.*I3Units.m,
                   SphereRadius = 800.*I3Units.m
                   #DiskRadius = 1200 *I3Units.m,
                   #SphereRadius = 1700 *I3Units.m
                   )

    tray.AddSegment(PropagateStaus, "PropagateStaus",
                    RandomService = randomService, 
                    debug = args.debug)
                   
    outfilename = args.out
    if args.weightshow:
        outfilename = f'{(args.out).split(".i")[0]}_r{args.minloge:.2f}-{args.maxloge:.2f}_inverseEweighted_run{str(args.run).zfill(8)}.i3'
    
    tray.AddModule("I3Writer","writer",Filename = outfilename)
    tray.Execute(args.nevents+3)
    
    return outfilename
    
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

    tray.Add('i3ppc','ppc',
            MCTree = "I3MCTree",
            infoName = "PPCInfoDict", # create a dict
            verbose=args.debug,
            photons=True,
            pseries=True)

    tray.Add("Rename", "PPCRename", Keys=["MCPESeriesMap", "I3MCPESeriesMap"])

    outfilename = infilename.split('.i')[0] + '_ppc.i3'

    tray.Add('I3Writer', filename=outfilename, streams=[icetray.I3Frame.DAQ])
    tray.Execute()

if __name__=='__main__':
    outfile = stau_generation(args)
    #stautrack_ppc(args, outfile)
