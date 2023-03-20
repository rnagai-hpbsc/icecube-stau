#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_staus.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default=expandvars("$I3_TESTDATA/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=1,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("--no-apply-mmc", action="store_false", default=True,
                  dest="APPLYMMC", help="do not apply MMC to the I3MCTree after generating the muons")
parser.add_option("-e", "--energy", type="float", default=10,
                  dest="ENERGY", help="input energy in TeV")
parser.add_option("--loge", type="float", default=None, 
                  dest="LOGE", help="if log is better, it can be used (priority)")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

from I3Tray import *
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, sim_services

from PropagateStaus import PropagateStaus

#
import math
import numpy

import validation

class mySimpleStau(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.STauMinus)
        self.AddParameter("EnergyMin", "", 10*I3Units.TeV)
        self.AddParameter("EnergyMax", "", 10*I3Units.TeV)
        self.AddParameter("ZenithMin", "", 0*I3Units.degree)
        self.AddParameter("ZenithMax", "", 90*I3Units.degree)
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

        nlogE, nflux = norm_icflux(np.log10(self.energyMin), np.log10(self.energyMax), *validation.makeicflux(0))
        self.normflux = nflux
        self.normLogE = nlogE

    def DAQ(self, frame):

        validFlag = False
        while not validFlag:
            azi = self.rs.uniform(self.azimuthMin,self.azimuthMax)

            cos_zen_low = math.cos(self.zenithMin / I3Units.radian)
            cos_zen_high = math.cos(self.zenithMax / I3Units.radian )
            zen = math.acos(self.rs.uniform(cos_zen_low,cos_zen_high))

            r = self.ConstructPerpVector(zen,azi) * math.sqrt(self.rs.uniform(0,self.diskRadius**2))

            diskCenter = self.sphereRadius * numpy.array([math.sin(zen) * math.cos(azi),\
                                                          math.sin(zen) * math.sin(azi),
                                                          math.cos(zen)])

            pos = diskCenter + r
    
            # set the particle's energy
            energy = self.rs.uniform(self.energyMin,self.energyMax) * I3Units.GeV        
            validFlag = validation.valid_event(np.log10(energy),zen,self.normLogE,self.normflux)

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

        self.PushFrame(frame)

    def ConstructPerpVector(self, zenith,azimuth):
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

    def isValid():
        return True




tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

if options.LOGE is not None:
    primary_energy_in_TeV = 10**options.LOGE
else:
    primary_energy_in_TeV = options.ENERGY

tray.AddModule(mySimpleStau, "injectStau",
               I3RandomService = randomService,
               Type = dataclasses.I3Particle.ParticleType.STauMinus,
               NEvents = options.NUMEVENTS,
               EnergyMin = energyMin_in_TeV*I3Units.TeV,
               EnergyMax = energyMax_in_TeV*I3Units.TeV,
               DiskRadius = 200.*I3Units.m,
               SphereRadius = 800.*I3Units.m
               #DiskRadius = 1200 *I3Units.m,
               #SphereRadius = 1700 *I3Units.m
               )


if options.APPLYMMC:
    tray.AddSegment(PropagateStaus, "PropagateStaus",
            RandomService = randomService)
               

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)



tray.Execute(options.NUMEVENTS+3)




