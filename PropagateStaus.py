#!/usr/bin/env python3

from os.path import expandvars
from icecube import icetray, dataclasses
from I3Tray import I3Units

import math

def MakePropagator(mediadef=None):
    """
    Create a muon propagator service.

    """
    from icecube import sim_services, PROPOSAL
    # in PROPOSAL everything can be defined in the configuration file
    if mediadef is None:
        mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/config_icesim.json')
    return PROPOSAL.I3PropagatorServicePROPOSAL(config_file=mediadef)

def DefinePropagators():
    from icecube import sim_services, phys_services, simclasses, dataclasses
    from icecube import cmc
    propagator      = MakePropagator(mediadef='copy_config_icesim.json')
    cascadePropagator = cmc.I3CascadeMCService(phys_services.I3GSLRandomService(1)) # dummy RNG
    
    propagators = sim_services.I3ParticleTypePropagatorServiceMap()
    for pt in 'MuMinus', 'MuPlus', 'TauMinus', 'TauPlus', 'STauMinus', 'STauPlus',\
              'STauMinus100', 'STauPlus100',\
              'STauMinus200', 'STauPlus200',\
              'STauMinus300', 'STauPlus300',\
              'STauMinus400', 'STauPlus400',\
              'STauMinus500', 'STauPlus500',\
              'STauMinus600', 'STauPlus600',\
              'STauMinus700', 'STauPlus700',\
              'STauMinus800', 'STauPlus800',\
              'STauMinus900', 'STauPlus900',:
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = propagator
    for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cascadePropagator
    
    return propagators
 
@icetray.traysegment
def PropagateStaus(tray, name, particle_type=None,
    RandomService = None, debug=False):

    from I3Tray import I3Units

    from icecube import icetray, dataclasses, phys_services, sim_services, simclasses
    from icecube import cmc

    if particle_type is None:
        particle_type = dataclasses.I3Particle.STauMinus
    propagator      = MakePropagator(mediadef='copy_config_icesim.json')
    cascadePropagator = cmc.I3CascadeMCService(phys_services.I3GSLRandomService(1)) # dummy RNG

    print(propagator)
    stau = dataclasses.I3Particle()
    stau.type = particle_type
    stau.pos = dataclasses.I3Position(0,0,0)
    stau.dir = dataclasses.I3Direction(0,0)
    stau.energy = 100 * I3Units.TeV
    stau.time = 0*I3Units.ns
    stau.location_type = dataclasses.I3Particle.InIce
    propagator.register_particletype(stau.type)
    if debug:
        print(stau)
        print(propagator.Propagate(stau))

    # set up propagators
    propagators = sim_services.I3ParticleTypePropagatorServiceMap()
    for pt in 'MuMinus', 'MuPlus', 'TauMinus', 'TauPlus', 'STauMinus', 'STauPlus',\
              'STauMinus100', 'STauPlus100',\
              'STauMinus200', 'STauPlus200',\
              'STauMinus300', 'STauPlus300',\
              'STauMinus400', 'STauPlus400',\
              'STauMinus500', 'STauPlus500',\
              'STauMinus600', 'STauPlus600',\
              'STauMinus700', 'STauPlus700',\
              'STauMinus800', 'STauPlus800',\
              'STauMinus900', 'STauPlus900',:
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = propagator
    for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cascadePropagator

    tray.AddModule('I3PropagatorModule', name+'_propagator',
        PropagatorServices=propagators,
        RandomService=RandomService,
        RNGStateName="RNGState",
        InputMCTreeName="I3MCTree_preStauProp",
        OutputMCTreeName="I3MCTree")

    # add empty MMCTrackList objects for events that have none
    def addEmptyMMCTrackList(frame):
        if "MMCTrackList" not in frame:
            frame["MMCTrackList"] = simclasses.I3MMCTrackList()
    tray.AddModule(addEmptyMMCTrackList, name+'_addEmptyMMCTrackList',
        Streams=[icetray.I3Frame.DAQ])


if __name__=="__main__":
    from optparse import OptionParser
    from os.path import expandvars

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile",default="test_muons_propagated.i3",
                      dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
    parser.add_option("-s", "--seed",type="int",default=12345,
                      dest="SEED", help="Initial seed for the random number generator")
    parser.add_option("-r", "--runnumber", type="int", default=1,
                      dest="RUNNUMBER", help="The run number for this simulation")

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

    from icecube import icetray, dataclasses, dataio, phys_services

    tray = I3Tray()

    # set up a random number generator
    randomService = phys_services.I3SPRNGRandomService(
        seed = options.SEED*2,
        nstreams = 10000,
        streamnum = options.RUNNUMBER)


    # re-use the same RNG for modules that need it on the context
    tray.context['I3RandomService'] = randomService

    tray.AddSegment(PropagateStaus, "PropagateStaus",
        RandomService = randomService, debug=True)

    tray.AddModule("I3Writer","writer",
        Filename = options.OUTFILE)

    tray.Execute(10)

    print(tray)


