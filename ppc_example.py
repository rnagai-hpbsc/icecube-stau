import os
import argparse
import numpy as np
from I3Tray import I3Tray, load
from icecube import icetray, clsim, phys_services, simclasses, photonics_service


from PropagateStaus import DefinePropagators

icetray.set_log_level(icetray.I3LogLevel.LOG_TRACE)


def main():
    """ resim for monopod
    """
    parser = argparse.ArgumentParser(
        description='Resimulate benchmark i3s (e.g. under different ice)')
    parser.add_argument('infiles', nargs='+')
    parser.add_argument('-o', '--out', default='out.i3.zst',
                        help='output file')
    parser.add_argument('-s', '--seed', default=8888, type=int,
                        help='seed for sprng')
    parser.add_argument('-r', '--run', default=0, type=int,
                        help='run for sprng')
    parser.add_argument('--nframes', type=int, default=None, help='number of frames to process')
    parser.add_argument('--usecpus', default=False, action='store_true')
    parser.add_argument('--domos', default=5, type=float)
    parser.add_argument('--domeff', default=1., type=float)
    parser.add_argument('--gcd', default=os.path.os.path.expandvars('$I3_DATA/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz'),
                        type=str,
                        help='gcd file')
    parser.add_argument('--icemodel', choices=('SpiceMie', 'SpiceLea', 'Spice3.2', 'Spice3.2.1',
                                               'Spicebfr-v1', 'Spicebfr-v2', 'Spicebfr-v2_flat',
                                               'Spicebfr-v2_tilt'),
                        default='Spicebfr-v2', help='The ice model to use for photon propagation.')
    parser.add_argument("--holeice",default='as.h2-50cm', choices=('as.flasher_p1_0.30_p2_0',
                                                                   'as.flasher_p1_0.30_p2_-1',
                                                                   'as.flasher_p1_0.30_p2_+1',
                                                                   'as.flasher_p1_0.30_p2_-2',
                                                                   'as.flasher_p1_0.30_p2_-3',
                                                                   'as.h2-50cm',
                                                                   'as.nominal',
                                                                   'as.set0_p0=0.0_p1=0.0'), help="Holeice file")
    parser.add_argument('--nocascadeextension', default=False, action='store_true')
    parser.add_argument('--notilt', default=False, action='store_true')
    parser.add_argument('--unweighted', default=False, action='store_true')
    parser.add_argument('--scale', default=None, type=float)
    parser.add_argument('--all', default=False, action='store_true')
    parser.add_argument('--prescale', default=0.01, type=float)
    parser.add_argument('--history', default=0, type=int)
    parser.add_argument('--prop',default=None, type=str)
    parser.add_argument('--usegpus',type=int,default=-1)
    args = parser.parse_args()

    print("RNG being initiated...")
    # set up a random number generator
    randomService = phys_services.I3SPRNGRandomService(
        seed = args.seed,
        nstreams = 100000000,
        streamnum = args.run)

    load('ppc')
    tray = I3Tray()

    tray.context['I3RandomService'] = randomService
    tray.Add('I3Reader', Filenamelist=[args.gcd]+args.infiles)

    tray.Add('i3ppc','ppc',
                     #gpu = args.usegpus,
                     MCTree = "I3MCTree",
                     infoName = "PPCInfoDict", # create a dict
                     verbose=True,
                     photons=True,
                     pseries=True)
        
    tray.Add("Rename", "PPCRename", Keys=["MCPESeriesMap", "I3MCPESeriesMap"])

    tray.Add('I3Writer', filename=args.out,
             streams=[icetray.I3Frame.DAQ])
    if args.nframes is None:
        tray.Execute()
    else:
        tray.Execute(args.nframes)


if __name__=='__main__':
    main()
