import os, sys, collections
from I3Tray import *
from icecube import icetray, dataclasses, simclasses, dataio, phys_services, sim_services

import click

import numpy as np
import matplotlib.pyplot as plt

filename = 'data/RunTest1/IncTest3_r5.00-6.00_w1.0371670564530074e-25_run00000000_ppc.i3'

@click.group()
def cli():
    pass

@cli.command()
def read_one():
    i3f = dataio.I3File(filename)

    while i3f.more():
        frame = i3f.pop_frame()
        mcpe = frame['I3MCPESeriesMap']
        for key in mcpe.keys():
            print(mcpe[key].value)
    
    i3f.close()


@cli.command()
def read_tray():
    tray = I3Tray()
    
    tray.Add("I3Reader", filename=filename)
    
    def read(frame):
        print('do something')
        for key in frame.keys():
            fr = frame[key]
            #print(fr)
    
    tray.AddModule(read,"read",streams=[icetray.I3Frame.DAQ])
    
    tray.Execute()

@cli.command()
def mcpe_read():
    i3f = dataio.I3File(filename)

    while i3f.more():
        frame = i3f.pop_frame()
        mcpes = frame.Get('I3MCPESeriesMap')
        for key in mcpes.keys():
            for mcpe in mcpes[key]:
                print(f'{key}, npe: {mcpe.npe}, time: {mcpe.time}, id: {mcpe.ID}')
            

if __name__=='__main__':
    cli()
