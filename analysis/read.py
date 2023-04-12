import os, sys, collections
from I3Tray import *
from icecube import icetray, dataclasses, simclasses, dataio, phys_services, sim_services

import click

import numpy as np
import matplotlib.pyplot as plt

DEFAULTFILE = 'data/RunTest1/IncTest3_r5.00-6.00_w1.0371670564530074e-25_run00000000_ppc.i3'

@click.group()
def cli():
    pass

@cli.command()
def read_one():
    filename = DEFAULTFILE
    i3f = dataio.I3File(filename)

    while i3f.more():
        frame = i3f.pop_frame()
        mcpe = frame['I3MCPESeriesMap']
        for key in mcpe.keys():
            print(mcpe[key].value)
    
    i3f.close()


@cli.command()
def read_tray():
    filename = DEFAULTFILE
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
    filename = DEFAULTFILE
    i3f = dataio.I3File(filename)

    while i3f.more():
        frame = i3f.pop_frame()
        mcpes = frame.Get('I3MCPESeriesMap')
        for key in mcpes.keys():
            for mcpe in mcpes[key]:
                print(f'{key}, npe: {mcpe.npe}, time: {mcpe.time}, id: {mcpe.ID}')
            
@cli.command()
def mctree_read():
    filename = DEFAULTFILE
    i3f = dataio.I3File(filename)

    while i3f.more():
        frame = i3f.pop_frame()
        mctree = frame.Get('I3MCTree')
        for mcp in mctree:
            print(mcp.major_id, mcp.minor_id, mcp.dir, mcp.pos, mcp.time, mcp.energy, mcp.speed, mcp.length, mcp.type, mcp.pdg_encoding)
            try:
                mass = mcp.mass
            except: 
                print('no mass')
            else:
                print(mass)

@cli.command()
@click.option('--filename',default=DEFAULTFILE)
def energy_plot(filename):
    i3f = dataio.I3File(filename)

    try: 
        weight = float(filename.split('_w')[1].split('_')[0])
    except: 
        weight = 1

    energy_info = {'ienergy':[], 'eloss':[], 'logEi': [], 'logEloss':[]}
    while i3f.more():
        frame = i3f.pop_frame()
        mmc = frame.Get('MMCTrackList')
        mmc0 = mmc[0]
        mmc0particle = mmc0.GetI3Particle()
        init_energy = mmc0particle.energy
        elost = mmc0.GetElost()
        energy_info['ienergy'].append(init_energy)
        energy_info['eloss'].append(elost)
        energy_info['logEi'].append(np.log10(init_energy))
        energy_info['logEloss'].append(np.log10(elost))
    
    print(energy_info)
    h_logEi, edge = np.histogram(energy_info['logEi'],bins=10,range=(5,7))
    h_logEi_weighted = h_logEi * weight / len(h_logEi)
    bincenter = [(edge[i]+edge[i+1])/2 for i in range(len(edge)-1)]
    plt.plot(bincenter, h_logEi_weighted,ds='steps-mid')
    plt.yscale('log')
    plt.savefig('hist.pdf')
    plt.show()

@cli.command()
@click.option('--filename',default=DEFAULTFILE)
@click.option('--weighted','-w',is_flag=True,default=False)
def initial_energy_plot(filename,weighted):
    i3f = dataio.I3File(filename)

    try: 
        rangeinfo = filename.split('_r')[1].split('_')[0]
    except: 
        logEmin = 5
        logEmax = 12
    else: 
        logEmin = float(rangeinfo.split('-')[0])
        logEmax = float(rangeinfo.split('-')[1])

    energies = []
    weights = []
    while i3f.more():
        frame = i3f.pop_frame()
        try:
            mctree = frame['I3MCTree']
        except KeyError: 
            print("I3MCTree not found!")
            continue
        ienergy = mctree[0].energy
        energies.append(np.log10(ienergy))
        weight = 1
        if weighted: 
            try: 
                #weight = float(f"{frame['EventWeight']}".split('I3Double(')[-1].split(')')[0])
                weight = float(f"{frame['EventWeight']}".split('I3String("')[-1].split('")')[0])
            except KeyError:
                weight = 1
        weights.append(weight)
    print(energies, weights)
    weights = np.array(weights)/len(weights)*((logEmax-logEmin)*10)

    h_energies, edge = np.histogram(energies,bins=70,range=(5,12),weights=weights)
    bincenter = [(edge[i]+edge[i+1])/2 for i in range(len(edge)-1)]
    plt.plot(bincenter, h_energies, ds='steps-mid')
    plt.yscale('log')
    plt.savefig('analysis/eslope_test_init_energy.pdf')
    plt.show()

if __name__=='__main__':
    cli()
