import os, sys
from I3Tray import *
from icecube import icetray, dataclasses, simclasses, dataio, phys_services, sim_services

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

import glob

import click

from tqdm import tqdm

@click.command()
@click.option('--endn',type=int,default=None)
@click.option('--files',type=str,default='../../cobalt_data/Stau150_dataset-1_202303/*0_ppc.i3')
@click.option('--datadir',type=str,default=None)
@click.option('--alpha',type=float,default=1.5)
def main(endn,files,datadir,alpha):
    if datadir is None:
        filelist = glob.glob(files)
    else:
        filelist = glob.glob(f'{datadir}/*.i3')
    eidists = {}
    eldists = {}
    zeniths = {}
    coszens = {}
    azimuths = {}
    weights = {}
    eventweights = {}

    ET = {}

    i = 0
    for fname in tqdm(filelist):
        i += 1
        tqdm.write(fname)
        i3f = dataio.I3File(fname)

        try: 
            weight = float(fname.split('_w')[1].split('_')[0])
        except:
            weight = 1

        try: 
            hrange = float(fname.split('_r')[1].split('_')[0].split('-')[0])
            lrange = float(fname.split('_r')[1].split('_')[0].split('-')[1])
        except:
            hrange = None
            lrange = None

        eidists.setdefault(f'{hrange}-{lrange}',[])
        eldists.setdefault(f'{hrange}-{lrange}',[])
        zeniths.setdefault(f'{hrange}-{lrange}',[])
        coszens.setdefault(f'{hrange}-{lrange}',[])
        azimuths.setdefault(f'{hrange}-{lrange}',[])
        weights.setdefault(f'{hrange}-{lrange}',[])
        eventweights.setdefault(f'{hrange}-{lrange}',[])
        

        energy_info = {'ienergy':[], 'eloss':[], 'logEi': [], 'logEloss':[], 'zenith':[], 'coszen':[], 'azimuth':[], 'eventweight':[]}
        while i3f.more():
            try: 
                frame = i3f.pop_daq()
            except RuntimeError:
                continue
            mmc = frame.Get('MMCTrackList')
            mmc0 = mmc[0]
            mmc0particle = mmc0.GetI3Particle()
            init_energy = mmc0particle.energy
            elost = mmc0.GetElost()
            energy_info['ienergy'].append(init_energy)
            energy_info['eloss'].append(elost)
            energy_info['logEi'].append(np.log10(init_energy))
            energy_info['logEloss'].append(np.log10(elost))
            energy_info['zenith'].append(mmc0particle.dir.zenith/np.pi*180)
            energy_info['azimuth'].append(mmc0particle.dir.azimuth/np.pi*180)
            energy_info['coszen'].append(np.cos(mmc0particle.dir.zenith))
            try: 
                event_weight = float(f"{frame['EventWeight']}".split('I3String("')[-1].split('")')[0])
            except: 
                event_weight = 1
            energy_info['eventweight'].append(event_weight)

            try:
                mctree = frame['I3MCTree']
            except KeyError:
                print('I3MCTree not found')
            else:
                iLogE = int(np.log10(mctree[0].energy)*2)/2
                ET.setdefault(f'{iLogE}',[])
                try:
                    t0 = float(f"{frame['Time_0']}".split('I3String("')[-1].split('")')[0])
                    t1 = float(f"{frame['Time_1']}".split('I3String("')[-1].split('")')[0])
                except:
                    t0 = 0
                    t1 = 999
                timediff = t1 - t0
                ET[f'{iLogE}'].append(timediff)

        i3f.close()

        thisweight = (-1)*np.array(energy_info['ienergy'])**alpha*(10**(hrange*(1-alpha))-10**(lrange*(1-alpha)))/(1-alpha)
        eidists[f'{hrange}-{lrange}'].extend(energy_info['logEi'])
        eldists[f'{hrange}-{lrange}'].extend(energy_info['logEloss'])
        zeniths[f'{hrange}-{lrange}'].extend(energy_info['zenith'])
        coszens[f'{hrange}-{lrange}'].extend(energy_info['coszen'])
        azimuths[f'{hrange}-{lrange}'].extend(energy_info['azimuth'])
        eventweights[f'{hrange}-{lrange}'].extend(energy_info['eventweight'])
        weights[f'{hrange}-{lrange}'].extend(thisweight)
        
        if i==endn:
            break

    ienergies = [float(key) for key in ET]
    tmean = [np.mean(np.array(ET[key])) for key in ET]
    tsigma = [np.std(np.array(ET[key])) for key in ET]
    
    if endn is None:
        pdf = PdfPages('basic_plots_all.pdf')
    else:
        pdf = PdfPages(f'basic_plots_{endn}.pdf')

    plt.errorbar(ienergies, tmean, yerr=tsigma,capsize=5,fmt='o')
    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel('Mean time for single event generation [s]')
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    weighting_array = {}
    for key in eidists.keys():
        weighting_array[key] = np.array(weights[key])/len(eidists[key])*np.array(eventweights[key])

    #print(eidists)
    for key in eidists.keys():
        plt.hist(eidists[key],bins=50,range=(5,12),weights=np.ones(len(eidists[key])),stacked=True,label=f'{key}: unweighted')

    #plt.yscale('log')
    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()
    
    for key in eidists.keys():
        plt.hist(eidists[key],bins=64,range=(5,11.4),weights=np.array(weights[key])/len(eidists[key])*np.array(eventweights[key]),stacked=True,label=f'{key}: weighted')
        print(len(eidists[key]))

    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel('$\mathrm{cm^{-2}\ s^{-1}\ sr^{-2}\ GeV^{-1}}$')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in eldists.keys():
        plt.hist(eldists[key],bins=64,range=(0,11.4),weights=np.ones(len(eldists[key])),stacked=True,label=f'{key}: unweighted')

    plt.xlabel('$\log_{10}(E_\mathrm{loss}\ [\mathrm{GeV}])$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in eldists.keys():
        plt.hist(eldists[key],bins=64,range=(0,11.4),weights=np.array(weights[key])/len(eidists[key])*np.array(eventweights[key]),stacked=True,label=f'{key}: weighted')

    plt.xlabel('$\log_{10}(E_\mathrm{loss}\ [\mathrm{GeV}])$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in zeniths.keys():
        plt.hist(zeniths[key],bins=50,range=(0,180),weights=np.ones(len(zeniths[key]))*weights[key],stacked=True,label=f'{key}: unweighted')

    plt.xlabel('Zenith [deg]')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in zeniths.keys():
        plt.hist(zeniths[key],bins=50,range=(0,180),weights=weighting_array[key],stacked=True,label=f'{key}: weighted')

    plt.xlabel('Zenith [deg]')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in azimuths.keys():
        plt.hist(azimuths[key],bins=50,range=(0,360),weights=np.ones(len(azimuths[key])),stacked=True,label=f'{key}: unweighted')

    plt.xlabel('Azimuth [deg]')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in azimuths.keys():
        plt.hist(azimuths[key],bins=50,range=(0,360),weights=weighting_array[key],stacked=True,label=f'{key}: weighted')

    plt.xlabel('Azimuth [deg]')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in coszens.keys():
        plt.hist(coszens[key],bins=50,range=(-1,1),weights=np.ones(len(coszens[key])),stacked=True,label=f'{key}: unweighted')

    plt.xlabel('$\cos\\theta_{zen}$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in coszens.keys():
        plt.hist(coszens[key],bins=50,range=(-1,1),weights=weighting_array[key],stacked=True,label=f'{key}: weighted')

    plt.xlabel('$\cos\\theta_{zen}$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    pdf.close()

if __name__=='__main__':
    main()
