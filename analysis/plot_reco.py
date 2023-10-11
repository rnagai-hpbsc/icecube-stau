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
@click.option('--files',type=str,default='outdata/Stau150_Muon_20230518_Gen_Events/Stau/*.i3')
@click.option('--datadir',type=str,default=None)
@click.option('--alpha',type=float,default=1.5)
@click.option('--out',type=str,default=None)
def main(endn,files,datadir,alpha,out):
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

    reco_energy = {}
    reco_zenith = {}
    reco_azimuth = {}
    reco_coszen = {}

    ET = {}

    truth_keys = ['ienergy','eloss','logEi','logEloss','zenith','coszen','azimuth','eventweight','weight']
    reco_keys = ['energy','zenith','coszen','azimuth']

    truth_info = {}
    reco_info = {}

    for key in truth_keys:
        truth_info.setdefault(key,[])
    for key in reco_keys:
        reco_info.setdefault(key,[])

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
            hrange = 5
            lrange = 12

        eidists.setdefault(f'{hrange}-{lrange}',[])
        eldists.setdefault(f'{hrange}-{lrange}',[])
        zeniths.setdefault(f'{hrange}-{lrange}',[])
        coszens.setdefault(f'{hrange}-{lrange}',[])
        azimuths.setdefault(f'{hrange}-{lrange}',[])
        weights.setdefault(f'{hrange}-{lrange}',[])
        eventweights.setdefault(f'{hrange}-{lrange}',[])
        
        reco_energy.setdefault(f'{hrange}-{lrange}',[])
        reco_zenith.setdefault(f'{hrange}-{lrange}',[])
        reco_azimuth.setdefault(f'{hrange}-{lrange}',[])
        reco_coszen.setdefault(f'{hrange}-{lrange}',[])

        truth_energy_info = {}
        reco_energy_info = {}

        for key in truth_keys:
            truth_energy_info.setdefault(key,[])
        for key in reco_keys:
            reco_energy_info.setdefault(key,[])
        #truth_energy_info = {'ienergy':[], 'eloss':[], 'logEi': [], 'logEloss':[], 'zenith':[], 'coszen':[], 'azimuth':[], 'eventweight':[]}
        #reco_energy_info = {'energy':[], 'zenith':[], 'coszen':[], 'azimuth':[]}

        while i3f.more():
            frame = i3f.pop_physics()

            try:
                mmc = frame.Get('MMCTrackList')
                mctree = frame.Get('I3MCTree')
            except KeyError:
                continue
            try:
                sp = frame.Get('OnlineL2_SplineMPE_MuEx')
            except KeyError:
                sp = None
            mmc0 = mmc[0]
            mmc0particle = mmc0.GetI3Particle()
            init_energy = mmc0particle.energy
            elost = mmc0.GetElost()
            truth_energy_info['ienergy'].append(init_energy)
            truth_energy_info['eloss'].append(elost)
            truth_energy_info['logEi'].append(np.log10(init_energy))
            truth_energy_info['logEloss'].append(np.log10(elost))
            truth_energy_info['zenith'].append(mmc0particle.dir.zenith/np.pi*180)
            truth_energy_info['azimuth'].append(mmc0particle.dir.azimuth/np.pi*180)
            truth_energy_info['coszen'].append(np.cos(mmc0particle.dir.zenith))
            try: 
                event_weight = float(f"{frame['EventWeight']}".split('I3String("')[-1].split('")')[0])
            except: 
                event_weight = 1
            truth_energy_info['eventweight'].append(event_weight)

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

            if sp is not None:
                reco_energy_info['energy'].append(sp.energy)
                reco_energy_info['zenith'].append(sp.dir.zenith/np.pi*180)
                reco_energy_info['azimuth'].append(sp.dir.azimuth/np.pi*180)
                reco_energy_info['coszen'].append(np.cos(sp.dir.zenith))
            else:
                reco_energy_info['energy'].append(np.nan)
                reco_energy_info['zenith'].append(np.nan)
                reco_energy_info['azimuth'].append(np.nan)
                reco_energy_info['coszen'].append(np.nan)

        i3f.close()

        for key in truth_keys:
            truth_info[key].extend(truth_energy_info[key])

        thisweight = (-1)*np.array(truth_energy_info['ienergy'])**alpha*(10**(hrange*(1-alpha))-10**(lrange*(1-alpha)))/(1-alpha)
        truth_info['weight'].extend(thisweight)

        for key in reco_keys:
            reco_info[key].extend(reco_energy_info[key])

        #eidists[f'{hrange}-{lrange}'].extend(truth_energy_info['logEi'])
        #eldists[f'{hrange}-{lrange}'].extend(truth_energy_info['logEloss'])
        #zeniths[f'{hrange}-{lrange}'].extend(truth_energy_info['zenith'])
        #coszens[f'{hrange}-{lrange}'].extend(truth_energy_info['coszen'])
        #azimuths[f'{hrange}-{lrange}'].extend(truth_energy_info['azimuth'])
        #eventweights[f'{hrange}-{lrange}'].extend(truth_energy_info['eventweight'])
        #weights[f'{hrange}-{lrange}'].extend(thisweight)

        #reco_energy[f'{hrange}-{lrange}'].extend(reco_energy_info['energy'])
        #reco_zenith[f'{hrange}-{lrange}'].extend(reco_energy_info['zenith'])
        #reco_azimuth[f'{hrange}-{lrange}'].extend(reco_energy_info['azimuth'])
        #reco_coszen[f'{hrange}-{lrange}'].extend(reco_energy_info['coszen'])
        
        if i==endn:
            break

    ienergies = [float(key) for key in ET]
    tmean = [np.mean(np.array(ET[key])) for key in ET]
    tsigma = [np.std(np.array(ET[key])) for key in ET]
    
    if out is not None:
        suffix = f'_{out}'
    else:
        suffix = ''

    prefix = 'plots/'
    if endn is None:
        pdf = PdfPages(f'{prefix}reco_plots{suffix}_all.pdf')
    else:
        pdf = PdfPages(f'{prefix}reco_plots{suffix}_{endn}.pdf')

    plt.errorbar(ienergies, tmean, yerr=tsigma,capsize=5,fmt='o')
    plt.title('Mean time for a single event generation')
    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel('Time [s]')
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    flux_str = 'Flux [$\mathrm{cm^{-2}\ s^{-1}\ sr^{-2}\ GeV^{-1}}$]'

    #weighting_array = {}
    #for key in eidists.keys():
    #    weighting_array[key] = np.array(weights[key])/len(eidists[key])*np.array(eventweights[key])
    weighting_array = np.array(truth_info['weight'])/len(truth_info['ienergy'])*np.array(truth_info['eventweight'])

    #for key in eidists.keys():
    #    plt.hist(eidists[key],bins=120,range=(0,12),weights=np.ones(len(eidists[key])),stacked=True,label=f'{key}')
    plt.hist(truth_info['logEi'],bins=120,range=(0,12),label='Truth')

    plt.title('Injection Energy: unweighted')
    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()
    
    #for key in eidists.keys():
    #    plt.hist(eidists[key],bins=120,range=(0,12),weights=np.array(weights[key])/len(eidists[key])*np.array(eventweights[key]),stacked=True,label=f'{key}')
    #    print(len(eidists[key]))
    plt.hist(truth_info['logEi'],bins=120,range=(0,12),weights=weighting_array,label='Truth')

    plt.title('Injection Energy: weighted')
    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel(flux_str)
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in eldists.keys():
    #    plt.hist(eldists[key],bins=64,range=(0,11.4),weights=np.ones(len(eldists[key])),stacked=True,label=f'{key}')
    plt.hist(truth_info['logEloss'],bins=64,range=(0,11.4),label='Truth')

    plt.title('Total energy loss: unweighted')
    plt.xlabel('$\log_{10}(E_\mathrm{loss}\ [\mathrm{GeV}])$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in eldists.keys():
    #    plt.hist(eldists[key],bins=64,range=(0,11.4),weights=np.array(weights[key])/len(eidists[key])*np.array(eventweights[key]),stacked=True,label=f'{key}')
    plt.hist(truth_info['logEloss'],bins=64,range=(0,11.4),weights=weighting_array,label='Truth')

    plt.title('Total energy loss: weighted')
    plt.xlabel('$\log_{10}(E_\mathrm{loss}\ [\mathrm{GeV}])$')
    plt.ylabel(flux_str)
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_zenith.keys():
    #    plt.hist(reco_zenith[key],bins=50,range=(0,180),weights=np.ones(len(reco_zenith[key])),label=f'Reco',histtype='step')
    #for key in zeniths.keys():
    #    plt.hist(zeniths[key],bins=50,range=(0,180),weights=np.ones(len(zeniths[key])),label=f'Truth',histtype='step')
    plt.hist(reco_info['zenith'],bins=50,range=(0,180),label='Reco',histtype='step')
    plt.hist(truth_info['zenith'],bins=50,range=(0,180),label='Truth',histtype='step')

    plt.title('Zenith: unweighted')
    plt.xlabel('Zenith [deg]')
    plt.ylabel('#Events')
    plt.legend(loc='upper left')
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_zenith.keys():
    #    plt.hist(reco_zenith[key],bins=50,range=(0,180),weights=weighting_array[key],label=f'Reco',histtype='step')
    #for key in zeniths.keys():
    #    plt.hist(zeniths[key],bins=50,range=(0,180),weights=weighting_array[key],label=f'Truth',histtype='step')
    plt.hist(reco_info['zenith'],bins=50,range=(0,180),weights=weighting_array,label='Reco',histtype='step')
    plt.hist(truth_info['zenith'],bins=50,range=(0,180),weights=weighting_array,label='Truth',histtype='step')

    plt.title('Zenith: weighted')
    plt.xlabel('Zenith [deg]')
    plt.ylabel(flux_str)
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_azimuth.keys():
    #    plt.hist(reco_azimuth[key],bins=50,range=(0,360),weights=np.ones(len(reco_azimuth[key])),label=f'Reco',histtype='step')
    #for key in azimuths.keys():
    #    plt.hist(azimuths[key],bins=50,range=(0,360),weights=np.ones(len(azimuths[key])),label=f'Truth',histtype='step')
    plt.hist(reco_info['azimuth'],bins=50,range=(0,360),label='Reco',histtype='step')
    plt.hist(truth_info['azimuth'],bins=50,range=(0,360),label='Truth',histtype='step')

    plt.title('Azimuth: unweighted')
    plt.xlabel('Azimuth [deg]')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_azimuth.keys():
    #    plt.hist(reco_azimuth[key],bins=50,range=(0,360),weights=weighting_array[key],label=f'Reco',histtype='step')
    #for key in azimuths.keys():
    #    plt.hist(azimuths[key],bins=50,range=(0,360),weights=weighting_array[key],label=f'Truth',histtype='step')
    plt.hist(reco_info['azimuth'],bins=50,range=(0,360),weights=weighting_array,label='Reco',histtype='step')
    plt.hist(truth_info['azimuth'],bins=50,range=(0,360),weights=weighting_array,label='Truth',histtype='step')

    plt.title('Azimuth: weighted')
    plt.xlabel('Azimuth [deg]')
    plt.ylabel(flux_str)
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_coszen.keys():
    #    plt.hist(reco_coszen[key],bins=50,range=(-1,1),weights=np.ones(len(reco_coszen[key])),label=f'Reco',histtype='step')
    #for key in coszens.keys():
    #    plt.hist(coszens[key],bins=50,range=(-1,1),weights=np.ones(len(coszens[key])),label=f'Truth',histtype='step')
    plt.hist(reco_info['coszen'],bins=50,range=(-1,1),label='Reco',histtype='step')
    plt.hist(truth_info['coszen'],bins=50,range=(-1,1),label='Truth',histtype='step')

    plt.title('Cosine of Zenith: unweighted')
    plt.xlabel('$\cos\\theta_{zen}$')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_coszen.keys():
    #    plt.hist(reco_coszen[key],bins=50,range=(-1,1),weights=weighting_array[key],label=f'Reco',histtype='step')
    #for key in coszens.keys():
    #    plt.hist(coszens[key],bins=50,range=(-1,1),weights=weighting_array[key],label=f'Truth',histtype='step')
    plt.hist(reco_info['coszen'],bins=50,range=(-1,1),weights=weighting_array,label='Reco',histtype='step')
    plt.hist(truth_info['coszen'],bins=50,range=(-1,1),weights=weighting_array,label='Truth',histtype='step')

    plt.title('Cosine of Zenith: weighted')
    plt.xlabel('$\cos\\theta_{zen}$')
    plt.ylabel(flux_str)
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_energy.keys():
    #    plt.hist(np.log10(np.array(reco_energy[key])),bins=120,range=(0,12),weights=np.ones(len(reco_energy[key])),label=f'Reco',histtype='step')
    #for key in eidists.keys():
    #    plt.hist(eidists[key],bins=120,range=(0,12),weights=np.ones(len(eidists[key])),label=f'Truth',histtype='step')
    plt.hist(np.log10(np.array(reco_info['energy'])),bins=120,range=(0,12),label='Reco',histtype='step')
    plt.hist(truth_info['logEi'],bins=120,range=(0,12),label='Truth',histtype='step')

    plt.title('Track energy: unweighted')
    plt.xlabel('$\log_{10}$(Energy [GeV])')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    #for key in reco_energy.keys():
    #    y1, x1, _ = plt.hist(np.log10(np.array(reco_energy[key])),bins=120,range=(0,12),weights=weighting_array[key],label=f'Reco',histtype='step')
    #for key in eidists.keys():
    #    y2, x2, _ = plt.hist(eidists[key],bins=120,range=(0,12),weights=weighting_array[key],label=f'Truth',histtype='step')
    y1, x1, _ = plt.hist(np.log10(np.array(reco_info['energy'])),bins=120,range=(0,12),weights=weighting_array,label='Reco',histtype='step')
    y2, x2, _ = plt.hist(truth_info['logEi'],bins=120,range=(0,12),weights=weighting_array,label='Truth',histtype='step')

    print(y1,x1,y2,x2)

    plt.title('Track energy: weighted')
    plt.xlabel('$\log_{10}$(Energy [GeV])')
    plt.ylabel(flux_str)
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    if (min(y1[y1!=0]) < 1e-30) * (min(y2[y2!=0]) < 1e-30):
        maxy = max(max(y1),max(y2))
        plt.ylim(1e-30, maxy*10)
    pdf.savefig()
    plt.clf()

    #for key in eidists.keys():
    #    truth_init_energy = np.array(eidists[key])
    #    y1, x1, _ = plt.hist(truth_init_energy[np.isnan(np.array(reco_energy[key]))],bins=120,range=(0,12),label=f'Truth: untriggered',histtype='step')
    #    y2, x2, _ = plt.hist(truth_init_energy[~np.isnan(np.array(reco_energy[key]))],bins=120,range=(0,12),label=f'Truth: triggered',histtype='step')
    truth_init_energy = np.array(truth_info['logEi'])
    reco_nan_array = np.isnan(np.array(reco_info['energy']))
    y1, x1, _ = plt.hist(truth_init_energy[reco_nan_array],bins=120,range=(0,12),label='Truth: untriggered',histtype='step')
    y2, x2, _ = plt.hist(truth_init_energy[~reco_nan_array],bins=120,range=(0,12),label='Truth: triggered',histtype='step')

    plt.title('Untriggered Events: unweighted')
    plt.xlabel('$\log_{10}$(Energy [GeV])')
    plt.ylabel('#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    if (min(y1[y1!=0]) < 1e-30) * (min(y2[y2!=0]) < 1e-30):
        maxy = max(max(y1),max(y2))
        plt.ylim(1e-30, maxy*10)
    pdf.savefig()
    plt.clf()

    pdf.close()

if __name__=='__main__':
    main()
