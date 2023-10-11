import os, sys
from I3Tray import *
from icecube import icetray, dataclasses, simclasses, dataio, phys_services, sim_services

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

import glob
from natsort import natsorted
import click
from tqdm import tqdm

import copy

flux_str = 'Flux [$\mathrm{cm^{-2}\ s^{-1}\ sr^{-2}\ GeV^{-1}}$]'
truth_info = {}
reco_info = {}

FILENAME = 'reco2_plots'

configs = {'Energy'   : [None, "Energy [GeV]"],
           'Eloss'    : [None, "Energy Loss [GeV]"],
           'logEi'    : [[120,(0,12)],"$\log_{10}$ ($E$ [GeV])"],
           'logEloss' : [[120,(0,12)],"$\log_{10}$ ($E_\mathrm{loss}$ [GeV])"],
           'Zenith'   : [[50,(0,180)],"Zenith [deg]"],
           'cosZen'   : [[50,(-1,1)],"$\cos\\theta_\mathrm{zen}$"],
           'Azimuth'  : [[50,(0,360)],"Azimuth [deg]"],
           'Eventweight' : [None,"Event Weight"],
           'Weight'   : [None, ""]
           }

#ranges = {'Energy':None, 'Eloss':None, 'logEi':[120,(0,12)], 'logEloss':[120,(0,12)], 'Zenith':[50,(0,180)], 'cosZen':[50,(-1,1)], 'Azimuth':[50,(0,360)], 'Eventweight':None, 'Weight':None}
#labels = {'Energy':"Energy [GeV]", 'Eloss':"Energy Loss [GeV]", 'logEi':"$\log_{10} E$ [GeV]", 'logEloss':"$\log_{10} E_\mathrm{loss}$ [GeV]", 
#          'Zenith':"Zenith [deg]", 'cosZen':"$\cos\\theta_\mathrm{zen}$", 'Azimuth':"Azimuth [deg]", 'Eventweight':"Event Weight", 'Weight':[]}

def plot_truth_reco(pdf, key, truth, reco, weight=False):
    h_ranges = configs[key][0]
    xlabel = configs[key][1]
    if h_ranges is None:
        return 
    print(key,len(truth[key]),len(reco[key]))
    weights = np.array(truth['Weight'])/len(truth['Energy'])*np.array(truth['Eventweight'])*1e6 if weight else np.ones(len(truth[key]))
    
    if len(reco[key])!=0:
        y1, x1, _ = plt.hist(reco[key],bins=h_ranges[0],range=h_ranges[1],weights=weights,label='Reco',histtype='step')
    else:
        y1 = np.ones(10)
    y2, x2, _ = plt.hist(truth[key],bins=h_ranges[0],range=h_ranges[1],weights=weights,label='Truth',histtype='step')

    plt.title(f'{key}: weighted' if weight else f'{key}: unweighted')
    plt.xlabel(xlabel)
    plt.ylabel(flux_str if weight else '#Events')
    plt.legend()
    pdf.savefig()
    plt.yscale('log')
    if (min(y1[y1!=0]) < 1e-30) + (min(y2[y2!=0]) < 1e-30):
        maxy = max(max(y1),max(y2))
        plt.ylim(1e-30, maxy*10)
    pdf.savefig()
    plt.clf()

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
        filelist = natsorted(glob.glob(f'{datadir}/*.i3'))

    ET = {}

    for key in configs.keys():
        truth_info.setdefault(key,[])
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

        truth_energy_info = {}
        reco_energy_info = {}

        for key in configs.keys():
            truth_energy_info.setdefault(key,[])
            reco_energy_info.setdefault(key,[])

        while i3f.more():
            try:
                frame = i3f.pop_physics()
            except:
                break

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
            truth_energy_info['Energy'].append(init_energy)
            truth_energy_info['Eloss'].append(elost)
            truth_energy_info['logEi'].append(np.log10(init_energy))
            truth_energy_info['logEloss'].append(np.log10(elost))
            truth_energy_info['Zenith'].append(mmc0particle.dir.zenith/np.pi*180)
            truth_energy_info['Azimuth'].append(mmc0particle.dir.azimuth/np.pi*180)
            truth_energy_info['cosZen'].append(np.cos(mmc0particle.dir.zenith))
            try: 
                event_weight = float(f"{frame['EventWeight']}".split('I3String("')[-1].split('")')[0])
            except: 
                event_weight = 1
            truth_energy_info['Eventweight'].append(event_weight)

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
                reco_energy_info['Energy'].append(sp.energy)
                reco_energy_info['logEi'].append(np.log10(sp.energy))
                reco_energy_info['Zenith'].append(sp.dir.zenith/np.pi*180)
                reco_energy_info['Azimuth'].append(sp.dir.azimuth/np.pi*180)
                reco_energy_info['cosZen'].append(np.cos(sp.dir.zenith))
            else:
                reco_energy_info['Energy'].append(np.nan)
                reco_energy_info['logEi'].append(np.nan)
                reco_energy_info['Zenith'].append(np.nan)
                reco_energy_info['Azimuth'].append(np.nan)
                reco_energy_info['cosZen'].append(np.nan)

        i3f.close()

        for key in configs.keys():
            truth_info[key].extend(truth_energy_info[key])
            reco_info[key].extend(reco_energy_info[key])

        thisweight = (-1)*np.array(truth_energy_info['Energy'])**alpha*(10**(hrange*(1-alpha))-10**(lrange*(1-alpha)))/(1-alpha)
        truth_info['Weight'].extend(thisweight)

        if i==endn:
            break

    ienergies = [float(key) for key in ET]
    tmean = [np.mean(np.array(ET[key])) for key in ET]
    tsigma = [np.std(np.array(ET[key])) for key in ET]
    
    if out is not None:
        suffix = f'_{out}'
    else:
        suffix = ''

    prefix = 'plots'
    os.makedirs(prefix,exist_ok=True)
    if endn is None:
        pdf = PdfPages(f'{prefix}/{FILENAME}{suffix}_all.pdf')
    else:
        pdf = PdfPages(f'{prefix}/{FILENAME}{suffix}_{endn}.pdf')

    plt.errorbar(ienergies, tmean, yerr=tsigma,capsize=5,fmt='o')
    plt.title('Mean time for a single event generation')
    plt.xlabel('$\log_{10}(E\ [\mathrm{GeV}])$')
    plt.ylabel('Time [s]')
    pdf.savefig()
    plt.yscale('log')
    pdf.savefig()
    plt.clf()

    for key in configs.keys():
        plot_truth_reco(pdf, key, truth_info, reco_info)
        plot_truth_reco(pdf, key, truth_info, reco_info, weight=True)

    truth_init_energy = np.array(truth_info['logEi'])
    reco_nan_array = np.isnan(np.array(reco_info['Energy']))
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

    ## 2d hist
    cmap = copy.copy(plt.cm.jet)
    #cmap.set_under('w',1) # less than 1 will be filled in white

    import matplotlib
    h, _, _, im = plt.hist2d(reco_info['cosZen'],reco_info['logEi'],bins=(50,60),range=(configs['cosZen'][0][1],configs['logEi'][0][1]),weights=truth_info['Eventweight'],cmap=cmap,norm=matplotlib.colors.LogNorm())
    plt.xlabel(configs['cosZen'][1])
    plt.ylabel(configs['logEi'][1])
    plt.colorbar(im,label=flux_str)
    if (min(h[h!=0]) < 1e-30):
        maxy = h.max()
        plt.clim(1e-30, maxy*10)
    pdf.savefig()
    plt.clf()

    pdf.close()

if __name__=='__main__':
    main()
