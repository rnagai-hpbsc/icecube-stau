import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import click

import os, sys

from modules import formMMCTrackData

from tqdm import tqdm

@click.command()
@click.option('--muon',type=str,default=None)
@click.option('--stau',type=str,default=None)
@click.option('--nmax',type=int,default=None)
@click.option('--out',type=str,default='')
@click.option('--debug',is_flag=True,default=False)
def main(muon,stau,nmax,out,debug):
    if muon is None: 
        print('Set the muon directory with "--muon"')
        sys.exit(1)
    if stau is None:
        print('Set the stau directory with "--stau"')
        sys.exit(1)

    muon_mmc = formMMCTrackData(f'{muon}/*.i3',num=nmax,debug=debug)
    stau_mmc = formMMCTrackData(f'{stau}/*.i3',num=nmax,debug=debug)

    outfilename = f'{out}' if nmax is None else f'{out}_{nmax}'
    with PdfPages(f'plots/compare_{outfilename}.pdf') as pdf:
        for key in muon_mmc:
            muon_data = [data for data in muon_mmc[key]['Muon']]
            stau_data = [data for data in stau_mmc[key]['Stau']]
            stau_weight = [data for data in stau_mmc['evt_weight']['Stau']]
            while len(stau_data) != len(muon_data):
                if len(stau_data) != len(muon_data):
                    corr_muon_data = []
                    for i in tqdm(range(len(stau_data))):
                        zen_stau = stau_mmc['zenith']['Stau'][i]
                        for j in range(len(muon_data)):
                            zen_muon = muon_mmc['zenith']['Muon'][j]
                            if zen_stau == zen_muon:
                                corr_muon_data.append(muon_mmc[key]['Muon'][j])
                                break
                    muon_data = corr_muon_data
                if len(muon_data) != len(stau_data):
                    corr_stau_data =  []
                    corr_stau_weight = []
                    for i in tqdm(range(len(muon_data))):
                        zen_muon = muon_mmc['zenith']['Muon'][i]
                        for j in range(len(stau_data)):
                            zen_stau = stau_mmc['zenith']['Stau'][j]
                            if zen_muon == zen_stau:
                                corr_stau_data.append(stau_mmc[key]['Stau'][j])
                                corr_stau_weight.append(stau_mmc['evt_weight']['Stau'][j])
                                break
                    stau_data = corr_stau_data
                    stau_weight = corr_stau_weight
            print(f'stau size: {len(stau_data)}')
            print(f'muon size: {len(muon_data)}')
            print(f'stau weight size: {len(stau_weight)}')
            print(f'\r {key} size: {len(stau_data)}',end='')
            stau_weight = np.array(stau_weight)/len(stau_data)
            gmax = max(max(muon_data), max(stau_data))
            gmin = min(min(muon_data), min(stau_data))
            gx = np.linspace(gmin, gmax, 100)
            gy = np.linspace(gmin, gmax, 100)
            if key=='Einit':
                gy = np.linspace(gmin*105.66/150e3,gmax*105.66/150e3,100)
            if key=='logEi':
                gy = np.linspace(gmin+np.log10(105.66/150e3),gmax+np.log10(105.66/150e3),100)

            fig = plt.figure(figsize=(7.5,6))
            #plt.scatter(stau_data,muon_data)
            #heatmap, xedges, yedges = np.histogram2d(stau_data, muon_data, bins=50)
            #extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
            #plt.imshow(heatmap,extent=extent)
            if key=='logEloss':
                plt.hist2d(stau_data,muon_data,bins=100, norm=colors.LogNorm(),weights=stau_weight)
                plt.clim(1e-40,1e-22)
            else:
                plt.hist2d(stau_data,muon_data,bins=100, norm=colors.LogNorm(),weights=stau_weight)
            plt.colorbar()
            plt.plot(gx,gy,linestyle=':',color='magenta')
            plt.title(key)
            plt.xlabel('stau')
            plt.ylabel('muon')
            pdf.savefig()
            plt.clf()
            
            plt.hist(stau_data,label='stau',bins=100,range=(gmin,gmax),weights=stau_weight,histtype='step')
            plt.hist(muon_data,label='muon',bins=100,range=(gmin,gmax),weights=stau_weight,histtype='step')
            plt.title(key)
            plt.legend()
            pdf.savefig()
            plt.clf()
        
    

if __name__ == '__main__':
    main()
