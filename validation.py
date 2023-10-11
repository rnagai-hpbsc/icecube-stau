import numpy as np
from juliet_read import getPropMtx
from tqdm import tqdm
import random

def valid_event(logE, zenith, logE_array, normflux):
    valid_min = logE_array >  logE - 0.1
    valid_max = logE_array <= logE
    energy_bin = logE_array[valid_min*valid_max]
    value = normflux[logE_array==energy_bin]
    validator = random.random()
    if validator <= value[0]:
        return True
    else:
        return False

def valid_event_energy_slope(logE, logEnergyMin, factor=1):
    validator = random.random()
    if validator <= 10**(factor*(logEnergyMin-logE)):
        return True
    else:
        return False

def getTheWeight(logE, zenith, logE_array, fluxes):
    index = (logE_array > logE - 0.1) * (logE_array <= logE)
    energy_bin = logE_array[index]
    flux = fluxes[int(zenith)]
    value = flux[logE_array==energy_bin]
    return value

def makeicflux1d(zenith,
                 inkey='stau',
                 outkey='stau',
                 prodflux_path='inputs/mass150/prodflux/totalhist.npy',
                 propMtx_path=None
                 ):

    if propMtx_path is None:
        inputMatrixFilename = f'inputs/mass150/propMtx/stau150GeV_inf_{int(zenith)}deg'
    else:
        inputMatrixFilename = propMtx_path
    propMtx = getPropMtx(inputMatrixFilename,show=False)
    #earthfluxdata = np.load('inputs/mass150/prodflux/totalhist.npy')
    earthfluxdata = np.load(prodflux_path)

    E = [10**(5+0.01*i) for i in range(700)]
    E10 = [10**(5+0.1*i) for i in range(70)]
    logE10 = [5+0.1*i for i in range(70)]
    totalflux = np.zeros(700)
    totalflux10 = np.zeros(70)
    for iloge in range(70):
        energy = iloge*10
        flux = earthfluxdata[50+iloge]
        totalflux += propMtx[f'{inkey}2{outkey}'][energy]*flux
        for i in range(10):
            for j, element in enumerate(propMtx[f'{inkey}2{outkey}'][iloge*10+i]*flux):
                totalflux10[int(j/10)] += element/10
    totalnumber = 0
    for i in range(69):
        totalnumber += (E10[i+1]-E10[i])*(totalflux10[i+1]+totalflux10[i])/2 
    totalnumber = totalnumber * (60*60*24*365*10) * 2*np.pi * (1000*100)**2
    #print(totalnumber)
    return np.array(logE10), np.array(totalflux10)

def normicflux1d(logE_min, logE_max, logE_array, flux_array):
    valid_min = logE_array >= logE_min
    valid_max = logE_array <  logE_max
    valid = valid_min * valid_max
    logEcut = logE_array[valid]
    flux_energycut = flux_array[valid]
    norm_flux_energycut = flux_energycut / np.max(flux_energycut)
    return logEcut, norm_flux_energycut

def makeicflux2d():
    fluxes = []
    for i in tqdm(range(90)):
        logE, flux = makeicflux1d(i)
        fluxes.append(flux)
    return np.array(logE), np.array(fluxes)

def normicflux2d(logE_min, logE_max, logE_array, fluxes, debug=False):
    valid = (logE_array >= logE_min) * (logE_array < logE_max)
    logEcut = logE_array[valid]
    maxfluxvalue = 0
    fluxes_Ecut = []
    for flux_ in tqdm(fluxes):
        flux_energycut = flux_[valid]
        thismax = np.max(flux_energycut)
        if thismax > maxfluxvalue:
            maxfluxvalue = thismax
        fluxes_Ecut.append(flux_energycut)
    norm_fluxes_Ecut = [np.array(flux_) / maxfluxvalue for flux_ in fluxes_Ecut]
    if debug:
        print(f'The value: {np.max(norm_fluxes_Ecut)} should be 1.')
    return logEcut, norm_fluxes_Ecut

def getWeight2d(logE_min, logE_max, logE_array, fluxes):
    valid = (logE_array >= logE_min) * (logE_array < logE_max)
    logEcut = logE_array[valid]
    maxfluxvalue = 0
    for flux_ in tqdm(fluxes):
        flux_energycut = flux_[valid]
        thismax = np.max(flux_energycut)
        if thismax > maxfluxvalue:
            maxfluxvalue = thismax
    return maxfluxvalue

if __name__=='__main__':
    import sys
    import matplotlib.pyplot as plt
    
    nevents = 1000000
    logE_min = 5
    logE_max = 8

    print(sys.argv)

    if len(sys.argv) > 1:
        if sys.argv[1] == '1':
            nlogE, nflux = normicflux1d(logE_min, logE_max, *makeicflux1d(89))
            plt.plot(10**nlogE, nflux/np.sum(nflux)*nevents)

            energies = []
            for i in tqdm(range(nevents)):
                breakflag = False
                while not breakflag:
                    logE = random.uniform(logE_min, logE_max)
                    if valid_event(logE,89,nlogE,nflux):
                        energies.append(logE)
                        breakflag = True

            hist, bins = np.histogram(energies, bins=int((logE_max-logE_min)*10), range=(logE_min, logE_max))
            x = np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])

            plt.errorbar(10**x, hist, yerr=np.sqrt(hist), capsize=2, fmt='o', ecolor='black', color='black', markerfacecolor='black')

            plt.xscale('log')
            plt.yscale('log')
            plt.savefig('testflux.pdf')
        elif sys.argv[1] == '2':
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages('testflux2d.pdf')

            logEs, fluxes = makeicflux2d()
            nlogE, nfluxes = normicflux2d(logE_min, logE_max, logEs, fluxes)

            energies = [[] for i in range(90)]
            for i in tqdm(range(nevents)):
                breakflag = False
                while not breakflag:
                    zenith = random.random()*90
                    logE = random.uniform(logE_min, logE_max)
                    if valid_event(logE,zenith,nlogE,nfluxes[int(zenith)]):
                        energies[int(zenith)].append(logE)
                        breakflag = True

            for i in tqdm(range(90)):
                plt.plot(10**nlogE, nfluxes[i]/np.sum(nfluxes)*nevents,label='Flux')
                hist, bins = np.histogram(energies[i], bins=int((logE_max-logE_min)*10), range=(logE_min, logE_max))
                x = np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])
                plt.errorbar(10**x, hist, yerr=np.sqrt(hist), capsize=2, fmt='o', ecolor='black', color='black', markerfacecolor='black',label='MC')
                plt.title(f'Zenith:{i}deg')
                plt.xscale('log')
                plt.yscale('log')
                plt.legend()
                pdf.savefig()
                plt.close()
            pdf.close()

