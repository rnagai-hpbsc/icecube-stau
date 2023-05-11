import numpy as np

from I3Tray import *
from icecube import icetray, dataclasses, simclasses, dataio, phys_services, sim_services

import glob 

from tqdm import tqdm

def formMMCTrackData(files,alpha=1.5):
    filelist = glob.glob(files)
    if len(filelist)<1:
        print('Invalid filename!')
        raise Exception

    ofilename = filelist[0].split('/')[-1].split('_')[0]

    mmc_dict = {'Einit': {},
                'Eloss': {},
                'logEi': {},
                'logEloss':{},
                'zenith':{},
                'azimuth':{},
                'coszen':{},
                'weight':{},
                'evt_weight':{}
               }

    i = 0
    for fname in tqdm(filelist):
        i += 1

        try: 
            category = fname.split('/')[-1].split('_')[0].split('-')[1]
        except:
            category = 'None'

        energy_info = {}
        for key in mmc_dict:
            mmc_dict[key].setdefault(category,[])
            energy_info.setdefault(key,[])

        with dataio.I3File(fname) as i3f:
            while i3f.more():
                try:
                    frame = i3f.pop_daq()
                except RuntimeError:
                    continue
                mmc = frame.Get('MMCTrackList')
                mmc0 = mmc[0]
                mmc0particle = mmc0.GetI3Particle()
                init_energy = mmc0particle.energy
                eloss = mmc0.GetElost()
                energy_info['Einit'].append(init_energy)
                energy_info['Eloss'].append(eloss)
                energy_info['logEi'].append(np.log10(init_energy))
                energy_info['logEloss'].append(np.log10(eloss))
                energy_info['zenith'].append(mmc0particle.dir.zenith/np.pi*180)
                energy_info['azimuth'].append(mmc0particle.dir.azimuth/np.pi*180)
                energy_info['coszen'].append(np.cos(mmc0particle.dir.zenith))
                try:
                    event_weight = float(f"{frame['EventWeight']}".split('I3String("')[-1].split('")')[0])
                except:
                    event_weight = 1
                energy_info['evt_weight'].append(event_weight)
                energy_info['weight'].append(1)

        for key in mmc_dict:
            mmc_dict[key][category].extend(energy_info[key])

    return mmc_dict

if __name__ == '__main__':
    ex_dict = formMMCTrackData('../data/RunTest5/Stau_ppc/*0001_ppc.i3')
    print(ex_dict)
