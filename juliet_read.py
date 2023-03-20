import javaobj
import numpy as np
from data_input_stream import DataInputStream

def getIntMtx(javafile,show=True):
    if show:
        print(f'Read {javafile}...')
    with open(javafile,"rb") as f:
        sobj = f.read()

    intMtx = javaobj.loads(sobj)
    return intMtx

def getIntermedFile(javadat,show=True):
    if show:
        print(f'Read JULIeT intermediate file: {javadat}...')
    logYList = []
    logEList = []
    logYDsigmaDyList = []
    DsDyList = []
    with open(javadat,'rb') as f:
        IntTable = DataInputStream(f)
        for jLogY in range(35):
            logYList.append(IntTable.read_double())
        for iLogE in range(140):
            rawE = IntTable.read_double()
            logEList.append(rawE)
            yDsDyList_tmp = []
            DsDyList_tmp = []
            for jLogY in range(35):
                raw = IntTable.read_double()
                yDsDyList_tmp.append(raw)
                DsDyList_tmp.append(10**raw * 10**rawE * 1e-32 / 10**logYList[jLogY])
            logYDsigmaDyList.append(yDsDyList_tmp)
            DsDyList.append(DsDyList_tmp)
    logYArray = np.array(logYList)
    logEArray = np.array(logEList)
    logYDsigmaDyArray = np.array(logYDsigmaDyList)
    DsDyArray = np.array(DsDyList)
    return logYArray, logEArray, logYDsigmaDyArray, DsDyArray

def getPropMtx(javadat,show=True): 
    if show:
        print(f'Read JULIeT Propagation Matrix: {javadat}...')
    propMtx = {}
    key_1 = ['nue','numu','nutau','mu','tau','stau'] #input particle
    key_2 = ['nue','numu','nutau','e','mu','tau','hadron'] #out particle

    for pin in key_1:
        for pout in key_2:
            propMtx[f'{pin}2{pout}'] = np.zeros((700,700))
    propMtx['stau2stau'] = np.zeros((700,700))

    with open(javadat,'rb') as f:
        jPropMtx = DataInputStream(f)
        for iLogE in range(700):
            for jLogE in range(iLogE+1):
                for pin in key_1:
                    for pout in key_2:
                        propMtx[f'{pin}2{pout}'][iLogE][jLogE] = jPropMtx.read_double()
                propMtx['stau2stau'][iLogE][jLogE] = jPropMtx.read_double()
    return propMtx 

def getLaTeXStr(number):
    n = f'{number:.3e}'.split('e')
    if len(n)>1:
        num = f'${n[0]}\\times' +  '10^{' + f'{int(n[1])}' + '}$'
    else:
        num = n[0]
    return num
    
