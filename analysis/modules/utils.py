
def getLaTeX(inchar):
    outchar = ""
    if inchar == 'nue':
        outchar = '\\nu_{e}'
    elif inchar == 'numu':
        outchar = '\\nu_{\mu}'
    elif inchar == 'nutau':
        outchar = '\\nu_{\\tau}'
    elif inchar == 'mu':
        outchar = '\mu'
    elif inchar == 'tau':
        outchar = '\\tau'
    elif inchar == 'hadron':
        outchar = '\mathrm{had}'
    elif inchar == 'stau':
        outchar = '\\tilde{\\tau}'
    else:
        outchar = inchar
    return outchar 
