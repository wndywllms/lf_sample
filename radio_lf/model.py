import os
from astropy.table import Table
import numpy as np

MODPATH = os.path.dirname(__file__)

##### LF's from the literature #####

def get_novak_lf_model(z=0, scalef=150.):
    dlogLrange = 0.1
    lLrange = np.arange(20, 28.5, dlogLrange)
    Lrange = 10**lLrange
    #lLrange150 = np.log10(Lrange*(150./1400)**-0.7)
    #Lrange150 = 10**lLrange150
    
    # at 150 MHz
    lLrange150 = np.log10(Lrange*(scalef/1400)**-0.7)
    Lrange150 = 10**lLrange150
    
    #lLrange150 = lLrange
    #Lrange150 = Lrange
    
    #local phi    
    def phi0(L, phiS=3.55e-3, LS=1.85e21, alpha=1.22, sigma=0.63):
        return phiS*((L/LS)**(1-alpha))*np.exp(-1.*(np.log10(1+L/LS))**2./(2*sigma**2.))
    
    alphaL = 3.16
    betaL = -0.32
    rho = phi0( Lrange / (1+z)**(alphaL+z*betaL) )
    
    
            
    
    rho = rho*lLrange/lLrange150
    
    return lLrange150, rho



def get_BH(ttype='all', f=150.):
    '''
    f in MHz
    '''
    # load Best & Heckman LF
    BHLF = Table.read(os.path.join(MODPATH,'data/LFs/bestheckmanLF.fits'))
    logPlow = BHLF['Plow']
    logPhigh = BHLF['Phigh']
    logp_BH = (logPlow+logPhigh)/2.

    # scale to f
    if f != 1400.:
        alpha = -0.7
        logp_BH = logp_BH + alpha*np.log10(f/1400)
        

    if ttype == 'all':
        log_rho_BH = BHLF['A_log_rho']
        log_rho_BH_erru = BHLF['A_err_up']
        log_rho_BH_errl = BHLF['A_err_low']
    elif ttype == 'lerg':
        log_rho_BH = BHLF['L_log_rho']
        log_rho_BH_erru = BHLF['L_err_up']
        log_rho_BH_errl = BHLF['L_err_low']
    elif ttype == 'herg':
        log_rho_BH = BHLF['H_log_rho']
        log_rho_BH_erru = BHLF['H_err_up']
        log_rho_BH_errl = BHLF['H_err_low']
    elif ttype == 'agn':
        log_rho_BH = BHLF['AGN_log_rho']
        log_rho_BH_erru = BHLF['AGN_err_up']
        log_rho_BH_errl = BHLF['AGN_err_low']
    elif ttype == 'sf':
        log_rho_BH = BHLF['SF_log_rho']
        log_rho_BH_erru = BHLF['SF_err_up']
        log_rho_BH_errl = BHLF['SF_err_low']
    else:
        raise Exception('ttype not handled')

    x = logp_BH
    y = 10**log_rho_BH
    yerru = 10**(log_rho_BH+log_rho_BH_erru) - 10**(log_rho_BH)
    yerrl = 10**(log_rho_BH) - 10**(log_rho_BH-log_rho_BH_errl)
    
    return x, y, np.array([yerru,yerrl])

def get_MS(ttype='agn', f=150.):
    '''
    f in MHz
    '''
    if ttype == 'agn':
        ff = os.path.join(MODPATH,'data/LFs/ms-agn.txt') 
    elif ttype =='sf':
        ff = os.path.join(MODPATH,'data/LFs/ms-sf.txt') 
        
    elif ttype == 'all':
        xA, yA, yerrA = get_MS(ttype='agn', f=f)
        xS, yS, yerrS = get_MS(ttype='sf', f=f)
        
        x = np.unique(np.hstack((xA,xS)))
        y = np.nan*np.ones_like(x)
        yerr = np.nan*np.ones((2,len(x)))
        for i,xx in enumerate(x):
            iA = np.where(xA == xx)[0]
            iS = np.where(xS == xx)[0]
            # must be both!!
            if len(iA) ==1 and len(iS) == 1:
                y[i] = yA[iA[0]] + yS[iS[0]]
                yerr[0][i] = yerrA[0][iA[0]] + yerrS[0][iS[0]]
                yerr[1][i] = yerrA[1][iA[0]] + yerrS[1][iS[0]]
        return x, y, yerr
                
        
    t = np.genfromtxt(ff, dtype=[('x','f8'),('y','e'),('yerru','e'),('yerrl','e')])
    
    
    
    x = t['x']
    
    
    # scale to f MHz
    if f != 1400.:
        alpha = -0.7
        x = x + alpha*np.log10(f/1400)
    
    y = 2.5*10**t['y']
    yerru = 2.5*(10**(t['y']+t['yerru']) - 10**(t['y']))
    yerrl = 2.5*(10**(t['y']) - 10**(t['y']-t['yerrl']))
    
    return x, y, np.array([yerru,yerrl])


def get_P(ttype='agn', f=150.):
    '''
    f in MHz
    '''
    if ttype == 'agn':
        ff = os.path.join(MODPATH,'data/LFs/prescott-agn.txt')
    elif ttype =='sf':
        ff = os.path.join(MODPATH,'data/LFs/prescott-sf.txt')
    elif ttype == 'all':
        xA, yA, yerrA = get_P(ttype='agn', f=f)
        xS, yS, yerrS = get_P(ttype='sf', f=f)
        
        x = np.unique(np.hstack((xA,xS)))
        y = np.nan*np.ones_like(x)
        yerr = np.nan*np.ones((2,len(x)))
        for i,xx in enumerate(x):
            iA = np.where(xA == xx)[0]
            iS = np.where(xS == xx)[0]
            # must be both!!
            if len(iA) ==1 and len(iS) == 1:
                y[i] = yA[iA[0]] + yS[iS[0]]
                yerr[0][i] = yerrA[0][iA[0]] + yerrS[0][iS[0]]
                yerr[1][i] = yerrA[1][iA[0]] + yerrS[1][iS[0]]
        return x, y, yerr
    
    t = np.genfromtxt(ff, dtype=[('x','f8'),('y','e'),('yerrl','e'),('yerru','e')])
    
    x = t['x']
    
    
    # scale to f MHz
    if f != 1400.:
        alpha = -0.7
        x = x + alpha*np.log10(f/325)
    
    y = 2.5*t['y']
    yerru = 2.5*t['yerru']
    yerrl = 2.5*t['yerrl']
    
    return x, y, np.array([yerru,yerrl])


def get_mjh(ttype='agn', zbin=None):
    if ttype == 'agn':
        f = os.path.join(MODPATH,'data/LFs/lofar-agn.txt')
    elif ttype =='sf':
        if zbin is not None:
            f = os.path.join(MODPATH,'data/LFs/lofar-sf-zbin{i:d}.txt'.format(i=zbin))
        else:
            f = os.path.join(MODPATH,'data/LFs/lofar-sf.txt')
    elif ttype == 'all':
        xA, yA, yerrA = get_mjh(ttype='agn')
        xS, yS, yerrS = get_mjh(ttype='sf')
        
        x = np.unique(np.hstack((xA,xS)))
        y = np.nan*np.ones_like(x)
        yerr = np.nan*np.ones((2,len(x)))
        for i,xx in enumerate(x):
            iA = np.where(xA == xx)[0]
            iS = np.where(xS == xx)[0]
            # must be both!!
            if len(iA) ==1 and len(iS) == 1:
                y[i] = yA[iA[0]] + yS[iS[0]]
                yerr[0][i] = yerrA[0][iA[0]] + yerrS[0][iS[0]]
                yerr[1][i] = yerrA[1][iA[0]] + yerrS[1][iS[0]]
        return x, y, yerr
    
    t = np.genfromtxt(f, dtype=[('x','f8'),('y','e'),('yerr','e')])
    
    
    x = t['x']
    y = 2.5*t['y']
    yerru = 2.5*t['yerr']
    yerrl = 2.5*t['yerr']
    
    return x, y, np.array([yerru,yerrl])


def get_pracy_LF(ttype='all', f=150.):
    '''
    f in MHz
    '''
    # load Pracy
    PLF = Table.read(os.path.join(MODPATH,'data/LFs/pracy.fits'))
    logp_P = PLF['P']

    # scale to f MHz
    if f != 1400.:
        alpha = -0.7
        logp_P = logp_P + alpha*np.log10(f/1400)
        

    if ttype == 'all':
        log_rho_P = PLF['All_rho']
        log_rho_P_erru = PLF['All_rhoerrup']
        log_rho_P_errl = PLF['All_rhoerrlow']
    elif ttype == 'lerg':
        log_rho_P = PLF['LERG_rho']
        log_rho_P_erru = PLF['LERG_rhoerrup']
        log_rho_P_errl = PLF['LERG_rhoerrlow']
    elif ttype == 'herg':
        log_rho_P = PLF['HERG_rho']
        log_rho_P_erru = PLF['HERG_rhoerrup']
        log_rho_P_errl = PLF['HERG_rhoerrlow']
    elif ttype == 'agn':
        log_rho_P = PLF['AGN_rho']
        log_rho_P_erru = PLF['AGN_rhoerrup']
        log_rho_P_errl = PLF['AGN_rhoerrlow']
    elif ttype == 'sf':
        log_rho_P = PLF['SF_rho']
        log_rho_P_erru = PLF['SF_rhoerrup']
        log_rho_P_errl = PLF['SF_rhoerrlow']
    else:
        raise Exception('ttype not handled')

    x = logp_P
    y = 2.5*10**log_rho_P
    yerru = 10**(log_rho_P+log_rho_P_erru) - 10**(log_rho_P)
    yerrl = 10**(log_rho_P) - 10**(log_rho_P-log_rho_P_errl)
    
    return x, y, np.array([yerru,yerrl])


def get_best_lf_model(z=0, model='1a', scalef=150.):
    '''from 2014MNRAS.445..955B'''
    dlogLrange = 0.1
    lLrange = np.arange(22.8, 28.5, dlogLrange)
    Lrange = 10**lLrange
    #lLrange150 = np.log10(Lrange*(150./1400)**-0.7)
    #Lrange150 = 10**lLrange150
    
    # at 150 MHz
    lLrange150 = np.log10(Lrange*(scalef/1400)**-0.7)
    Lrange150 = 10**lLrange150
    
    #lLrange150 = lLrange
    #Lrange150 = Lrange
    
    # z=0 jetmode LF
    lL0 = 24.81
    lrho0 = -5.3
    beta = 0.39
    gamma = 1.61
    
    L0 = 10**lL0
    rho0 = 10**lrho0
    
    fp = 0.
    fL = 0.
    
    if model=='':
        hostevol=False
        lumevol=False
        rhoevol=False
        delay=False
    elif model=='1a':
        #as potential hosts
        hostevol=True
        lumevol=False
        rhoevol=False
        delay=False
    elif model=='1b':
        #as potential hosts with delay
        hostevol=True
        lumevol=False
        rhoevol=False
        delay=True
        tau = 2.0  # Gyr
    elif model=='2a':
        #lum evol
        #as potential hosts
        hostevol=True
        lumevol=True
        rhoevol=False
        delay=False
        delta = 1.6
    elif model=='2b':
        #lum evol
        #as potential hosts with delay
        hostevol=True
        lumevol=True
        rhoevol=False
        delay=True
        tau = 1.5  # Gyr
        delta = 2.8
    elif model=='2c':
        #lum evol
        #as potential hosts with dens evol
        hostevol=False
        lumevol=True
        rhoevol=True
        delay=False
        eta = -1.6
        delta = 2.8
    elif model=='3a':
        #as potential hosts
        hostevol=True
        lumevol=False
        rhoevol=False
        delay=False
        fp = 1.2
        fL = 0.18
        pass
    elif model=='3b':
        #lum evol
        #as potential hosts with delay
        hostevol=True
        lumevol=False
        rhoevol=False
        delay=True
        delta = 1.6
        tau = 1.4  # Gyr
        #delta = 2.8
        fp = 2.0
        fL = 0.14
        pass
    elif model=='3c':
        #lum evol
        #as potential hosts with dens evol
        hostevol=False
        lumevol=False
        rhoevol=True
        delay=False
        eta = -0.9
        #delta = 2.8
        fp = 2.0
        fL = 0.14
        pass
    else:
        print(model + ' not supported')
        return [np.nan], [np.nan]
    
    print(model)
        
        
        
    if lumevol:
        L0 = L0*(1+z)**delta 
    
        
    
    lrat = Lrange/L0
    #lrat = lLrange/lL0
    rho = rho0 / (lrat**beta + lrat**gamma)  # beta is faint end, gama is bright end
    rho1 = rho0 * (lrat**gamma)
    rho2 = rho0 * (lrat**beta)
    
    if delay:

        cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        zt = z_at_value(acosmo.age, acosmo.age(z) - tau*u.Gyr)
        
    else:
        zt = z
    
    print(z, zt)
    
    if rhoevol:
        rho = rho*(1+zt)**eta
        
    if hostevol:
        if zt < 0.8:
            rho = rho*(1+zt)**-0.1
        else:
            rho = rho*((1.8**-0.1)/(1.8**-6.5))*((1+zt)**-6.5)

    def radLF(L, z):
        L0 = 10**26.62
        rho0 = 10**-7.32
        beta = 0.35
        gamma = 1.7
        lrat = L/L0
        rho = rho0 / (lrat**beta + lrat**gamma)
        
        #z0 = 0.  # sevenfold increase between z=0 and z=0.75 gives rho ~ rho0*(1+z)^3.5
        #rhoz = rho * 10.*(z-z0)
        rhoz = rho*(1+z)**3.5
        
        return rhoz
            
    if fp > 0:
        rho = rho + fp*radLF(Lrange/fL, z)
        #rho = radLF(Lrange, z)
    
    rho = rho*lLrange/lLrange150
    #rho = rho  #/Lrange
    #return lLrange, rho
    return lLrange150, rho


def get_best_lf(mode='local-all'):
    '''from 2014MNRAS.445..955B'''
    lLrange = np.arange(23, 29, 0.01)
    Lrange = 10**lLrange
    lLrange150 = np.log10(Lrange*(150./1400)**-0.7)
    Lrange150 = 10**lLrange150
    if mode == 'local-all':
        lL0 = 24.95
        lrho0 = -5.33
        beta = 0.42
        gamma = 1.66
    elif mode == 'local-jet':
        lL0 = 24.81
        lrho0 = -5.3
        beta = 0.39
        gamma = 1.61
    elif mode == 'local-radiative':
        lL0 = 26.62
        lrho0 = -7.32
        beta = 0.35
        gamma = 1.70
    elif mode == '0.5-1-jet':
        lL0 = 25.5
        lrho0 = -5.62
        beta = 0.35
        gamma = 1.7
    elif mode == '0.5-1-radiative':
        lL0 = 26.45
        lrho0 = -6.37
        beta = 0.30
        gamma = 1.70
    
    lrat = 10**lLrange/10**lL0
    lrho = lrho0 - np.log10(lrat**beta+lrat**gamma)
    rho = 10**lrho
    rho = rho*lLrange/lLrange150
    return lLrange150, rho

