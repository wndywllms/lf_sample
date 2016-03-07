import scipy as sp
import pylab as pl
import numpy as np
import cosmolopy as cosmo
import os

# LRT+ default
#If no file is provided a cosmology of (Omega_Matter, Omega_Lambda, Omega_k, H0) = (0.3, 0.7, 0.0, 70.0) is assumed.
default_cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.70}

def RadioPower(Flux, z, alpha=0.8):
    """RadioPower(Flux, z, alpha=0.8)
args
    Flux - in Jy
    z - redshift 
kwargs
    alpha  - spectral index (default = 0.8)
    """
    # flux in Jy to Power in W/Hz
    DL = cosmo.distance.luminosity_distance(z, **default_cosmo)
    power = Flux*(4.*np.pi*DL**2.)*1.e-26*(3.08667758e22**2.)*(1. + z)**(-1.+alpha)
    return power


def RadioFlux (power, z, alpha=0.8): 
    # function to calculate flux density given in Jy some radio power
    DL = cosmo.distance.luminosity_distance(z, **default_cosmo)
    S = 1e26* power / ( 4.*np.pi*DL**2*(3.08667758e22**2.)*(1. + z)**(-1.+alpha) )
    return S

def OpticalLuminosity(Flux, z):
    # flux in W/Hz/m^2 to luminosity in W/Hz
    DL = cosmo.distance.luminosity_distance(z, **default_cosmo)
    luminosity = Flux*(4.*np.pi*DL**2.)*(3.08667758e22**2.)*(1. + z)
    return luminosity

def OpticalLuminosity2(Flux, z, alpha):
    # flux in W/Hz/m^2 to luminosity in W/Hz
    DL = cosmo.distance.luminosity_distance(z, **default_cosmo)
    luminosity = Flux*(4.*np.pi*DL**2.)*(3.08667758e22**2.)*(1. + z)**(-1.+alpha)
    return luminosity

def OpticalFlux (luminosity, z): 
    # function to calculate flux density given some optical luminosity
    DL = cosmo.distance.luminosity_distance(z, **default_cosmo)
    S = luminosity / ( 4.*np.pi*DL**2*(3.08667758e22**2.)*(1. + z) )
    return S

def OpticalMag(mag, z):
    dm = cosmo.magnitudes.distance_modulus(z, **default_cosmo)
    Mag = mag - dm
    return Mag

def XrayLuminosity(Flux, z):
    # flux in W/m^2 to luminosity in W
    DL = cosmo.distance.luminosity_distance(z, **default_cosmo)
    luminosity = Flux*(4.*np.pi*DL**2.)*(3.08667758e22**2.)  #*(1. + z) these are X-ray bolometric
    luminosity = luminosity*1e7  # W -> erg/s
    return luminosity


def match_indices(x1,x2):
    xind1 = []
    xind2 = []
    for i1,x1i in enumerate(x1):
        i2 = np.where(x2==x1i)[0]
        if len(i2) > 0:
            xind1.append(i1)
            xind2.append(i2[0])
    return xind1,xind2
    
def zlim_func(x,m,z,mlim):
    '''
    function whose zero is the maximum z at which target with magnitude m can be observed given the magnitude limit mlim
    f(x) = mlim - DM(x) - m + DM(z)
    '''
    f = mlim -m + cosmo.magnitudes.distance_modulus(z, **default_cosmo) - cosmo.magnitudes.distance_modulus(x, **default_cosmo)
    #print mlim, m, cosmo.magnitudes.distance_modulus(z, **default_cosmo), cosmo.magnitudes.distance_modulus(x, **default_cosmo), x
    return f
    
def vmax(m,z,mlim,area):
    '''
    input - m : observed magnitude of source
            z : redshift of source
         mlim : magnitude limit of survey
         area : in sq degrees
   output - vmax : maximum volume in which source could be observed
    '''
    # find zero - where M = Mlim
    f0 = zlim_func(0.,m,z,mlim)
    f10 = zlim_func(10.,m,z,mlim)
    if f0*f10 < 0:
        try:
            zlim = sp.optimize.brentq(zlim_func, 0., 10., args=(m,z,mlim))
        except RuntimeError:
            print 'solve not converging %.3f, %.3f, %.3f' %(m,z,mlim)
            zlim = np.nan
            return np.nan
    else:
        zlim = np.nan
        return np.nan
    
    ## for checking
    #M = m - cosmo.magnitudes.distance_modulus(z, **default_cosmo)
    #Mlim = mlim - cosmo.magnitudes.distance_modulus(zlim, **default_cosmo)
    #print 'M(z=%.3f) = %.2f' %(z,M)
    #print 'Mlim(zlim=%.3f) = %.2f' %(zlim,Mlim)
    
    degtost = 4.*180.**2./np.pi  # sq degrees to steradians
    domega = area/degtost
    DL = cosmo.distance.luminosity_distance(zlim, **default_cosmo)
    vc = (4.*np.pi/3.)*(DL/(1+zlim))**3.  # volume at zlim
    vmax = domega*vc
    return vmax

def get_zmax(z, L, fluxlim2, stype='Radio',filename='zmax.sav.npy', clobber=False):
    if os.path.isfile(filename) and (not clobber):
        print 'read zmax from '+filename
        zmax = np.load(filename)
    else:
        if os.path.isfile(filename): os.system('rm -rf '+filename)
        zmax = np.zeros(0)
    Nsrc2 = len(z)
    if len(zmax) != len(z):
        print 'calculating zmax for '+filename
        zmax   = np.zeros(Nsrc2)
        for i in range( 0, Nsrc2):
            zt = np.arange(5,-0.001,-0.01) + z[i]
            if stype == 'Radio':
                ft = RadioFlux(L[i], zt, alpha=0.8)
            elif stype == 'Optical':
                ft = OpticalFlux(L[i], zt)
            zm = np.interp(fluxlim2, ft, zt)
            zmax[i] = zm
        np.save(filename, (zmax))
    return zmax

def get_zmin(z, L, fluxlim2, stype='Radio',filename='zmin.sav.npy', clobber=False):
    if os.path.isfile(filename) and (not clobber):
        print 'read zmin from '+filename
        zmin = np.load(filename)
    else:
        if os.path.isfile(filename): os.system('rm -rf '+filename)
        zmin = np.zeros(0)
    Nsrc2 = len(z)
    if len(zmin) != len(z):
        print 'calculating zmin for '+filename
        zmin   = np.zeros(Nsrc2)
        for i in range( 0, Nsrc2):
            zt = np.arange(z[i],-0.001,-0.001)
            if stype == 'Radio':
                ft = RadioFlux(L[i], zt, alpha=0.8)
            elif stype == 'Optical':
                ft = OpticalFlux(L[i], zt)
            zm = np.interp(fluxlim2, ft, zt)
            zmin[i] = zm
        np.save(filename, (zmin))
    return zmin

def get_LF(pbins, power, zmin, zmax, area, ind=None, verbose=True):
    '''
    pbins - bins of power (giving bin boundaries)
    power - log10(power) values
    zmin - for each source minimum z it could be observed (taking into account all selections)
    zmax - for each source maximum z it could be observed
    area - of survey in steradians
    ind - optically specify a selection index
    
    The differential comoving volume element dV_c/dz/dSolidAngle.
    Dimensions are volume per unit redshift per unit solid angle.
    Units are Mpc**3 Steradians^-1.
    '''
    print "Calculating LF: {n} sources".format(n=len(power))
    if ind is not None:
        power = power[ind]
        zmin = zmin[ind]
        zmax = zmax[ind]
        print "sub-selected {n} sources".format(n=len(power))
    
    # number of bins
    N_P_bins = len(pbins) - 1
    
    # bin widths
    dp = pbins[1:] - pbins[:-1]
    
    rho  = np.nan*np.ones((N_P_bins))
    rhoerr  = np.nan*np.ones((N_P_bins))
    num   = np.nan*np.ones((N_P_bins))
    for P_i in range( N_P_bins ): #loop over log P
        # count  #
        sel_ind = np.where( (power > pbins[P_i])  & (power <= pbins[P_i + 1]) )[0]
        #powerhigh = pbins[P_i + 1]
        
        if len(sel_ind) > 0:
            zmax_h = zmax[sel_ind]
            zmin_h = zmin[sel_ind]
            
            vzmax = cosmo.distance.comoving_volume(zmax_h, **default_cosmo)
            vzmin = cosmo.distance.comoving_volume(zmin_h, **default_cosmo)
            
            vi = (vzmax - vzmin)*(area/(4.*np.pi))
            rho[P_i]  =  np.sum(1./vi)      #sum of 1/v_max for each source
            #rhoerr[P_i] =  rho[P_i]*np.sqrt(float(len(sel_ind)))/float(len(sel_ind))
            rhoerr[P_i] = np.sqrt(np.sum((1./vi)**2.))
            num[P_i]  =  float(len(sel_ind))
        if verbose:
            print "{p1:7.2f} < P <= {p2:6.2f} ({n:.0f}) : {rho:6.2e} +/- {rhoerr:6.2e}".format(p1=pbins[P_i], p2=pbins[P_i + 1], rho=rho[P_i], rhoerr=rhoerr[P_i], n=num[P_i])
            
    # per P bin
    rho = rho / dp
    rhoerr = rhoerr / dp
            
    return rho, rhoerr, num


def get_rho_z(zbins, pbins, power, zmin, zmax, area, ind=None, verbose=True):
    '''
    pbins - bin of power (giving bin boundaries)
    power - log10(power) values
    zmin - for each source minimum z it could be observed (taking into account all selections)
    zmax - for each source maximum z it could be observed
    area - of survey in steradians
    ind - optically specify a selection index
    
    The differential comoving volume element dV_c/dz/dSolidAngle.
    Dimensions are volume per unit redshift per unit solid angle.
    Units are Mpc**3 Steradians^-1.
    '''
    print "Calculating LF: {n} sources".format(n=len(power))
    if ind is not None:
        power = power[ind]
        zmin = zmin[ind]
        zmax = zmax[ind]
        print "sub-selected {n} sources".format(n=len(power))
    
    # number of bins
    #N_P_bins = len(pbins) - 1
    N_z_bins = len(zbins) - 1
    
    # bin widths
    dp = pbins[1:] - pbins[:-1]
    
    rho  = np.nan*np.ones((N_z_bins))
    rhoerr  = np.nan*np.ones((N_z_bins))
    num   = np.nan*np.ones((N_z_bins))
    for P_i in range( N_z_bins ): #loop over log P
        # count  #
        sel_ind = np.where( (power > pbins[P_i])  & (power <= pbins[P_i + 1]) )[0]
        #powerhigh = pbins[P_i + 1]
        
        if len(sel_ind) > 0:
            zmax_h = zmax[sel_ind]
            zmin_h = zmin[sel_ind]
            
            vzmax = cosmo.distance.comoving_volume(zmax_h, **default_cosmo)
            vzmin = cosmo.distance.comoving_volume(zmin_h, **default_cosmo)
            
            vi = (vzmax - vzmin)*(area/(4.*np.pi))
            rho[P_i]  =  np.sum(1./vi)      #sum of 1/v_max for each source
            #rhoerr[P_i] =  rho[P_i]*np.sqrt(float(len(sel_ind)))/float(len(sel_ind))
            rhoerr[P_i] = np.sqrt(np.sum((1./vi)**2.))
            num[P_i]  =  float(len(sel_ind))
        if verbose:
            print "{p1:7.2f} < P <= {p2:6.2f} ({n:.0f}) : {rho:6.2e} +/- {rhoerr:6.2e}".format(p1=pbins[P_i], p2=pbins[P_i + 1], rho=rho[P_i], rhoerr=rhoerr[P_i], n=num[P_i])
            
    # per P bin
    rho = rho / dp
    rhoerr = rhoerr / dp
            
    return rho, rhoerr, num


def get_CLF(pbins, power, zmin, zmax, area, ind=None, verbose=True):
    '''
    pbins - bins of power (giving bin boundaries)
    power - log10(power) values
    zmin - for each source minimum z it could be observed (taking into account all selections)
    zmax - for each source maximum z it could be observed
    area - of survey in steradians
    ind - optically specify a selection index
    
    The differential comoving volume element dV_c/dz/dSolidAngle.
    Dimensions are volume per unit redshift per unit solid angle.
    Units are Mpc**3 Steradians^-1.
    '''
    print "Calculating CLF: {n} sources".format(n=len(power))
    if ind is not None:
        power = power[ind]
        zmin = zmin[ind]
        zmax = zmax[ind]
        print "sub-selected {n} sources".format(n=len(power))
    
    # number of bins
    N_P_bins = len(pbins) - 1
    
    # bin widths
    dp = pbins[1:] - pbins[:-1]
    
    rho  = np.nan*np.ones((N_P_bins))
    rhoerr  = np.nan*np.ones((N_P_bins))
    num   = np.nan*np.ones((N_P_bins))
    for P_i in range( N_P_bins ): #loop over log P
        # count  #
        sel_ind = np.where( (power > pbins[P_i])  )[0]
        #powerhigh = pbins[P_i + 1]
        
        if len(sel_ind) > 0:
            zmax_h = zmax[sel_ind]
            zmin_h = zmin[sel_ind]
            
            vzmax = cosmo.distance.comoving_volume(zmax_h, **default_cosmo)
            vzmin = cosmo.distance.comoving_volume(zmin_h, **default_cosmo)
            
            vi = (vzmax - vzmin)*(area/(4.*np.pi))
            rho[P_i]  =  np.sum(1./vi)      #sum of 1/v_max for each source
            #rhoerr[P_i] =  rho[P_i]*np.sqrt(float(len(sel_ind)))/float(len(sel_ind))
            rhoerr[P_i] = np.sqrt(np.sum((1./vi)**2.))
            num[P_i]  =  float(len(sel_ind))
        if verbose:
            print "{p1:7.2f} < P ({n:.0f}) : {rho:6.2e} +/- {rhoerr:6.2e}".format(p1=pbins[P_i], rho=rho[P_i], rhoerr=rhoerr[P_i], n=num[P_i])
            
    # per P bin
    rho = rho 
    rhoerr = rhoerr 
            
    return rho, rhoerr, num

def get_LF_f_areal(pbins_in, power, zmin, zmax, fcor, areal, area, ind=None, verbose=True, xstr="P", ignoreMinPower=False):
    '''
    pbins - bins of power (giving bin boundaries)
    power - log10(power) values
    zmin - for each source minimum z it could be observed (taking into account all selections)
    zmax - for each source maximum z it could be observed
    fcor - correction factor for completeness
    areal - the fractional areal coverage per source
    area - of survey in steradians
    ind - optically specify a selection index
    
    The differential comoving volume element dV_c/dz/dSolidAngle.
    Dimensions are volume per unit redshift per unit solid angle.
    Units are Mpc**3 Steradians^-1.
    '''
    pbins = pbins_in.copy()  # copy the pbins cos we're going to mess with them #
    #print pbins
    #print type(pbins)
    #print pbins_in
    #print type(pbins_in)
    print "Calculating {s}F (f_areal): {n} sources".format(n=len(power),s=xstr)
    if ind is not None:
        power = power[ind]
        zmin = zmin[ind]
        zmax = zmax[ind]
        fcor = fcor[ind]
        areal = areal[ind]
        print "sub-selected {n} sources".format(n=len(power))
    
    minpower = np.min(power)
    
    # number of bins
    N_P_bins = len(pbins) - 1
    
    # bin widths
    dpnom = pbins[1:] - pbins[:-1]
    
    dp  = np.ones((N_P_bins)) * dpnom
    rho  = np.nan*np.ones((N_P_bins))
    rhoerr  = np.nan*np.ones((N_P_bins))
    num   = np.nan*np.ones((N_P_bins))
    for P_i in range( N_P_bins ): #loop over log P
        
    #for P_i in [5]: #loop over log P
        # count  #
        sel_ind = np.where( (power >= pbins[P_i])  & (power < pbins[P_i + 1]) )[0]
        if len(sel_ind) > 0:
        #powerhigh = pbins[P_i + 1]
        
        
            # discard the bin if the minimum value lies in this bin  - i.e. we are certainly not complete here #
            if not ignoreMinPower:
                if minpower > pbins[P_i]:
                    continue
                
            zmax_h = zmax[sel_ind]
            zmin_h = zmin[sel_ind]
            areal_h = areal[sel_ind]
            fcor_h = fcor[sel_ind]
            
            vzmax = cosmo.distance.comoving_volume(zmax_h, **default_cosmo) ## Mpc^3
            vzmin = cosmo.distance.comoving_volume(zmin_h, **default_cosmo)
            
            vi = (vzmax - vzmin)*(area/(4.*np.pi))  #area/4pi gives per sterad
            rho[P_i]  =  np.sum(fcor_h/(areal_h*vi))      #sum of 1/v_max for each source
            #rhoerr[P_i] =  rho[P_i]*np.sqrt(float(len(sel_ind)))/float(len(sel_ind))
            rhoerr[P_i] = np.sqrt(np.sum((fcor_h/(areal_h*vi))**2.))
            num[P_i]  =  float(len(sel_ind))
            #print vzmax.min(), vzmax.max()#, vzmax
            #print vzmin.min(), vzmin.max()#, vzmin
            #print vi.min(), vi.max()#, vzmin
            #print np.sum(1./vi)#, vzmin
        if verbose:
            print "{p1:7.2f} < {x} <= {p2:6.2f} ({n:.0f}) : {rho:6.2e} +/- {rhoerr:6.2e}".format(x=xstr, p1=pbins[P_i], p2=pbins[P_i + 1], rho=rho[P_i]/dp[P_i], rhoerr=rhoerr[P_i]/dp[P_i], n=num[P_i])
            
    # per P bin
    rho = rho / dp
    rhoerr = rhoerr / dp
    
            
    return rho, rhoerr, num


def get_CLF_f_areal(pbins, power, zmin, zmax, fcor, areal, area, ind=None, verbose=True, xstr="P"):
    '''
    pbins - bins of power (giving bin boundaries)
    power - log10(power) values
    zmin - for each source minimum z it could be observed (taking into account all selections)
    zmax - for each source maximum z it could be observed
    fcor - correction factor for completeness
    areal - the fractional areal coverage per source
    area - of survey in steradians
    ind - optically specify a selection index
    
    The differential comoving volume element dV_c/dz/dSolidAngle.
    Dimensions are volume per unit redshift per unit solid angle.
    Units are Mpc**3 Steradians^-1.
    '''
    print "Calculating CLF (f_areal): {n} sources".format(n=len(power))
    if ind is not None:
        power = power[ind]
        zmin = zmin[ind]
        zmax = zmax[ind]
        fcor = fcor[ind]
        areal = areal[ind]
        print "sub-selected {n} sources".format(n=len(power))
    
    minpower = np.min(power)
    
    # number of bins
    N_P_bins = len(pbins) - 1
    
    # bin widths
    dp = pbins[1:] - pbins[:-1]
    
    rho  = np.nan*np.ones((N_P_bins))
    rhoerr  = np.nan*np.ones((N_P_bins))
    num   = np.nan*np.ones((N_P_bins))
    for P_i in range( N_P_bins ): #loop over log P
        
        # discard the bin if the minimum value lies in this bin  - i.e. we are certainly not complete here #
        if minpower > pbins[P_i]:
            continue
        # count  #
        sel_ind = np.where( (power >= pbins[P_i])  )[0]
        #powerhigh = pbins[P_i + 1]
        
        if len(sel_ind) > 0:
            zmax_h = zmax[sel_ind]
            zmin_h = zmin[sel_ind]
            areal_h = areal[sel_ind]
            fcor_h = fcor[sel_ind]
            
            vzmax = cosmo.distance.comoving_volume(zmax_h, **default_cosmo)
            vzmin = cosmo.distance.comoving_volume(zmin_h, **default_cosmo)
            
            vi = (vzmax - vzmin)*(area/(4.*np.pi))  #area/4pi gives per sterad
            rho[P_i]  =  np.sum(fcor_h/(areal_h*vi))      #sum of 1/v_max for each source
            #rhoerr[P_i] =  rho[P_i]*np.sqrt(float(len(sel_ind)))/float(len(sel_ind))
            rhoerr[P_i] = np.sqrt(np.sum((fcor_h/(areal_h*vi))**2.))
            num[P_i]  =  float(len(sel_ind))
        if verbose:
            print "{x} > {p1:7.2f} ({n:.0f}) : {rho:6.2e} +/- {rhoerr:6.2e}".format(x=xstr, p1=pbins[P_i], rho=rho[P_i], rhoerr=rhoerr[P_i], n=num[P_i])
            
    # per P bin - no, this is CLF
    rho = rho #/ dp
    rhoerr = rhoerr #/ dp
            
    return rho, rhoerr, num


def get_rho_Plim_f_areal(plimit, power, zmin, zmax, fcor, areal, area, ind=None, verbose=True, xstr="P"):
    '''
    plimit - bins of power (giving bin boundaries)
    power - log10(power) values
    zmin - for each source minimum z it could be observed (taking into account all selections)
    zmax - for each source maximum z it could be observed
    fcor - correction factor for completeness
    areal - the fractional areal coverage per source
    area - of survey in steradians
    ind - optically specify a selection index
    
    The differential comoving volume element dV_c/dz/dSolidAngle.
    Dimensions are volume per unit redshift per unit solid angle.
    Units are Mpc**3 Steradians^-1.
    '''
    print "Calculating density above power limit (f_areal): {n} sources".format(n=len(power))
    if ind is not None:
        power = power[ind]
        zmin = zmin[ind]
        zmax = zmax[ind]
        fcor = fcor[ind]
        areal = areal[ind]
        print "sub-selected {n} sources".format(n=len(power))
    
    
    rho  = np.nan
    rhoerr  = np.nan
    num   = np.nan
    
    minpower = np.min(power)
    if minpower > plimit:
        print " missing sources between {p2:.2f} and {p1:.2f}".format(p1=minpower, p2=plimit)
        return rho, rhoerr, num
        
    sel_ind = np.where( (power >= plimit)  )[0]
    #powerhigh = pbins[P_i + 1]
    
    if len(sel_ind) > 0:
        zmax_h = zmax[sel_ind]
        zmin_h = zmin[sel_ind]
        areal_h = areal[sel_ind]
        fcor_h = fcor[sel_ind]
        
        vzmax = cosmo.distance.comoving_volume(zmax_h, **default_cosmo)
        vzmin = cosmo.distance.comoving_volume(zmin_h, **default_cosmo)
        
        vi = (vzmax - vzmin)*(area/(4.*np.pi))  #area/4pi gives per sterad
        rho  =  np.sum(fcor_h/(areal_h*vi))      #sum of 1/v_max for each source
        #rhoerr[P_i] =  rho[P_i]*np.sqrt(float(len(sel_ind)))/float(len(sel_ind))
        rhoerr = np.sqrt(np.sum((fcor_h/(areal_h*vi))**2.))
        num  =  float(len(sel_ind))
    if verbose:
        print "{x} > {p1:7.2f} ({n:.0f}) : {rho:6.2e} +/- {rhoerr:6.2e}".format(x=xstr, p1=plimit, rho=rho, rhoerr=rhoerr, n=num)
            
    ## per P bin - nope this is > plimit
    #rho = rho #/ dp
    #rhoerr = rhoerr #/ dp
            
    return rho, rhoerr, num


def vmax_arr(m,z,mlim,area):
    '''
    vmax for array input
    '''
    N = len(m)
    if len(z) != N:
        print 'mismatch in input array sizes, m and z must be same length'
        return
    Avmax = np.nan*np.zeros(N)
    for i in range(N):
        Avmax[i] = vmax(m[i],z[i],mlim,area)
    return Avmax
    
def count_in_bins(xbins,xdata, norm=False):
    '''
    count_in_bins : histogram with errors
     input - xbins : array of x bins
             xdata : data to bin
    output - xmidrange : midpoints of bins
                nrange : counts per bin
               enrange : errors per bin [Poissonian] 
    '''
    nbins = len(xbins)-1
    nrange = np.zeros(nbins)
    enrange = np.zeros(nbins)
    xmidrange = np.zeros(nbins)
    for xi in range(nbins):
        xmidrange[xi] = 0.5*(xbins[xi] + xbins[xi+1])
        # count xdata between xi and xi+1
        nrange[xi] = np.sum((xdata>=xbins[xi])*(xdata<xbins[xi+1]) )
        enrange[xi] = np.sqrt(nrange[xi])
    # normalise:
    # note things outside the range covered will NOT be counted in the normalisation
    if norm:
        nF = 1./np.sum(nrange)
        nrange = nF*nrange
        enrange = nF*enrange
    
    return xmidrange, nrange, enrange
    
    
def sum_in_bins(xbins,xdata, ydata, norm=False):
    '''
    sum_in_bins : histogram with errors
     input - xbins : array of x bins
             xdata : data for bins
             ydata : data to sum in each x bin
    output - xmidrange : midpoints of bins
                nrange : counts per bin
               enrange : errors per bin [Poissonian] 
    '''
    nbins = len(xbins)-1
    Srange = np.zeros(nbins)
    eSrange = np.zeros(nbins)
    xmidrange = np.zeros(nbins)
    for xi in range(nbins):
        xmidrange[xi] = 0.5*(xbins[xi] + xbins[xi+1])
        # count xdata between xi and xi+1
        N = np.sum((xdata>=xbins[xi])*(xdata<xbins[xi+1]) )
        Srange[xi] = np.sum(ydata*(xdata>xbins[xi])*(xdata<xbins[xi+1]) )
        eSrange[xi] = Srange[xi] * np.sqrt(N)
    # normalise:
    # note things outside the range covered will NOT be counted in the normalisation
    if norm:
        nF = 1./np.sum(Srange)
        Srange = nF*Srange
        eSrange = nF*eSrange
    
    return xmidrange, Srange, eSrange
