from LF_util import RadioFlux
import scipy.optimize as so
from astropy.io import fits
import numpy as np
default_cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.70}
import time
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def get_zmax(L,rms,factor,alpha=0.8):
    #find z such that radio flux = factor * rms by root-finding
    return so.brentq(lambda z: RadioFlux(L,z,alpha)-factor*rms,0,10)

def get_vmax(L,rms,factor,area,alpha=0.8):
    # directly copied from vmax in Wendy's code
    # return units are Mpc^3

    zlim=get_zmax(L,rms,factor,alpha)

    degtost = 4.*180.**2./np.pi  # sq degrees to steradians
    domega = area/degtost
    DL = cosmo.luminosity_distance(zlim)
    vc = (4.*np.pi/3.)*(DL/(1.+zlim))**3.  # volume at zlim
    vmax = domega*vc
    return vmax

class rmsmapz(object):
    def __init__(self,map,sampling=100):
        # map is a 2D nan-blanked FITS array as produced by PyBDSM
        hdu=fits.open(map)
        self.data=hdu[0].data
        self.sampling=sampling
        self.min=np.nanmin(self.data)
        self.max=np.nanmax(self.data)
        print self.min,self.max
        self.bins=np.logspace(np.log10(self.min),np.log10(self.max),sampling+1)
        self.centres=0.5*(self.bins[:-1]+self.bins[1:])
        cdata=self.data[~np.isnan(self.data)]
        self.hist=np.asfarray(np.histogram(cdata,self.bins)[0])/len(cdata)
        self.area=len(cdata)*hdu[0].header['CDELT1']**2.0 # in sq. deg
        print 'Read',map,'area is',self.area,'deg^2'

    def vmax(self,L):
        # compute the Smolcic Ak*V_max(L,rms_k) sum
        vm=self.get_vmax(L,self.centres)
        return np.sum(self.hist*vm)

    def interp_setup(self,Lmin,Lmax,factor,alpha=0.8,sampling=100):
        # we solve in terms of r = L/(factor*rms)
        self.factor=factor
        self.alpha=alpha
        rmin=Lmin/(factor*self.max)
        rmax=Lmax/(factor*self.min)
        print 'Using r range',rmin,'to',rmax
        rvals=np.linspace(np.log10(rmin),np.log10(rmax),sampling)
        zvals=np.zeros_like(rvals)
        for i in range(len(rvals)):
            print i,rvals[i]
            zvals[i]=so.brentq(lambda z: RadioFlux(10**rvals[i],z,alpha)-1,0,20)
        self.get_zmax_interp=interp1d(rvals,zvals,kind='cubic')

    def get_vmax(self,L,rms):
        # uses interpolation so can work with array-like rms
        zlim=self.get_zmax_interp(np.log10(L/(rms*self.factor)))

        degtost = 4.*180.**2./np.pi  # sq degrees to steradians
        domega = self.area/degtost
        DL = cosmo.luminosity_distance(zlim)/u.Mpc
        vc = (4.*np.pi/3.)*(DL/(1.+zlim))**3.  # volume at zlim
        vmax = domega*vc
        return vmax


#if 1:
    #import compute_vmax_rms as vm
    #CATPATH = '/local/wwilliams/phd/multi_bootes/'
    #rms=vm.rmsmapz(CATPATH+'/LOFAR150/rms.radio_opt_masked.fits')
    #print vm.get_vmax(1e24,rms.min,5.0,rms.area)
    #rms.interp_setup(1e21,1e27,5.0)
    #print rms.vmax(1e24)
    
    #steps=100
    #t0 = time.time()
    #for i in range(steps):
        #result=rms.vmax(1e24)
    #t1 = time.time()
    #print 'Each call took',(t1-t0)/steps,'seconds'