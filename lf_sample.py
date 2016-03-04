import pyfits as pf
from fits_util import *
#from LF_util import *
import LF_util
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

import cosmolopy as cosmo


#from cosmos_sample_util import *
#from sdss_sample_util import *
#from lf_by_mass_util import *


#If no file is provided a cosmology of (Omega_Matter, Omega_Lambda, Omega_k, H0) = (0.3, 0.7, 0.0, 70.0) is assumed.
#default_cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.70}
#default_cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':0.70}
#cosmo.distance.set_omega_k_0(default_cosmo)

#print "NB: have you run sdss_BH_rlf.py?"

#name2 = 'COSMOS'


### selections on full sample ###

#define samples
#sampleA = {'name': "1", 'zlim_low': , 'zlim_high':  }







#class lf:
    #pass
    
    

#class lf_subsample(lf_sample):
    #pass

class lf_sample:
  
  
    def __init__(self, name, cat, zlow=np.nan, zhigh=np.nan, radio_fluxlim_faint = np.nan, opt_fluxlim_faint = np.nan, opt_fluxlim_bright = np.nan, area=np.nan):
        
        
        self.name = name
        self.zlim_low = zlow
        self.zlim_high = zhigh
        self.cat = cat
        self.Nsrcs = len(cat)
        
        
        self.area = area
        
        self.radio_fluxlim_faint = radio_fluxlim_faint
        self.opt_fluxlim_faint = opt_fluxlim_faint
        self.opt_fluxlim_bright = opt_fluxlim_bright
    
        self.set_power_z_limit()
        
        
        self.LF_x = None
        self.LF_xerr = None
        self.LF_rho = None
        self.LF_rhoerrup = None
        self.LF_rhoerrlow = None
        self.LF_num = None
        
        self.CLF_x = None
        self.CLF_xerr = None
        self.CLF_rho = None
        self.CLF_rhoerrup = None
        self.CLF_rhoerrlow = None
        self.CLF_num = None
        
        self.rhoPlim_x = None
        self.rhoPlim_xerr = None
        self.rhoPlim_rho = None
        self.rhoPlim_rhoerrup = None
        self.rhoPlim_rhoerrlow = None
        self.rhoPlim_num = None
        
        return
    
    def copy(self):
        return lf_sample(self.name, self.cat, zlow=self.zlim_low, zhigh=self.zlim_high, radio_fluxlim_faint = self.radio_fluxlim_faint, opt_fluxlim_faint = self.opt_fluxlim_faint, opt_fluxlim_bright = self.opt_fluxlim_bright, area=self.area)
    
    def copy_subcat(self, name, cat):
        return lf_sample(name, cat, zlow=self.zlim_low, zhigh=self.zlim_high, radio_fluxlim_faint = self.radio_fluxlim_faint, opt_fluxlim_faint = self.opt_fluxlim_faint, opt_fluxlim_bright = self.opt_fluxlim_bright, area=self.area)
    
    
    def sub_z_sample(self, name, zlow, zhigh):
        ''' make a new subsample with name 'name' from the z range provided'''
        
        #new_self = self.copy()
        #new_self.name = name
        thisname = self.name
        
        
        ind_z = np.where((self.cat['z'] > zlow) & (self.cat['z'] <= zhigh))[0]
        # handle no sources in selection
        if len(ind_z) == 0:
            return None
        print 'subsample z : {n1} out of {n2} sources selected '.format(n1=len(ind_z),n2=len(self.cat))
        new_self = self.sub_sample_ind(thisname+name, ind_z)
        
        new_self.zlim_low = zlow
        new_self.zlim_high = zhigh
        
        new_self.set_power_z_limit()
        
        new_self.calc_zmin_zmax()
        
        return new_self
    
    
    def sub_sample_by_field(self, name, field, fieldlimlow, fieldlimhigh, req_new_volumes=False):
        ''' make a new subsample with name 'name' from the z range provided'''
        
        #new_self = self.copy()
        #new_self.name = name
                
        ind = np.where((self.cat[field] > fieldlimlow) & (self.cat[field] <= fieldlimhigh))[0]
        # handle no sources in selection
        if len(ind) == 0:
            return None
        print 'subsample {n} : {n1} out of {n2} sources selected on field {f}'.format(n=name,n1=len(ind),n2=len(self.cat),f=field)
        new_self = self.sub_sample_ind(name, ind)
        
        # calculate new zmin, zmax if needed
        if req_new_volumes:
            new_self.calc_zmin_zmax()
        
        return new_self
    
    def sub_sample_ind(self, name, ind):
        ''' make a new subsample with name 'name' from the catalogue indicies'''
        
        # handle no sources in selection
        if len(ind) == 0:
            return None
        new_self = self.copy_subcat(name, self.cat[ind])
        #new_self.name = name
        #new_self.cat = new_self.cat[ind]
        #new_self.Nsrcs = len(new_self.cat)
        
        return new_self
    
    #def load_cat(self, cat, fluxunit='mJy'):
        
        ## subsample #

        #cat = cat.copy()



        #self.zz = cat.z
        
        #fKmag = cat.Ks_tot   # in fluxes in cat
        ##fKmag2 = 3631.*10.**(-0.4*Kmag2)    # AB mag
        #self.Klum = OpticalLuminosity(fKmag, zz)
        
        #scaleflux = 1.
        #if fluxunit == 'mJy':
            #scaleflux = 1e-3
        #elif fluxunit == 'Jy':
            #scaleflux = 1.
            
        #self.ff = (cat.Si_1_4_)*scaleflux   #in Jy
        #self.power = np.log10(RadioPower(self.ff, self.zz, alpha=0.8))  # assumed spec ind

        #self.Fcor = np.ones(len(ff))
        #self.areal = get_area_for_flux(ff)

        ## stellar masses
        #self.smass = cat_vlacosmos.lmass


        #self.Nsrc = len(zz)
        
        #vals = func()
        
        #self.z = vals['z']
        #self.optical_flux = vals['optical_flux']
        #self.optical_flux = vals['optical_flux']
        #self.Nsrc =  
        
        #return
    
    def plot_zmin_zmax(self):
        
        f=plt.figure()
        if 'zmax' in self.cat.dtype.names:
            plt.scatter(self.cat['power'], self.cat['zmax'], c=self.cat['z'], marker='v', s=20, edgecolor='none')
        if 'Pzmax' in self.cat.dtype.names:
            plt.scatter(self.cat['power'], self.cat['Pzmax'], c=self.cat['z'], marker='v', s=10, edgecolor='none')
        if 'Optzmax' in self.cat.dtype.names:
            plt.scatter(self.cat['power'], self.cat['Optzmax'], c=self.cat['z'], marker='v', s=10, edgecolor='none')
        if 'Optzmin' in self.cat.dtype.names:
            plt.scatter(self.cat['power'], self.cat['Optzmin'], c=self.cat['z'], marker='^', s=10, edgecolor='none')
        if 'zmin' in self.cat.dtype.names:
            plt.scatter(self.cat['power'], self.cat['zmin'], c=self.cat['z'], marker='^', s=20, edgecolor='none')
        plt.minorticks_on()
        plt.xlabel("power")
        plt.ylabel("z min/max")
        plt.savefig('sample_zmax_zmin_{name}.png'.format(name=self.name))
        plt.close(f)
        
        return
    
    
    def calc_zmin_zmax(self, plot=False):
        
        from numpy.lib.recfunctions import rec_append_fields
        
        
        haspower = False
        if 'power' in self.cat.dtype.names:
            haspower = True
        
        if haspower:
            print "getting zmin zmax for radio-optical sample {n}".format(n=self.name)
        else:
            print "getting zmin zmax for optical sample {n}".format(n=self.name)
        
        #at what redshift does each source fall below the flux density limit?
        if haspower:
            Pzmax = LF_util.get_zmax(self.cat['z'], 10.**self.cat['power'], self.radio_fluxlim_faint, stype='Radio',filename='zmax.radio.sav.%s.npy' %(self.name), clobber=0)
        Optzmax = LF_util.get_zmax(self.cat['z'], self.cat['opt_lum'], self.opt_fluxlim_faint, stype='Optical',filename='zmax.optical.sav.%s.npy' %(self.name), clobber=0)
        Optzmin = LF_util.get_zmin(self.cat['z'], self.cat['opt_lum'], self.opt_fluxlim_bright, stype='Optical',filename='zmin.optical.sav.%s.npy' %(self.name), clobber=0)

        if haspower:
            if 'Pzmax' not in self.cat.dtype.names:
                self.cat = rec_append_fields(self.cat, ('Pzmax'),  (Pzmax) )
            else:
                self.cat['Pzmax'] = Pzmax
            
        if 'Optzmax' not in self.cat.dtype.names:
            self.cat = rec_append_fields(self.cat, ('Optzmax', 'Optzmin'),  (Optzmax, Optzmin) )
        else:
            self.cat['Optzmax'] = Optzmax
            self.cat['Optzmin'] = Optzmin

        #np.savez('sample_{name}.npz'.format(name=self.name), z=self.cat['z'], sm=self.smass, P=self.power )
        
        if plot:
            self.plot_zmin_zmax()
        
        if haspower:
            ### Combine zmax's from radio, optical and z selections
            zmax = self.zlim_high*np.ones(len(Optzmax))
            t1 = np.minimum(Optzmax, Pzmax)
            t2 = np.minimum(t1, zmax)
            zmax = t2
        else:
            ### Combine zmax's from radio, optical and z selections
            zmax = self.zlim_high*np.ones(len(Optzmax))
            t2 = np.minimum(Optzmax, zmax)
            zmax = t2
                
                
        ### Combine zmins's from optical and z selections
        zmin = self.zlim_low*np.ones(len(Optzmin))
        t1 = np.maximum(Optzmin, zmin)
        zmin = t1
        
        if 'zmin' not in self.cat.dtype.names:
            self.cat = rec_append_fields(self.cat, ('zmin'),  (zmin) )
        else:
            self.cat['zmin'] = zmin
            
        if 'zmax' not in self.cat.dtype.names:
            self.cat = rec_append_fields(self.cat, ('zmax'),  (zmax) )
        else:
            self.cat['zmax'] = zmax
            
            
        if 'Fcor' not in self.cat.dtype.names:
            self.cat = rec_append_fields(self.cat, ('Fcor'),  ( np.ones(len(zmax)) ) )
        else:
            self.cat['Fcor'] = np.ones(len(zmax))
            
        if 'areal' not in self.cat.dtype.names:
            self.cat = rec_append_fields(self.cat, ('areal'),  (np.ones(len(zmax))) )
        else:
            self.cat['areal'] = np.ones(len(zmax))

        
        return
    
    
    def set_power_z_limit(self):
        zz = np.linspace(self.zlim_low, self.zlim_high, 10)
        Plim = np.log10(LF_util.RadioPower(self.radio_fluxlim_faint, zz, alpha=0.8)) 
        self.power_z_limit = np.vstack((zz, Plim))
        return
    
    
    
    def compute_LF(self, pgrid, maskbins=None, CV_f=None, ignoreMinPower=False):
        
        
        logp_radio_lf = (pgrid[:-1]+pgrid[1:])/2.
        dlogp_radio_lf = (pgrid[1:]-pgrid[:-1])/2.
        
        #print self.cat['power']
        rho, rhoerr, num = LF_util.get_LF_f_areal(pgrid, self.cat['power'], self.cat['zmin'], self.cat['zmax'], self.cat['Fcor'], self.cat['areal'], self.area, ignoreMinPower=ignoreMinPower)
        #print self.cat['power']
        #print pgrid
        #print type(pgrid)
        
        rhoerrlow = rhoerr
        rhoerrup = rhoerr
        
        if maskbins is not None:
            if isinstance(maskbins, np.ndarray):
                if len(maskbins) == 2:
                    if np.isfinite(maskbins[0]): 
                        rho[logp_radio_lf<=maskbins[0]] *= np.nan
                        rhoerrlow[logp_radio_lf<=maskbins[0]] *= np.nan
                        rhoerrup[logp_radio_lf<=maskbins[0]] *= np.nan
                    if np.isfinite(maskbins[1]): 
                        rho[logp_radio_lf>=maskbins[1]] *= np.nan
                        rhoerrlow[logp_radio_lf>=maskbins[1]] *= np.nan
                        rhoerrup[logp_radio_lf>=maskbins[1]] *= np.nan
                else:
                    print "maskbins invalid"
            else:
                print "maskbins invalid"
            
        #print '###', rhoerrup
        # add cosmic variance errors #
        if CV_f is not None:
            #self.LF_rhoerrup = np.sqrt(self.LF_rhoerrup**2. + (num*CV_f)**2.)
            #self.LF_rhoerrlow = np.sqrt(self.LF_rhoerrlow**2. + (num*CV_f)**2.)
            #rhoerrup = np.sqrt(rhoerrup**2. + (CV_f*rho)**2.)
            #rhoerrlow = np.sqrt(rhoerrlow**2. + (CV_f*rho)**2.)
            rhoerrup = rhoerrup*np.sqrt(1 + (CV_f)**2.)
            rhoerrlow = rhoerrlow*np.sqrt(1 + (CV_f)**2.)
            
        #print '###', rhoerrup
            
        # fix lower errors to be rho - for plotting logscale #
        rhoerrlow[rhoerrlow>=rho] = rho[rhoerrlow>=rho]*0.9999
        
        self.LF_x = logp_radio_lf
        self.LF_xerr = dlogp_radio_lf
        self.LF_rho = rho
        self.LF_rhoerrup = rhoerrup
        self.LF_rhoerrlow = rhoerrlow
        self.LF_num = num
        
        
        #print '* xact',logp_radio_lf
        #print '* x',self.LF_x
        #print '* xerract',dlogp_radio_lf
        #print '* xerr',self.LF_xerr
        #print '* rhoact',rho
        #print '* rho',self.LF_rho
        
        return rho, rhoerrlow, rhoerrup, num

    def compute_SMF(self, smgrid, masscompleteness):
        '''
masscompleteness - mass(z) function describing the completeness envelope
        '''
        
        logsm = (smgrid[:-1]+smgrid[1:])/2.
        dlogsm = (smgrid[1:]-smgrid[:-1])/2.
        
        
        # select mass-complete #
        smenv_ind =  np.where(self.cat['smass'] >= masscompleteness(z=self.cat['z']))[0]  # Stellar mass cut #
        print "{n1} of {n2} sources selected on SM envelope".format(n2=len(self.cat),n1= len(smenv_ind))
        sm_complete_sample = self.sub_sample_ind('m', smenv_ind )
        
        # zmax determined from mass-completeness also #
        SMzmax = masscompleteness(m=sm_complete_sample.cat['smass'])
        sm_complete_sample.cat['zmax']= np.minimum(sm_complete_sample.cat['zmax'], SMzmax)
    
        
        #print self.cat['power']
        rho, rhoerr, num = LF_util.get_LF_f_areal(smgrid, sm_complete_sample.cat['smass'], sm_complete_sample.cat['zmin'], sm_complete_sample.cat['zmax'], sm_complete_sample.cat['Fcor'], sm_complete_sample.cat['areal'], sm_complete_sample.area, xstr="SM")
        #print self.cat['power']
        #print pgrid
        #print type(pgrid)
        
        # fix lower errors to be rho #
        rhoerrlow = rhoerr
        rhoerrup = rhoerr
        rhoerrlow[rhoerrlow>=rho] = rho[rhoerrlow>=rho]*0.9999
        
        self.SMF_x = logsm
        self.SMF_xerr = dlogsm
        self.SMF_rho = rho
        self.SMF_rhoerrup = rhoerrup
        self.SMF_rhoerrlow = rhoerrlow
        self.SMF_num = num
        
        #print '* xact',logp_radio_lf
        #print '* x',self.LF_x
        #print '* xerract',dlogp_radio_lf
        #print '* xerr',self.LF_xerr
        #print '* rhoact',rho
        #print '* rho',self.LF_rho
        
        return rho, rhoerrlow, rhoerrup, num
    
    
    def compute_CLF(self,  pgrid):
        
        logp_radio_lf = (pgrid[:-1]+pgrid[1:])/2.
        dlogp_radio_lf = (pgrid[1:]-pgrid[:-1])/2.
        
        rho, rhoerr, num = LF_util.get_CLF_f_areal(pgrid, self.cat['power'], self.cat['zmin'], self.cat['zmax'], self.cat['Fcor'], self.cat['areal'], self.area)
        
        # fix lower errors to be rho #
        rhoerrlow = rhoerr
        rhoerrup = rhoerr
        rhoerrlow[rhoerrlow>=rho] = rho[rhoerrlow>=rho]*0.9999
        
        self.CLF_x = logp_radio_lf
        self.CLF_xerr = dlogp_radio_lf
        self.CLF_rho = rho
        self.CLF_rhoerrup = rhoerrup
        self.CLF_rhoerrlow = rhoerrlow
        self.CLF_num = num
        
        return rho, rhoerrlow, rhoerrup, num
    
    
    def compute_rhoPlim(self,  Plimit):
        
        
        rho, rhoerr, num = LF_util.get_rho_Plim_f_areal(Plimit, self.cat['power'], self.cat['zmin'], self.cat['zmax'], self.cat['Fcor'], self.cat['areal'], self.area, xstr="SM")
        
        # fix lower errors to be rho #
        rhoerrlow = rhoerr
        rhoerrup = rhoerr
        if rhoerrlow>=rho:
            rhoerrlow = rho*0.9999
        
        self.rhoPlim_plim = Plimit
        self.rhoPlim_rho = rho
        self.rhoPlim_rhoerrup = rhoerrup
        self.rhoPlim_rhoerrlow = rhoerrlow
        self.rhoPlim_num = num
        
        return rho, rhoerrlow, rhoerrup, num
    
    
        
    
    
    
    
    
    
    
    
    
    
    
  
