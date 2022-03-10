import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import astropy.io.fits as pf
#from utils.fits_util import *
from astropy.table import Table, Column
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
acosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from . import util as LF_util


class lf_sample:
  
  
    def __init__(self, name, cat, zlow=0, zhigh=20, radio_fluxlim_faint=np.nan, radio_fluxlim_bright=np.nan, opt_fluxlim_faint=np.nan, opt_fluxlim_bright=np.nan, domega=np.nan, area=np.nan, rmsmap=None, completeness=None, savedir=None):
        '''
        define the survey area either as:
            area - in square degrees!
            domega - is the faction of the full sky
        '''
        
        if savedir is None:
            savedir='lfdata_'+name
        self.savedir = savedir
        if not os.path.isdir(self.savedir):
            print('creating sample directory ',self.savedir)
            os.mkdir(self.savedir)
        #add trailing /
        if self.savedir[-1] != '/':
            self.savedir = self.savedir + '/'
        
        self.name = name
        self.zlim_low = zlow
        self.zlim_high = zhigh
        
        self.cat = cat
        self.Nsrcs = len(cat)
        
        if np.isfinite(domega):
            self.domega = domega
            self.area = domega *4*np.pi * (180./np.pi) **2.
        elif np.isfinite(area):
            self.area = area 
            self.domega = area / (4*np.pi * (180./np.pi) **2.)
        
        self.rmsmap = None
        if rmsmap is not None:
            if isinstance(rmsmap, LF_util.rmsmapz):
                self.rmsmap = rmsmap
                self.rmsmap.interp_setup(1e30,1e20,5.0)
                # update area
                if self.domega != self.rmsmap.domega:
                    self.domega = self.rmsmap.domega
                    self.area = self.rmsmap.area
                    print('WARNING updating area, using {A:.3f} deg^2 from the rms map ({p:.3f})'.format(A=self.area, p=self.domega)) 
                
            elif isinstance(rmsmap, str):
                if os.path.exists(rmsmap):
                    self.rmsmap = LF_util.rmsmapz(rmsmap)
                    self.rmsmap.interp_setup(1e30,1e20,5.0)
                    # update area
                    if self.domega != self.rmsmap.domega:
                        self.domega = self.rmsmap.domega
                        self.area = self.rmsmap.area
                        print('WARNING updating area, using {A:.3f} deg^2 from the rms map ({p:.3f})'.format(A=self.area, p=self.domega)) 
                    
                else:
                    print("WARNING given rmsmap {map:s} does not exist, continuing without".format(map=rmsmap))
            else:
                print("WARNING given rmsmap type not understood, continuing without")
            
            
        Vzlow = self.domega*acosmo.comoving_volume(self.zlim_low).value
        Vzhigh = self.domega*acosmo.comoving_volume(self.zlim_high).value
        self.Vzlim_low = Vzlow
        self.Vzlim_high = Vzhigh
            
        self.completeness = None
        if completeness is not None:
            if isinstance(completeness, LF_util.completenessf):
                self.completeness = completeness
            elif isinstance(completeness, str):
                if os.path.exists(completeness):
                    self.completeness = LF_util.completenessf(completeness)
                else:
                    print("WARNING given completeness {map:s} does not exist, continuing without".format(map=completeness))
            else:
                    print("WARNING given completeness type not understoof, continuing without".format(map=completeness))
        
        
        self.radio_fluxlim_faint = radio_fluxlim_faint
        self.radio_fluxlim_bright = radio_fluxlim_bright
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
        return lf_sample(self.name, self.cat, zlow=self.zlim_low, zhigh=self.zlim_high, radio_fluxlim_faint=self.radio_fluxlim_faint, radio_fluxlim_bright=self.radio_fluxlim_bright, opt_fluxlim_faint=self.opt_fluxlim_faint, opt_fluxlim_bright=self.opt_fluxlim_bright, domega=self.domega, rmsmap=self.rmsmap, completeness=self.completeness, savedir=self.savedir)
    
    def copy_subcat(self, name, cat):
        return lf_sample(name, cat, zlow=self.zlim_low, zhigh=self.zlim_high, radio_fluxlim_faint = self.radio_fluxlim_faint, radio_fluxlim_bright=self.radio_fluxlim_bright, opt_fluxlim_faint=self.opt_fluxlim_faint, opt_fluxlim_bright=self.opt_fluxlim_bright, domega=self.domega, rmsmap=self.rmsmap, completeness=self.completeness, savedir=self.savedir)
    
    
    def make_z_samples(self, zbins, dolf=True, pbins=None, savelf=True, plot=False, forcecalc=False, ignoreMinPower=False, catcol='power'):
        '''make_z_samples(self, zbins, dolf=True, savelf=True, plot=False, forcecalc=False)
        make a set of samples over redshift bins
        '''
        z_samples = []
        for i in range(len(zbins)):
            print(zbins[i][0], zbins[i][1])
            z_sample_i = self.sub_z_sample('zbin{i:02d}_{z1:.2f}_{z2:.2f}'.format(i=i, z1=zbins[i][0], z2=zbins[i][1]), zbins[i][0], zbins[i][1], forcecalc=forcecalc, savefiles=True, ignoreMinPower=ignoreMinPower)

            if plot:
                z_sample_i.plot_Vzmin_Vzmax(keep=True, catcol=catcol)
                z_sample_i.plot_Vzmin_Vzmax(keep=True, ccol='opt_mag', cbarlabel='$r$', saveext='rmag', catcol=catcol)
                try:
                    z_sample_i.plot_Vzmin_Vzmax(keep=True, ccol='opt_col', cbarlabel='$g-r$', saveext='col', catcol=catcol)
                except:
                    print('cant plot colour: opt_col not present')

            if z_sample_i is not None:
                if dolf and pbins is not None:
                    z_sample_i.compute_LF(pbins, catcol=catcol)
                    if savelf:
                        np.savez('{ddir}{n:s}.npz'.format(ddir=self.savedir, n=z_sample_i.name),x=z_sample_i.LF_x, xerr=z_sample_i.LF_xerr, N=z_sample_i.LF_num,  y=z_sample_i.LF_rho, yerrup=z_sample_i.LF_rhoerrup, yerrlow=z_sample_i.LF_rhoerrlow)
                elif dolf:
                    print('pbins is not set and dolf is on, please specify pbins')

                
                
            z_samples.append(z_sample_i)
        
        return z_samples
    
    
    def sub_z_sample(self, name, zlow, zhigh,forcecalc=False, savefiles=True, plot=False, ignoreMinPower=False):
        ''' make a new subsample with name 'name' from the z range provided'''
        
        #new_self = self.copy()
        #new_self.name = name
        thisname = self.name
        
        
        ind_z = np.where((self.cat['z'] >= zlow) & (self.cat['z'] < zhigh))[0]
        # handle no sources in selection
        if len(ind_z) == 0:
            return None
        print('subsample z : {n1} out of {n2} sources selected '.format(n1=len(ind_z),n2=len(self.cat)))
        new_self = self.sub_sample_ind(thisname+'_'+name, ind_z)
        
        new_self.zlim_low = zlow
        new_self.zlim_high = zhigh
        new_self.Vzlim_low = self.domega * acosmo.comoving_volume(zlow).value
        new_self.Vzlim_high = self.domega * acosmo.comoving_volume(zhigh).value
        
        new_self.set_power_z_limit()
        
        #new_self.calc_zmin_zmax()
        new_self.calc_Vzmin_Vzmax(forcecalc=forcecalc, savefiles=savefiles, plot=plot)
        
        return new_self
    
    
    def sub_sample_by_field(self, name, field, fieldlimlow, fieldlimhigh, req_new_volumes=False, plot=False):
        ''' make a new subsample with name 'name' from the z range provided'''
        
        #new_self = self.copy()
        #new_self.name = name
                
        ind = np.where((self.cat[field] > fieldlimlow) & (self.cat[field] <= fieldlimhigh))[0]
        # handle no sources in selection
        if len(ind) == 0:
            return None
        print('subsample {n} : {n1} out of {n2} sources selected on field {f}'.format(n=name,n1=len(ind),n2=len(self.cat),f=field))
        new_self = self.sub_sample_ind(name, ind)
        
        # calculate new zmin, zmax if needed
        if req_new_volumes:
            #new_self.calc_zmin_zmax()
            new_self.calc_Vzmin_Vzmax(plot=plot)
        
        return new_self
    
    
    def sub_sample_by_field_value(self, name, field, value, req_new_volumes=False, plot=False):
        ''' make a new subsample with name 'name' from the z range provided'''
        
        #new_self = self.copy()
        #new_self.name = name
                
        ind = np.where((self.cat[field] == value ))[0]
        # handle no sources in selection
        if len(ind) == 0:
            return None
        print('subsample {n} : {n1} out of {n2} sources selected on field {f}'.format(n=name,n1=len(ind),n2=len(self.cat),f=field))
        new_self = self.sub_sample_ind(name, ind)
        
        # calculate new zmin, zmax if needed
        if req_new_volumes:
            #new_self.calc_zmin_zmax()
            new_self.calc_Vzmin_Vzmax(plot=plot)
        
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
        plt.savefig('inspect_zmax_zmin_{name}.png'.format(name=self.name))
        plt.close(f)
        
        return
    
    def plot_Vzmin_Vzmax(self, logV=True, vmin=0, vmax=1, keep=False, ccol='z', cbarlabel='$z$', saveext='', catcol='power'):
        
        fig=plt.figure()
        ax = fig.add_subplot(111)
        if 'Vzmax' in self.cat.dtype.names:
            c = ax.scatter(self.cat[catcol], self.cat['Vzmax'], c=self.cat[ccol], marker='v', s=20, edgecolor='none',label='Max')
        if 'PVzmax' in self.cat.dtype.names:
            c = ax.scatter(self.cat[catcol], self.cat['PVzmax'], c=self.cat[ccol], marker='v', s=10, edgecolor='none',label='Pmax')
        if 'PVzmin' in self.cat.dtype.names:
            c = ax.scatter(self.cat[catcol], self.cat['PVzmin'], c=self.cat[ccol], marker='^', s=10, edgecolor='none',label='Pmin')
        if 'OptVzmax' in self.cat.dtype.names:
            c = ax.scatter(self.cat[catcol], self.cat['OptVzmax'], c=self.cat[ccol], marker='1', s=10, edgecolor='none',label='Optmax')
        if 'OptVzmin' in self.cat.dtype.names:
            c = ax.scatter(self.cat[catcol], self.cat['OptVzmin'], c=self.cat[ccol], marker='^', s=10, edgecolor='none',label='Optmin')
        if 'Vzmin' in self.cat.dtype.names:
            c = ax.scatter(self.cat[catcol], self.cat['Vzmin'], c=self.cat[ccol], marker='2', s=20, edgecolor='none',label='Min')
        if logV:
            ax.semilogy()
        cbar=fig.colorbar(c)
        cbar.set_label(cbarlabel)
        ax.legend(loc='upper left')
        plt.minorticks_on()
        ax.set_xlabel("$\log L$  [W/Hz]")
        x1,x2 = ax.get_xlim()
        ax.hlines(self.Vzlim_low, x1, x2, 'gray', alpha=0.5)
        ax.hlines(self.Vzlim_high, x1, x2, 'gray', alpha=0.5)
        ax.text(x1, self.Vzlim_low, str(self.zlim_low))
        ax.text(x1, self.Vzlim_high, str(self.zlim_high))
        ax.set_xlim(x1,x2)
        ax.set_ylabel("$V_{\mathrm{zmin/zmax}}$  [Mpc$^3$]")
        if saveext == '':
            fig.savefig('{ddir}{name}.inspect_Vzmax_Vzmin.png'.format(ddir=self.savedir, name=self.name))
        else:
            fig.savefig('{ddir}{name}.inspect_Vzmax_Vzmin_{s}.png'.format(ddir=self.savedir, name=self.name, s=saveext))
        if not keep:
            plt.close(fig)
        
        return
    
    def calc_zmin_zmax(self, plot=False):
        
        haspower = False
        if 'power' in self.cat.dtype.names:
            haspower = True
        
        if haspower:
            print("getting zmin zmax for radio-optical sample {n}".format(n=self.name))
        else:
            print("getting zmin zmax for optical sample {n}".format(n=self.name))
        
        #at what redshift does each source fall below the flux density limit?
        if haspower:
            Pzmax = LF_util.get_zmax(self.cat['z'], 10.**self.cat['power'], self.radio_fluxlim_faint, stype='Radio',filename='{ddir}{name}.zmax.radio.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=0)
        Optzmax = LF_util.get_zmax(self.cat['z'], self.cat['opt_lum'], self.opt_fluxlim_faint, stype='Optical',filename='{ddir}{name}.zmax.optical.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=0)
        Optzmin = LF_util.get_zmin(self.cat['z'], self.cat['opt_lum'], self.opt_fluxlim_bright, stype='Optical',filename='{ddir}{name}.zmin.optical.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=0)

        if haspower:
            if 'Pzmax' not in self.cat.dtype.names:
                self.cat.add_column(Column(Pzmax, 'Pzmax') )
            else:
                self.cat['Pzmax'] = Pzmax
            
        if 'Optzmax' not in self.cat.dtype.names:
            self.cat.add_column(Column(Optzmax, 'Optzmax') )
            self.cat.add_column(Column(Optzmin, 'Optzmin') )
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
            self.cat.add_column(Column(zmin, 'zmin'))
        else:
            self.cat['zmin'] = zmin
            
        if 'zmax' not in self.cat.dtype.names:
            self.cat.add_column(Column(zmax, 'zmax'))
        else:
            self.cat['zmax'] = zmax
            
            
        if 'Fcor' not in self.cat.dtype.names:
            self.cat.add_column(Column(np.ones(len(zmax)), 'Fcor'))
        else:
            self.cat['Fcor'] = np.ones(len(zmax))
            
        if 'areal' not in self.cat.dtype.names:
            self.cat.add_column(Column(np.ones(len(zmax)), 'areal'))
        else:
            self.cat['areal'] = np.ones(len(zmax))

        
        return
    
    
    def calc_Vzmin_Vzmax(self, plot=True, verbose=True, forcecalc=False, savefiles=True):
        
        haspower = False
        if 'power' in self.cat.dtype.names:
            haspower = True
        
        if haspower:
            print("getting Vzmin Vzmax for radio-optical sample {n}".format(n=self.name))
        else:
            print("getting Vzmin Vzmax for optical sample {n}".format(n=self.name))
        
        #at what redshift does each source fall below the flux density limit?
        if haspower:
            PVzmax = LF_util.get_Vzmax(self.cat['z'], 10.**self.cat['power'], self.radio_fluxlim_faint, self.domega, zmax=self.zlim_high, rmsmap=self.rmsmap, completeness=self.completeness, stype='Radio',filename='{ddir}{name}.Vzmax.radio.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=forcecalc, savefile=savefiles)
            ## argh, why???
            PVzmin = LF_util.get_Vzmin(self.cat['z'], 10.**self.cat['power'], self.radio_fluxlim_bright, self.domega, zmin=self.zlim_low, rmsmap=self.rmsmap, completeness=self.completeness, stype='Radio',filename='{ddir}{name}.Vzmin.radio.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=forcecalc, savefile=savefiles)
            #PVzmin = np.ones(len(self.cat)) * self.domega*acosmo.comoving_volume(self.zlim_low).value
        OptVzmax = LF_util.get_Vzmax(self.cat['z'], self.cat['opt_lum'], self.opt_fluxlim_faint, self.domega, stype='Optical',filename='{ddir}{name}.Vzmax.optical.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=forcecalc, savefile=savefiles)
        OptVzmin = LF_util.get_Vzmin(self.cat['z'], self.cat['opt_lum'], self.opt_fluxlim_bright, self.domega, stype='Optical',filename='{ddir}{name}.Vzmin.optical.sav.npy'.format(ddir=self.savedir, name=self.name), clobber=forcecalc, savefile=savefiles)
        #import ipdb
        #if haspower:
            #if np.any(PVzmin >= PVzmax):
                #ipdb.set_trace()

        if haspower:
            if 'PVzmax' not in self.cat.dtype.names:
                self.cat.add_column(Column(PVzmax, 'PVzmax'))
            else:
                self.cat['PVzmax'] = PVzmax
            if 'PVzmin' not in self.cat.dtype.names:
                self.cat.add_column(Column(PVzmin, 'PVzmin'))
            else:
                self.cat['PVzmin'] = PVzmin
            
        if 'OptVzmax' not in self.cat.dtype.names:
            self.cat.add_column(Column(OptVzmax, 'OptVzmax'))
            self.cat.add_column(Column(OptVzmin, 'OptVzmin'))
        else:
            self.cat['OptVzmax'] = OptVzmax
            self.cat['OptVzmin'] = OptVzmin

        #np.savez('sample_{name}.npz'.format(name=self.name), z=self.cat['z'], sm=self.smass, P=self.power )
        
        
        Vzlim_high = self.domega*acosmo.comoving_volume(self.zlim_high).value
        if haspower:
            ### Combine Vzmax's from radio, optical and z selections
            Vzmax = Vzlim_high*np.ones(len(OptVzmax))
            t1 = np.minimum(OptVzmax, PVzmax)
            t2 = np.minimum(t1, Vzmax)
            Vzmax = t2
        else:
            ### Combine Vzmax's from radio, optical and z selections
            Vzmax = Vzlim_high*np.ones(len(OptVzmax))
            t2 = np.minimum(OptVzmax, Vzmax)
            Vzmax = t2
                
                
        ### Combine Vzmins's from optical and z selections
        Vzlim_low = self.domega*acosmo.comoving_volume(self.zlim_low).value
        Vzmin = Vzlim_low*np.ones(len(OptVzmin))
        if haspower:
            ### Combine Vzmin's from radio, optical and z selections
            #Vzmin = Vzlim_low*np.ones(len(OptVzmax))
            #t1 = np.maximum(Vzmin, PVzmin)
            #t2 = np.maximum(t1, OptVzmin)
            #Vzmin = Vzlim_low*np.ones(len(OptVzmax))
            #t1 = np.maximum(Vzmin, PVzmin)
            t2 = np.maximum(PVzmin, OptVzmin)
            Vzmin = t2
            
        else:
            t1 = np.maximum(OptVzmin, Vzmin)
            Vzmin = t1
        
        if 'Vzmin' not in self.cat.dtype.names:
            self.cat.add_column(Column(Vzmin, 'Vzmin'))
        else:
            self.cat['Vzmin'] = Vzmin
            
        if 'Vzmax' not in self.cat.dtype.names:
            self.cat.add_column(Column(Vzmax, 'Vzmax'))
        else:
            self.cat['Vzmax'] = Vzmax
            
            
        if 'Fcor' not in self.cat.dtype.names:
            self.cat.add_column(Column(np.ones(len(Vzmax)), 'Fcor'))
        else:
            self.cat['Fcor'] = np.ones(len(Vzmax))
            
        if 'areal' not in self.cat.dtype.names:
            self.cat.add_column(Column(np.ones(len(Vzmax)), 'areal'))
        else:
            self.cat['areal'] = np.ones(len(Vzmax))

        
        if plot and haspower:
            self.plot_Vzmin_Vzmax()
        
        return
    
    
    def set_power_z_limit(self):
        zz = np.linspace(self.zlim_low, self.zlim_high, 10)
        Plim = np.log10(LF_util.RadioPower(self.radio_fluxlim_faint, zz, alpha=-0.7)) 
        self.power_z_limit = np.vstack((zz, Plim))
        return
    
    
    
    def compute_LF(self, grid, maskbins=None, CV_f=None, ignoreMinPower=False, catcol='power'):
        
        
        logp_radio_lf = (grid[:-1]+grid[1:])/2.
        dlogp_radio_lf = (grid[1:]-grid[:-1])/2.
        
        rho, rhoerr, num = LF_util.get_LF_rms_f_areal(grid, self.cat[catcol], self.cat['Vzmin'], self.cat['Vzmax'], self.cat['Fcor'], self.cat['areal'], ignoreMinPower=ignoreMinPower)
        
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
                    print("maskbins invalid")
            else:
                print("maskbins invalid")
            
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
        print("{n1} of {n2} sources selected on SM envelope".format(n2=len(self.cat),n1= len(smenv_ind)))
        sm_complete_sample = self.sub_sample_ind('m', smenv_ind )
        
        # zmax determined from mass-completeness also #
        SMzmax = masscompleteness(m=sm_complete_sample.cat['smass'])
        VSMzmax = acosmo.comoving_volume(SMzmax).value
        #sm_complete_sample.cat['Vzmax']=
        newVSMzmax = np.minimum(sm_complete_sample.cat['Vzmax'], VSMzmax)
        
        #import ipdb
        #ipdb.set_trace()
    
        
        #print self.cat['power']
        #rho, rhoerr, num = LF_util.get_LF_rms_f_areal(smgrid, sm_complete_sample.cat['smass'], sm_complete_sample.cat['Vzmin'], sm_complete_sample.cat['Vzmax'], sm_complete_sample.cat['Fcor'], sm_complete_sample.cat['areal'], sm_complete_sample.area, xstr="SM")
        rho, rhoerr, num = LF_util.get_LF_rms_f_areal(smgrid, sm_complete_sample.cat['smass'], sm_complete_sample.cat['Vzmin'], newVSMzmax, sm_complete_sample.cat['Fcor'], sm_complete_sample.cat['areal'], sm_complete_sample.domega, xstr="SM")
        
        
        #print self.cat['power']
        #print grid
        #print type(grid)
        
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
        
        rho, rhoerr, num = LF_util.get_CLF_f_areal(pgrid, self.cat['power'], self.cat['zmin'], self.cat['zmax'], self.cat['Fcor'], self.cat['areal'], self.domega)
        
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
        
        
        rho, rhoerr, num = LF_util.get_rho_Plim_f_areal(Plimit, self.cat['power'], self.cat['zmin'], self.cat['zmax'], self.cat['Fcor'], self.cat['areal'], self.domega, xstr="SM")
        
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
    
    
        
    
