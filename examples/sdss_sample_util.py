import LF_util
import numpy as np
from astropy.table import Table, MaskedColumn



maglim_bright = 14.5
maglim_faint = 17.77
fmaglim_bright =  3631.*10**(-0.4*maglim_bright)  # AB mag
fmaglim_faint =  3631.*10**(-0.4*maglim_faint)    # AB mag

radiofluxlim = 5e-3               ## in Jy


CATPATH = './'


def SMenvelope(z=None, m=None):
    if z is not None:
        #zenv,lgMenv = np.loadtxt("Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
        zenv = np.linspace(0,0.7, 100)
        lgMenv = np.linspace(9, 13., 100)
        #lgMenv = np.linspace(10, 14.5, 100)
        #lgMenv = np.linspace(10, 13, 100)
        #lgMenv = np.linspace(0, 0, 100)
        zenv = np.append(zenv, 100.)
        lgMenv = np.append(lgMenv, 10.991006)
        return np.interp(z, zenv, lgMenv)

    elif m is not None:
        #zenv,lgMenv = np.loadtxt("Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
        zenv = np.linspace(0,0.7, 100)
        lgMenv = np.linspace(9, 13., 100)
        #lgMenv = np.linspace(10, 14.5, 100)
        #lgMenv = np.linspace(0, 0, 100)
        zenv = np.append(zenv, 100.)
        lgMenv = np.append(lgMenv, 10.991006)
        return np.interp(m, lgMenv, zenv)
    else:
        return

def SMenvelope_sdss(z):
    #zenv,lgMenv = np.loadtxt("Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
    zenv = np.linspace(0,0.7, 100)
    lgMenv = np.linspace(9, 13., 100)
    #lgMenv = np.linspace(10, 14.5, 100)
    #lgMenv = np.linspace(10, 13, 100)
    #lgMenv = np.linspace(0, 0, 100)
    zenv = np.append(zenv, 100.)
    lgMenv = np.append(lgMenv, 10.991006)
    return np.interp(z, zenv, lgMenv)


def SMenvelope_sdss_z(M):
    #zenv,lgMenv = np.loadtxt("Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
    zenv = np.linspace(0,0.7, 100)
    lgMenv = np.linspace(9, 13., 100)
    #lgMenv = np.linspace(10, 14.5, 100)
    #lgMenv = np.linspace(0, 0, 100)
    zenv = np.append(zenv, 100.)
    lgMenv = np.append(lgMenv, 10.991006)
    return np.interp(M, lgMenv, zenv)





def select_good_sample(good=True):
    
    global maglim_bright
    global maglim_faint
    # cutout a block
    ralim1_full = 140.
    ralim2_full = 220.
    declim1_full = 0.
    declim2_full = 60.
    
    tt = np.load(CATPATH+'projects/sdss/non_dup_ind.npz')
    non_dup_ind = tt['ind']

    print "read full sdss catalogue"
    cat_sdss = Table.read(CATPATH+'projects/sdss/gal_fnal_dr7_v5_2.fit')
    #cat_sdss = Table.read(PATH+'projects/sdss/gal_fnal_dr7_v5_2.copy.fit')
    # cut out the duplicates
    cat_sdss = cat_sdss[non_dup_ind]
    
    # cut out a given block with known area
    selind_block = (cat_sdss.RA >= ralim1_full)*(cat_sdss.RA <= ralim2_full)*(cat_sdss.DEC >= declim1_full)*(cat_sdss.DEC <= declim2_full)
    cat_sdss = cat_sdss[selind_block]

    # make a copy (to be the full sample)
    #cat_sdss1 = cat_sdss.copy()


    # main_samp =1 included in main sample
    #cat_sdss = cat_sdss[np.where(cat_sdss.main_samp==1)]

    #### MAGNITUDE LIMITS ###
    #cat_sdss = cat_sdss[np.where((cat_sdss.mag_r >= maglim_bright) & (cat_sdss.mag_r <= maglim_faint))]   # magnitude cut


    # The mass file
    print "read full sdss sm catalogue"
    cat_sdss_sm = Table.read(CATPATH+'projects/sdss/totlgm_dr7_v5_2b.fit')
    #cat_sdss_sm = Table.read(PATH+'projects/sdss/totlgm_dr7_v5_2.fit')
    cat_sdss_sm = cat_sdss_sm[non_dup_ind]
    cat_sdss_sm = cat_sdss_sm[selind_block]

    # The sfr file
    print "read full sdss sfr catalogue"
    cat_sdss_sfr = Table.read(CATPATH+'projects/sdss/gal_totsfr_dr7_v5_2.fits')
    #cat_sdss_sfr = Table.read(PATH+'projects/sdss/totlgm_dr7_v5_2.fit')
    cat_sdss_sfr = cat_sdss_sfr[non_dup_ind]
    cat_sdss_sfr = cat_sdss_sfr[selind_block]
    
    # make a copy (to be the full sample)
    #cat_sdss_sm1 = cat_sdss_sm.copy()

    # the redshift file
    print "read full sdss z catalogue"
    cat_sdss_z = Table.read(CATPATH+'projects/sdss/gal_info_dr7_v5_2.fit')
    cat_sdss_z = cat_sdss_z[non_dup_ind]
    cat_sdss_z = cat_sdss_z[selind_block]

    # make a copy (to be the full sample)
    #cat_sdss_z1 = cat_sdss_z.copy()
    # 927552 - VASC - http://www.mpa-garching.mpg.de/SDSS/

    # rad_agn = 1 for AGN, 0 for SF
    # lerg = 1 for lerg
    # herg = 1 for herg

    ## main_samp =1 included in main sample
    #cat_sdss = cat_sdss[np.where(cat_sdss.main_samp==1)]

    #import sys
    #sys.exit()

    if good:
        # main_samp =1 included in main sample
        #cat_sdss = cat_sdss[np.where(cat_sdss.main_samp==1)]

        #### MAGNITUDE LIMITS ###
        ind = np.where((cat_sdss.MAG[:,2] >= maglim_bright) & (cat_sdss.MAG[:,2] <= maglim_faint))
        
    
    
        cat_sdss = cat_sdss[ind]   # magnitude cut
        cat_sdss_sm = cat_sdss_sm[ind] 
        cat_sdss_sfr = cat_sdss_sfr[ind]
        cat_sdss_z = cat_sdss_z[ind] 
        
        
        # select z !=0 
        ind = np.where(cat_sdss_z.Z > 0)
        
        cat_sdss = cat_sdss[ind]   # magnitude cut
        cat_sdss_sm = cat_sdss_sm[ind] 
        cat_sdss_sfr = cat_sdss_sfr[ind]
        cat_sdss_z = cat_sdss_z[ind] 

    ##### REDSHIFT LIMITS ###
    #ind = np.where((cat_sdss_z.Z > zlim_low_sdss) & (cat_sdss_z.Z < zlim_high_sdss))

    zz = cat_sdss_z.Z
    smass = cat_sdss_sm.MEDIAN
    sfr = cat_sdss_sfr.MEDIAN

    # Magnitudes
    rmag = cat_sdss.MAG[:,2]  
    frmag = 3631.*10.**(-0.4*rmag)    # AB mag
    rlum = LF_util.OpticalLuminosity(frmag, zz)


    lf_cat = np.core.records.fromarrays((zz, rmag, rlum, smass, sfr), names='z, opt_mag, opt_lum, smass, sfr', formats = 'f8, f8, f8, f8, f8')

    return cat_sdss, lf_cat

def select_good_radio_sample(good=True):
    
    global maglim_bright
    global maglim_faint
    global radiofluxlim
    
    cat_sdss1 = Table.read(CATPATH+'/sdss_dr7_radiosources_with_wise_matches_allsources_vasc_sdss.fits')
    cat_sdss = cat_sdss1.copy()
    # rad_agn = 1 for AGN, 0 for SF
    # lerg = 1 for lerg
    # herg = 1 for herg

    if good:
        # main_samp =1 included in main sample
        cat_sdss = cat_sdss[np.where(cat_sdss['main_samp']==1)]

        #### FLUX LIMIT ###
        cat_sdss = cat_sdss[np.where(cat_sdss['S_NVSS'] > radiofluxlim)]   # radio flux limit

        #### MAGNITUDE LIMITS ###
        cat_sdss = cat_sdss[np.where((cat_sdss['mag_r'] >= maglim_bright) & (cat_sdss['mag_r'] <= maglim_faint))]   # magnitude cut


    # Redshift
    zz = cat_sdss['z']

    # Magnitudes
    rmag = cat_sdss['mag_r']   
    frmag = 3631.*10.**(-0.4*rmag)    # AB mag
    rlum = LF_util.OpticalLuminosity(frmag, zz)

    # Radio flux
    ff = cat_sdss['S_NVSS']   #in Jy
    power = np.log10(LF_util.RadioPower(ff, zz, alpha=0.8))  # assumed spec ind

    smass = cat_sdss['Mass']
    sfr = cat_sdss['SFR']
    

    #lf_cat = np.core.records.fromarrays((zz, ff, power, rmag, rlum, smass, sfr, cat_sdss['rad_agn'], cat_sdss['herg'], cat_sdss['lerg']), names='z, radio_flux, power, opt_mag, opt_lum, smass, sfr, agn, herg, lerg', formats = 'f8, f8, f8, f8, f8, f8, f8, i8, i8, i8')
    
    cat_sdss.rename_column('rad_agn','agn')

    cat_sdss.add_column(MaskedColumn(ff ,'radio_flux'))
    cat_sdss.add_column(MaskedColumn(power ,'power'))
    cat_sdss.add_column(MaskedColumn(rmag ,'opt_mag'))
    cat_sdss.add_column(MaskedColumn(rlum ,'opt_lum'))
    cat_sdss.add_column(MaskedColumn(smass ,'smass'))
    cat_sdss.add_column(MaskedColumn(sfr ,'sfr'))

    
    keeplist = ['z', 'radio_flux', 'power', 'opt_mag', 'opt_lum', 'smass', 'sfr', 'agn', 'herg', 'lerg']
    cat_sdss.keep_columns(keeplist)

    return cat_sdss




areaSDSS = 2.17    ## 2.17 st calculated by Best and Heckman
areaSDSS_full = 2.4    ## area of full sdss sample
areaSDSS_full = 1.2    ## test
areaSDSS_full = 0.6    ## test
areaSDSS_full = 4.8    ## test


ralim1_full = 140.
ralim2_full = 220.
declim1_full = 0.
declim2_full = 60.

#zlim2_low_sdss = 0.01
#zlim2_high_sdss = 0.3
zlim_low_sdss = 0.01
zlim_high_sdss = 0.3
