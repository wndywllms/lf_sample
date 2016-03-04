



def select_good_radio_sample(maglim_bright, fmaglim_faint, radiofluxlim, good=True):
    '''
    read sdss catalogue and output a lf_sample array
    '''
    
    CATPATH = '/home/wwilliams/'

    
    
    cat_sdss1 = load_fits(CATPATH+'projects/sdss/sdss_dr7_radiosources_with_wise_matches_allsources_vasc_sdss.fits')
    cat_sdss = cat_sdss1.copy()
    # rad_agn = 1 for AGN, 0 for SF
    # lerg = 1 for lerg
    # herg = 1 for herg

    if good:
        # main_samp =1 included in main sample
        cat_sdss = cat_sdss[np.where(cat_sdss.main_samp==1)]

        #### FLUX LIMIT ###
        cat_sdss = cat_sdss[np.where(cat_sdss.S_NVSS > radiofluxlim)]   # radio flux limit

        #### MAGNITUDE LIMITS ###
        cat_sdss = cat_sdss[np.where((cat_sdss.mag_r >= maglim_bright) & (cat_sdss.mag_r <= maglim_faint))]   # magnitude cut


    # Redshift
    zz = cat_sdss.z

    # Magnitudes
    rmag = cat_sdss.mag_r   
    frmag = 3631.*10.**(-0.4*rmag)    # AB mag
    rlum = LF_util.OpticalLuminosity(frmag, zz)

    # Radio flux
    ff = cat_sdss.S_NVSS   #in Jy
    power = np.log10(LF_util.RadioPower(ff, zz, alpha=0.8))  # assumed spec ind

    smass = cat_sdss.Mass
    sfr = cat_sdss.SFR
    

    lf_cat = np.core.records.fromarrays((zz, ff, power, rmag, rlum, smass, sfr, cat_sdss.rad_agn, cat_sdss.herg, cat_sdss.lerg), names='z, radio_flux, power, opt_mag, opt_lum, smass, sfr, agn, herg, lerg', formats = 'f8, f8, f8, f8, f8, f8, f8, i8, i8, i8')
    

    return cat_sdss, lf_cat


areaSDSS = 2.17    ## 2.17 st calculated by Best and Heckman
maglim_bright = 14.5
maglim_faint = 17.77
fmaglim_bright =  3631.*10**(-0.4*maglim_bright)  # AB mag
fmaglim_faint =  3631.*10**(-0.4*maglim_faint)    # AB mag

radiofluxlim = 5e-3               ## in Jy

# load sdss catalogue
# load sdss #
cat_nvsssdss_all, lfcat_nvsssdss_all = sdss_sample_util.select_good_radio_sample()

# initialise lf_sample instance
sdss_radio_sample = lf_sample.lf_sample("NVSS_SDSS", lfcat_nvsssdss_all, zlow=0.01, zhigh=0.3, radio_fluxlim_faint = sdss_sample_util.radiofluxlim, opt_fluxlim_faint = sdss_sample_util.fmaglim_faint, opt_fluxlim_bright=sdss_sample_util.fmaglim_bright, area=sdss_sample_util.areaSDSS)

# calc zmin/zmax
sdss_radio_sample.calc_zmin_zmax()

# subselection on redshift
low_z_bins =  np.array([[0.01,0.3]])
low_z_sample = sdss_radio_sample.sub_z_sample('zbin01', low_z_bins[0][0], low_z_bins[0][1])
# compute LF
dp = 0.5
pgrid_radio_lf = np.arange(22.,26.01,dp)
low_z_sample.compute_LF(pgrid_radio_lf)

# subselection on a different field in the lf_sample array
# hergs/lergs for low z
low_z_herg_sample = low_z_sample.sub_sample_by_field('herg','herg',0.9,1.1)
low_z_lerg_sample = low_z_sample.sub_sample_by_field('lerg','lerg',0.9,1.1)
low_z_herg_sample.compute_LF(pgrid_radio_lf)
low_z_lerg_sample.compute_LF(pgrid_radio_lf)





## PLOT ##

f, ax = pp.paper_single_ax()
pl.minorticks_on()

ax.set_yscale('log')
ax.xaxis.set_major_locator( pl.MaxNLocator(nbins=6, prune='lower') )
ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator( mpl.ticker.LogLocator(numticks=5))


ax.errorbar(low_z_sample.LF_x, low_z_sample.LF_rho, [low_z_sample.LF_rhoerrup, low_z_sample.LF_rhoerrlow], low_z_sample.LF_xerr, label='All' )
ax.errorbar(low_z_herg_sample.LF_x, low_z_herg_sample.LF_rho, [low_z_herg_sample.LF_rhoerrup, low_z_herg_sample.LF_rhoerrlow], low_z_herg_sample.LF_xerr, label='HERG' )
ax.errorbar(low_z_lerg_sample.LF_x, low_z_lerg_sample.LF_rho, [low_z_lerg_sample.LF_rhoerrup, low_z_lerg_sample.LF_rhoerrlow], low_z_lerg_sample.LF_xerr, label='LERG' )

ax.set_ylabel(r'$\rho \,[\mathrm{Mpc}^{-3}~\log P ^{-1}]$')
    
ax.set_xlabel(r'$\log P \, [\log \mathrm{W~Hz}^{-1}]$')
ax.set_ylim(10**-8.5, 10**-2.5)
ax.set_xlim(22.0, 26.5)

ax.legend(loc='upper left')
pl.minorticks_on()
pp.fig_save_many(f, 'RLF_local_HERG_LERG')
