import sdss_sample_util
import lf_sample
import LF_util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# load sdss catalogue
# load sdss #
lfcat_nvsssdss_all = sdss_sample_util.select_good_radio_sample()

# initialise lf_sample instance
sdss_radio_sample = lf_sample.lf_sample("NVSS_SDSS", lfcat_nvsssdss_all, zlow=0.01, zhigh=0.3, radio_fluxlim_faint = sdss_sample_util.radiofluxlim, opt_fluxlim_faint = sdss_sample_util.fmaglim_faint, opt_fluxlim_bright=sdss_sample_util.fmaglim_bright, area=sdss_sample_util.areaSDSS)

# calc zmin/zmax
sdss_radio_sample.calc_zmin_zmax()

# subselection on redshift
low_z_bins =  np.array([[0.01,0.3]])
low_z_sample = sdss_radio_sample.sub_z_sample('zbin01', low_z_bins[0][0], low_z_bins[0][1])

# subselection on a different field in the lf_sample array
# hergs/lergs for low z
low_z_sf_sample = low_z_sample.sub_sample_by_field('sf','agn',-0.1,0.1)
low_z_agn_sample = low_z_sample.sub_sample_by_field('agn','agn',0.9,1.1)
low_z_lerg_sample = low_z_sample.sub_sample_by_field('lerg','lerg',0.9,1.1)
low_z_herg_sample = low_z_sample.sub_sample_by_field('herg','herg',0.9,1.1)


# compute LF
dp = 0.3
pgrid_radio_lf = np.arange(22.,26.61,dp)
low_z_sample.compute_LF(pgrid_radio_lf)
low_z_sf_sample.compute_LF(pgrid_radio_lf)
low_z_agn_sample.compute_LF(pgrid_radio_lf)
low_z_lerg_sample.compute_LF(pgrid_radio_lf)
low_z_herg_sample.compute_LF(pgrid_radio_lf)



## PLOT ##

#f, ax = pp.paper_single_ax()
mpl.rc('text', usetex=True)


f = plt.figure()
ax = f.add_subplot(111)
plt.minorticks_on()

ax.set_yscale('log')
ax.xaxis.set_major_locator( plt.MaxNLocator(nbins=6, prune='lower') )
ax.xaxis.set_minor_locator( mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_major_locator( mpl.ticker.LogLocator(numticks=5))

ee = []
l = ax.errorbar(low_z_sample.LF_x, low_z_sample.LF_rho, [low_z_sample.LF_rhoerrup, low_z_sample.LF_rhoerrlow], low_z_sample.LF_xerr, label='All' )
ee.append(l)
l = ax.errorbar(low_z_sf_sample.LF_x, low_z_sf_sample.LF_rho, [low_z_sf_sample.LF_rhoerrup, low_z_sf_sample.LF_rhoerrlow], low_z_sf_sample.LF_xerr, label='SF' )
ee.append(l)
l = ax.errorbar(low_z_agn_sample.LF_x, low_z_agn_sample.LF_rho, [low_z_agn_sample.LF_rhoerrup, low_z_agn_sample.LF_rhoerrlow], low_z_agn_sample.LF_xerr, label='AGN' )
ee.append(l)

ee2 = []
x,y, yerr =  LF_util.get_BH(ttype='all', f=1400.)
l = ax.errorbar(x,y,yerr, label='BH ALL') 
ee2.append(l)
x,y, yerr =  LF_util.get_BH(ttype='sf', f=1400.)
l = ax.errorbar(x,y,yerr, label='BH SF') 
ee2.append(l)
x,y, yerr =  LF_util.get_BH(ttype='agn', f=1400.)
l = ax.errorbar(x,y,yerr, label='BH AGN') 
ee2.append(l)


ax.set_ylabel(r'$\rho \,[\mathrm{Mpc}^{-3}~\log P ^{-1}]$')
    
ax.set_xlabel(r'$\log P \, [\log \mathrm{W~Hz}^{-1}]$')
ax.set_ylim(10**-8.5, 10**-2.5)
ax.set_xlim(22.0, 26.5)

l1 = ax.legend(handles=ee, loc='lower left')
ax.add_artist(l1)
ax.legend(handles=ee2, loc='upper right')
plt.minorticks_on()
plt.savefig('RLF_local_AGN_SF')




f = plt.figure()
ax = f.add_subplot(111)
plt.minorticks_on()

ax.set_yscale('log')
ax.xaxis.set_major_locator( plt.MaxNLocator(nbins=6, prune='lower') )
ax.xaxis.set_minor_locator( mpl.ticker.MultipleLocator(0.25))
ax.yaxis.set_major_locator( mpl.ticker.LogLocator(numticks=5))

ee = []
l = ax.errorbar(low_z_agn_sample.LF_x, low_z_agn_sample.LF_rho, [low_z_agn_sample.LF_rhoerrup, low_z_agn_sample.LF_rhoerrlow], low_z_agn_sample.LF_xerr, label='AGN' )
ee.append(l)
l = ax.errorbar(low_z_herg_sample.LF_x, low_z_herg_sample.LF_rho, [low_z_herg_sample.LF_rhoerrup, low_z_herg_sample.LF_rhoerrlow], low_z_herg_sample.LF_xerr, label='HERG' )
ee.append(l)
l = ax.errorbar(low_z_lerg_sample.LF_x, low_z_lerg_sample.LF_rho, [low_z_lerg_sample.LF_rhoerrup, low_z_lerg_sample.LF_rhoerrlow], low_z_lerg_sample.LF_xerr, label='LERG' )
ee.append(l)

ee2 = []
x,y, yerr =  LF_util.get_BH(ttype='agn', f=1400.)
l = ax.errorbar(x,y,yerr, label='BH AGN') 
ee2.append(l)
x,y, yerr =  LF_util.get_BH(ttype='herg', f=1400.)
l = ax.errorbar(x,y,yerr, label='BH HERG') 
ee2.append(l)
x,y, yerr =  LF_util.get_BH(ttype='lerg', f=1400.)
l = ax.errorbar(x,y,yerr, label='BH LERG') 
ee2.append(l)


ax.set_ylabel(r'$\rho \,[\mathrm{Mpc}^{-3}~\log P ^{-1}]$')
    
ax.set_xlabel(r'$\log P \, [\log \mathrm{W~Hz}^{-1}]$')
ax.set_ylim(10**-8.5, 10**-2.5)
ax.set_xlim(22.0, 26.5)

l1 = ax.legend(handles=ee, loc='lower left')
ax.add_artist(l1)
ax.legend(handles=ee2, loc='upper right')
plt.minorticks_on()
plt.savefig('RLF_local_HERG_LERG')


#pp.fig_save_many(f, 'RLF_local_HERG_LERG')
