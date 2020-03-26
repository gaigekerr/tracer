#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Conduct data analysis and plot quantities related to GEOSChem Rn–Pb–Be 
tracer simulations 

Revision History
    15102019 -- initial version created
    27022020 -- data analysis to examine the zonal mean jet/tracer relationship
"""

def nhmap(lat, lng, field, title, cbar_label, clevs, cmap, fstr, 
    oceanon=False, extend='both', hatch=None, ebar=None): 
    """Plot map (Robinson projection) of desired field. 
    
    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, n.b., the first or last longitude value might need to be 
        fudged to wrap contours around lng = 0 deg, [lng,]          
    field : numpy.ndarray
        Desired field/quantity, [lat, lng]
    title : str
        Title for plot        
    cbar_label : str
        Label for the colorbar (field and units)
    clevs : numpy.ndarray
        Filled contour levels values
    cmap : str
        Colormap name
    fstr : str
        Output filename suffix, should specify field, levels, and time period
    oceanon : bool, optional 
        If True, map adds ocean shapefiles; delfault is False
    extend : str, optional
        Extend settings for matplotlib colormap/colorbar; i.e., 'neither', 
        'both', 'min', 'max'; default is 'both'
    hatch : numpy.ndarray, optional
        Hatching to juxtapose on filled contour; grid cells equal to 1 are 
        hatched and values equal to NaN are left empty, [lat, lng]
    ebar : numpy.ndarray, option
        Errorbars to indicate mean location and variability of the eddy-driven
        jet, units of degrees, [time, lng]

    Returns
    -------
    None              
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    fig = plt.figure(figsize=(9,3.5))
    ax=plt.subplot2grid((1,2), (0,0), colspan=2,
                        projection=ccrs.PlateCarree(central_longitude=0.))
    ax.set_title(title, fontsize=16, x=0.02, ha='left')
    if oceanon == 'yes':
        ax.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='lightgrey')
    ax.coastlines(lw=0.25, color='k', zorder=3)
    # Plot filled contours
    cmap = plt.get_cmap(cmap)
    # For wrapping around the Prime Meridian
    if lng[-1] != 360:
        lng[-1] = 360.
    mb = ax.contourf(lng, lat, field, clevs, cmap=cmap, extend=extend,
        transform=ccrs.PlateCarree(), zorder=1)
    # If specified, add hatching (e.g., for significance)
    if hatch is None: pass
    else:
        # Hatching for significance
        ax.contourf(lng, lat, hatch, hatches=['//////'], colors='none', 
            transform=ccrs.PlateCarree())
    # If specified add errorbars (e.g., for location of eddy-driven jet)
    if ebar is None: pass
    else:
        # Plot only every X number of longitude values
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(ebar, axis=0)[::skiplng], 
            yerr=np.nanstd(ebar, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
    # Label parallels and meridians
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax.xaxis.set_major_formatter(lng_formatter)      
    ax.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)  
    ax.set_extent([lng.min()-180., lng.max()-180., 0., 90.])  
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.08, right=0.82)
    colorbar_axes = plt.gcf().add_axes([0.85, ax.get_position().y0, 
        0.02, (ax.get_position().y1-ax.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=clevs, extend=extend)
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=16)
    ax.outline_patch.set_zorder(20)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
        'nhmap_%s.png'%(fstr), dpi = 350)
    return

import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate
sys.path.append('/Users/ghkerr/phd/tracer/')
import tracer_open, tracer_calculate
years = [2008, 2009, 2010]
months = ['jun', 'jul', 'aug']
latmin = 0.
latmax = 90.
lngmin = 0.
lngmax = 360.
pmin = 800.
pmax = 1005.
# import matplotlib.pyplot as plt
# # Load GEOSChem tracers
# TRAC_0_10, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_0_10', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_10_20, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_10_20', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_20_30, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_20_30', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_30_40, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_30_40', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_40_50, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_40_50', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_50_60, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_50_60', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_60_70, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_60_70', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_70_80, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_70_80', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_80_90, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, months, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_80_90', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# # Load MERRA-2 column meridional wind with the same vertical levels as 
# # GEOSChem 
# Vcolumn, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, months, 'V',
#     lngmin, latmax, lngmax, latmin, pmin, pmax)
# Vcolumn = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, Vcolumn)
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
# edj, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, months, 'U', 
#     lngmin, latmax, lngmax, latmin, 487., 526., operation='mean')
# # Degrade to resolution of GEOSChem    
# edj = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, lat_edj, 
#     lng_edj, edj)
# # Subset fields in mid-latitudes
# edj, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj, lat_gc, lng_gc, 
#     0., 360., 20., 70.)
# # Determine jet latitude
# edj = tracer_calculate.find_jetlatitude(edj, lat_edj, lng_edj)
# # Find distance from the jet where positive distances imply that the jet is
# # north of a particular location
# edj_dist = tracer_calculate.determine_jetdist(edj, lat_gc, lng_gc)
# # Load GEOSCHem PM2.5 and O3 
# O3, lat_gc_c, lng_gc_c, lev_gc_c = tracer_open.open_geoschem([2008], months, 
#     'merra2_4x5_tropchem', 'SpeciesConc', 'SpeciesConc_O3', latmin, latmax, 
#     lngmin, lngmax, 985., 1000.)
# PM25, lat_gc_c, lng_gc_c, lev_gc_c = tracer_open.open_geoschem([2008], months, 
#     'merra2_4x5_tropchem', 'AerosolMass', 'PM25', latmin, latmax, lngmin, 
#     lngmax, 985., 1000.)
# O3 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_gc_c, lng_gc_c, O3[:,0]*1e9)
# PM25 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_gc_c, lng_gc_c, PM25[:,0])

"""O3, NOX, AND CO-JET RELATIONSHIPS"""
# # Load Northern Hemisphere HindcastMR2 GMI CTM O3, CO, and NOx
# lat_gmi, lng_gmi, times_gmi, o3_gmi = \
#     globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
#     ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'O3', 'HindcastMR2')
# o3_gmi = o3_gmi*1e9
# lat_gmi, lng_gmi, times_gmi, co_gmi = \
#     globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
#     ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'CO', 'HindcastMR2')
# co_gmi = co_gmi*1e9
# lat_gmi, lng_gmi, times_gmi, no_gmi = \
#     globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
#     ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'NO', 'HindcastMR2')
# no_gmi = no_gmi*1e9
# lat_gmi, lng_gmi, times_gmi, no2_gmi = \
#     globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
#     ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'NO2', 'HindcastMR2')
# no2_gmi = no2_gmi*1e9
# nox_gmi = no_gmi+no2_gmi
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
# edj, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain([2008, 2009, 2010], 
#     ['jun', 'jul', 'aug'], 'U', 0., 90., 360., -1., 487., 526., 
#     operation='mean')
# # Degrade to resolution of GEOSChem    
# edj = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
#     lat_edj, lng_edj, edj)
# # Subset in mid-latitudes
# edj, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj, lat_gmi, 
#     lng_gmi, 0., 360., 20., 70.)
# # Determine jet latitude
# edj = tracer_calculate.find_jetlatitude(edj, lat_edj, lng_edj)
# edj_dist = tracer_calculate.determine_jetdist(edj, lat_gmi, lng_gmi)
# labels = ['nox', 'o3', 'co']
# titles = ['NO$_{x,\:PW}$ - NO$_{x,\:EW}$', 'O$_{3,\:PW}$ - O$_{3,\:EW}$', 
#     'CO$_{PW}$ - CO$_{EW}$']
# clevs = [np.linspace(-0.5, 0.5, 9), np.linspace(-8, 8, 9), 
#     np.linspace()]
# for i, species in enumerate([nox_gmi, o3_gmi, co_gmi]):
#     r_speciesjet = globalo3_calculate.calculate_r(species, edj_dist, lat_gmi, 
#         lng_gmi)
#     significance_r_speciesjet = \
#         globalo3_calculate.calculate_r_significance(species, edj_dist, 
#         r_speciesjet, lat_gmi, lng_gmi)
#     (eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_species, 
#         eqjet_species) = globalo3_calculate.segregate_field_bylat(species, 
#         lng_gmi, edj, np.arange(0, len(species), 1))
#     nhmap(lat_gmi, 
#         lng_gmi, 
#         pwjet_species-eqjet_species,
#         titles[i],
#         '[ppbv]', 
#         np.linspace(-30, 30, 11),
#         'coolwarm', 
#         'pwjet-eqjet_%s_gmi_sfc'%labels[i],
#         extend='both', 
#         hatch=significance_r_speciesjet,
#         ebar=edj)

"""O3-JET RELATIONSHIPS (MAPS AND ZONAL-MEAN FLUX) ON DAILY, MONTHLY, AND 
   SEASONAL TIMESCALES"""
# import pandas as pd
# import sys
# sys.path.append('/Users/ghkerr/phd/globalo3/')
# import globalo3_open
# months_str = ['jun', 'jul', 'aug']
# years = np.arange(2000, 2011, 1)
# # Load Northern Hemisphere HindcastMR2 GMI CTM O3 
# lat_gmi, lng_gmi, times_gmi, o3_gmi = \
#     globalo3_open.open_overpass2_specifieddomain(years, 
#     ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'O3', 'HindcastMR2')
# o3_gmi = o3_gmi*1e9
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
# edj, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, 
#     ['jun', 'jul', 'aug'], 'U', 0., 90., 360., -1., 487., 526., 
#     operation='mean')
# # Degrade to resolution of GEOSChem    
# edj = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
#     lat_edj, lng_edj, edj)
# # Subset in mid-latitudes
# edj, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj, lat_gmi, 
#     lng_gmi, 0., 360., 20., 70.)
# # Determine jet latitude
# edj = tracer_calculate.find_jetlatitude(edj, lat_edj, lng_edj)
# edj_dist = tracer_calculate.determine_jetdist(edj, lat_gmi, lng_gmi)
# # Load MERRA-2 near-surface meridional wind at the lowest reanalysis level
# Vcolumn, lat_merra, lng_merra, lev_merra_full = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, ['jun', 
#     'jul', 'aug'], 'V', 0., 90., 360., -1., 800., 1005.)
# Vcolumn = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#     lng_gmi, lat_merra, lng_merra, Vcolumn)
# # Load MERRA-2 temperature at the surface 
# Tsfc, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years,
#     ['jun', 'jul', 'aug'], 'T', 0., 90., 360., -1., 985., 1005.)
# Tsfc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#     lng_gmi, lat_merra, lng_merra, Tsfc)
# # Find summer days
# o3_days = pd.date_range(start='2000-06-01', end='2010-08-31')
# where_jja = np.where((o3_days.month.isin([6,7,8])==True) & 
#                       (o3_days.year.isin(years)))[0]
# o3_days = o3_days[where_jja] 
# # Reshape into on monthly, seasonal chunks
# o3_yearly, o3_monthly = [], []
# edj_dist_yearly, edj_dist_monthly = [], []
# edj_yearly, edj_monthly = [], []
# Vcolumn_yearly, Vcolumn_monthly = [], []
# Tsfc_yearly, Tsfc_monthly = [], []
# for year in years:
#     where_yr = np.where(o3_days.year==year)[0]
#     # Calculate seasonal average for year of interest
#     o3_yearly.append(np.nanmean(o3_gmi[where_yr], axis=0))
#     edj_dist_yearly.append(np.nanmean(edj_dist[where_yr], axis=0))
#     edj_yearly.append(np.nanmean(edj[where_yr], axis=0))
#     Vcolumn_yearly.append(np.nanmean(Vcolumn[where_yr], axis=0))
#     Tsfc_yearly.append(np.nanmean(Tsfc[where_yr], axis=0))    
#     for month in np.arange(6, 9, 1):
#         where_mon = np.where((o3_days.year==year) & 
#                               (o3_days.month==month))[0]
#         o3_monthly.append(np.nanmean(o3_gmi[where_mon], axis=0))
#         edj_dist_monthly.append(np.nanmean(edj_dist[where_mon], axis=0))
#         edj_monthly.append(np.nanmean(edj[where_mon], axis=0))
#         Vcolumn_monthly.append(np.nanmean(Vcolumn[where_mon], axis=0))
#         Tsfc_monthly.append(np.nanmean(Tsfc[where_mon], axis=0))        
# # Daily, monthly, and seasonal jet-O3, jet-V, O3-T correlation and PW-EW composites
# # Daily jet-O3
# r_o3jet_daily = globalo3_calculate.calculate_r(o3_gmi, edj_dist, lat_gmi, 
#     lng_gmi)
# significance_r_o3jet_daily = \
#     globalo3_calculate.calculate_r_significance(o3_gmi, edj_dist, 
#     r_o3jet_daily, lat_gmi, lng_gmi)
# (eqjet_lat_daily, eqjet_lat_var_daily, pwjet_lat_daily,
# pwjet_lat_var_daily, pwjet_o3_daily, eqjet_o3_daily) = \
#     globalo3_calculate.segregate_field_bylat(o3_gmi, lng_gmi, edj, 
#     np.arange(0, len(o3_gmi), 1))
# nhmap(lat_gmi, 
#     lng_gmi, 
#     pwjet_o3_daily-eqjet_o3_daily,
#     'O$_{3,\:PW}$ - O$_{3,\:EW}$',
#     '[ppbv]', 
#     np.linspace(-8, 8, 9), 
#     'coolwarm', 
#     'pwjet-eqjet_dailyo3_gmi_2000-2010',
#     extend='both', 
#     hatch=significance_r_o3jet_daily,
#     ebar=edj)
# # Daily jet-V
# r_Vjet_daily = globalo3_calculate.calculate_r(np.nanmean(Vcolumn, axis=1), 
#     edj_dist, lat_gmi, lng_gmi)
# significance_r_Vjet_daily = \
#     globalo3_calculate.calculate_r_significance(np.nanmean(Vcolumn, axis=1),  
#     edj_dist, r_Vjet_daily, lat_gmi, lng_gmi)
# (eqjet_lat_daily, eqjet_lat_var_daily, pwjet_lat_daily,
# pwjet_lat_var_daily, pwjet_V_daily, eqjet_V_daily) = \
#     globalo3_calculate.segregate_field_bylat(np.nanmean(Vcolumn, axis=1),   
#     lng_gmi, edj, np.arange(0, len(Vcolumn), 1))
# nhmap(lat_gmi, 
#     lng_gmi, 
#     pwjet_V_daily-eqjet_V_daily,
#     'V$_{PW}$ - V$_{EW}$',
#     '[m s$^{-1}$]', 
#     np.linspace(-5, 5, 11), 
#     'coolwarm', 
#     'pwjet-eqjet_daily_v_mean1000-800hPa_2000-2010',
#     extend='both', 
#     hatch=significance_r_Vjet_daily,
#     ebar=edj)
# # Daily O3-T
# r_to3_daily = globalo3_calculate.calculate_r(np.nanmean(Tsfc, axis=1), 
#     o3_gmi, lat_gmi, lng_gmi)
# significance_r_to3_daily = \
#     globalo3_calculate.calculate_r_significance(np.nanmean(Tsfc, axis=1),  
#     o3_gmi, r_to3_daily, lat_gmi, lng_gmi)
# nhmap(lat_gmi, 
#     lng_gmi, 
#     r_to3_daily,
#     'r(T, O$_{3}$)',
#     '[$\cdot$]', 
#     np.linspace(-1., 1., 11), 
#     'coolwarm', 
#     'r_t_985hPa_o3_gmi_daily_2000-2010',
#     extend='neither', 
#     hatch=significance_r_to3_daily,
#     ebar=edj)      
# Monthly jet-O3
# r_o3jet_monthly = globalo3_calculate.calculate_r(np.array(o3_monthly),
#     np.array(edj_dist_monthly), lat_gmi, lng_gmi)
# significance_r_o3jet_monthly = \
#     globalo3_calculate.calculate_r_significance(np.array(o3_monthly),
#     np.array(edj_dist_monthly), r_o3jet_monthly, lat_gmi, lng_gmi)
# (eqjet_lat_monthly, eqjet_lat_var_monthly, pwjet_lat_monthly,
# pwjet_lat_var_monthly, pwjet_o3_monthly, eqjet_o3_monthly) = \
#     globalo3_calculate.segregate_field_bylat(np.array(o3_monthly), lng_gmi, 
#     np.array(edj_monthly), np.arange(0, len(o3_monthly), 1))
# nhmap(lat_gmi,
#     lng_gmi,
#     pwjet_o3_monthly-eqjet_o3_monthly,
#     'O$_{3,\:PW}$ - O$_{3,\:EW}$',
#     '[ppbv]',
#     np.linspace(-8, 8, 9),
#     'coolwarm',
#     'pwjet-eqjet_monthlyo3_gmi_2000-2010',
#     extend='both',
#     hatch=significance_r_o3jet_monthly,
#     ebar=np.array(edj_monthly))
# # Monthly jet-V
# r_Vjet_monthly = globalo3_calculate.calculate_r(
#     np.nanmean(Vcolumn_monthly, axis=1), np.array(edj_dist_monthly), lat_gmi, 
#     lng_gmi)
# significance_r_Vjet_monthly = \
#     globalo3_calculate.calculate_r_significance(
#     np.nanmean(Vcolumn_monthly, axis=1), np.array(edj_dist_monthly), 
#     r_Vjet_monthly, lat_gmi, lng_gmi)
# (eqjet_lat_monthly, eqjet_lat_var_monthly, pwjet_lat_monthly,
# pwjet_lat_var_monthly, pwjet_V_monthly, eqjet_V_monthly) = \
#     globalo3_calculate.segregate_field_bylat(
#     np.nanmean(Vcolumn_monthly, axis=1), lng_gmi, np.array(edj_monthly),
#     np.arange(0, len(Vcolumn_monthly), 1))
# nhmap(lat_gmi, 
#     lng_gmi, 
#     pwjet_V_monthly-eqjet_V_monthly,
#     'V$_{PW}$ - V$_{EW}$',
#     '[m s$^{-1}$]', 
#     np.linspace(-5, 5, 11), 
#     'coolwarm', 
#     'pwjet-eqjet_monthly_v_mean1000-800hPa_2000-2010',
#     extend='both', 
#     hatch=significance_r_Vjet_monthly,
#     ebar=np.array(edj_monthly))
# # Monthly O3-T
# r_to3_monthly = globalo3_calculate.calculate_r(
#     np.nanmean(Tsfc_monthly, axis=1), np.array(o3_monthly), lat_gmi, lng_gmi)
# significance_r_to3_monthly = \
#     globalo3_calculate.calculate_r_significance(
#     np.nanmean(Tsfc_monthly, axis=1),  np.array(o3_monthly), r_to3_monthly, 
#     lat_gmi, lng_gmi)
# nhmap(lat_gmi, 
#     lng_gmi, 
#     r_to3_monthly,
#     'r(T, O$_{3}$)',
#     '[$\cdot$]', 
#     np.linspace(-1., 1., 11), 
#     'coolwarm', 
#     'r_t_985hPa_o3_gmi_monthly_2000-2010',
#     extend='neither', 
#     hatch=significance_r_to3_monthly,
#     ebar=np.array(edj_monthly))
# Yearly jet-O3
# r_o3jet_yearly = globalo3_calculate.calculate_r(np.array(o3_yearly), 
#     np.array(edj_dist_yearly), lat_gmi, lng_gmi)
# significance_r_o3jet_yearly = \
#     globalo3_calculate.calculate_r_significance(np.array(o3_yearly), 
#     np.array(edj_dist_yearly), r_o3jet_yearly, lat_gmi, lng_gmi)
# (eqjet_lat_yearly, eqjet_lat_var_yearly, pwjet_lat_yearly,
# pwjet_lat_var_yearly, pwjet_o3_yearly, eqjet_o3_yearly) = \
#     globalo3_calculate.segregate_field_bylat(np.array(o3_yearly), lng_gmi, 
#     np.array(edj_yearly), np.arange(0, len(o3_yearly), 1))
# nhmap(lat_gmi, 
#     lng_gmi, 
#     pwjet_o3_yearly-eqjet_o3_yearly,
#     'O$_{3,\:PW}$ - O$_{3,\:EW}$',
#     '[ppbv]', 
#     np.linspace(-8, 8, 9), 
#     'coolwarm', 
#     'pwjet-eqjet_yearlyo3_gmi_2000-2010',
#     extend='both', 
#     hatch=significance_r_o3jet_yearly,
#     ebar=np.array(edj_yearly))
# # Yearly jet-V
# r_Vjet_yearly = globalo3_calculate.calculate_r(
#     np.nanmean(Vcolumn_yearly, axis=1), np.array(edj_dist_yearly), lat_gmi, 
#     lng_gmi)
# significance_r_Vjet_yearly = \
#     globalo3_calculate.calculate_r_significance(
#     np.nanmean(Vcolumn_yearly, axis=1), np.array(edj_dist_yearly), 
#     r_Vjet_yearly, lat_gmi, lng_gmi)
# (eqjet_lat_yearly, eqjet_lat_var_yearly, pwjet_lat_yearly,
# pwjet_lat_var_yearly, pwjet_V_yearly, eqjet_V_yearly) = \
#     globalo3_calculate.segregate_field_bylat(
#     np.nanmean(Vcolumn_yearly, axis=1), lng_gmi, np.array(edj_yearly),
#     np.arange(0, len(Vcolumn_yearly), 1))
# nhmap(lat_gmi, 
#     lng_gmi, 
#     pwjet_V_yearly-eqjet_V_yearly,
#     'V$_{PW}$ - V$_{EW}$',
#     '[m s$^{-1}$]', 
#     np.linspace(-5, 5, 11), 
#     'coolwarm', 
#     'pwjet-eqjet_yearly_v_mean1000-800hPa_2000-2010',
#     extend='both', 
#     hatch=significance_r_Vjet_yearly,
#     ebar=np.array(edj_yearly))
# # Yearly O3-T
# r_to3_yearly = globalo3_calculate.calculate_r(
#     np.nanmean(Tsfc_yearly, axis=1), np.array(o3_yearly), lat_gmi, lng_gmi)
# significance_r_to3_yearly = \
#     globalo3_calculate.calculate_r_significance(
#     np.nanmean(Tsfc_yearly, axis=1),  np.array(o3_yearly), r_to3_yearly, 
#     lat_gmi, lng_gmi)
# nhmap(lat_gmi, 
#     lng_gmi, 
#     r_to3_yearly,
#     'r(T, O$_{3}$)',
#     '[$\cdot$]', 
#     np.linspace(-1., 1., 11), 
#     'coolwarm', 
#     'r_t_985hPa_o3_gmi_yearly_2000-2010',
#     extend='neither', 
#     hatch=significance_r_to3_yearly,
#     ebar=np.array(edj_yearly))
# Calculate zonal-mean and eddy flux for all, PW, and EW days for different 
# averaging lengths (i.e., daily, monthly, and yearly)
# Daily 
# Vsfc_daily_eqjet, Vsfc_daily_pwjet, o3_gmi_daily_eqjet, o3_gmi_daily_pwjet = \
#     globalo3_calculate.sortfield_byjetlat_column(
#     np.reshape(Vsfc, (1012,1,92,288)), np.reshape(o3_gmi, (1012,1,92,288)),
#     edj, lng_gmi, lat_gmi, np.array([985.]), psize=30)
# o3_daily_mean, o3_daily_stationary, o3_daily_transient, o3_daily_total = \
#     globalo3_calculate.meridional_flux(Vsfc[:,0], o3_gmi, o3_days, lat_gmi, 
#     lng_gmi)
# o3_daily_mean_pwjet, o3_daily_stationary_pwjet, o3_daily_transient_pwjet, \
#     o3_daily_total_pwjet = globalo3_calculate.meridional_flux(
#     Vsfc_daily_pwjet[:,0], o3_gmi_daily_pwjet[:,0], 
#     o3_days[:len(Vsfc_daily_pwjet)], lat_gmi, lng_gmi)
# o3_daily_mean_eqjet, o3_daily_stationary_eqjet, o3_daily_transient_eqjet, \
#     o3_daily_total_eqjet = globalo3_calculate.meridional_flux(
#     Vsfc_daily_eqjet[:,0], o3_gmi_daily_eqjet[:,0], 
#     o3_days[:len(Vsfc_daily_eqjet)], lat_gmi, lng_gmi)
# # Monthly 
# Vsfc_monthly_eqjet, Vsfc_monthly_pwjet, o3_gmi_monthly_eqjet, \
#     o3_gmi_monthly_pwjet = globalo3_calculate.sortfield_byjetlat_column(
#     np.reshape(Vsfc_monthly, (33,1,92,288)), np.reshape(o3_monthly, (33,1,92,288)),
#     np.array(edj_monthly), lng_gmi, lat_gmi, np.array([985.]), psize=30)
# o3_monthly_mean, o3_monthly_stationary, o3_monthly_transient, \
#     o3_monthly_total = globalo3_calculate.meridional_flux(
#     np.array(Vsfc_monthly)[:,0], np.array(o3_monthly), 
#     np.arange(0, len(o3_monthly), 1), lat_gmi, lng_gmi)
# o3_monthly_mean_pwjet, o3_monthly_stationary_pwjet, o3_monthly_transient_pwjet, \
#     o3_monthly_total_pwjet = globalo3_calculate.meridional_flux(
#     np.array(Vsfc_monthly_pwjet)[:,0], np.array(o3_gmi_monthly_pwjet)[:,0], 
#     o3_days[:len(Vsfc_monthly_pwjet)], lat_gmi, lng_gmi)
# o3_monthly_mean_eqjet, o3_monthly_stationary_eqjet, o3_monthly_transient_eqjet, \
#     o3_monthly_total_eqjet = globalo3_calculate.meridional_flux(
#     np.array(Vsfc_monthly_eqjet)[:,0], np.array(o3_gmi_monthly_eqjet)[:,0], 
#     o3_days[:len(Vsfc_monthly_eqjet)], lat_gmi, lng_gmi)
# # Yearly 
# Vsfc_yearly_eqjet, Vsfc_yearly_pwjet, o3_gmi_yearly_eqjet, \
#     o3_gmi_yearly_pwjet = globalo3_calculate.sortfield_byjetlat_column(
#     np.reshape(Vsfc_yearly, (11,1,92,288)), np.reshape(o3_yearly, (11,1,92,288)),
#     np.array(edj_yearly), lng_gmi, lat_gmi, np.array([985.]), psize=30)
# o3_yearly_mean, o3_yearly_stationary, o3_yearly_transient, \
#     o3_yearly_total = globalo3_calculate.meridional_flux(
#     np.array(Vsfc_yearly)[:,0], np.array(o3_yearly), 
#     np.arange(0, len(o3_yearly), 1), lat_gmi, lng_gmi)
# o3_yearly_mean_pwjet, o3_yearly_stationary_pwjet, o3_yearly_transient_pwjet, \
#     o3_yearly_total_pwjet = globalo3_calculate.meridional_flux(
#     np.array(Vsfc_yearly_pwjet)[:,0], np.array(o3_gmi_yearly_pwjet)[:,0], 
#     o3_days[:len(Vsfc_yearly_pwjet)], lat_gmi, lng_gmi)
# o3_yearly_mean_eqjet, o3_yearly_stationary_eqjet, o3_yearly_transient_eqjet, \
#     o3_yearly_total_eqjet = globalo3_calculate.meridional_flux(
#     np.array(Vsfc_yearly_eqjet)[:,0], np.array(o3_gmi_yearly_eqjet)[:,0], 
#     o3_days[:len(Vsfc_yearly_eqjet)], lat_gmi, lng_gmi)
# # Plot
# total = o3_daily_total
# mean = o3_daily_mean
# eddy = o3_daily_stationary+o3_daily_transient_pwjet
# pwjet_total = o3_daily_total_pwjet
# pwjet_mean = o3_daily_mean_pwjet
# pwjet_eddy = o3_daily_stationary_pwjet+o3_daily_transient_pwjet
# eqjet_total = o3_daily_total_eqjet
# eqjet_mean = o3_daily_mean_eqjet
# eqjet_eddy = o3_daily_stationary_eqjet+o3_daily_transient_eqjet
# f, axes = plt.subplots(1, 3, sharey=True, figsize=(8,4))
# # All days
# axes[0].plot(lat_gmi, total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[0].plot(lat_gmi, mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[0].plot(lat_gmi, eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)
# # PW days  
# axes[1].plot(lat_gmi, pwjet_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[1].plot(lat_gmi, pwjet_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[1].plot(lat_gmi, pwjet_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)    
# # EW days  
# axes[2].plot(lat_gmi, eqjet_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[2].plot(lat_gmi, eqjet_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[2].plot(lat_gmi, eqjet_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)        
# for ax in [axes[0], axes[1], axes[2]]:
#     ax.set_xlim([0,90])
#     ax.set_xticks([0,30,60,90])
#     ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
#     ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
#         linestyles='--', linewidths=0.75)
# axes[0].set_ylabel('O$_{3}$ flux [ppbv m s$^{\mathregular{-1}}$]', fontsize=14)
# axes[0].set_title('All days', fontsize=14)
# axes[1].set_title('PW days', fontsize=14)
# axes[2].set_title('EW days', fontsize=14)
# plt.subplots_adjust(left=0.2, bottom=0.15)    
# axes[2].legend(ncol=1, frameon=False)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'zonalavg_o3flux_daily_2000-2010.png', dpi=300)

"""GEOSCHEM TRACER-JET RELATIONSHIP MAPS"""
# labels = ['$\chi_{0-10^{\circ}}$', '$\chi_{10-20^{\circ}}$', 
#     '$\chi_{20-30^{\circ}}$', '$\chi_{30-40^{\circ}}$', 
#     '$\chi_{40-50^{\circ}}$', '$\chi_{50-60^{\circ}}$',
#     '$\chi_{60-70^{\circ}}$', '$\chi_{70-80^{\circ}}$',
#     '$\chi_{80-90^{\circ}}$', '$\chi_{0-90^{\circ}}$']
# filename = ['trac_0_10', 'trac_10_20', 'trac_20_30', 'trac_30_40',
#     'trac_40_50', 'trac_50_60', 'trac_60_70', 'trac_70_80', 'trac_80_90', 
#     'global']
# # Loop through tracers 
# for i, tracer in enumerate([TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
#     TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90, GLOBAL]):
#     r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
#         edj_dist, lat_gc, lng_gc)
#     significance_r_tracerjet = \
#         globalo3_calculate.calculate_r_significance(np.nanmean(tracer, axis=1), 
#         edj_dist, r_tracerjet, lat_gc, lng_gc)
#     # Find (PW - EW) jet composites
#     eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
#         eqjet_tracer = globalo3_calculate.segregate_field_bylat(
#         np.nanmean(tracer, axis=1), lng_gc, edj, np.arange(0, len(tracer), 1))
#     # Plot (PW-EW) jet tracer difference
#     nhmap(lat_gc, 
#         lng_gc, 
#         (pwjet_tracer-eqjet_tracer)*1e9,
#         labels[i],
#         '[ppbv]', 
#         np.linspace(-200, 200, 11), 
#         'coolwarm', 
#         'pwjet-eqjet_'+filename[i]+'_geoschem_985-800_DJF',
#         extend='both', 
#         hatch=significance_r_tracerjet, 
#         ebar = edj)
#     # Calculate components of zonal- and time-mean transport for tracer from 
#     # GEOSChem for all days 
#     tracer_total, tracer_mean, tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer, Vcolumn_gc, np.arange(0, len(tracer), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))
#     # Separate into PW and EW days
#     Vcolumn_gc_eqjet, Vcolumn_gc_pwjet, tracer_eqjet, tracer_pwjet = \
#         globalo3_calculate.sortfield_byjetlat_column(Vcolumn_gc, tracer, 
#         edj, lng_gc, lat_gc, lev_gc, psize=30) 
#     # Zonal- and time-mean transport for PW jet 
#     pwjet_tracer_total, pwjet_tracer_mean, pwjet_tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer_pwjet, Vcolumn_gc_pwjet, np.arange(0, len(tracer_pwjet), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))        
#     # Zonal- and time-mean transport for EW jet 
#     eqjet_tracer_total, eqjet_tracer_mean, eqjet_tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer_eqjet, Vcolumn_gc_eqjet, np.arange(0, len(tracer_eqjet), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))          
#     # Plot 
#     f, axes = plt.subplots(1, 3, sharey=True, figsize=(8,4))
#     # All days
#     axes[0].plot(lat_gc, tracer_total, ls='-', color='#A2C3E0', lw=2, 
#         label='Total', zorder=2)
#     axes[0].plot(lat_gc, tracer_mean, ls='-', color='#EF9802', lw=2, 
#         label='Mean', zorder=3)
#     axes[0].plot(lat_gc, tracer_eddy, ls='-', color='#3F79B7', lw=2, 
#         label='Eddy', zorder=4)
#     # PW days  
#     axes[1].plot(lat_gc, pwjet_tracer_total, ls='-', color='#A2C3E0', lw=2, 
#         label='Total', zorder=2)
#     axes[1].plot(lat_gc, pwjet_tracer_mean, ls='-', color='#EF9802', lw=2, 
#         label='Mean', zorder=3)
#     axes[1].plot(lat_gc, pwjet_tracer_eddy, ls='-', color='#3F79B7', lw=2, 
#         label='Eddy', zorder=4)    
#     # EW days  
#     axes[2].plot(lat_gc, eqjet_tracer_total, ls='-', color='#A2C3E0', lw=2, 
#         label='Total', zorder=2)
#     axes[2].plot(lat_gc, eqjet_tracer_mean, ls='-', color='#EF9802', lw=2, 
#         label='Mean', zorder=3)
#     axes[2].plot(lat_gc, eqjet_tracer_eddy, ls='-', color='#3F79B7', lw=2, 
#         label='Eddy', zorder=4)        
#     for ax in [axes[0], axes[1], axes[2]]:
#         ax.set_xlim([0,90])
#         ax.set_xticks([0,30,60,90])
#         ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
#         ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
#             linestyles='--', linewidths=0.75)
#     axes[0].set_ylabel('GEOSChem %s\n800-955 hPa Flux '%labels[i]+
#         '[kg s$^{\mathregular{-1}}$]', fontsize=14)
#     axes[0].set_title('All days', fontsize=14)
#     axes[1].set_title('PW days', fontsize=14)
#     axes[2].set_title('EW days', fontsize=14)
#     plt.subplots_adjust(left=0.2, bottom=0.15)    
#     axes[2].legend(ncol=1, frameon=False)
#     plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#         'zonalavg_verticallyintegratedflux_%s_geoschem_955-800_DJF.png'
#         %filename[i], dpi=300)

"""ZONALLY-AVERAGED TRACER-JET CORRELATIONS"""
# labels = ['$\chi_{0-10^{\circ}}$', '$\chi_{10-20^{\circ}}$', 
#     '$\chi_{20-30^{\circ}}$', '$\chi_{30-40^{\circ}}$', 
#     '$\chi_{40-50^{\circ}}$', '$\chi_{50-60^{\circ}}$',
#     '$\chi_{60-70^{\circ}}$', '$\chi_{70-80^{\circ}}$',
#     '$\chi_{80-90^{\circ}}$']
# colors = ['#081d58','#253494','#225ea8','#1d91c0','#41b6c4','#7fcdbb',
#           '#c7e9b4','#edf8b1','#ffffd9']
# fig = plt.figure()
# ax = plt.subplot2grid((1,1),(0,0))
# # Loop through tracers 
# for i, tracer in enumerate([TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
#     TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90]):
#     r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
#         edj_dist, lat_gc, lng_gc)
#     r_tracerjet_zm = np.nanmean(r_tracerjet, axis=-1)
#     ax.plot(lat_gc, r_tracerjet_zm, color=colors[i], ls='-', lw=2, label=labels[i])
# plt.legend(ncol=3, loc=2, fontsize=6)
# ax.set_xlim([0, 90])
# ax.set_xlabel('Latitude [$^{\circ}$N]')
# ax.set_ylabel('r($\chi$, $\phi_{jet} - \phi$)')
# plt.savefig('/Users/ghkerr/Desktop/zonalavg_r_tracerjet_JJA.png', dpi=300)

"""COMPUTE ~950-800 HPA VERTICALLY INTEGRATED NH50 AND CO50 FLUX FROM 
GEOSCHEM AND CO50 FLUX FROM REPLAY AND DETERMINE RELATIONSHIP WITH THE JET 
STREAM"""
# # Calculate components of zonal- and time-mean transport for CO50 from 
# # GEOSChem
# co_50_gc_total, co_50_gc_mean, co_50_gc_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     co_50_gc, Vcolumn_gc, np.arange(0, len(co_50_gc), 1), 
#     lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))
# fig = plt.figure()
# ax = plt.subplot2grid((1,1),(0,0))
# ax.plot(lat_gc, co_50_gc_total, ls='-', color='#A2C3E0', lw=2, label='Total', 
#     zorder=2)
# ax.plot(lat_gc, co_50_gc_mean, ls='-', color='#EF9802', lw=2, label='Mean', 
#     zorder=3)
# ax.plot(lat_gc, co_50_gc_eddy, ls='-', color='#3F79B7', lw=2, label='Eddy', 
#     zorder=4)        
# ax.set_xlim([0,90])
# ax.set_xticks([0,30,60,90])
# ax.set_ylim([-1000, 3000])
# ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
# ax.set_ylabel('800-955 hPa Flux [kg s$^{\mathregular{-1}}$]', fontsize=14)
# ax.set_title('GEOSChem CO$_{\mathregular{50}}$', fontsize=14)
# ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
#     linestyles='--', linewidths=0.75)
# plt.legend(ncol=3, frameon=False)
# plt.subplots_adjust(left=0.2)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'zonalavg_verticallyintegratedflux_co50_geoschem_955-800.png', dpi=300)
# plt.show()
# # Calculate components of zonal- and time-mean transport for CO50 from 
# # Replay
# co_50_replay_total, co_50_replay_mean, co_50_replay_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     co_50_replay/1e9, Vcolumn_replay, np.arange(0, len(co_50_replay), 1), 
#     lat_gc, lng_gc, lev_replay, 950., 800., (28/28.97))
# # Separate into PW and EW days
# Vcolumn_replay_eqjet, Vcolumn_replay_pwjet, co_50_replay_eqjet, \
#     co_50_replay_pwjet = globalo3_calculate.sortfield_byjetlat_column(
#     Vcolumn_replay, co_50_replay/1e9, edj, lng_gc, lat_gc, lev_replay, psize=30)
# # Zonal- and time-mean transport of CO50 for PW jet 
# pwjet_co_50_total, pwjet_co_50_mean, pwjet_co_50_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     co_50_replay_pwjet, Vcolumn_replay_pwjet, 
#     np.arange(0, len(co_50_replay_pwjet), 1), lat_gc, lng_gc, lev_replay, 
#     950., 800., (28/28.97))        
# # Zonal- and time-mean transport of CO50 for EW jet 
# eqjet_co_50_total, eqjet_co_50_mean, eqjet_co_50_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     co_50_replay_eqjet, Vcolumn_replay_eqjet, 
#     np.arange(0, len(co_50_replay_eqjet), 1), lat_gc, lng_gc, lev_replay, 
#     950., 800., (28/28.97))
# # Plot
# f, axes = plt.subplots(1, 3, sharey=True, figsize=(8,4))
# # All days
# axes[0].plot(lat_gc, co_50_replay_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[0].plot(lat_gc, co_50_replay_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[0].plot(lat_gc, co_50_replay_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)
# # PW days  
# axes[1].plot(lat_gc, pwjet_co_50_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[1].plot(lat_gc, pwjet_co_50_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[1].plot(lat_gc, pwjet_co_50_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)    
# # EW days  
# axes[2].plot(lat_gc, eqjet_co_50_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[2].plot(lat_gc, eqjet_co_50_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[2].plot(lat_gc, eqjet_co_50_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)        
# for ax in [axes[0], axes[1], axes[2]]:
#     ax.set_xlim([0,90])
#     ax.set_xticks([0,30,60,90])
#     ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
#     ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
#         linestyles='--', linewidths=0.75)
# axes[0].set_ylabel('Replay CO$_{\mathregular{50}}$\n800-955 hPa Flux '+
#     '[kg s$^{\mathregular{-1}}$]', fontsize=14)
# axes[0].set_title('All days', fontsize=14)
# axes[1].set_title('PW days', fontsize=14)
# axes[2].set_title('EW days', fontsize=14)
# plt.subplots_adjust(left=0.2, bottom=0.15)    
# axes[2].legend(ncol=1, frameon=False)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'zonalavg_verticallyintegratedflux_co_50_replay_950-800.png',
#     dpi=300)
# # Zonal- and time-mean transport of NH50 for all days 
# nh_50_replay_total, nh_50_replay_mean, nh_50_replay_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     nh_50_replay/1e9, Vcolumn_replay, np.arange(0, len(nh_50_replay), 1), 
#     lat_gc, lng_gc, lev_replay, 950., 800., (28/28.97))
# # Separate into PW and EW days
# Vcolumn_replay_eqjet, Vcolumn_replay_pwjet, nh_50_replay_eqjet, \
#     nh_50_replay_pwjet = globalo3_calculate.sortfield_byjetlat_column(
#     Vcolumn_replay, nh_50_replay/1e9, edj, lng_gc, lat_gc, lev_replay, psize=30)
# # Zonal- and time-mean transport of NH50 for PW jet 
# pwjet_nh_50_total, pwjet_nh_50_mean, pwjet_nh_50_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     nh_50_replay_pwjet, Vcolumn_replay_pwjet, 
#     np.arange(0, len(nh_50_replay_pwjet), 1), lat_gc, lng_gc, lev_replay, 
#     950., 800., (28/28.97))        
# # Zonal- and time-mean transport of NH50 for EW jet 
# eqjet_nh_50_total, eqjet_nh_50_mean, eqjet_nh_50_eddy = \
#     globalo3_calculate.verticallyintegrated_meridional_flux(
#     nh_50_replay_eqjet, Vcolumn_replay_eqjet, 
#     np.arange(0, len(nh_50_replay_eqjet), 1), lat_gc, lng_gc, lev_replay, 
#     950., 800., (28/28.97))  
# # Plot 
# f, axes = plt.subplots(1, 3, sharey=True, figsize=(8,4))
# # All days
# axes[0].plot(lat_gc, nh_50_replay_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[0].plot(lat_gc, nh_50_replay_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[0].plot(lat_gc, nh_50_replay_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)
# # PW days  
# axes[1].plot(lat_gc, pwjet_nh_50_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[1].plot(lat_gc, pwjet_nh_50_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[1].plot(lat_gc, pwjet_nh_50_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)    
# # EW days  
# axes[2].plot(lat_gc, eqjet_nh_50_total, ls='-', color='#A2C3E0', lw=2, 
#     label='Total', zorder=2)
# axes[2].plot(lat_gc, eqjet_nh_50_mean, ls='-', color='#EF9802', lw=2, 
#     label='Mean', zorder=3)
# axes[2].plot(lat_gc, eqjet_nh_50_eddy, ls='-', color='#3F79B7', lw=2, 
#     label='Eddy', zorder=4)        
# for ax in [axes[0], axes[1], axes[2]]:
#     ax.set_xlim([0,90])
#     ax.set_xticks([0,30,60,90])
#     ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
#     ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
#         linestyles='--', linewidths=0.75)
# axes[0].set_ylabel('Replay NH$_{\mathregular{50}}$\n800-955 hPa Flux '+
#     '[kg s$^{\mathregular{-1}}$]', fontsize=14)
# axes[0].set_title('All days', fontsize=14)
# axes[1].set_title('PW days', fontsize=14)
# axes[2].set_title('EW days', fontsize=14)
# plt.subplots_adjust(left=0.2, bottom=0.15)    
# axes[2].legend(ncol=1, frameon=False)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'zonalavg_verticallyintegratedflux_nh_50_replay_950-800.png',
#     dpi=300)
# plt.show()
# # Calculate and plot GEOSChem CO50-jet relationship
# r_co_50jet = globalo3_calculate.calculate_r(np.nanmean(co_50_gc*1e9,axis=1),
#     edj_dist, lat_gc, lng_gc)
# significance_r_co_50jet = globalo3_calculate.calculate_r_significance(
#     np.nanmean(co_50_gc*1e9,axis=1), edj_dist, r_co_50jet, lat_gc, lng_gc)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_co_50, \
#     eqjet_co_50 = globalo3_calculate.segregate_field_bylat(
#     np.nanmean(co_50_gc*1e9,axis=1), lng_gc, edj, 
#     np.arange(0, len(co_50_gc), 1))        
# nhmap(lat_gc, 
#     lng_gc, 
#     (pwjet_co_50-eqjet_co_50),
#     'GEOSChem CO$_{\mathregular{50}}$',
#     '[ppbv]', 
#     np.linspace(-10, 10, 11), 
#     'coolwarm', 
#     'pwjet-ewjet_co50_geoschem_1000-800',
#     extend='both', 
#     hatch=significance_r_co_50jet, 
#     ebar = edj)
# # Calculate and plot Replay CO50-jet relationship
# r_co_50jet = globalo3_calculate.calculate_r(np.nanmean(co_50_replay,axis=1),
#     edj_dist, lat_gc, lng_gc)
# significance_r_co_50jet = globalo3_calculate.calculate_r_significance(
#     np.nanmean(co_50_replay,axis=1), edj_dist, r_co_50jet, lat_gc, lng_gc)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_co_50, \
#     eqjet_co_50 = globalo3_calculate.segregate_field_bylat(
#     np.nanmean(co_50_replay,axis=1), lng_gc, edj, 
#     np.arange(0, len(co_50_replay), 1))        
# nhmap(lat_gc, 
#     lng_gc, 
#     (pwjet_co_50-eqjet_co_50),
#     'Replay CO$_{\mathregular{50}}$',
#     '[ppbv]', 
#     np.linspace(-10, 10, 11), 
#     'coolwarm', 
#     'pwjet-eqjet_co50_replay_1000-800',
#     extend='both', 
#     hatch=significance_r_co_50jet, 
#     ebar = edj)
# # Calculate and plot Replay NH50-jet relationship
# r_nh_50jet = globalo3_calculate.calculate_r(np.nanmean(nh_50_replay,axis=1), 
#     edj_dist, lat_gc, lng_gc)
# significance_r_nh_50jet = globalo3_calculate.calculate_r_significance(
#     np.nanmean(nh_50_replay,axis=1), edj_dist, r_nh_50jet, lat_gc, lng_gc)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_nh_50, \
#     eqjet_nh_50 = globalo3_calculate.segregate_field_bylat(
#     np.nanmean(nh_50_replay,axis=1), lng_gc, edj, 
#     np.arange(0, len(nh_50_replay), 1))
# nhmap(lat_gc, 
#     lng_gc, 
#     (pwjet_nh_50-eqjet_nh_50),
#     'Replay NH$_{\mathregular{50}}$',
#     '[ppbv]', 
#     np.linspace(-15000, 15000, 11), 
#     'coolwarm', 
#     'pwjet-eqjet_nh50_replay_1000-800',
#     extend='both', 
#     hatch=significance_r_nh_50jet, 
#     ebar = edj)
    
"""ZONALLY-AVERAGED TRACER CONCENTRATION, LATITUDINAL GRADIENT, AND
JET-TRACER RELATIONSHIP"""
# labels = ['$\chi_{0-10^{\circ}}$', '$\chi_{10-20^{\circ}}$', 
#     '$\chi_{20-30^{\circ}}$', '$\chi_{30-40^{\circ}}$', 
#     '$\chi_{40-50^{\circ}}$', '$\chi_{50-60^{\circ}}$',
#     '$\chi_{60-70^{\circ}}$', '$\chi_{70-80^{\circ}}$',
#     '$\chi_{80-90^{\circ}}$', '$\chi_{0-90^{\circ}}$', 'CO$_{50}$']
# filename = ['trac_0_10', 'trac_10_20', 'trac_20_30', 'trac_30_40',
#     'trac_40_50', 'trac_50_60', 'trac_60_70', 'trac_70_80', 'trac_80_90', 
#     'trac_0_90', 'co_50_geoschem']
# # Determine near-surface meridional wind-jet correlation
# r_Vjet = globalo3_calculate.calculate_r(V, edj_dist, lat_gc, lng_gc)
# # # Loop through tracers 
# for i, tracer in enumerate([TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
#     TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90, GLOBAL,
#     co_50_gc]):
#     # Find column-averaged mixing ratio and convert from mol mol-1 to ppm 
#     tracer = np.nanmean(tracer, axis=1)*1e6
#     V = np.nanmean(Vcolumn_gc, axis=1)
#     # Determine tracer-jet correlation
#     r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist, lat_gc, lng_gc)
#     lat = lat_gc
#     lng = lng_gc
#     # Zonal mean 
#     tracer_zm = np.nanmean(tracer, 
#         axis=np.where(np.array(tracer.shape)==lng.shape[0])[0][0])
#     # Latitudinal gradient
#     grad_tracer_zm = np.gradient(np.nanmean(tracer_zm, axis=0))
#     # Plotting
#     fig = plt.figure(figsize=(10,5))
#     ax1 = plt.subplot2grid((1,2),(0,0))
#     ax2 = plt.subplot2grid((1,2),(0,1))
#     # (Left) Zonal mean tracer concentration and (right) latitudinal gradient 
#     ax1.plot(lat, np.nanmean(tracer_zm, axis=0), '-k')
#     ax1.set_ylabel(r'%s [ppm]'%labels[i], fontsize=14)
#     ax1t = ax1.twinx()
#     ax1t.plot(lat, grad_tracer_zm, '-b')
#     ax1t.set_ylabel('d%s/d$\phi$ [ppm/$^{\circ}$]'%labels[i], fontsize=14, 
#         color='b')
#     # (Left) Jet-tracer correlation and (right) Correlation * gradient
#     ax2.plot(lat, np.nanmean(r_tracerjet, axis=1), '-k')
#     ax2.set_ylabel('r(%s, $\phi_{jet} - \phi$) [$\cdot$]'%labels[i], 
#         fontsize=14, color='k')
#     ax2t = ax2.twinx()
#     ax2t.plot(lat, -1.*(np.nanmean(r_Vjet, axis=1)*grad_tracer_zm), '-b')
#     ax2t.set_ylabel('$-$r(%s, $\phi_{jet} - \phi$) $\cdot$ '%labels[i]+\
#         'd%s/d$\phi$ [ppm/$^{\circ}$]'%labels[i], fontsize=14, color='b')
#     for ax in [ax1, ax1t, ax2, ax2t]:
#         ax.set_xlim([0, 90])
#         ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
#     plt.subplots_adjust(wspace=0.7, right=0.85)
#     plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#         'zonalavg_%s_grad%s_rjet%s.png'%(filename[i],filename[i],filename[i]), 
#         dpi=350)

"""1000-200 HPA V ANOMALY AND U WIND FOR NORTH AMERICA"""
# Ucolumnfull_gc, lat_merra_gc, lng_merra_gc, lev_merrafull_gc = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, months, 'U',
#     lngmin, latmax, lngmax, latmin, 150., 1005.)
# Ucolumnfull_gc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra_gc, lng_merra_gc, Ucolumnfull_gc)
# Vcolumnfull_gc, lat_merra_gc, lng_merra_gc, lev_merrafull_gc = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, months, 'V',
#     lngmin, latmax, lngmax, latmin, 150., 1005.)
# Vcolumnfull_gc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra_gc, lng_merra_gc, Vcolumnfull_gc)
# # Select fields over North America (125-60 deg W) (i.e., 235-300 deg E)
# Vna = Vcolumnfull_gc[:,:,:,94:121]
# Vnaavg = np.nanmean(Vna, axis=tuple((0,3)))
# Una = Ucolumnfull_gc[:,:,:,94:121]
# Unaavg = np.nanmean(Una, axis=tuple((0,3)))
# edjna = edj[:,94:121]
# edjnaavg = np.nanmean(edjna, axis=-1)
# # Generate dates to label plots
# import datetime
# start_date = datetime.date(2008, 6, 1)
# end_date = datetime.date(2010, 8, 31)
# dates_all = [start_date + datetime.timedelta(n) for n in 
#     range(int((end_date-start_date).days))]
# dates = []
# for date in dates_all:
#     if (date.month == 6) or (date.month == 7) or (date.month == 8):
#         dates.append(date)
# for day in np.arange(0, 276, 1):
#     fig = plt.figure()
#     ax = plt.subplot2grid((1,1),(0,0))
#     ax.set_title('%s\nZonal average over North America (235-300$^{\circ}$E)'
#         %dates[day])
#     clevs = np.linspace(-5, 5, 11)
#     # V
#     cf = ax.contourf(lat_gc, lev_merrafull_gc, 
#         np.nanmean(Vna[day], axis=-1)-Vnaavg, clevs, cmap=plt.get_cmap('bwr'), 
#         extend='both')
#     # Eddy-driven jet
#     ax.errorbar(np.nanmean(edjna[day]), 500., 
#         xerr=np.nanstd(edjna[day]), zorder=10, color='k', markersize=3, 
#         elinewidth=1.25, ecolor='k', fmt='o')
#     CS = ax.contour(lat_gc, lev_merrafull_gc, 
#         np.nanmean(Una[day], axis=-1),
#         [0, 5, 10, 15, 20, 25], colors='k')
#     ax.clabel(CS, fontsize=9, fmt='%1d', inline=1)    
#     ax.set_xlim([25, 75])
#     ax.set_xlabel('Latitude [$^{\circ}$N]')
#     ax.set_ylabel('Pressure [hPa]')
#     plt.colorbar(cf, label='$V - \overline{V}$ [m s$^{-1}$]')
#     plt.gca().invert_yaxis()
#     plt.savefig('/Users/ghkerr/Desktop/tempjet/zana_%s.jpg'%str(day).zfill(3), 
#         dpi=200)
#     plt.show()
# # JJA average
# fig = plt.figure()
# ax = plt.subplot2grid((1,1),(0,0))
# clevs = np.linspace(-2, 2, 9)
# cf = ax.contourf(lat_gc, lev_merrafull_gc, Vnaavg, clevs, 
#     cmap=plt.get_cmap('bwr'), extend='both')
# CS = ax.contour(lat_gc, lev_merrafull_gc, Unaavg, [0, 5, 10, 15, 20, 25], 
#     colors='k')
# ax.clabel(CS, fontsize=9, fmt='%1d', inline=1)    
# ax.set_xlim([25, 75])
# ax.set_xlabel('Latitude [$^{\circ}$N]')
# ax.set_ylabel('Pressure [hPa]')
# plt.colorbar(cf, label='$\overline{V}$ [m s$^{-1}$]')
# plt.gca().invert_yaxis()
# plt.savefig('/Users/ghkerr/Desktop/trial222.png', dpi=350)

"""DJF TRACERS AND JET"""
# years = [2008, 2009, 2010]
# djf = ['jan', 'feb', 'dec']
# # Load column GEOSChem CO50 and MERRA-2 column meridional wind with the 
# # same vertical levels as CO50 in GEOSChem
# co_50_gc_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_CO_50', latmin, latmax, lngmin, lngmax, 800., 1000.)
# Vcolumn_gc_djf, lat_merra_gc, lng_merra_gc, lev_merra_gc = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf, 'V',
#     lngmin, latmax, lngmax, latmin, pmin, pmax)
# Vcolumn_gc_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra_gc, lng_merra_gc, Vcolumn_gc_djf)
# TRAC_0_10_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_0_10', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_10_20_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_10_20', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_20_30_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_20_30', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_30_40_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_30_40', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_40_50_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_40_50', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_50_60_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_50_60', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_60_70_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_60_70', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_70_80_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_70_80', latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_80_90_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_TRAC50_80_90', latmin, latmax, lngmin, lngmax, pmin, pmax)
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
# edj_djf, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf, 'U', 
#     lngmin, latmax, lngmax, latmin, 487., 526., operation='mean')
# # Degrade to resolution of GEOSChem    
# edj_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, 
#     lat_edj, lng_edj, edj_djf)
# # Subset fields in mid-latitudes
# edj_djf, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj_djf, lat_gc, 
#     lng_gc, 0., 360., 20., 70.)
# # Determine jet latitude
# edj_djf = tracer_calculate.find_jetlatitude(edj_djf, lat_edj, lng_edj)
# # Find distance from the jet where positive distances imply that the jet is
# # north of a particular location
# edj_dist_djf = tracer_calculate.determine_jetdist(edj_djf, lat_gc, lng_gc)

"""JJA AND DJF EDDY AND MEAN MERIDIONAL FLUXES FOR COMPARISON WITH HARTMANN 
   (2007)"""
# import numpy as np
# import sys
# sys.path.append('/Users/ghkerr/phd/globalo3/')
# import globalo3_open, globalo3_calculate
# sys.path.append('/Users/ghkerr/phd/tracer/')
# import tracer_open, tracer_calculate
# years = [2008, 2009, 2010]
# djf = ['jan', 'feb', 'dec']
# jja = ['jun', 'jul', 'aug']
# latmin = 0.
# latmax = 90.
# lngmin = 0.
# lngmax = 360.
# pmin = 150.
# pmax = 1005.
# import numpy as np
# import matplotlib.pyplot as plt
# # Load column GEOSChem CO50 to obtain GEOSChem coordinates 
# co_50_djf, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, djf, 
#     'SpeciesConc_CO_50', latmin, latmax, lngmin, lngmax, pmin, pmax)
# # For DJF
# Tcol_djf, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf,
#     'T', lngmin, latmax, lngmax, latmin, pmin, pmax)
# Tcol_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, Tcol_djf)
# Vcol_djf, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf, 'V', 
#     lngmin, latmax, lngmax, latmin, pmin, pmax)
# Vcol_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, Vcol_djf)
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
# edj_djf, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf, 'U', 
#     lngmin, latmax, lngmax, latmin, 487., 526., operation='mean')
# # Degrade to resolution of GEOSChem    
# edj_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, 
#     lat_edj, lng_edj, edj_djf)
# # Subset fields in mid-latitudes
# edj_djf, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj_djf, lat_gc, 
#     lng_gc, 0., 360., 20., 70.)
# # Determine jet latitude
# edj_djf = tracer_calculate.find_jetlatitude(edj_djf, lat_edj, lng_edj)
# # Find distance from the jet where positive distances imply that the jet is
# # north of a particular location
# edj_dist_djf = tracer_calculate.determine_jetdist(edj_djf, lat_gc, lng_gc)
# # For JJA
# co_50_jja, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, jja, 
#     'SpeciesConc_CO_50', latmin, latmax, lngmin, lngmax, pmin, pmax)
# Tcol_jja, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, jja,
#     'T', lngmin, latmax, lngmax, latmin, pmin, pmax)
# Tcol_jja = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, Tcol_jja)
# Vcol_jja, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, jja, 'V', 
#     lngmin, latmax, lngmax, latmin, pmin, pmax)
# Vcol_jja = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, Vcol_jja)
# edj_jja, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, jja, 'U', 
#     lngmin, latmax, lngmax, latmin, 487., 526., operation='mean')
# # Degrade to resolution of GEOSChem    
# edj_jja = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, 
#     lat_edj, lng_edj, edj_jja)
# # Subset fields in mid-latitudes
# edj_jja, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj_jja, lat_gc, 
#     lng_gc, 0., 360., 20., 70.)
# # Determine jet latitude
# edj_jja = tracer_calculate.find_jetlatitude(edj_jja, lat_edj, lng_edj)
# # Find distance from the jet where positive distances imply that the jet is
# # north of a particular location
# edj_dist_jja = tracer_calculate.determine_jetdist(edj_jja, lat_gc, lng_gc)
# meancol_jja, stationarycol_jja, transientcol_jja, totalcol_jja = [], [], [], []
# meancol_djf, stationarycol_djf, transientcol_djf, totalcol_djf = [], [], [], []
# # Loop through levels
# for lev in np.arange(lev_merra.shape[0]):
#     print('Determining %d hPa fluxes...'%lev_merra[lev])
#     # Single level fluxes for JJA
#     mean_jja, stationary_jja, transient_jja, total_jja = \
#         globalo3_calculate.meridional_flux(Vcol_jja[:,lev], Tcol_jja[:,lev], 
#         np.arange(0, len(Tcol_jja), 1), lat_gc, lng_gc)
#     # Append to multi-level lists
#     meancol_jja.append(mean_jja)
#     stationarycol_jja.append(stationary_jja)
#     transientcol_jja.append(transient_jja)
#     totalcol_jja.append(total_jja)
#     # For DJF
#     mean_djf, stationary_djf, transient_djf, total_djf = \
#         globalo3_calculate.meridional_flux(Vcol_djf[:,lev], Tcol_djf[:,lev], 
#         np.arange(0, len(Tcol_djf), 1), lat_gc, lng_gc)
#     meancol_djf.append(mean_djf)
#     stationarycol_djf.append(stationary_djf)
#     transientcol_djf.append(transient_djf)
#     totalcol_djf.append(total_djf)
# # Stack and combine transient and stationary eddies into single term
# eddy_jja = np.vstack(transientcol_jja) + np.vstack(stationarycol_jja)
# eddy_djf = np.vstack(transientcol_djf) + np.vstack(stationarycol_djf)
# mean_jja = np.vstack(meancol_jja)
# mean_djf = np.vstack(meancol_djf)
# # Plotting
# fig = plt.figure()
# ax1 = plt.subplot2grid((2,2),(0,0))
# ax2 = plt.subplot2grid((2,2),(0,1))
# ax3 = plt.subplot2grid((2,2),(1,0))
# ax4 = plt.subplot2grid((2,2),(1,1))
# clevs_eddy = np.linspace(-20, 20, 11)
# clevs_mean = np.linspace(-700, 700, 11)
# # DJF eddy 
# cf = ax1.contourf(lat_gc, lev_merra, eddy_djf, clevs_eddy, 
#     cmap=plt.get_cmap('bwr'), extend='both')
# ax1.invert_yaxis()
# ax1.set_title('DJF Eddy')
# ax1.set_ylabel('Pressure [hPa]')
# fig.colorbar(cf, ax=ax1)
# # DJF mean 
# cf = ax2.contourf(lat_gc, lev_merra, mean_djf, clevs_mean, 
#     cmap=plt.get_cmap('bwr'), extend='both')
# ax2.invert_yaxis()
# ax2.set_title('DJF Mean')
# fig.colorbar(cf, ax=ax2)
# # JJA eddy 
# cf = ax3.contourf(lat_gc, lev_merra, eddy_jja, clevs_eddy, 
#     cmap=plt.get_cmap('bwr'), extend='both')                  
# ax3.invert_yaxis()
# ax3.set_title('JJA Eddy')
# ax3.set_xlabel('Latitude [$^{\circ}$N]')
# ax3.set_ylabel('Pressure [hPa]')
# fig.colorbar(cf, ax=ax3)
# # JJA mean 
# cf = ax4.contourf(lat_gc, lev_merra, mean_jja, clevs_mean, 
#     cmap=plt.get_cmap('bwr'), extend='both')                  
# ax4.invert_yaxis()
# ax4.set_title('JJA Mean')
# ax4.set_xlabel('Latitude [$^{\circ}$N]')
# fig.colorbar(cf, ax=ax4)
# plt.subplots_adjust(hspace=0.4, wspace=0.25)
# plt.savefig('/Users/ghkerr/Desktop/hartmann_flux.png', dpi=300)

"""xxx """
# # Indices corresponding to the "mid-latitudes"
# mls = 0 # lat_gc[15] = 30 deg N
# mln = 45 # lat_gc[35] = 70 deg N
# # For JJA fields
# eddy_ml, mean_ml = [], []
# pweddy_ml, pwmean_ml = [], []
# eqeddy_ml, eqmean_ml = [], []
# for tracer in [TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
#     TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90]:
#     # Calculate components of zonal- and time-mean transport for tracer for 
#     # all days
#     tracer_total, tracer_mean, tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer, Vcolumn_gc, np.arange(0, len(tracer), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))
#     eddy_ml.append(np.sum(np.array(tracer_eddy)[mls:mln]))
#     mean_ml.append(np.sum(np.array(tracer_mean)[mls:mln]))
#     # Find (PW - EW) jet composites
#     eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
#         eqjet_tracer = globalo3_calculate.segregate_field_bylat(
#         np.nanmean(tracer, axis=1), lng_gc, edj, np.arange(0, len(tracer), 1))
#     # Separate into PW and EW days
#     Vcolumn_gc_eqjet, Vcolumn_gc_pwjet, tracer_eqjet, tracer_pwjet = \
#         globalo3_calculate.sortfield_byjetlat_column(Vcolumn_gc, tracer, 
#         edj, lng_gc, lat_gc, lev_gc, psize=30) 
#     # Zonal- and time-mean transport for PW jet 
#     pwjet_tracer_total, pwjet_tracer_mean, pwjet_tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer_pwjet, Vcolumn_gc_pwjet, np.arange(0, len(tracer_pwjet), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))        
#     # Zonal- and time-mean transport for EW jet 
#     eqjet_tracer_total, eqjet_tracer_mean, eqjet_tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer_eqjet, Vcolumn_gc_eqjet, np.arange(0, len(tracer_eqjet), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))    
#     pweddy_ml.append(np.sum(np.array(pwjet_tracer_eddy)[mls:mln]))
#     pwmean_ml.append(np.sum(np.array(pwjet_tracer_mean)[mls:mln]))
#     eqeddy_ml.append(np.sum(np.array(eqjet_tracer_eddy)[mls:mln]))
#     eqmean_ml.append(np.sum(np.array(eqjet_tracer_mean)[mls:mln]))
# # For DJF fields
# eddy_ml_djf, mean_ml_djf = [], []
# pweddy_ml_djf, pwmean_ml_djf = [], []
# eqeddy_ml_djf, eqmean_ml_djf = [], []
# for tracer in [TRAC_0_10_djf, TRAC_10_20_djf, TRAC_20_30_djf, TRAC_30_40_djf, 
#     TRAC_40_50_djf, TRAC_50_60_djf, TRAC_60_70_djf, TRAC_70_80_djf, 
#     TRAC_80_90_djf]:
#     # Calculate components of zonal- and time-mean transport for tracer for 
#     # all days    
#     tracer_total, tracer_mean, tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer, Vcolumn_gc_djf, np.arange(0, len(tracer), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))
#     eddy_ml_djf.append(np.sum(np.array(tracer_eddy)[mls:mln]))
#     mean_ml_djf.append(np.sum(np.array(tracer_mean)[mls:mln]))
#     # Find (PW - EW) jet composites
#     eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
#         eqjet_tracer = globalo3_calculate.segregate_field_bylat(
#         np.nanmean(tracer, axis=1), lng_gc, edj_djf, 
#         np.arange(0, len(tracer), 1))
#     # Separate into PW and EW days
#     Vcolumn_gc_eqjet, Vcolumn_gc_pwjet, tracer_eqjet, tracer_pwjet = \
#         globalo3_calculate.sortfield_byjetlat_column(Vcolumn_gc_djf, tracer, 
#         edj, lng_gc, lat_gc, lev_gc, psize=30) 
#     # Zonal- and time-mean transport for PW jet 
#     pwjet_tracer_total, pwjet_tracer_mean, pwjet_tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer_pwjet, Vcolumn_gc_pwjet, np.arange(0, len(tracer_pwjet), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))        
#     # Zonal- and time-mean transport for EW jet 
#     eqjet_tracer_total, eqjet_tracer_mean, eqjet_tracer_eddy = \
#         globalo3_calculate.verticallyintegrated_meridional_flux(
#         tracer_eqjet, Vcolumn_gc_eqjet, np.arange(0, len(tracer_eqjet), 1), 
#         lat_gc, lng_gc, lev_gc, 955., 800., (28/28.97))    
#     pweddy_ml_djf.append(np.sum(np.array(pwjet_tracer_eddy)[mls:mln]))
#     pwmean_ml_djf.append(np.sum(np.array(pwjet_tracer_mean)[mls:mln]))
#     eqeddy_ml_djf.append(np.sum(np.array(eqjet_tracer_eddy)[mls:mln]))
#     eqmean_ml_djf.append(np.sum(np.array(eqjet_tracer_mean)[mls:mln]))
# labels = ['$\chi_{0-10^{\circ}}$', '$\chi_{10-20^{\circ}}$', 
#     '$\chi_{20-30^{\circ}}$', '$\chi_{30-40^{\circ}}$', 
#     '$\chi_{40-50^{\circ}}$', '$\chi_{50-60^{\circ}}$',
#     '$\chi_{60-70^{\circ}}$', '$\chi_{70-80^{\circ}}$',
#     '$\chi_{80-90^{\circ}}$']
# x = np.arange(len(labels))  # label locations
# width = 0.35  # width of the bars
# # Plotting 
# fig = plt.figure(figsize=(15,5))
# ax1 = plt.subplot2grid((2,3),(0,0))
# ax2 = plt.subplot2grid((2,3),(0,1))
# ax3 = plt.subplot2grid((2,3),(0,2))
# ax4 = plt.subplot2grid((2,3),(1,0))
# ax5 = plt.subplot2grid((2,3),(1,1))
# ax6 = plt.subplot2grid((2,3),(1,2))
# # JJA (left) all, (middle) PW, and (right) EW eddy 
# rects1 = ax1.bar(x-width/2, np.array(mean_ml), width, label='Mean')
# rects2 = ax1.bar(x+width/2, np.array(eddy_ml), width, label='Eddy')
# ax2.bar(x-width/2, np.array(pwmean_ml), width)
# ax2.bar(x+width/2, np.array(pweddy_ml), width)
# ax3.bar(x-width/2, np.array(eqmean_ml), width)
# ax3.bar(x+width/2, np.array(eqeddy_ml), width)
# # Aesthetics
# ax1.set_title('All')
# ax2.set_title('PW')
# ax3.set_title('EW')
# ax2.set_yticklabels([''])
# ax3.set_yticklabels([''])
# ax1.set_ylabel('JJA 0$-$90$^{\circ}$N flux [kg s$^{-1}$]')
# for ax in [ax1, ax2, ax3]:
#     ax.set_ylim([-900000,900000])
#     ax.set_xticks(np.arange(0, 9, 1))
#     ax.set_xticklabels([''])
# # DJF (left) all, (middle) PW, and (right) EW eddy 
# ax4.bar(x-width/2, np.array(mean_ml_djf), width, label='Mean')
# ax4.bar(x+width/2, np.array(eddy_ml_djf), width, label='Eddy')
# ax5.bar(x-width/2, np.array(pwmean_ml_djf), width)
# ax5.bar(x+width/2, np.array(pweddy_ml_djf), width)
# ax6.bar(x-width/2, np.array(eqmean_ml_djf), width)
# ax6.bar(x+width/2, np.array(eqeddy_ml_djf), width)
# ax5.set_yticklabels([''])
# ax6.set_yticklabels([''])
# ax4.set_ylabel('DJF 0$-$9$^{\circ}$N flux [kg s$^{-1}$]')
# for ax in [ax4, ax5, ax6]:
#     ax.set_ylim([-700000,700000])    
#     ax.set_xticks(np.arange(0, 9, 1))
#     ax.set_xticklabels(labels, rotation=35)
# ax1.legend()
# plt.subplots_adjust(wspace=0.05)
# plt.savefig('/Users/ghkerr/Desktop/bar_meaneddyall_nh.png', dpi=300)

"""ZONALLY-AVERAGED JET-TRACER CORRELATIONS AT DIFFERENT PRESSURE LEVELS"""
# import numpy as np
# import sys
# sys.path.append('/Users/ghkerr/phd/globalo3/')
# import globalo3_open, globalo3_calculate
# sys.path.append('/Users/ghkerr/phd/tracer/')
# import tracer_open, tracer_calculate
# years = [2008, 2009, 2010]
# # months = ['jun', 'jul', 'aug']
# months = ['jan', 'feb', 'dec']
# latmin = 0.
# latmax = 90.
# lngmin = 0.
# lngmax = 360.
# pmin = 800.
# pmax = 1005.
# tracer, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, months, 
#     'SpeciesConc_TRAC50_70_80', latmin, latmax, lngmin, lngmax, pmin, pmax)
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
# edj, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, months, 'U', 
#     lngmin, latmax, lngmax, latmin, 487., 526., operation='mean')
# # Degrade to resolution of GEOSChem    
# edj = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, lat_edj, 
#     lng_edj, edj)
# # Subset fields in mid-latitudes
# edj, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj, lat_gc, lng_gc, 
#     0., 360., 20., 70.)
# # Determine jet latitude
# edj = tracer_calculate.find_jetlatitude(edj, lat_edj, lng_edj)
# # Find distance from the jet where positive distances imply that the jet is
# # north of a particular location
# edj_dist = tracer_calculate.determine_jetdist(edj, lat_gc, lng_gc)    
# label = '$\chi_{70\mathregular{-}80^{\circ}\mathregular{N}}$'
# filename = 'trac_70_80'
# fig = plt.figure()
# ax = plt.subplot2grid((1,1),(0,0))
# ax.set_title(label, fontsize=16)
# rcol_tracerjet = [] 
# # Levels of interest in GEOSChem pressure coordinates, units of hPa
# lev_interest = [985., 895., 775., 675., 563., 488., 375., 245., 208.]
# lev_interest_i = [np.where(lev_gc == lev)[0][0] for lev in lev_interest]           
# # Colors corresponding to the levels of interest; note the length of this 
# # list should be the same lenght as lev_interest_i 
# colors = ['#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8',
#     '#253494','#081d58', 'k']
# # Loop through vertical pressure levels of interest and calculate the 
# # jet-tracer correlation at each level
# for i,lev in enumerate(lev_interest_i):
#     r_tracerjet = globalo3_calculate.calculate_r(tracer[:,lev], edj_dist, 
#             lat_gc, lng_gc)
#     ax.plot(lat_gc, np.nanmean(r_tracerjet, axis=-1), color=colors[i],
#         label='%d hPa'%lev_interest[i])
# plt.legend(ncol=2, fontsize=8)
# ax.set_xlim([0,90])    
# ax.set_xticks([0,30,60,90])
# ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
# ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
#     linestyles='--', linewidths=0.75)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/zonalavg_r_%sjet_presscoords.png'
#     %filename, dpi=300)

"""r(TRACER, JET) OVER LAND VS OCEAN FOR DIFFERENT SOURCE REGIONS """
# import numpy as np
# import sys
# sys.path.append('/Users/ghkerr/phd/globalo3/')
# import globalo3_open, globalo3_calculate, observations_open
# datapath = '/Users/ghkerr/phd/globalo3/data/parsed/'      
# import pandas as pd
# import netCDF4 as nc
# o3_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['O3_control'][:].data
# lat_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['lat'][:].data
# lng_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['lng'][:].data
# # Interpolate to the resolution of GEOSChem 
# o3_gmi = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_gmi, lng_gmi, o3_gmi, checkplot='yes')
# # FIGURE OUT WHAT OCEAN DOESN"T WORK FOR GEOSCHEMM!!
# land = globalo3_calculate.find_grid_overland(lat_gmi[:-5], lng_gmi)
# wherenan = np.where(land != 1)
# ocean = np.empty(shape=land.shape)
# ocean[:] = np.nan
# ocean[wherenan] = 1.
# land_gc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_gmi[:-5], lng_gmi, land, checkplot='no')
# ocean_gc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_gmi[:-5], lng_gmi, ocean, checkplot='no')
# import numpy as np
# import matplotlib as mpl
# mpl.rcParams['hatch.linewidth'] = 0.3     
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
# o3_zm_land = np.nanmean((o3_gmi*land_gc), axis=tuple((0,-1)))
# o3_zm_ocean = np.nanmean((o3_gmi*ocean_gc), axis=tuple((0,-1)))
# grad_o3_zm_land = np.gradient(o3_zm_land)
# grad_o3_zm_ocean = np.gradient(o3_zm_ocean)
# r_o3jet = globalo3_calculate.calculate_r(o3_gmi, edj_dist, lat_gc, 
#     lng_gc)
# significance_r_o3jet = \
#     globalo3_calculate.calculate_r_significance(o3_gmi, edj_dist, 
#     r_o3jet, lat_gc, lng_gc)
# TRAC_40_50, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, months, 
#     'SpeciesConc_TRAC50_40_50', latmin, latmax, lngmin, lngmax, 800., 956.)
# TRAC_40_50 = np.nanmean(TRAC_40_50, axis=1)*1e6
# TRAC_40_50zm_ocean = np.nanmean((TRAC_40_50*ocean_gc), axis=tuple((0,-1)))
# grad_TRAC_40_50zm_ocean = np.gradient(TRAC_40_50zm_ocean)
# r_TRAC_40_50jet = globalo3_calculate.calculate_r(TRAC_40_50, edj_dist, lat_gc, 
#     lng_gc)
# significance_r_TRAC_40_50jet = \
#     globalo3_calculate.calculate_r_significance(TRAC_40_50, edj_dist, 
#     r_TRAC_40_50jet, lat_gc, lng_gc)
# TRAC_30_40, lat_gc, lng_gc, lev_gc = \
#     tracer_open.open_geoschem_merra2_2x25_RnPbBe(years, months, 
#     'SpeciesConc_TRAC50_30_40', latmin, latmax, lngmin, lngmax, 800., 956.)
# TRAC_30_40 = np.nanmean(TRAC_30_40, axis=1)*1e6
# TRAC_30_40zm_land = np.nanmean((TRAC_30_40*land_gc), axis=tuple((0,-1)))
# grad_TRAC_30_40zm_land = np.gradient(TRAC_30_40zm_land)
# r_TRAC_30_40jet = globalo3_calculate.calculate_r(TRAC_30_40, edj_dist, lat_gc, 
#     lng_gc)
# significance_r_TRAC_30_40jet = \
#     globalo3_calculate.calculate_r_significance(TRAC_30_40, edj_dist, 
#     r_TRAC_30_40jet, lat_gc, lng_gc)
# # Plotting
# fig = plt.figure(figsize=(10,2.5))
# ax1 = plt.subplot2grid((2,3),(0,0)) 
# ax2 = plt.subplot2grid((2,3),(0,1), 
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax3 = plt.subplot2grid((2,3),(0,2), 
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax4 = plt.subplot2grid((2,3),(1,0)) 
# ax5 = plt.subplot2grid((2,3),(1,1), 
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax6 = plt.subplot2grid((2,3),(1,2), 
#     projection=ccrs.PlateCarree(central_longitude=0.))
# cmap = plt.get_cmap('coolwarm')
# lng_gc[-1] = 360. # for wrapping around the Prime Meridian 
# # r(TRAC_40_50, jet) over oceans
# ax2.set_title('$\chi_{40\mathregular{-}50^{\circ}}$', fontsize=16, x=0.02, 
#     ha='left')
# mb = ax2.contourf(lng_gc, lat_gc, r_TRAC_40_50jet,
#     np.linspace(-1, 1., 11), cmap=cmap, extend='neither',
#     transform=ccrs.PlateCarree(), zorder=1)
# ax2.coastlines(lw=0.25, color='k', zorder=3)
# ax2.add_feature(cfeature.LAND, zorder=2, lw=0.0, color='w')
# ax2.contourf(lng_gc, lat_gc, significance_r_TRAC_40_50jet, hatches=['//////'], 
#     colors='none', transform=ccrs.PlateCarree())
# # r(O3, jet) over land
# ax3.set_title('O$_{\mathregular{3}}$', fontsize=10, x=0.02, 
#     ha='left')
# mb = ax3.contourf(lng_gc, lat_gc, r_o3jet,
#     np.linspace(-1, 1, 11), cmap=cmap, extend='neither',
#     transform=ccrs.PlateCarree(), zorder=1)
# ax3.contourf(lng_gc, lat_gc, significance_r_o3jet, hatches=['//////'], 
#     colors='none', transform=ccrs.PlateCarree())
# ax3.coastlines(lw=0.25, color='k', zorder=3)
# ax3.add_feature(cfeature.LAND, zorder=2, lw=0.0, color='w')
# # Latitudinal gradient of O3 and TRAC_40_50 over oceans
# ax1.plot(grad_o3_zm_ocean, lat_gc, '-k')
# [t.set_fontsize(9) for t in ax1.xaxis.get_ticklabels()]
# # ax1.set_xlabel('dO$_{3}$/d$\phi$', fontsize=10)
# ax1b = ax1.twiny()
# ax1b.plot(grad_TRAC_40_50zm_ocean, lat_gc, '-r')
# [t.set_color('red') for t in ax1b.xaxis.get_ticklabels()]
# [t.set_fontsize(9) for t in ax1b.xaxis.get_ticklabels()]
# pos1 = ax2.get_position() # get the original position 
# pos1 = [pos1.x0-0.15, pos1.y0,  pos1.width/2, pos1.height] 
# ax1.set_position(pos1) # set a new position
# ax1.set_ylim([0, 90])
# ax1.set_yticks([0, 30, 60, 90])
# ax1.set_ylabel('Latitude [$^\circ$N]')
# # r(TRAC_30-40, jet) over land
# ax5.set_title('$\chi_{30\mathregular{-}40^{\circ}}$', fontsize=16, x=0.02, 
#     ha='left')
# mb = ax5.contourf(lng_gc, lat_gc, r_TRAC_30_40jet,
#     np.linspace(-1, 1, 11), cmap=cmap, extend='neither',
#     transform=ccrs.PlateCarree(), zorder=1)
# ax5.contourf(lng_gc, lat_gc, significance_r_TRAC_30_40jet, hatches=['//////'], 
#     colors='none', transform=ccrs.PlateCarree())
# ax5.coastlines(lw=0.25, color='k', zorder=3)
# ax5.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='w')
# # r(O3, jet) over land
# ax6.set_title('O$_{\mathregular{3}}$', fontsize=10, x=0.02, 
#     ha='left')
# mb = ax6.contourf(lng_gc, lat_gc, r_o3jet,
#     np.linspace(-1, 1, 11), cmap=cmap, extend='neither',
#     transform=ccrs.PlateCarree(), zorder=1)
# ax6.contourf(lng_gc, lat_gc, significance_r_o3jet, hatches=['//////'], 
#     colors='none', transform=ccrs.PlateCarree())
# ax6.coastlines(lw=0.25, color='k', zorder=3)
# ax6.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='w')
# colorbar_axes = plt.gcf().add_axes([0.92, ax6.get_position().y0, 
#     0.02, (ax3.get_position().y1-ax6.get_position().y0)]) 
# colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
# colorbar.ax.tick_params(labelsize=12)
# colorbar.set_label('[$\cdot$]', fontsize=10)
# ax4.plot(grad_o3_zm_land, lat_gc, '-k')
# [t.set_fontsize(9) for t in ax4.xaxis.get_ticklabels()]
# # ax4.set_xlabel('dO$_{3}$/d$\phi$', fontsize=10)
# ax4b = ax4.twiny()
# ax4b.plot(grad_TRAC_30_40zm_land, lat_gc, '-r')
# [t.set_color('red') for t in ax4b.xaxis.get_ticklabels()]
# [t.set_fontsize(9) for t in ax4b.xaxis.get_ticklabels()]
# ax4b.xaxis.set_label_coords(1.4, 1.35)
# pos1 = ax5.get_position() # get the original position 
# pos1 = [pos1.x0-0.15, pos1.y0,  pos1.width/2, pos1.height] 
# ax4.set_position(pos1) # set a new position
# ax4.set_ylim([0, 90])
# ax4.set_yticks([0, 30, 60, 90])
# ax4.set_ylabel('Latitude [$^\circ$N]')
# plt.savefig('/Users/ghkerr/Desktop/hihigaige.png', dpi=300)

""" CO50, NH50, AND TRAC_30_40 OVER LAND VERSUS OCEANS"""
# # Tracer with 30-40 deg source most similar to CO_50 (at least regarding)
# # where concentrations peak; i.e., uncomment for visual inspection
# plt.plot(lat_gc, np.nanmean(co_50_gc, axis=tuple((0,1,3))));
# plt.twinx();plt.plot(lat_gc, np.nanmean(TRAC_30_40, axis=tuple((0,1,3))))
# # Create land and ocean masks
# land = globalo3_calculate.find_grid_overland(lat_gc[:-3], lng_gc)
# # Since urcrnrlat must be between -90. and 90. degrees in function 
# # find_grid_overland, and this raises a ValueError without stripping off 
# # the final few latitude circles from the latitude grid. Add NaNs back
# # to preserve dimensionality
# add_nan = np.empty(shape=(lat_gc.shape[0]-land.shape[0], land.shape[1]))
# add_nan[:] = np.nan
# land = np.vstack([land, add_nan])
# wherenan = np.where(land != 1)
# ocean = np.empty(shape=land.shape)
# ocean[:] = np.nan
# ocean[wherenan] = 1.
# # Jet-tracer relationship  
# r_TRAC_30_40jet = globalo3_calculate.calculate_r(np.nanmean(TRAC_30_40, axis=1), 
#     edj_dist, lat_gc, lng_gc)
# r_co_50jet = globalo3_calculate.calculate_r(np.nanmean(co_50_gc, axis=1), 
#     edj_dist, lat_gc, lng_gc)
# r_nh_50jet = globalo3_calculate.calculate_r(np.nanmean(nh_50_replay, axis=1), 
#     edj_dist, lat_gc, lng_gc)
# # Plotting
# fig = plt.figure()
# ax1 = plt.subplot2grid((2,1),(0,0))
# ax2 = plt.subplot2grid((2,1),(1,0))
# for ax in [ax1, ax2]:
#     ax.set_xlim([0, 89])
#     ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
#     ax.set_xticklabels(['0', '', '30', '', '60', '', '90'])
#     ax.set_ylabel('r($\chi$, $\phi_{jet} - \phi$)')
# # For ocean 
# ax1.set_title('(a) Ocean')
# ax1.plot(lat_gc, np.nanmean(r_TRAC_30_40jet*ocean, axis=-1), color='#1b9e77')
# ax1.plot(lat_gc, np.nanmean(r_co_50jet*ocean, axis=-1), color='#d95f02')
# ax1.plot(lat_gc, np.nanmean(r_nh_50jet*ocean, axis=-1), color='#7570b3')
# # For land
# ax2.set_title('(b) Land')
# ax2.plot(lat_gc, np.nanmean(r_TRAC_30_40jet*land, axis=-1), 
#     color='#1b9e77', label='$\chi_{30\mathregular{-}40^{\circ}}$')
# ax2.plot(lat_gc, np.nanmean(r_co_50jet*land, axis=-1), color='#d95f02',
#     label='CO$_{50}$')
# ax2.plot(lat_gc, np.nanmean(r_nh_50jet*land, axis=-1), color='#7570b3', 
#     label='NO$_{50}$')
# ax2.set_xlabel('Latitude [$^{\circ}$]')
# plt.subplots_adjust(hspace=0.45)
# ax2.legend(loc=3, ncol=3, frameon=False)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'zonalavg_r_trac_30_40_co_50_nh_50jet.png', dpi=300)

"""STATIONARY VERSUS TRANSIENT EDDIES FROM DIFFERENT TRACERS"""
# n.b., the original plots weren't mass streamfunctions; they were vertically
# integrated ppm m s-1 from 985-800hPa
# # Separate into PW and EW days
# Vcolumn_gc_eqjet, Vcolumn_gc_pwjet, tracer_eqjet, tracer_pwjet = \
#     globalo3_calculate.sortfield_byjetlat_column(Vcolumn_gc, TRAC_10_20, 
#     edj, lng_gc, lat_gc, lev_gc, psize=30) 
# labels = ['$\chi_{0-10^{\circ}}$', '$\chi_{10-20^{\circ}}$', 
#     '$\chi_{20-30^{\circ}}$', '$\chi_{30-40^{\circ}}$', 
#     '$\chi_{40-50^{\circ}}$', '$\chi_{50-60^{\circ}}$',
#     '$\chi_{60-70^{\circ}}$', '$\chi_{70-80^{\circ}}$',
#     '$\chi_{80-90^{\circ}}$']
# colors = ['#081d58','#253494','#225ea8','#1d91c0','#41b6c4','#7fcdbb',
#           '#c7e9b4','#edf8b1','#ffffd9']
# fig = plt.figure()
# ax = plt.subplot2grid((1,1),(0,0))
# # Loop through tracers 
# for i, tracer in enumerate([TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
#     TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90]):
#     # 
#     Vcolumn_gc_eqjet, Vcolumn_gc_pwjet, tracer_eqjet, tracer_pwjet = \
#         globalo3_calculate.sortfield_byjetlat_column(Vcolumn_gc, tracer, 
#         edj, lng_gc, lat_gc, lev_gc, psize=30)     
#     print(i)
#     stationary_all = []
#     stationary_pw = []
#     stationary_ew = []
#     trans_all = []
#     trans_pw = []
#     trans_ew = []
#     for lev in np.arange(0, len(lev_gc), 1):
#         species_mean, species_stationary, species_transient, species_total = \
#             globalo3_calculate.meridional_flux(Vcolumn_gc[:,lev], tracer[:,lev], 
#             np.arange(0, len(tracer), 1), lat_gc, lng_gc)
#         (species_mean_pwjet, species_stationary_pwjet, species_transient_pwjet, 
#             species_total_pwjet) = globalo3_calculate.meridional_flux(
#             Vcolumn_gc_pwjet[:,lev], tracer_pwjet[:,lev], 
#             np.arange(0, len(tracer_pwjet), 1), lat_gc, lng_gc)
#         (species_mean_eqjet, species_stationary_eqjet, species_transient_eqjet, 
#             species_total_eqjet) = globalo3_calculate.meridional_flux(
#             Vcolumn_gc_eqjet[:,lev], tracer_eqjet[:,lev], 
#             np.arange(0, len(tracer_eqjet), 1), lat_gc, lng_gc)
#         stationary_all.append(species_stationary)
#         stationary_pw.append(species_stationary_pwjet)
#         stationary_ew.append(species_stationary_eqjet)        
#         trans_all.append(species_transient)
#         trans_pw.append(species_transient_pwjet)
#         trans_ew.append(species_transient_eqjet)
#     stationary_all = np.nansum(np.vstack(stationary_all), axis=0)
#     stationary_pw = np.nansum(np.vstack(stationary_pw), axis=0)
#     stationary_ew = np.nansum(np.vstack(stationary_ew), axis=0)
#     trans_all = np.nansum(np.vstack(trans_all), axis=0)
#     trans_pw = np.nansum(np.vstack(trans_pw), axis=0)
#     trans_ew = np.nansum(np.vstack(trans_ew), axis=0)
#     # Plotting    
#     ax.plot(lat_gc, stationary_all*1e6,  color=colors[i], ls='-', lw=1., 
#             label=labels[i])
#     ax.set_title('Stationary Eddies')
#     ax.set_xlim([0, 90])
#     ax.set_xlabel('Latitude [$^{\circ}$N]')
#     ax.set_ylabel('Flux [ppm m s$^{-1}$]')
# ax.legend(ncol=3, loc=4)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'zonalavg_stationaryeddyflux_tracers_geoschem_985-800.png', dpi=300)
# plt.show()


lng = lng_gc
lat = lat_gc
V = Vcolumn_gc


# To generate colormap, I found these three "midpoint" colors from 
# https://learnui.design/tools/data-color-picker.html#palette: 
#003f5c, #bc5090, #ffa600
# Then I created a nine color palette using 
# https://vis4.net/palettes/#/9|s|00429d,96ffea,ffffe0|ffffe0,ff005e,93003a|1|1
COLORS = ['#003f5c', '#434766', '#6a4f6d', '#8b586f', '#a9636c', 
    '#c47064', '#db7f56', '#ef913e', '#ffa600']
LABELS = ['$\chi_{0-10^{\circ}}$', '$\chi_{10-20^{\circ}}$', 
    '$\chi_{20-30^{\circ}}$', '$\chi_{30-40^{\circ}}$', 
    '$\chi_{40-50^{\circ}}$', '$\chi_{50-60^{\circ}}$',
    '$\chi_{60-70^{\circ}}$', '$\chi_{70-80^{\circ}}$',
    '$\chi_{80-90^{\circ}}$']



import numpy as np
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.3  
mpl.rcParams['xtick.major.width'] = 1 
mpl.rcParams['ytick.major.width'] = 1   
mpl.rcParams['axes.spines.bottom'] = 1
mpl.rcParams['axes.spines.left'] = 1
mpl.rcParams['axes.spines.right'] = 1
mpl.rcParams['axes.spines.top'] = 1
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter  


def fig2():
    """
    """    
    if lng[-1] != 360:
        lng[-1] = 360.    
    fig = plt.figure(figsize=(10,6))
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=1, rowspan=3)
    ax2 = plt.subplot2grid((3,3), (0,1), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax3 = plt.subplot2grid((3,3), (1,1), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax4 = plt.subplot2grid((3,3), (2,1), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    # Zonally-averaged tracer concentrations
    ax1.set_title('(a)', fontsize=12, x=0.02, ha='left')
    for i, tracer in enumerate([TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
        TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90]):
        # Zonal, time, and vertical average; convert to ppm
        tracer = np.nanmean(tracer, axis=tuple((0,1,3)))*1e6
        ax1.plot(tracer, lat, color=COLORS[i], label=LABELS[i], lw=2)
    ax1.set_xlim([0, 1.6])
    ax1.set_xticks([0, 0.4, 0.8, 1.2, 1.6])
    ax1.set_xticklabels(['0.0', '', '0.8', '', '1.6'])
    ax1.set_ylim([0, lat.max()])
    ax1.set_yticks([0, 15, 30, 45, 60, 75, lat.max()])
    ax1.set_yticklabels(['0$^{\circ}$', '15$^{\circ}$N', '30$^{\circ}$N',
        '45$^{\circ}$N', '60$^{\circ}$N', '75$^{\circ}$N', '90$^{\circ}$N'])
    ax1.set_xlabel('[ppm]', fontsize=12)
    ax1.tick_params(which='major', labelsize=9)
    # Add Legend 
    ax1.legend(bbox_to_anchor=[3.6, -0.1, 0, 0], ncol=5, frameon=False, 
        fontsize=12)
    # # Distribution of X10-20
    ax2.set_title('(b) $\chi_{10-20^{\circ}}$', fontsize=12, x=0.02, ha='left')
    mb = ax2.contourf(lng, lat, np.nanmean(TRAC_10_20, axis=tuple((0,1)))*1e6, 
        np.linspace(0., 1.6, 9), cmap=plt.get_cmap('pink_r'), extend='max',
        transform=ccrs.PlateCarree(), zorder=2)
    # Distribution of X40-50
    ax3.set_title('(c) $\chi_{40-50^{\circ}}$', fontsize=12, x=0.02, ha='left')
    mb = ax3.contourf(lng, lat, np.nanmean(TRAC_40_50, axis=tuple((0,1)))*1e6, 
        np.linspace(0., 1.6, 9), cmap=plt.get_cmap('pink_r'), extend='max',
        transform=ccrs.PlateCarree(), zorder=2)
    # Distribution of X70-80
    ax4.set_title('(d) $\chi_{70-80^{\circ}}$', fontsize=12, x=0.02, ha='left')
    mb = ax4.contourf(lng, lat, np.nanmean(TRAC_70_80, axis=tuple((0,1)))*1e6, 
        np.linspace(0., 1.6, 9), cmap=plt.get_cmap('pink_r'), extend='max',
        transform=ccrs.PlateCarree(), zorder=2)
    # Aesthetics
    for ax in [ax2, ax3, ax4]:
        # Add eddy-driven jet
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(edj, axis=0)[::skiplng], 
            yerr=np.nanstd(edj, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        ax.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        ax.set_extent([lng.min()-180., lng.max()-180., lat.min(), lat.max()])
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)         
        ax.get_xaxis().set_ticklabels([])
        ax.set_yticks([0, 30, 60, 90], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()    
        ax.yaxis.set_major_formatter(lat_formatter)    
        ax.tick_params(which='major', labelsize=8)
    ax4.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
    lng_formatter = LongitudeFormatter()
    ax4.xaxis.set_major_formatter(lng_formatter)       
    ax4.tick_params(which='major', labelsize=9)
    plt.subplots_adjust(top=0.95, right=0.85, bottom=0.18)
    # Add colorbar
    cbaxes = fig.add_axes([ax4.get_position().x1+0.03, ax4.get_position().y0, 
        0.02, ax2.get_position().y1-ax4.get_position().y0]) 
    cb = plt.colorbar(mb, cax=cbaxes, orientation='vertical')
    cb.set_label(label='[ppm]', size=12)
    cb.ax.tick_params(labelsize=9)
    # Move ax1 such that its top/bottom line up with ax2 and ax4
    pos1 = ax1.get_position() # get the original position 
    pos1 = [pos1.x0-0.04, ax4.get_position().y0, pos1.width, 
        (ax2.get_position().y1-ax4.get_position().y0)] 
    ax1.set_position(pos1) # set a new position
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig2.png', dpi=500)
    return 

def fig3(): 
    """
    """  
    if lng[-1] != 360:
        lng[-1] = 360.    
    fig = plt.figure(figsize=(7,6))
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2 = plt.subplot2grid((3,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax3 = plt.subplot2grid((3,2), (2,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title('(a) $\chi_{10-20^{\circ}}$', fontsize=12, x=0.02, ha='left')
    ax2.set_title('(b) $\chi_{40-50^{\circ}}$', fontsize=12, x=0.02, ha='left')
    ax3.set_title('(c) $\chi_{70-80^{\circ}}$', fontsize=12, x=0.02, ha='left')
    # Loop through axes and plot
    axes = [ax1, ax2, ax3]
    for i, tracer in enumerate([TRAC_10_20, TRAC_40_50, TRAC_70_80]):
        # Calculate tracer-jet correlation and significance 
        r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
            edj_dist, lat, lng)
        significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
            np.nanmean(tracer, axis=1), edj_dist, r_tracerjet, lat, lng)
        # Find (PW - EW) jet composites
        eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
            eqjet_tracer = globalo3_calculate.segregate_field_bylat(
            np.nanmean(tracer, axis=1), lng_gc, edj, np.arange(0, len(tracer), 1))
        mb = axes[i].contourf(lng, lat, (pwjet_tracer-eqjet_tracer)*1e6,
            np.linspace(-0.2, 0.2, 11), cmap=plt.get_cmap('coolwarm'), 
            extend='both', transform=ccrs.PlateCarree(), zorder=2)
        # Add hatching for significance
        axes[i].contourf(lng, lat, significance_r_tracerjet, 
            hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
            zorder= 4)
        # Add eddy-driven jet
        skiplng = 5
        axes[i].errorbar(lng[::skiplng], np.nanmean(edj, axis=0)[::skiplng], 
            yerr=np.nanstd(edj, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., lat.min(), lat.max()])
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        axes[i].xaxis.set_major_formatter(lng_formatter)         
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].set_yticks([0, 30, 60, 90], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()    
        axes[i].yaxis.set_major_formatter(lat_formatter)    
        axes[i].tick_params(which='major', labelsize=8)
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)       
    ax3.tick_params(which='major', labelsize=9)
    plt.subplots_adjust(left=0.05, right=0.85, hspace=0.3)
    # Add colorbar
    cbaxes = fig.add_axes([ax1.get_position().x1+0.03, ax3.get_position().y0, 
        0.02, ax1.get_position().y1-ax3.get_position().y0]) 
    cb = plt.colorbar(mb, cax=cbaxes, orientation='vertical')
    cb.set_label(label='[ppm]', size=12)
    cb.set_ticks(np.linspace(-0.2, 0.2, 11))
    cb.ax.tick_params(labelsize=9)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig3.png', dpi=500)
    return

def fig4():
    """
    """
    import scipy.interpolate
    fig = plt.figure(figsize=(7,6))
    ax1 = plt.subplot2grid((2,1), (0,0))
    ax2 = plt.subplot2grid((2,1), (1,0))
    ax1.set_title('(a) r(V, $\phi_{jet}$)', fontsize=12, x=0.02, ha='left')
    ax2.set_title('(b) r($\chi$, $\phi_{jet}$)', fontsize=12, x=0.02, ha='left')
    # Draw y = 0 (no correlation) line
    ax1.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='grey')
    ax2.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='grey')
    # Determine near-surface meridional wind-jet correlation
    r_Vjet = globalo3_calculate.calculate_r(np.nanmean(V, axis=1), edj_dist, 
        lat, lng)
    ax1.plot(lat, np.nanmean(r_Vjet, axis=1), '-k', lw=2)
    # Indicate jet; errorbars represent the zonally-averaged standard deviation 
    # of the daily variations in the position of the jet
    ax1.errorbar(np.nanmean(edj, axis=tuple((0,1))), 0., 
        xerr=np.nanmean(np.nanstd(edj, axis=0)), zorder=10, color='k', 
        markersize=5, elinewidth=2, ecolor='k', fmt='o')
    # Loop through tracers 
    for i, tracer in enumerate([TRAC_0_10, TRAC_10_20, TRAC_20_30, TRAC_30_40, 
        TRAC_40_50, TRAC_50_60, TRAC_60_70, TRAC_70_80, TRAC_80_90]):
        tracer = np.nanmean(tracer, axis=1)
        # Calculate tracer-jet correlation
        r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist, lat, 
            lng)
        ax2.plot(lat, np.nanmean(r_tracerjet, axis=1), color=COLORS[i], 
            label=LABELS[i], lw=2)
    # Aesthetics
    for ax in [ax1, ax2]:
        ax.set_xlim([0, lat.max()])
        ax.set_xticks([0, 15, 30, 45, 60, 75, lat.max()])
        ax.set_xticklabels([])
        ax.set_ylabel('[$\cdot$]', fontsize=12) 
        ax.tick_params(which='major', labelsize=9)
    ax1.set_ylim([-0.16, 0.16])
    ax1.set_yticks(np.linspace(-0.16, 0.16, 5))
    ax2.set_ylim([-0.28, 0.28])
    ax2.set_yticks(np.linspace(-0.28, 0.28, 5))
    ax2.set_xticklabels(['0$^{\circ}$', '15$^{\circ}$N', '30$^{\circ}$N',
        '45$^{\circ}$N', '60$^{\circ}$N', '75$^{\circ}$N', '90$^{\circ}$N']) 
    ax2.legend(bbox_to_anchor=[1.13, -0.2, 0, 0], ncol=5, frameon=False, 
        fontsize=12)
    plt.subplots_adjust(hspace=0.3, bottom=0.2)
    # # Interpolate the zonally-averaged r(V, jet) and find where the 
    # # r(V, jet) = 0.; this is kind of kludgey, so the values are hard-coded in
    # interp = scipy.interpolate.interp1d(lat, np.nanmean(r_Vjet, axis=-1))
    # r_Vjet_zm_interp = []
    # lat_interp = np.linspace(0, lat.max(), 500)
    # for x in lat_interp:
    #     r_Vjet_zm_interp.append(interp(x))
    # r_Vjet_zm_interp = np.abs(r_Vjet_zm_interp)
    # # Find index of 4 smallest values since r(V, jet) = 0. four times
    # r_Vjet_zm_interp.argsort()[:4] # however this has some repeats, so increase
    # # to greater values to pull off values 
    # lat_where_0 = lat_interp[[114, 211, 372, 494]]
    lat_where_0 = np.array([20.44689379, 37.84468938, 66.72144289, 88.60320641])
    for x in lat_where_0: 
        ax1.axvline(x=x, c='grey', lw=0.75, ls='--', zorder=0)
        ax2.axvline(x=x, c='grey', lw=0.75, ls='--', zorder=0)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig4.png', dpi=500)
    return