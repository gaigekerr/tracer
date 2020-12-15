#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Conduct data analysis and plot quantities related to GEOSChem Rn–Pb–Be 
tracer simulations 

Revision History
    15102019 -- initial version created
    27022020 -- data analysis to examine the zonal mean jet/tracer relationship
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# To generate colormap, I found these three "midpoint" colors from 
# https://learnui.design/tools/data-color-picker.html#palette: 
#003f5c, #bc5090, #ffa600
# Then I created a nine color palette using 
# https://vis4.net/palettes/#/9|s|00429d,96ffea,ffffe0|ffffe0,ff005e,93003a|1|1
COLORS = ['#003f5c', '#514969', '#80556f', '#a9636c', '#cc7560', 
    '#e98b47', '#ffa600']
LABELS = ['$\chi_{10-20}$', '$\chi_{20-30}$', 
    '$\chi_{30-40}$', '$\chi_{40-50}$', 
    '$\chi_{50-60}$','$\chi_{60-70}$', 
    '$\chi_{70-80}$']
mpl.rcParams['hatch.linewidth'] = 0.3  
mpl.rcParams['xtick.major.width'] = 1 
mpl.rcParams['ytick.major.width'] = 1   
mpl.rcParams['axes.spines.bottom'] = 1
mpl.rcParams['axes.spines.left'] = 1
mpl.rcParams['axes.spines.right'] = 1
mpl.rcParams['axes.spines.top'] = 1

def fig1(lat, lng, TRAC_10_20_jja, TRAC_20_30_jja, TRAC_30_40_jja, 
    TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, TRAC_70_80_jja, edj_jja):
    """
    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]    
    TRAC_10_20_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 10-20 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_20_30_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 20-30 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_30_40_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 30-40 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_40_50_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_50_60_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 50-60 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_60_70_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 60-70 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_70_80_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 70-80 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    edj_jja : numpy.ndarray
        Latitude of the jet stream for JJA 2008-2010, units of degrees north,
        [time, lng]
        
    Returns
    -------
    None
    """    
    if lng[-1] != 360:
        lng[-1] = 360.    
    fig = plt.figure(figsize=(10,6))
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=1, rowspan=3)
    ax2 = plt.subplot2grid((3,3), (0,1), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax3 = plt.subplot2grid((3,3), (1,1), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax4 = plt.subplot2grid((3,3), (2,1), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    # Zonally-averaged tracer concentrations
    ax1.set_title('(a)', fontsize=12, x=0.02, ha='left')
    lw = [2.5, 1, 1, 2.5, 1, 1, 2.5]
    for i, tracer in enumerate([TRAC_10_20_jja, TRAC_20_30_jja, TRAC_30_40_jja, 
        TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, TRAC_70_80_jja]):
        # Zonal, time, and vertical average; convert to ppm
        tracer = np.nanmean(tracer, axis=tuple((0,1,3)))*1e6
        ax1.plot(tracer, lat, color=COLORS[i], label=LABELS[i], 
            lw=lw[i])
    ax1.set_xlim([0, 1.6])
    ax1.set_xticks([0, 0.4, 0.8, 1.2, 1.6])
    ax1.set_xticklabels(['0.0', '', '0.8', '', '1.6'])  
    ax1.set_ylim([10, 80])
    ax1.set_yticks([10, 20, 30, 40, 50, 60, 70, 80])
    ax1.set_yticklabels(['10$^{\circ}$N', '20$^{\circ}$N', '30$^{\circ}$N',
        '40$^{\circ}$N', '50$^{\circ}$N', '60$^{\circ}$N', '70$^{\circ}$N',
        '80$^{\circ}$N'])
    ax1.set_xlabel('[ppm]', fontsize=12)
    ax1.tick_params(which='major', labelsize=9)
    # Add Legend 
    ax1.legend(bbox_to_anchor=[3.3, -0.07, 0, 0], ncol=4, frameon=False, 
        fontsize=12)
    # Distribution of X70-80
    ax2.set_title('(b) $\chi_{70-80}$', fontsize=12, x=0.02, ha='left')
    mb = ax2.contourf(lng, lat, np.nanmean(TRAC_70_80_jja, 
        axis=tuple((0,1)))*1e6, np.linspace(0., 1.6, 9), 
        cmap=plt.get_cmap('pink_r'), extend='max', 
        transform=ccrs.PlateCarree(), zorder=2)
    # Distribution of X40-50
    ax3.set_title('(c) $\chi_{40-50}$', fontsize=12, x=0.02, ha='left')
    mb = ax3.contourf(lng, lat, 
        np.nanmean(TRAC_40_50_jja, axis=tuple((0,1)))*1e6, 
        np.linspace(0., 1.6, 9), cmap=plt.get_cmap('pink_r'), extend='max',
        transform=ccrs.PlateCarree(), zorder=2)  
    # Distribution of X10-20
    ax4.set_title('(d) $\chi_{10-20}$', fontsize=12, x=0.02, ha='left')
    mb = ax4.contourf(lng, lat, 
        np.nanmean(TRAC_10_20_jja, axis=tuple((0,1)))*1e6, 
        np.linspace(0., 1.6, 9), cmap=plt.get_cmap('pink_r'), extend='max',
        transform=ccrs.PlateCarree(), zorder=2)    
    # Aesthetics
    for ax in [ax2, ax3, ax4]:
        # Add eddy-driven jet
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(edj_jja, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_jja, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        ax.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        ax.set_extent([lng.min()-180., lng.max()-180., 10, 80])
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)         
        ax.get_xaxis().set_ticklabels([])
        ax.set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()    
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.tick_params(which='major', labelsize=9)
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
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig1.pdf', dpi=500)
    plt.show()    
    return 

def fig2(lat, lng, TRAC_10_20_jja, TRAC_40_50_jja, TRAC_70_80_jja, 
    TRAC_10_20_djf, TRAC_40_50_djf, TRAC_70_80_djf, edj_dist_jja, 
    edj_dist_djf, edj_jja, edj_djf): 
    """
    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]    
    TRAC_10_20_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 10-20 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_40_50_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_70_80_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 70-80 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_10_20_djf : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 10-20 
        deg north for days during DJF 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_40_50_djf : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during DJF 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_70_80_djf : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 70-80 
        deg north for days during DJF 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]        
    edj_dist_jja : numpy.ndarray
        Distance from the jet stream for JJA 2008-2010, units of degrees, 
        [time, lat, lng]
    edj_dist_djf : numpy.ndarray
        Distance from the jet stream for DJF 2008-2010, units of degrees, 
        [time, lat, lng]    
    edj_jja : numpy.ndarray
        Latitude of the jet stream for JJA 2008-2010, units of degrees north,
        [time, lng]
    edj_djf : numpy.ndarray
        Latitude of the jet stream for DJF 2008-2010, units of degrees north,
        [time, lng]

    Returns
    -------
    None    
    """  
    if lng[-1] != 360:
        lng[-1] = 360.    
    fig = plt.figure(figsize=(11,5.5))
    ax1 = plt.subplot2grid((3,4), (0,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax2 = plt.subplot2grid((3,4), (1,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax3 = plt.subplot2grid((3,4), (2,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax4 = plt.subplot2grid((3,4), (0,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax5 = plt.subplot2grid((3,4), (1,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax6 = plt.subplot2grid((3,4), (2,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))    
    ax1.set_title('(a) JJA $\chi_{70-80}$ (PW$\:$-$\:$EW)', fontsize=12, 
        x=0.02, ha='left')
    ax2.set_title('(c) JJA $\chi_{40-50}$ (PW$\:$-$\:$EW)', fontsize=12, 
        x=0.02, ha='left')
    ax3.set_title('(e) JJA $\chi_{10-20}$ (PW$\:$-$\:$EW)', fontsize=12, 
        x=0.02, ha='left')
    ax4.set_title('(b) DJF $\chi_{70-80}$ (PW$\:$-$\:$EW)', fontsize=12, 
        x=0.02, ha='left')
    ax5.set_title('(d) DJF $\chi_{40-50}$ (PW$\:$-$\:$EW)', fontsize=12, 
        x=0.02, ha='left')
    ax6.set_title('(f) DJF $\chi_{10-20}$ (PW$\:$-$\:$EW)', fontsize=12, 
        x=0.02, ha='left')
    # Loop through JJA axes and plot
    axes = [ax1, ax2, ax3]
    for i, tracer in enumerate([TRAC_70_80_jja, TRAC_40_50_jja, 
        TRAC_10_20_jja]):
        # Calculate tracer-jet correlation and significance 
        r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
            edj_dist_jja, lat, lng)
        significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
            np.nanmean(tracer, axis=1), edj_dist_jja, r_tracerjet, lat, lng)
        # Find (PW - EW) jet composites
        eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
            eqjet_tracer = globalo3_calculate.segregate_field_bylat(
            np.nanmean(tracer, axis=1), lng, edj_jja, 
            np.arange(0, len(tracer), 1))
        mb = axes[i].contourf(lng, lat, (pwjet_tracer-eqjet_tracer)*1e6,
            np.linspace(-0.2, 0.2, 11), cmap=plt.get_cmap('coolwarm'), 
            extend='both', transform=ccrs.PlateCarree(), zorder=2)
        # Add hatching for significance
        axes[i].contourf(lng, lat, significance_r_tracerjet, 
            hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
            zorder= 4)
        # Add eddy-driven jet
        skiplng = 5
        axes[i].errorbar(lng[::skiplng], np.nanmean(edj_jja, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_jja, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()    
        axes[i].yaxis.set_major_formatter(lat_formatter)
        axes[i].tick_params(which='major', labelsize=9)
    # For DJF
    axes = [ax4, ax5, ax6]
    for i, tracer in enumerate([TRAC_70_80_djf, TRAC_40_50_djf, 
        TRAC_10_20_djf]):
        r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
            edj_dist_djf, lat, lng)
        significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
            np.nanmean(tracer, axis=1), edj_dist_djf, r_tracerjet, lat, lng)
        eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
            eqjet_tracer = globalo3_calculate.segregate_field_bylat(
            np.nanmean(tracer, axis=1), lng, edj_djf, 
            np.arange(0, len(tracer), 1))
        mb = axes[i].contourf(lng, lat, (pwjet_tracer-eqjet_tracer)*1e6,
            np.linspace(-0.2, 0.2, 11), cmap=plt.get_cmap('coolwarm'), 
            extend='both', transform=ccrs.PlateCarree(), zorder=2)
        axes[i].contourf(lng, lat, significance_r_tracerjet, 
            hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
            zorder= 4)
        skiplng = 5
        axes[i].errorbar(lng[::skiplng], np.nanmean(edj_djf, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_djf, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        axes[i].xaxis.set_major_formatter(lng_formatter)         
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()    
        axes[i].get_yaxis().set_ticklabels([])
        axes[i].tick_params(which='major', labelsize=9)
    for ax in [ax3, ax6]:
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())    
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)       
        ax.tick_params(which='major', labelsize=9)
    plt.subplots_adjust(left=0.05, right=0.86)
    # Add colorbar
    cbaxes = fig.add_axes([ax4.get_position().x1+0.03, ax6.get_position().y0, 
        0.02, ax4.get_position().y1-ax6.get_position().y0]) 
    cb = plt.colorbar(mb, cax=cbaxes, orientation='vertical')
    cb.set_label(label='[ppm]', size=12)
    cb.set_ticks(np.linspace(-0.2, 0.2, 11))
    cb.ax.tick_params(labelsize=9)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig2.pdf', dpi=500)
    plt.show()    
    return

def fig3(lat, lng, TRAC_10_20_jja, TRAC_20_30_jja, TRAC_30_40_jja, 
    TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, TRAC_70_80_jja, 
    TRAC_10_20_djf, TRAC_20_30_djf, TRAC_30_40_djf, TRAC_40_50_djf, 
    TRAC_50_60_djf, TRAC_60_70_djf, TRAC_70_80_djf, V_jja, V_djf, edj_dist_jja, 
    edj_dist_djf, edj_jja, edj_djf):
    """
    """
    import scipy.interpolate
    fig = plt.figure(figsize=(7.5,8))
    ax1 = plt.subplot2grid((3,1), (0,0))
    ax2 = plt.subplot2grid((3,1), (1,0))
    ax3 = plt.subplot2grid((3,1), (2,0))
    ax1.set_title('(a) JJA r($\chi$, $\phi_{jet}$)', fontsize=12, x=0.02, 
        ha='left')
    ax2.set_title('(b) DJF r($\chi$, $\phi_{jet}$)', fontsize=12, x=0.02, 
        ha='left')
    ax3.set_title('(c) r(V, $\phi_{jet}$)', fontsize=12, x=0.02, 
        ha='left')    
    # Draw y = 0 (no correlation) line
    ax1.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='black')
    ax2.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='black')
    ax3.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='black')    
    # Determine near-surface meridional wind-jet correlation
    r_Vjet = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
        edj_dist_jja, lat, lng)
    ax3.plot(lat, np.nanmean(r_Vjet, axis=1), ls='solid', color='k', lw=2,
        label='JJA')
    r_Vjet_djf = globalo3_calculate.calculate_r(np.nanmean(V_djf, axis=1), 
        edj_dist_djf, lat, lng)
    ax3.plot(lat, np.nanmean(r_Vjet_djf, axis=1), ls='dashdot', color='k', 
        lw=2, label='DJF')
    # Loop through tracers for JJA
    for i, tracer in enumerate([TRAC_10_20_jja, TRAC_20_30_jja, 
        TRAC_30_40_jja, TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, 
        TRAC_70_80_jja]):
        tracer = np.nanmean(tracer, axis=1)
        # Calculate tracer-jet correlation
        r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_jja, lat, 
            lng)
        ax1.plot(lat, np.nanmean(r_tracerjet, axis=1), color=COLORS[i], 
            label=LABELS[i], lw=2)
    # Loop through tracers for DJF
    for i, tracer in enumerate([TRAC_10_20_djf, TRAC_20_30_djf, 
        TRAC_30_40_djf, TRAC_40_50_djf, TRAC_50_60_djf, TRAC_60_70_djf, 
        TRAC_70_80_djf]):
        tracer = np.nanmean(tracer, axis=1)
        # Calculate tracer-jet correlation
        r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_djf, lat, 
            lng)
        ax2.plot(lat, np.nanmean(r_tracerjet, axis=1), color=COLORS[i], 
            label=LABELS[i], lw=2)
    # Aesthetics
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim([10, 80])
        ax.set_xticks([10, 20, 30, 40, 50, 60, 70, 80])
        ax.set_xticklabels([])
        ax.set_ylabel('[$\cdot$]', fontsize=12) 
        ax.tick_params(which='major', labelsize=9)
    ax1.set_ylim([-0.28, 0.28])
    ax1.set_yticks(np.linspace(-0.28, 0.28, 5))  
    ax2.set_ylim([-0.36, 0.36])
    ax2.set_yticks(np.linspace(-0.36, 0.36, 5))
    # ax1.set_ylim([-0.36, 0.36])
    # ax1.set_yticks(np.linspace(-0.36, 0.36, 5))        
    # ax2.set_ylim([-0.44, 0.44])
    # ax2.set_yticks(np.linspace(-0.44, 0.44, 5))
    ax3.set_ylim([-0.2, 0.2])
    ax3.set_yticks(np.linspace(-0.2, 0.2, 5))
    ax3.set_xticklabels(['10$^{\circ}$N', '20$^{\circ}$N', '30$^{\circ}$N',
        '40$^{\circ}$N', '50$^{\circ}$N', '60$^{\circ}$N', '70$^{\circ}$N',
        '80$^{\circ}$N']) 
    ax3.legend(loc=1, ncol=1, frameon=False, fontsize=12, handlelength=3)
    ax2.legend(bbox_to_anchor=[0.95, -1.45, 0, 0], ncol=4, frameon=False, 
        fontsize=12)
    plt.subplots_adjust(hspace=0.3, top=0.95, bottom=0.15)
    # Interpolate the zonally-averaged r(V, jet) and find where the 
    # r(V, jet) = 0.; this is kind of kludgey, so the values are hard-coded in
    interp = scipy.interpolate.interp1d(lat, np.nanmean(r_Vjet, axis=-1))
    interp_djf = scipy.interpolate.interp1d(lat, np.nanmean(r_Vjet_djf, 
        axis=-1))
    r_Vjet_zm_interp = []
    r_Vjet_zm_interp_djf = []
    lat_interp = np.linspace(0, lat.max(), 500)
    for x in lat_interp:
        r_Vjet_zm_interp.append(interp(x))
    r_Vjet_zm_interp = np.abs(r_Vjet_zm_interp)
    for x in lat_interp:
        r_Vjet_zm_interp_djf.append(interp_djf(x))
    r_Vjet_zm_interp_djf = np.abs(r_Vjet_zm_interp_djf)
    # # Find index of 4 smallest values since r(V, jet) = 0. four times
    # r_Vjet_zm_interp.argsort()[:4] # however this has some repeats, so increase
    # # to greater values to pull off values 
    # lat_where_0 = lat_interp[[114, 211, 372, 494]]
    # Indicate jet; errorbars represent the zonally-averaged standard deviation 
    # of the daily variations in the position of the jet
    ax1.errorbar(np.nanmean(edj_jja, axis=tuple((0,1))), 0.0, 
        xerr=np.nanmean(np.nanstd(edj_jja, axis=0)), zorder=10, color='k', 
        markersize=8, elinewidth=3., ecolor='k', fmt='o')
    ax2.errorbar(np.nanmean(edj_djf, axis=tuple((0,1))), 0.0, 
        xerr=np.nanmean(np.nanstd(edj_djf, axis=0)), zorder=10, color='k', 
        markersize=8, elinewidth=3., ecolor='k', fmt='o')    
    lat_where_0 = np.array([20.44689379, 37.84468938, 66.72144289, 
        88.60320641])
    lat_where_0_djf = np.array([6.0981963, 20.8056, 64.2104])
    for x in lat_where_0: 
        ax1.axvline(x=x, c='k', lw=0.75, ls='--', zorder=0)
    for x in lat_where_0_djf: 
        ax2.axvline(x=x, c='k', lw=0.75, ls='--', zorder=0)        
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig3.pdf', dpi=500)
    plt.show()    
    return

def fig4(lat, lng, TRAC_40_50_jja, TRAC_40_50_djf, V_jja, V_djf, edj_dist_jja, 
    edj_dist_djf, edj_jja, edj_djf):
    """

    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]          
    TRAC_40_50_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_40_50_djf : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during DJF 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    V_jja : numpy.ndarray
        Meridional wind for JJA 2008-2010 (note that the levels should match 
        that of the GEOS-Chem tracers), units of m s-1, [time, lev, lat, lng]
    V_djf : numpy.ndarray
        Meridional wind for DJF 2008-2010 (note that the levels should match 
        that of the GEOS-Chem tracers), units of m s-1, [time, lev, lat, lng]
    edj_dist_jja : numpy.ndarray
        Distance from the jet stream for JJA 2008-2010, units of degrees, 
        [time, lat, lng]
    edj_dist_djf : numpy.ndarray
        Distance from the jet stream for DJF 2008-2010, units of degrees, 
        [time, lat, lng]    
    edj_jja : numpy.ndarray
        Latitude of the jet stream for JJA 2008-2010, units of degrees north,
        [time, lng]
    edj_djf : numpy.ndarray
        Latitude of the jet stream for DJF 2008-2010, units of degrees north,
        [time, lng]

    Returns
    -------
    None          
    """    
    from matplotlib.patches import Patch
    if lng[-1] != 360:
        lng[-1] = 360.
    fig = plt.figure(figsize=(11,3.85))
    ax1 = plt.subplot2grid((2,4), (0,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax2 = plt.subplot2grid((2,4), (0,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax3 = plt.subplot2grid((2,4), (1,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax4 = plt.subplot2grid((2,4), (1,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax1.set_title('(a) JJA V (PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax2.set_title('(b) DJF V (PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax3.set_title('(c) JJA r($\chi_{40-50}$, $\phi_{jet}$)', 
        fontsize=12, x=0.02, ha='left')
    ax4.set_title('(d) DJF r($\chi_{40-50}$, $\phi_{jet}$)', 
        fontsize=12, x=0.02, ha='left')
    cmap = plt.get_cmap('coolwarm')
    # Determine near-surface meridional V-jet correlation
    r_Vjet_jja = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
        edj_dist_jja, lat, lng)
    r_Vjet_djf = globalo3_calculate.calculate_r(np.nanmean(V_djf, axis=1), 
        edj_dist_djf, lat, lng)
    # Calculate V-jet significance, and PW/EW composites for 
    # JJA and DJF
    significance_r_Vjet_jja = \
        globalo3_calculate.calculate_r_significance(np.nanmean(V_jja, axis=1),  
        edj_dist_jja, r_Vjet_jja, lat, lng)
    significance_r_Vjet_djf = \
        globalo3_calculate.calculate_r_significance(np.nanmean(V_djf, axis=1),  
        edj_dist_djf, r_Vjet_djf, lat, lng)
    eqjet_lat_jja, eqjet_lat_var_jja, pwjet_lat_jja, pwjet_lat_var_jja, \
        pwjet_V_jja, eqjet_V_jja = globalo3_calculate.segregate_field_bylat(
        np.nanmean(V_jja, axis=1), lng, edj_jja, np.arange(0, len(V_jja), 1))
    eqjet_lat_djf, eqjet_lat_var_djf, pwjet_lat_djf, pwjet_lat_var_djf, \
        pwjet_V_djf, eqjet_V_djf = globalo3_calculate.segregate_field_bylat(
        np.nanmean(V_djf, axis=1), lng, edj_djf, np.arange(0, len(V_djf), 1))    
    # Plot (PW-EW) composites
    clevst = np.linspace(-5, 5, 13)
    mbt = ax1.contourf(lng, lat, (pwjet_V_jja-eqjet_V_jja),
        clevst, cmap=plt.get_cmap('coolwarm'), extend='both', 
        transform=ccrs.PlateCarree(), zorder=2)
    mbt = ax2.contourf(lng, lat, (pwjet_V_djf-eqjet_V_djf),
        clevst, cmap=plt.get_cmap('coolwarm'), extend='both', 
        transform=ccrs.PlateCarree(), zorder=2)
    # Hatching for significance
    ax1.contourf(lng, lat, significance_r_Vjet_jja, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=4)
    ax2.contourf(lng, lat, significance_r_Vjet_djf, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=4)
    # Contours indicating mean V 
    CS = ax1.contour(lng, lat, np.nanmean(V_jja, axis=tuple((0,1))), 
        [-5,5], linewidths=1., colors='k')
    CS = ax2.contour(lng, lat, np.nanmean(V_djf, axis=tuple((0,1))), 
        [-5,5], linewidths=1., colors='k')
    # JJA TRAC_40_50-jet correlation 
    clevsb = np.linspace(-0.5, 0.5, 11)
    r_tracerjet = globalo3_calculate.calculate_r(
        np.nanmean(TRAC_40_50_jja, axis=1), edj_dist_jja, lat, lng)
    # Column- and time-averaged tracer
    tracer = np.nanmean(np.nanmean(TRAC_40_50_jja, axis=1), axis=0)
    # Loop through each longitude and calculate the latitudinal gradient
    grad_tracer = np.empty(shape=tracer.shape)
    grad_tracer[:] = np.nan
    for spine_i in np.arange(0, len(lng), 1):
        # Tracer concentrations for given "spine" of latitudes at a given 
        # longitude
        spine = tracer[:,spine_i]
        grad_spine = np.gradient(spine)
        grad_tracer[:,spine_i] = grad_spine
    grad_tracer = np.array(grad_tracer)
    expected = -1.*r_Vjet_jja*grad_tracer
    # Find where the expected values are positive/negative
    poz = np.empty(grad_tracer.shape)
    poz[:] = np.nan
    where_poz = np.where(expected > 0.)
    poz[where_poz] = 1. 
    neg = np.empty(grad_tracer.shape)
    neg[:] = np.nan
    where_neg = np.where(expected < 0.)
    neg[where_neg] = 1.
    # Plotting
    mbb = ax3.contourf(lng, lat, r_tracerjet, clevsb, cmap=cmap, 
        extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng, lat, poz, hatches=['...'], colors='none', 
        transform=ccrs.PlateCarree())
    ax3.contourf(lng, lat, neg, hatches=['///'], colors='none', 
        transform=ccrs.PlateCarree())
    # DJF TRAC_40_50-jet correlation 
    r_tracerjet = globalo3_calculate.calculate_r(
        np.nanmean(TRAC_40_50_djf, axis=1), edj_dist_djf, lat, lng)
    # Column- and time-averaged tracer
    tracer = np.nanmean(np.nanmean(TRAC_40_50_djf, axis=1), axis=0)
    # Loop through each longitude and calculate the latitudinal gradient
    grad_tracer = np.empty(shape=tracer.shape)
    grad_tracer[:] = np.nan
    for spine_i in np.arange(0, len(lng), 1):
        # Tracer concentrations for given "spine" of latitudes at a given 
        # longitude
        spine = tracer[:,spine_i]
        grad_spine = np.gradient(spine)
        grad_tracer[:,spine_i] = grad_spine
    grad_tracer = np.array(grad_tracer)
    expected = -1.*r_Vjet_djf*grad_tracer
    # Find where the expected values are positive/negative
    poz = np.empty(grad_tracer.shape)
    poz[:] = np.nan
    where_poz = np.where(expected > 0.)
    poz[where_poz] = 1. 
    neg = np.empty(grad_tracer.shape)
    neg[:] = np.nan
    where_neg = np.where(expected < 0.)
    neg[where_neg] = 1.
    # Plotting
    mbb = ax4.contourf(lng, lat, r_tracerjet, clevsb, cmap=cmap, 
        extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax4.contourf(lng, lat, poz, hatches=['...'], colors='none', 
        transform=ccrs.PlateCarree())
    ax4.contourf(lng, lat, neg, hatches=['///'], colors='none', 
        transform=ccrs.PlateCarree())
    # Denote parallels and meridians
    axes = [ax1, ax2, ax3, ax4]
    for i in np.arange(0, len(axes), 1): 
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        axes[i].xaxis.set_major_formatter(lng_formatter)         
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        axes[i].yaxis.set_major_formatter(lat_formatter)
        axes[i].get_yaxis().set_ticklabels([])    
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
        axes[i].outline_patch.set_zorder(20)
        # Add eddy-driven jet
    for ax in [ax1, ax3]:
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(edj_jja, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_jja, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())    
    for ax in [ax2, ax4]:
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(edj_djf, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_djf, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())    
    # Label parallels and meridians, if appropriate
    axes = [ax1, ax3]
    for i in np.arange(0, len(axes), 1):
        axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        axes[i].yaxis.set_major_formatter(lat_formatter)
        axes[i].tick_params(which='major', labelsize=9)
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)       
    ax3.tick_params(which='major', labelsize=9)
    ax4.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
    lng_formatter = LongitudeFormatter()
    ax4.xaxis.set_major_formatter(lng_formatter)       
    ax4.tick_params(which='major', labelsize=9)
    plt.subplots_adjust(left=0.05, right=0.86)
    # Add custom legend for bottom row
    patch_1 = Patch(fill=False, 
        label='sgn[r($\chi_{40-50}$, $\phi_{jet}$)] > 0', hatch='...', 
        linewidth=0.5)
    patch_2 = Patch(fill=False, 
        label='sgn[r($\chi_{40-50}$, $\phi_{jet}$)] < 0', hatch='///', 
        linewidth=0.5)
    leg = ax3.legend(handles=[patch_1, patch_2], ncol=2, frameon=False, 
        bbox_to_anchor=[-0.04, -0.16, 0, 0], fontsize=12)
    for patch in leg.get_patches():
        patch.set_height(12)
        patch.set_width(28)
    # Add colorbar for V (PW - EW) plots
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.02, 
        ax2.get_position().y0, 0.012, 
        (ax2.get_position().y1-ax2.get_position().y0)]) 
    colorbar = plt.colorbar(mbt, colorbar_axes, orientation='vertical', 
        ticks=clevst[::3], extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[m s$^{-1}$]', fontsize=16)
    # Add colorbar for correlation plots 
    colorbar_axes = plt.gcf().add_axes([ax4.get_position().x1+0.02, 
        ax4.get_position().y0, 0.012, 
        (ax4.get_position().y1-ax4.get_position().y0)]) 
    colorbar = plt.colorbar(mbb, colorbar_axes, orientation='vertical', 
        ticks=clevsb[::2], extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\cdot$]', fontsize=16)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
        'fig4.pdf', dpi=500)
    plt.show()
    return

def figS2(lat, lng, TRAC_10_20_jja, TRAC_20_30_jja, TRAC_30_40_jja, 
    TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, TRAC_70_80_jja, 
    TRAC_10_20_djf, TRAC_20_30_djf, TRAC_30_40_djf, TRAC_40_50_djf, 
    TRAC_50_60_djf, TRAC_60_70_djf, TRAC_70_80_djf, V_jja, V_djf, edj_dist_jja, 
    edj_dist_djf, edj_jja, edj_djf):
    """
    """
    import scipy.interpolate
    fig = plt.figure(figsize=(7.5,8))
    ax1 = plt.subplot2grid((3,1), (0,0))
    ax2 = plt.subplot2grid((3,1), (1,0))
    ax3 = plt.subplot2grid((3,1), (2,0))
    ax1.set_title('(a) JJA r($\chi$, $\phi_{jet}$)', fontsize=12, x=0.02, 
        ha='left')
    ax2.set_title('(b) DJF r($\chi$, $\phi_{jet}$)', fontsize=12, x=0.02, 
        ha='left')
    ax3.set_title('(c) r(V, $\phi_{jet}$)', fontsize=12, x=0.02, 
        ha='left')    
    # Draw y = 0 (no correlation) line
    ax1.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='black')
    ax2.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='black')
    ax3.axhline(y=0, xmin=0, xmax=90, lw=0.75, ls='--', color='black')    
    # Determine near-surface meridional wind-jet correlation
    r_Vjet = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
        edj_dist_jja, lat, lng)
    ax3.plot(lat, np.nanmean(r_Vjet, axis=1), ls='solid', color='k', lw=2,
        label='JJA')
    r_Vjet_djf = globalo3_calculate.calculate_r(np.nanmean(V_djf, axis=1), 
        edj_dist_djf, lat, lng)
    ax3.plot(lat, np.nanmean(r_Vjet_djf, axis=1), ls='dashdot', color='k', 
        lw=2, label='DJF')
    # Loop through tracers for JJA
    for i, tracer in enumerate([TRAC_10_20_jja, TRAC_20_30_jja, 
        TRAC_30_40_jja, TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, 
        TRAC_70_80_jja]):
        tracer = np.nanmean(tracer, axis=1)
        # Calculate tracer-jet correlation
        r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_jja, lat, 
            lng)
        ax1.plot(lat, np.nanmean(r_tracerjet, axis=1), color=COLORS[i], 
            label=LABELS[i], lw=2)
    # Loop through tracers for DJF
    for i, tracer in enumerate([TRAC_10_20_djf, TRAC_20_30_djf, 
        TRAC_30_40_djf, TRAC_40_50_djf, TRAC_50_60_djf, TRAC_60_70_djf, 
        TRAC_70_80_djf]):
        tracer = np.nanmean(tracer, axis=1)
        # Calculate tracer-jet correlation
        r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_djf, lat, 
            lng)
        ax2.plot(lat, np.nanmean(r_tracerjet, axis=1), color=COLORS[i], 
            label=LABELS[i], lw=2)
    # Aesthetics
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim([10, 80])
        ax.set_xticks([10, 20, 30, 40, 50, 60, 70, 80])
        ax.set_xticklabels([])
        ax.set_ylabel('[$\cdot$]', fontsize=12) 
        ax.tick_params(which='major', labelsize=9)
    ax1.set_ylim([-0.36, 0.36])
    ax1.set_yticks(np.linspace(-0.36, 0.36, 5))        
    ax2.set_ylim([-0.44, 0.44])
    ax2.set_yticks(np.linspace(-0.46, 0.46, 5))
    ax3.set_ylim([-0.14, 0.14])
    ax3.set_yticks(np.linspace(-0.14, 0.14, 5))
    ax3.set_xticklabels(['10$^{\circ}$N', '20$^{\circ}$N', '30$^{\circ}$N',
        '40$^{\circ}$N', '50$^{\circ}$N', '60$^{\circ}$N', '70$^{\circ}$N',
        '80$^{\circ}$N']) 
    ax3.legend(loc=1, ncol=1, frameon=False, fontsize=12, handlelength=3)
    ax2.legend(bbox_to_anchor=[0.95, -1.45, 0, 0], ncol=4, frameon=False, 
        fontsize=12)
    plt.subplots_adjust(hspace=0.3, top=0.95, bottom=0.15)
    # Interpolate the zonally-averaged r(V, jet) and find where the 
    # r(V, jet) = 0.; this is kind of kludgey, so the values are hard-coded in
    interp = scipy.interpolate.interp1d(lat, np.nanmean(r_Vjet, axis=-1))
    interp_djf = scipy.interpolate.interp1d(lat, np.nanmean(r_Vjet_djf, 
        axis=-1))
    r_Vjet_zm_interp = []
    r_Vjet_zm_interp_djf = []
    lat_interp = np.linspace(0, lat.max(), 500)
    for x in lat_interp:
        r_Vjet_zm_interp.append(interp(x))
    r_Vjet_zm_interp = np.abs(r_Vjet_zm_interp)
    for x in lat_interp:
        r_Vjet_zm_interp_djf.append(interp_djf(x))
    r_Vjet_zm_interp_djf = np.abs(r_Vjet_zm_interp_djf)
    # # Find index of 4 smallest values since r(V, jet) = 0. four times
    r_Vjet_zm_interp.argsort()[:4] # however this has some repeats, so increase
    # # to greater values to pull off values 
    # lat_where_0 = lat_interp[[114, 211, 372, 494]]
    # Indicate jet; errorbars represent the zonally-averaged standard deviation 
    # of the daily variations in the position of the jet
    ax1.errorbar(np.nanmean(edj_jja, axis=tuple((0,1))), 0.0, 
        xerr=np.nanmean(np.nanstd(edj_jja, axis=0)), zorder=10, color='k', 
        markersize=8, elinewidth=3., ecolor='k', fmt='o')
    ax2.errorbar(np.nanmean(edj_djf, axis=tuple((0,1))), 0.0, 
        xerr=np.nanmean(np.nanstd(edj_djf, axis=0)), zorder=10, color='k', 
        markersize=8, elinewidth=3., ecolor='k', fmt='o')    
    lat_where_0 = np.array([21.70240480961924, 43.763527054108216, 
        67.08016032064128,])
    lat_where_0_djf = np.array([32.6432865731463, 65.46593186372746])
    for x in lat_where_0: 
        ax1.axvline(x=x, c='k', lw=0.75, ls='--', zorder=0)
    for x in lat_where_0_djf: 
        ax2.axvline(x=x, c='k', lw=0.75, ls='--', zorder=0)        
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'figS2.pdf', dpi=500)
    plt.show()    
    return

def figS3(lat, lng, TRAC_40_50_jja, TRAC_40_50_djf, V_jja, V_djf, edj_dist_jja, 
    edj_dist_djf, edj_jja, edj_djf):
    """

    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]          
    TRAC_40_50_jja : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during JJA 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    TRAC_40_50_djf : numpy.ndarray
        GEOS-Chem CO-like tracer for given pressure levels (n.b., column 
        average is performed in function) with prescribed emissions from 40-50 
        deg north for days during DJF 2008-2010, units of volume mixing ratio, 
        [time, lev, lat, lng]
    V_jja : numpy.ndarray
        Meridional wind for JJA 2008-2010 (note that the levels should match 
        that of the GEOS-Chem tracers), units of m s-1, [time, lev, lat, lng]
    V_djf : numpy.ndarray
        Meridional wind for DJF 2008-2010 (note that the levels should match 
        that of the GEOS-Chem tracers), units of m s-1, [time, lev, lat, lng]
    edj_dist_jja : numpy.ndarray
        Distance from the jet stream for JJA 2008-2010, units of degrees, 
        [time, lat, lng]
    edj_dist_djf : numpy.ndarray
        Distance from the jet stream for DJF 2008-2010, units of degrees, 
        [time, lat, lng]    
    edj_jja : numpy.ndarray
        Latitude of the jet stream for JJA 2008-2010, units of degrees north,
        [time, lng]
    edj_djf : numpy.ndarray
        Latitude of the jet stream for DJF 2008-2010, units of degrees north,
        [time, lng]

    Returns
    -------
    None          
    """    
    from matplotlib.patches import Patch
    if lng[-1] != 360:
        lng[-1] = 360.
    fig = plt.figure(figsize=(11,3.85))
    ax1 = plt.subplot2grid((2,4), (0,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax2 = plt.subplot2grid((2,4), (0,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax3 = plt.subplot2grid((2,4), (1,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax4 = plt.subplot2grid((2,4), (1,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax1.set_title('(a) JJA V (PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax2.set_title('(b) DJF V (PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax3.set_title('(c) JJA r($\chi_{40-50}$, $\phi_{jet}$)', 
        fontsize=12, x=0.02, ha='left')
    ax4.set_title('(d) DJF r($\chi_{40-50}$, $\phi_{jet}$)', 
        fontsize=12, x=0.02, ha='left')
    cmap = plt.get_cmap('coolwarm')
    # Determine near-surface meridional V-jet correlation
    r_Vjet_jja = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
        edj_dist_jja, lat, lng)
    r_Vjet_djf = globalo3_calculate.calculate_r(np.nanmean(V_djf, axis=1), 
        edj_dist_djf, lat, lng)
    # Calculate V-jet significance, and PW/EW composites for 
    # JJA and DJF
    significance_r_Vjet_jja = \
        globalo3_calculate.calculate_r_significance(np.nanmean(V_jja, axis=1),  
        edj_dist_jja, r_Vjet_jja, lat, lng)
    significance_r_Vjet_djf = \
        globalo3_calculate.calculate_r_significance(np.nanmean(V_djf, axis=1),  
        edj_dist_djf, r_Vjet_djf, lat, lng)
    eqjet_lat_jja, eqjet_lat_var_jja, pwjet_lat_jja, pwjet_lat_var_jja, \
        pwjet_V_jja, eqjet_V_jja = globalo3_calculate.segregate_field_bylat(
        np.nanmean(V_jja, axis=1), lng, edj_jja, np.arange(0, len(V_jja), 1))
    eqjet_lat_djf, eqjet_lat_var_djf, pwjet_lat_djf, pwjet_lat_var_djf, \
        pwjet_V_djf, eqjet_V_djf = globalo3_calculate.segregate_field_bylat(
        np.nanmean(V_djf, axis=1), lng, edj_djf, np.arange(0, len(V_djf), 1))    
    # Plot (PW-EW) composites
    clevst = np.linspace(-5, 5, 13)
    mbt = ax1.contourf(lng, lat, (pwjet_V_jja-eqjet_V_jja),
        clevst, cmap=plt.get_cmap('coolwarm'), extend='both', 
        transform=ccrs.PlateCarree(), zorder=2)
    mbt = ax2.contourf(lng, lat, (pwjet_V_djf-eqjet_V_djf),
        clevst, cmap=plt.get_cmap('coolwarm'), extend='both', 
        transform=ccrs.PlateCarree(), zorder=2)
    # Hatching for significance
    ax1.contourf(lng, lat, significance_r_Vjet_jja, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=4)
    ax2.contourf(lng, lat, significance_r_Vjet_djf, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=4)
    # Contours indicating mean V 
    CS = ax1.contour(lng, lat, np.nanmean(V_jja, axis=tuple((0,1))), 
        [-5,5], linewidths=1., colors='k')
    CS = ax2.contour(lng, lat, np.nanmean(V_djf, axis=tuple((0,1))), 
        [-5,5], linewidths=1., colors='k')
    # JJA TRAC_40_50-jet correlation 
    clevsb = np.linspace(-0.5, 0.5, 11)
    r_tracerjet = globalo3_calculate.calculate_r(
        np.nanmean(TRAC_40_50_jja, axis=1), edj_dist_jja, lat, lng)
    # Column- and time-averaged tracer
    tracer = np.nanmean(np.nanmean(TRAC_40_50_jja, axis=1), axis=0)
    # Loop through each longitude and calculate the latitudinal gradient
    grad_tracer = np.empty(shape=tracer.shape)
    grad_tracer[:] = np.nan
    for spine_i in np.arange(0, len(lng), 1):
        # Tracer concentrations for given "spine" of latitudes at a given 
        # longitude
        spine = tracer[:,spine_i]
        grad_spine = np.gradient(spine)
        grad_tracer[:,spine_i] = grad_spine
    grad_tracer = np.array(grad_tracer)
    expected = -1.*r_Vjet_jja*grad_tracer
    # Find where the expected values are positive/negative
    poz = np.empty(grad_tracer.shape)
    poz[:] = np.nan
    where_poz = np.where(expected > 0.)
    poz[where_poz] = 1. 
    neg = np.empty(grad_tracer.shape)
    neg[:] = np.nan
    where_neg = np.where(expected < 0.)
    neg[where_neg] = 1.
    # Plotting
    mbb = ax3.contourf(lng, lat, r_tracerjet, clevsb, cmap=cmap, 
        extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng, lat, poz, hatches=['...'], colors='none', 
        transform=ccrs.PlateCarree())
    ax3.contourf(lng, lat, neg, hatches=['///'], colors='none', 
        transform=ccrs.PlateCarree())
    # DJF TRAC_40_50-jet correlation 
    r_tracerjet = globalo3_calculate.calculate_r(
        np.nanmean(TRAC_40_50_djf, axis=1), edj_dist_djf, lat, lng)
    # Column- and time-averaged tracer
    tracer = np.nanmean(np.nanmean(TRAC_40_50_djf, axis=1), axis=0)
    # Loop through each longitude and calculate the latitudinal gradient
    grad_tracer = np.empty(shape=tracer.shape)
    grad_tracer[:] = np.nan
    for spine_i in np.arange(0, len(lng), 1):
        # Tracer concentrations for given "spine" of latitudes at a given 
        # longitude
        spine = tracer[:,spine_i]
        grad_spine = np.gradient(spine)
        grad_tracer[:,spine_i] = grad_spine
    grad_tracer = np.array(grad_tracer)
    expected = -1.*r_Vjet_djf*grad_tracer
    # Find where the expected values are positive/negative
    poz = np.empty(grad_tracer.shape)
    poz[:] = np.nan
    where_poz = np.where(expected > 0.)
    poz[where_poz] = 1. 
    neg = np.empty(grad_tracer.shape)
    neg[:] = np.nan
    where_neg = np.where(expected < 0.)
    neg[where_neg] = 1.
    # Plotting
    mbb = ax4.contourf(lng, lat, r_tracerjet, clevsb, cmap=cmap, 
        extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax4.contourf(lng, lat, poz, hatches=['...'], colors='none', 
        transform=ccrs.PlateCarree())
    ax4.contourf(lng, lat, neg, hatches=['///'], colors='none', 
        transform=ccrs.PlateCarree())
    # Denote parallels and meridians
    axes = [ax1, ax2, ax3, ax4]
    for i in np.arange(0, len(axes), 1): 
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        axes[i].xaxis.set_major_formatter(lng_formatter)         
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        axes[i].yaxis.set_major_formatter(lat_formatter)
        axes[i].get_yaxis().set_ticklabels([])    
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
        axes[i].outline_patch.set_zorder(20)
        # Add eddy-driven jet
    for ax in [ax1, ax3]:
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(edj_jja, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_jja, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())    
    for ax in [ax2, ax4]:
        skiplng = 5
        ax.errorbar(lng[::skiplng], np.nanmean(edj_djf, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_djf, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())    
    # Label parallels and meridians, if appropriate
    axes = [ax1, ax3]
    for i in np.arange(0, len(axes), 1):
        axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        axes[i].yaxis.set_major_formatter(lat_formatter)
        axes[i].tick_params(which='major', labelsize=9)
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)       
    ax3.tick_params(which='major', labelsize=9)
    ax4.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
    lng_formatter = LongitudeFormatter()
    ax4.xaxis.set_major_formatter(lng_formatter)       
    ax4.tick_params(which='major', labelsize=9)
    plt.subplots_adjust(left=0.05, right=0.86)
    # Add custom legend for bottom row
    patch_1 = Patch(fill=False, 
        label='sgn[r($\chi_{40-50}$, $\phi_{jet}$)] > 0', hatch='...', 
        linewidth=0.5)
    patch_2 = Patch(fill=False, 
        label='sgn[r($\chi_{40-50}$, $\phi_{jet}$)] < 0', hatch='///', 
        linewidth=0.5)
    leg = ax3.legend(handles=[patch_1, patch_2], ncol=2, frameon=False, 
        bbox_to_anchor=[-0.04, -0.16, 0, 0], fontsize=12)
    for patch in leg.get_patches():
        patch.set_height(12)
        patch.set_width(28)
    # Add colorbar for V (PW - EW) plots
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.02, 
        ax2.get_position().y0, 0.012, 
        (ax2.get_position().y1-ax2.get_position().y0)]) 
    colorbar = plt.colorbar(mbt, colorbar_axes, orientation='vertical', 
        ticks=clevst[::3], extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[m s$^{-1}$]', fontsize=16)
    # Add colorbar for correlation plots 
    colorbar_axes = plt.gcf().add_axes([ax4.get_position().x1+0.02, 
        ax4.get_position().y0, 0.012, 
        (ax4.get_position().y1-ax4.get_position().y0)]) 
    colorbar = plt.colorbar(mbb, colorbar_axes, orientation='vertical', 
        ticks=clevsb[::2], extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\cdot$]', fontsize=16)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
        'figS3.pdf', dpi=500)
    plt.show()
    return

def figS1():
    """

    Returns
    -------
    None.

    """    
    # 5 day lifetime
    TRAC5_40_50_jja, lat_gc, lng_gc, lev_gg = tracer_open.open_geoschem(years, 
        jja, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC5_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)
    TRAC5_40_50_djf, lat_gc, lng_gc, lev_gg = tracer_open.open_geoschem(years, 
        djf, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC5_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)
    # 25 day lifetime
    TRAC25_40_50_jja, lat_gc, lng_gc, lev_gg = tracer_open.open_geoschem(years, 
        jja, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC25_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)
    TRAC25_40_50_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, 
        djf, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC25_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)
    # 100 day lifetime 
    TRAC100_40_50_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, 
        jja, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC100_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)    
    TRAC100_40_50_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, 
        djf, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC100_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)    
    # 150 day lifetime
    TRAC150_40_50_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, 
        jja, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC150_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)
    TRAC150_40_50_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, 
        djf, 'merra2_2x25_RnPbBe_co50', 'avg24hr', 'SpeciesConc_TRAC150_40_50', 
        latmin, latmax, lngmin, lngmax, pmin, pmax)
    # For wrapping around the Prime Meridian
    if lng[-1] != 360:
        lng[-1] = 360.
    # Initialize figure, axes
    fig = plt.figure(figsize=(11,8))
    ax1 = plt.subplot2grid((4,4), (0,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax2 = plt.subplot2grid((4,4), (0,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax3 = plt.subplot2grid((4,4), (1,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax4 = plt.subplot2grid((4,4), (1,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax5 = plt.subplot2grid((4,4), (2,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax6 = plt.subplot2grid((4,4), (2,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax7 = plt.subplot2grid((4,4), (3,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax8 = plt.subplot2grid((4,4), (3,2), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    # Set titles for JJA
    ax1.set_title('(a) JJA $\chi_{5\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax2.set_title('(b) JJA $\chi_{25\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax3.set_title('(c) JJA $\chi_{100\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax4.set_title('(d) JJA $\chi_{150\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    # Loop through JJA axes and plot
    axes = [ax1, ax2, ax3, ax4]
    for i, tracer in enumerate([TRAC5_40_50_jja, TRAC25_40_50_jja, 
        TRAC100_40_50_jja, TRAC150_40_50_jja]):
        # Calculate tracer-jet correlation and significance 
        r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
            edj_dist_jja, lat, lng)
        significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
            np.nanmean(tracer, axis=1), edj_dist_jja, r_tracerjet, lat, lng)
        # Find (PW - EW) jet composites
        eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
            eqjet_tracer = globalo3_calculate.segregate_field_bylat(
            np.nanmean(tracer, axis=1), lng, edj_jja, 
            np.arange(0, len(tracer), 1))
        mb = axes[i].contourf(lng, lat, (pwjet_tracer-eqjet_tracer)*1e6,
            np.linspace(-0.2, 0.2, 11), cmap=plt.get_cmap('coolwarm'), 
            extend='both', transform=ccrs.PlateCarree(), zorder=2)
        # Add hatching for significance
        axes[i].contourf(lng, lat, significance_r_tracerjet, 
            hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
            zorder= 4)
        # Add eddy-driven jet
        skiplng = 5
        axes[i].errorbar(lng[::skiplng], np.nanmean(edj_jja, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_jja, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        axes[i].xaxis.set_major_formatter(lng_formatter)         
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].get_yaxis().set_ticklabels([])
    # Set titles for DJF
    ax5.set_title('(e) DJF $\chi_{5\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax6.set_title('(f) DJF $\chi_{25\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax7.set_title('(g) DJF $\chi_{100\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    ax8.set_title('(h) DJF $\chi_{150\:\mathregular{d}^{-1}; 40-50}$ '+
        '(PW$\:$-$\:$EW)', fontsize=12, x=0.02, ha='left')
    # Loop through DJF axes and plot
    axes = [ax5, ax6, ax7, ax8]
    for i, tracer in enumerate([TRAC5_40_50_djf, TRAC25_40_50_djf, 
        TRAC100_40_50_djf, TRAC150_40_50_djf]):
        r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
            edj_dist_djf, lat, lng)
        significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
            np.nanmean(tracer, axis=1), edj_dist_djf, r_tracerjet, lat, lng)
        eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_tracer, \
            eqjet_tracer = globalo3_calculate.segregate_field_bylat(
            np.nanmean(tracer, axis=1), lng, edj_djf, 
            np.arange(0, len(tracer), 1))
        mb = axes[i].contourf(lng, lat, (pwjet_tracer-eqjet_tracer)*1e6,
            np.linspace(-0.2, 0.2, 11), cmap=plt.get_cmap('coolwarm'), 
            extend='both', transform=ccrs.PlateCarree(), zorder=2)
        axes[i].contourf(lng, lat, significance_r_tracerjet, 
            hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
            zorder= 4)
        skiplng = 5
        axes[i].errorbar(lng[::skiplng], np.nanmean(edj_djf, axis=0)[::skiplng], 
            yerr=np.nanstd(edj_djf, axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
        axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
        axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        axes[i].xaxis.set_major_formatter(lng_formatter)         
        axes[i].get_xaxis().set_ticklabels([])
        axes[i].get_yaxis().set_ticklabels([])
    # Latitude/longitude labels
    for ax in [ax7, ax8]:
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], 
            crs=ccrs.PlateCarree())    
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)       
        ax.tick_params(which='major', labelsize=9)
    for ax in [ax1, ax3, ax5, ax7]:
        ax.set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()    
        ax.get_yaxis().set_ticklabels([])
        ax.tick_params(which='major', labelsize=9)
    # Add colorbar
    cbaxes = fig.add_axes([ax4.get_position().x1+0.01, 
        (ax8.get_position().y1+ax6.get_position().y0)/2., 
        0.02, (((ax4.get_position().y1+ax2.get_position().y0)/2.)-
        ((ax8.get_position().y1+ax6.get_position().y0)/2.))]) 
    cb = plt.colorbar(mb, cax=cbaxes, orientation='vertical')
    cb.set_label(label='[ppm]', size=12)
    cb.set_ticks(np.linspace(-0.2, 0.2, 11))
    cb.ax.tick_params(labelsize=9)
    plt.subplots_adjust(left=0.03, bottom=0.05, right=0.85, top=0.94)
    plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'figS1.pdf', dpi=500)
    plt.show()
    return

# import numpy as np
# import sys
# sys.path.append('/Users/ghkerr/phd/globalo3/')
# import globalo3_open, globalo3_calculate
# sys.path.append('/Users/ghkerr/phd/tracer/')
# import tracer_open, tracer_calculate
# years = [2008, 2009, 2010]
# jja = ['jun', 'jul', 'aug']
# djf = ['jan', 'feb', 'dec']
# latmin = 0.
# latmax = 90.
# lngmin = 0.
# lngmax = 360.
# pmin = 800.
# pmax = 1005.
# import matplotlib.pyplot as plt
# # Load GEOSChem tracers
# TRAC_10_20_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem([2008], jja,
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_10_20', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_10_20_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja,
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_10_20', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_20_30_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_20_30', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_30_40_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_30_40', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_40_50_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_40_50', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_50_60_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_50_60', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_60_70_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_60_70', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_70_80_jja, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, jja, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_70_80', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# # Load MERRA-2 column meridional wind with the same vertical levels as 
# # GEOSChem 
# V_jja, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, jja, 'V',
#     lngmin, latmax, lngmax, latmin, pmin, pmax)
# V_jja = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, V_jja)
# # Open MERRA-2 column zonal wind at ~500 hPa to identify the eddy-driven jet
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
# Load Northern Hemisphere HindcastMR2 GMI CTM O3 for JJA
# lat_gmi, lng_gmi, times_gmi, o3_jja = \
#     globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
#     ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'O3', 'HindcastMR2')
# o3_jja = o3_jja*1e9
# # Degrade GMI to resolution of GEOSChem    
# o3_jja = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, 
#     lat_gmi, lng_gmi, o3_jja)
# # Repeat above but for DJF
# TRAC_10_20_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_10_20', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_20_30_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_20_30', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_30_40_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_30_40', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_40_50_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_40_50', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_50_60_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_50_60', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_60_70_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_60_70', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# TRAC_70_80_djf, lat_gc, lng_gc, lev_gc = tracer_open.open_geoschem(years, djf, 
#     'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_70_80', 
#     latmin, latmax, lngmin, lngmax, pmin, pmax)
# V_djf, lat_merra_gc, lng_merra_gc, lev_merra_gc = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf, 'V',
#     lngmin, latmax, lngmax, latmin, pmin, pmax)
# V_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra_gc, lng_merra_gc, V_djf)    
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
# # lat_gmi, lng_gmi, times_gmi_djf, o3_djf = \
# #     globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
# #     ['jan', 'feb', 'dec'], -1., 90., 0., 360., 'O3', 'HindcastMR2')
# # o3_djf = o3_djf*1e9
# # o3_djf = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, lng_gc, 
# #     lat_gmi, lng_gmi, o3_djf)
# lng = lng_gc
# lat = lat_gc
# lev = lev_gc
    
# # Figures 
# fig1(lat, lng, TRAC_10_20_jja, TRAC_20_30_jja, TRAC_30_40_jja, 
#     TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, TRAC_70_80_jja, edj_jja)

# fig2(lat, lng, TRAC_10_20_jja, TRAC_40_50_jja, TRAC_70_80_jja, TRAC_10_20_djf, 
#     TRAC_40_50_djf, TRAC_70_80_djf, edj_dist_jja, edj_dist_djf, edj_jja, 
#     edj_djf)

# fig3(lat, lng, TRAC_10_20_jja, TRAC_20_30_jja, TRAC_30_40_jja, 
#     TRAC_40_50_jja, TRAC_50_60_jja, TRAC_60_70_jja, TRAC_70_80_jja, 
#     TRAC_10_20_djf, TRAC_20_30_djf, TRAC_30_40_djf, TRAC_40_50_djf, 
#     TRAC_50_60_djf, TRAC_60_70_djf, TRAC_70_80_djf, V_jja, V_djf, 
#     edj_dist_jja, edj_dist_djf, edj_jja, edj_djf)
    
# fig4(lat, lng, TRAC_40_50_jja, TRAC_40_50_djf, V_jja, V_djf, edj_dist_jja, 
#     edj_dist_djf, edj_jja, edj_djf)

# # JET-TRACER CORRELATIONS FOR MASS-WEIGHTED 1OOO-200 HPA TRACERS 
# # Load GEOSChem tracers for full column (FC)
# # JJA 
# TRAC_10_20_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_10_20', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_20_30_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_20_30', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_30_40_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_30_40', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_40_50_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_40_50', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_50_60_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_50_60', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_60_70_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_60_70', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_70_80_jja_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     jja, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_70_80', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# # DJF
# TRAC_10_20_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_10_20', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_20_30_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_20_30', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_30_40_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_30_40', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_40_50_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_40_50', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_50_60_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_50_60', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_60_70_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_60_70', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# TRAC_70_80_djf_fc, lat_gc, lng_gc, lev_gc_fc = tracer_open.open_geoschem(years, 
#     djf, 'merra2_2x25_RnPbBe_co50', 'SpeciesConc', 'SpeciesConc_TRAC50_70_80', 
#     latmin, latmax, lngmin, lngmax, 240, 1000)
# # Vertically-integrated tracer concentrations 
# # JJA
# TRAC_10_20_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_10_20_jja_fc, 
#     lev_gc_fc)
# TRAC_20_30_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_20_30_jja_fc, 
#     lev_gc_fc)
# TRAC_30_40_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_30_40_jja_fc, 
#     lev_gc_fc)
# TRAC_40_50_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_40_50_jja_fc, 
#     lev_gc_fc)
# TRAC_50_60_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_50_60_jja_fc, 
#     lev_gc_fc)
# TRAC_60_70_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_60_70_jja_fc, 
#     lev_gc_fc)
# TRAC_70_80_jja_vi = tracer_calculate.tracer_mass_weight(TRAC_70_80_jja_fc, 
#     lev_gc_fc)
# # DJF
# TRAC_10_20_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_10_20_djf_fc, 
#     lev_gc_fc)
# TRAC_20_30_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_20_30_djf_fc, 
#     lev_gc_fc)
# TRAC_30_40_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_30_40_djf_fc, 
#     lev_gc_fc)
# TRAC_40_50_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_40_50_djf_fc, 
#     lev_gc_fc)
# TRAC_50_60_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_50_60_djf_fc, 
#     lev_gc_fc)
# TRAC_60_70_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_60_70_djf_fc, 
#     lev_gc_fc)
# TRAC_70_80_djf_vi = tracer_calculate.tracer_mass_weight(TRAC_70_80_djf_fc, 
#     lev_gc_fc)
# # Load MERRA-2 column meridional wind with the same vertical levels as 
# # GEOSChem 
# V_jja_fc, lat_merra, lng_merra, lev_merra_fc = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, jja, 'V',
#     lngmin, latmax, lngmax, latmin, 208, 1000)
# V_jja_fc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, V_jja_fc)
# V_djf_fc, lat_merra, lng_merra, lev_merra_fc = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(years, djf, 'V',
#     lngmin, latmax, lngmax, latmin, 208, 1000)
# V_djf_fc = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#     lng_gc, lat_merra, lng_merra, V_djf_fc)

# figS1()
# figS2(lat, lng, TRAC_10_20_jja_vi, TRAC_20_30_jja_vi, TRAC_30_40_jja_vi, 
#     TRAC_40_50_jja_vi, TRAC_50_60_jja_vi, TRAC_60_70_jja_vi, TRAC_70_80_jja_vi, 
#     TRAC_10_20_djf_vi, TRAC_20_30_djf_vi, TRAC_30_40_djf_vi, TRAC_40_50_djf_vi, 
#     TRAC_50_60_djf_vi, TRAC_60_70_djf_vi, TRAC_70_80_djf_vi, V_jja_fc, 
#     V_djf_fc, edj_dist_jja, edj_dist_djf, edj_jja, edj_djf)
# figS3(lat, lng, TRAC_40_50_jja_vi, TRAC_40_50_djf_vi, V_jja_fc, V_djf_fc, 
#     edj_dist_jja, edj_dist_djf, edj_jja, edj_djf)

# """JUNK """
# def nhmap(lat, lng, field, title, cbar_label, clevs, cmap, fstr, 
#     oceanon=False, extend='both', hatch=None, ebar=None): 
#     """Plot map (Robinson projection) of desired field. 
    
#     Parameters
#     ----------
#     lat : numpy.ndarray
#         Latitude coordinates for the Northern Hemisphere, units of degrees 
#         north, [lat,]                
#     lng : numpy.ndarray
#         Longitude coordinates for the Northern Hemisphere, units of degrees 
#         east, n.b., the first or last longitude value might need to be 
#         fudged to wrap contours around lng = 0 deg, [lng,]          
#     field : numpy.ndarray
#         Desired field/quantity, [lat, lng]
#     title : str
#         Title for plot        
#     cbar_label : str
#         Label for the colorbar (field and units)
#     clevs : numpy.ndarray
#         Filled contour levels values
#     cmap : str
#         Colormap name
#     fstr : str
#         Output filename suffix, should specify field, levels, and time period
#     oceanon : bool, optional 
#         If True, map adds ocean shapefiles; delfault is False
#     extend : str, optional
#         Extend settings for matplotlib colormap/colorbar; i.e., 'neither', 
#         'both', 'min', 'max'; default is 'both'
#     hatch : numpy.ndarray, optional
#         Hatching to juxtapose on filled contour; grid cells equal to 1 are 
#         hatched and values equal to NaN are left empty, [lat, lng]
#     ebar : numpy.ndarray, option
#         Errorbars to indicate mean location and variability of the eddy-driven
#         jet, units of degrees, [time, lng]

#     Returns
#     -------
#     None              
#     """
#     import numpy as np
#     import matplotlib as mpl
#     mpl.rcParams['hatch.linewidth'] = 0.3     
#     import matplotlib.pyplot as plt
#     import cartopy.crs as ccrs
#     import cartopy.feature as cfeature
#     from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
#     fig = plt.figure(figsize=(9,3.5))
#     ax=plt.subplot2grid((1,2), (0,0), colspan=2,
#                         projection=ccrs.PlateCarree(central_longitude=0.))
#     ax.set_title(title, fontsize=16, x=0.02, ha='left')
#     if oceanon == 'yes':
#         ax.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='lightgrey')
#     ax.coastlines(lw=0.25, color='k', zorder=3)
#     # Plot filled contours
#     cmap = plt.get_cmap(cmap)
#     # For wrapping around the Prime Meridian
#     if lng[-1] != 360:
#         lng[-1] = 360.
#     mb = ax.contourf(lng, lat, field, clevs, cmap=cmap, extend=extend,
#         transform=ccrs.PlateCarree(), zorder=1)
#     # If specified, add hatching (e.g., for significance)
#     if hatch is None: pass
#     else:
#         # Hatching for significance
#         ax.contourf(lng, lat, hatch, hatches=['//////'], colors='none', 
#             transform=ccrs.PlateCarree())
#     # If specified add errorbars (e.g., for location of eddy-driven jet)
#     if ebar is None: pass
#     else:
#         # Plot only every X number of longitude values
#         skiplng = 5
#         ax.errorbar(lng[::skiplng], np.nanmean(ebar, axis=0)[::skiplng], 
#             yerr=np.nanstd(ebar, axis=0)[::skiplng], zorder=10, color='k', 
#             markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
#             transform=ccrs.PlateCarree())
#     # Label parallels and meridians
#     ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     ax.xaxis.set_major_formatter(lng_formatter)      
#     ax.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()
#     ax.yaxis.set_major_formatter(lat_formatter)  
#     ax.set_extent([lng.min()-180., lng.max()-180., 0., 90.])  
#     # Add colorbar
#     plt.gcf().subplots_adjust(left=0.08, right=0.82)
#     colorbar_axes = plt.gcf().add_axes([0.85, ax.get_position().y0, 
#         0.02, (ax.get_position().y1-ax.get_position().y0)]) 
#     colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#         ticks=clevs, extend=extend)
#     colorbar.ax.tick_params(labelsize=12)
#     colorbar.set_label(cbar_label, fontsize=16)
#     ax.outline_patch.set_zorder(20)
#     plt.savefig('/Users/ghkerr/Desktop/VEWPW_tracergrad_trial.png', dpi=350)    
#     # plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     #     'nhmap_%s.png'%(fstr), dpi=350)
#     return

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

"""CALCULATE EXPECTED JET-TRACER RELATIONSHIP WITH ZONAL MEAN JET-V 
RELATIONSHIP AND ZONAL MEAN BACKGROUND GRADIENT"""
# clevs = np.linspace(-12, 12, 13)
# # Plot filled contours
# cmap = plt.get_cmap('coolwarm')
# # For wrapping around the Prime Meridian
# if lng[-1] != 360:
#     lng[-1] = 360.
# # METHOD 1: local background gradient, local r(V, jet)    
# fig = plt.figure(figsize=(7,7.4))
# ax1 = plt.subplot2grid((4,2), (0,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax2 = plt.subplot2grid((4,2), (1,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax3 = plt.subplot2grid((4,2), (2,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax4 = plt.subplot2grid((4,2), (3,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax1.set_title('(a) E[r($\chi_{10-20^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax2.set_title('(b) E[r($\chi_{40-50^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax3.set_title('(c) E[r($\chi_{70-80^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax4.set_title('(d) E[r(O$_3$, $\phi_{jet}$)]*10', fontsize=12, x=0.02, ha='left')
# axes = [ax1, ax2, ax3, ax4]
# # Determine near-surface meridional wind-jet correlation
# r_Vjet = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
#     edj_dist_jja, lat, lng)
# # Loop through all tracers
# for i,tracer in enumerate([np.nanmean(TRAC_10_20_jja, axis=1), 
#     np.nanmean(TRAC_40_50_jja, axis=1), np.nanmean(TRAC_70_80_jja, axis=1), 
#     o3_jja/1e8]):
#     # Determine tracer-jet correlation
#     r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_jja, lat, 
#         lng)
#     # Time-averaged tracer
#     tracer = np.nanmean(tracer, axis=0)
#     # Loop through each longitude and calculate the latitudinal gradient
#     grad_tracer = np.empty(shape=tracer.shape)
#     grad_tracer[:] = np.nan
#     for spine_i in np.arange(0, len(lng), 1):
#         # Tracer concentrations for given "spine" of latitudes at a given 
#         # longitude
#         spine = tracer[:,spine_i]
#         grad_spine = np.gradient(spine)
#         grad_tracer[:,spine_i] = grad_spine
#     grad_tracer = np.array(grad_tracer)
#     expected = -1.*r_Vjet*grad_tracer
#     mb = axes[i].contourf(lng, lat, expected*1e9, clevs, cmap=cmap, 
#         extend='both', transform=ccrs.PlateCarree(), zorder=1)
#     skiplng = 5
#     axes[i].errorbar(lng[::skiplng], np.nanmean(edj, axis=0)[::skiplng], 
#         yerr=np.nanstd(edj, axis=0)[::skiplng], zorder=10, color='k', 
#         markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
#         transform=ccrs.PlateCarree())
#     axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#     axes[i].set_extent([lng.min()-180., lng.max()-180., lat.min(), 
#         lat.max()])
#     axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
#         crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     axes[i].xaxis.set_major_formatter(lng_formatter)         
#     axes[i].get_xaxis().set_ticklabels([])
#     axes[i].set_yticks([0, 30, 60, 90], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     axes[i].yaxis.set_major_formatter(lat_formatter)    
#     axes[i].tick_params(which='major', labelsize=9)
#     axes[i].outline_patch.set_zorder(20)    
# ax4.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
# lng_formatter = LongitudeFormatter()
# ax4.xaxis.set_major_formatter(lng_formatter)       
# ax4.tick_params(which='major', labelsize=9)
# plt.subplots_adjust(left=0.05, right=0.85, hspace=0.4)
# colorbar_axes = plt.gcf().add_axes([0.85, ax4.get_position().y0, 
#     0.02, (ax1.get_position().y1-ax4.get_position().y0)]) 
# colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#     ticks=clevs, extend='both')
# colorbar.ax.tick_params(labelsize=12)
# colorbar.set_label('[ppbv/$\circ$]', fontsize=16)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'expectedvalue_contourf.png', dpi=500)
# plt.show()

# # METHOD 2: local background gradient, zonal mean r(V, jet)
# fig = plt.figure(figsize=(7,7.4))
# ax1 = plt.subplot2grid((4,2), (0,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax2 = plt.subplot2grid((4,2), (1,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax3 = plt.subplot2grid((4,2), (2,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax4 = plt.subplot2grid((4,2), (3,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax1.set_title('(a) E[r($\chi_{10-20^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax2.set_title('(b) E[r($\chi_{40-50^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax3.set_title('(c) E[r($\chi_{70-80^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax4.set_title('(d) E[r(O$_3$, $\phi_{jet}$)]*10', fontsize=12, x=0.02, ha='left')
# axes = [ax1, ax2, ax3, ax4]
# # Determine near-surface meridional wind-jet correlation
# r_Vjet = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
#     edj_dist_jja, lat, lng)
# # Loop through all tracers
# for i,tracer in enumerate([np.nanmean(TRAC_10_20_jja, axis=1), 
#     np.nanmean(TRAC_40_50_jja, axis=1), np.nanmean(TRAC_70_80_jja, axis=1), 
#     o3_jja/1e8]):
#     # Determine tracer-jet correlation
#     r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_jja, lat, 
#         lng)
#     # Time-averaged tracer
#     tracer = np.nanmean(tracer, axis=0)
#     # Loop through each longitude and calculate the latitudinal gradient
#     grad_tracer = np.empty(shape=tracer.shape)
#     grad_tracer[:] = np.nan
#     r_Vjet = np.empty(shape=tracer.shape)
#     r_Vjet[:] = np.nan
#     for spine_i in np.arange(0, len(lng), 1):
#         spine = tracer[:,spine_i]
#         grad_spine = np.gradient(spine)
#         grad_tracer[:,spine_i] = grad_spine
#         r_Vjet[:,spine_i] = r_Vjet_zm
#     grad_tracer = np.array(grad_tracer)
#     expected = -1.*r_Vjet*grad_tracer
#     mb = axes[i].contourf(lng, lat, expected*1e9, clevs, cmap=cmap, 
#         extend='both', transform=ccrs.PlateCarree(), zorder=1)
#     skiplng = 5
#     axes[i].errorbar(lng[::skiplng], np.nanmean(edj, axis=0)[::skiplng], 
#         yerr=np.nanstd(edj, axis=0)[::skiplng], zorder=10, color='k', 
#         markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
#         transform=ccrs.PlateCarree())
#     axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#     axes[i].set_extent([lng.min()-180., lng.max()-180., lat.min(), 
#         lat.max()])
#     axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
#         crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     axes[i].xaxis.set_major_formatter(lng_formatter)         
#     axes[i].get_xaxis().set_ticklabels([])
#     axes[i].set_yticks([0, 30, 60, 90], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     axes[i].yaxis.set_major_formatter(lat_formatter)    
#     axes[i].tick_params(which='major', labelsize=9)
#     axes[i].outline_patch.set_zorder(20)    
# ax4.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
# lng_formatter = LongitudeFormatter()
# ax4.xaxis.set_major_formatter(lng_formatter)       
# ax4.tick_params(which='major', labelsize=9)
# plt.subplots_adjust(left=0.05, right=0.85, hspace=0.4)
# colorbar_axes = plt.gcf().add_axes([0.85, ax4.get_position().y0, 
#     0.02, (ax1.get_position().y1-ax4.get_position().y0)]) 
# colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#     ticks=clevs, extend='both')
# colorbar.ax.tick_params(labelsize=12)
# colorbar.set_label('[ppbv/$\circ$]', fontsize=16)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'expectedvalue_zonalmeanr_Vjet_contourf.png', dpi=500)
# plt.show()

# # METHOD 3: local r(V, jet), zonal mean background gradient
# fig = plt.figure(figsize=(7,7.4))
# ax1 = plt.subplot2grid((4,2), (0,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax2 = plt.subplot2grid((4,2), (1,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax3 = plt.subplot2grid((4,2), (2,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax4 = plt.subplot2grid((4,2), (3,0), colspan=2,
#     projection=ccrs.PlateCarree(central_longitude=0.))
# ax1.set_title('(a) E[r($\chi_{10-20^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax2.set_title('(b) E[r($\chi_{40-50^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax3.set_title('(c) E[r($\chi_{70-80^{\circ}}$, $\phi_{jet}$)]', fontsize=12, 
#     x=0.02, ha='left')
# ax4.set_title('(d) E[r(O$_3$, $\phi_{jet}$)]*10', fontsize=12, x=0.02, ha='left')
# axes = [ax1, ax2, ax3, ax4]
# # Determine near-surface meridional wind-jet correlation
# r_Vjet = globalo3_calculate.calculate_r(np.nanmean(V_jja, axis=1), 
#     edj_dist_jja, lat, lng)
# # Loop through all tracers
# for i,tracer in enumerate([np.nanmean(TRAC_10_20_jja, axis=1), 
#     np.nanmean(TRAC_40_50_jja, axis=1), np.nanmean(TRAC_70_80_jja, axis=1), 
#     o3_jja/1e8]):
#     # Determine tracer-jet correlation
#     r_tracerjet = globalo3_calculate.calculate_r(tracer, edj_dist_jja, lat, 
#         lng)
#     # Time-averaged tracer
#     tracer = np.nanmean(tracer, axis=0)
#     # Loop through each longitude and calculate the latitudinal gradient
#     grad_tracer = np.empty(shape=tracer.shape)
#     grad_tracer[:] = np.nan
#     grad_spine = np.gradient(np.nanmean(tracer, axis=-1))
#     for spine_i in np.arange(0, len(lng), 1):
#         spine = tracer[:,spine_i]
#         grad_tracer[:,spine_i] = grad_spine
#     grad_tracer = np.array(grad_tracer)
#     expected = -1.*r_Vjet*grad_tracer
#     mb = axes[i].contourf(lng, lat, expected*1e9, clevs, cmap=cmap, 
#         extend='both', transform=ccrs.PlateCarree(), zorder=1)
#     skiplng = 5
#     axes[i].errorbar(lng[::skiplng], np.nanmean(edj, axis=0)[::skiplng], 
#         yerr=np.nanstd(edj, axis=0)[::skiplng], zorder=10, color='k', 
#         markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
#         transform=ccrs.PlateCarree())
#     axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#     axes[i].set_extent([lng.min()-180., lng.max()-180., lat.min(), 
#         lat.max()])
#     axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
#         crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     axes[i].xaxis.set_major_formatter(lng_formatter)         
#     axes[i].get_xaxis().set_ticklabels([])
#     axes[i].set_yticks([0, 30, 60, 90], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     axes[i].yaxis.set_major_formatter(lat_formatter)    
#     axes[i].tick_params(which='major', labelsize=9)
#     axes[i].outline_patch.set_zorder(20)    
# ax4.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())    
# lng_formatter = LongitudeFormatter()
# ax4.xaxis.set_major_formatter(lng_formatter)       
# ax4.tick_params(which='major', labelsize=9)
# plt.subplots_adjust(left=0.05, right=0.85, hspace=0.4)
# colorbar_axes = plt.gcf().add_axes([0.85, ax4.get_position().y0, 
#     0.02, (ax1.get_position().y1-ax4.get_position().y0)]) 
# colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#     ticks=clevs, extend='both')
# colorbar.ax.tick_params(labelsize=12)
# colorbar.set_label('[ppbv/$\circ$]', fontsize=16)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+
#     'expectedvalue_zonalmeanlatgradient_contourf.png', dpi=500)
# plt.show()

"""COMPARE PW-EW TRACER COMPOSITES TO DAILY VARIABILITY"""
# if lng[-1] != 360:
#     lng[-1] = 360.    
# fig = plt.figure(figsize=(11,5.5))
# ax1 = plt.subplot2grid((3,4), (0,0), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax2 = plt.subplot2grid((3,4), (1,0), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax3 = plt.subplot2grid((3,4), (2,0), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax4 = plt.subplot2grid((3,4), (0,2), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax5 = plt.subplot2grid((3,4), (1,2), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax6 = plt.subplot2grid((3,4), (2,2), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))    
# ax1.set_title('(a) JJA $\chi_{70-80}$ (PW$\:$-$\:$EW)', fontsize=12, 
#     x=0.02, ha='left')
# ax2.set_title('(c) JJA $\chi_{40-50}$ (PW$\:$-$\:$EW)', fontsize=12, 
#     x=0.02, ha='left')
# ax3.set_title('(e) JJA $\chi_{10-20}$ (PW$\:$-$\:$EW)', fontsize=12, 
#     x=0.02, ha='left')
# ax4.set_title('(b) DJF $\chi_{70-80}$ (PW$\:$-$\:$EW)', fontsize=12, 
#     x=0.02, ha='left')
# ax5.set_title('(d) DJF $\chi_{40-50}$ (PW$\:$-$\:$EW)', fontsize=12, 
#     x=0.02, ha='left')
# ax6.set_title('(f) DJF $\chi_{10-20}$ (PW$\:$-$\:$EW)', fontsize=12, 
#     x=0.02, ha='left')
# # Loop through JJA axes and plot
# axes = [ax1, ax2, ax3]
# for i, tracer in enumerate([TRAC_70_80_jja, TRAC_40_50_jja, 
#     TRAC_10_20_jja]):
#     # Calculate tracer-jet correlation and significance 
#     r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
#         edj_dist_jja, lat, lng)
#     significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
#         np.nanmean(tracer, axis=1), edj_dist_jja, r_tracerjet, lat, lng)
#     mb = axes[i].contourf(lng, lat, np.nanstd(tracer,axis=tuple((0,1)))*1e6,
#         np.linspace(0, 0.6, 7), cmap=plt.get_cmap('Reds'), 
#         extend='max', transform=ccrs.PlateCarree(), zorder=2)
#     axes[i].contourf(lng, lat, significance_r_tracerjet, 
#         hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
#         zorder= 4)
#     # Add eddy-driven jet
#     skiplng = 5
#     axes[i].errorbar(lng[::skiplng], np.nanmean(edj_jja, axis=0)[::skiplng], 
#         yerr=np.nanstd(edj_jja, axis=0)[::skiplng], zorder=10, color='k', 
#         markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
#         transform=ccrs.PlateCarree())
#     axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#     axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
#     axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
#         crs=ccrs.PlateCarree())
#     axes[i].get_xaxis().set_ticklabels([])
#     axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     axes[i].yaxis.set_major_formatter(lat_formatter)
#     axes[i].tick_params(which='major', labelsize=9)
# # For DJF
# axes = [ax4, ax5, ax6]
# for i, tracer in enumerate([TRAC_70_80_djf, TRAC_40_50_djf, 
#     TRAC_10_20_djf]):
#     r_tracerjet = globalo3_calculate.calculate_r(np.nanmean(tracer, axis=1), 
#         edj_dist_djf, lat, lng)
#     significance_r_tracerjet = globalo3_calculate.calculate_r_significance(
#         np.nanmean(tracer, axis=1), edj_dist_djf, r_tracerjet, lat, lng)
#     mb = axes[i].contourf(lng, lat, np.nanstd(tracer, axis=tuple((0,1)))*1e6,
#         np.linspace(0, 0.6, 7), cmap=plt.get_cmap('Reds'), 
#         extend='max', transform=ccrs.PlateCarree(), zorder=2)
#     axes[i].contourf(lng, lat, significance_r_tracerjet, 
#         hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
#         zorder= 4)
#     skiplng = 5
#     axes[i].errorbar(lng[::skiplng], np.nanmean(edj_djf, axis=0)[::skiplng], 
#         yerr=np.nanstd(edj_djf, axis=0)[::skiplng], zorder=10, color='k', 
#         markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
#         transform=ccrs.PlateCarree())
#     axes[i].coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#     axes[i].set_extent([lng.min()-180., lng.max()-180., 10., 80.])
#     axes[i].set_xticks([-180, -120, -60, 0, 60, 120, 180], 
#         crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     axes[i].xaxis.set_major_formatter(lng_formatter)         
#     axes[i].get_xaxis().set_ticklabels([])
#     axes[i].set_yticks([15, 50, 75], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     axes[i].get_yaxis().set_ticklabels([])
#     axes[i].tick_params(which='major', labelsize=9)
# for ax in [ax3, ax6]:
#     ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], 
#         crs=ccrs.PlateCarree())    
#     lng_formatter = LongitudeFormatter()
#     ax.xaxis.set_major_formatter(lng_formatter)       
#     ax.tick_params(which='major', labelsize=9)
# plt.subplots_adjust(left=0.05, right=0.86)
# # Add colorbar
# cbaxes = fig.add_axes([ax4.get_position().x1+0.03, ax6.get_position().y0, 
#     0.02, ax4.get_position().y1-ax6.get_position().y0]) 
# cb = plt.colorbar(mb, cax=cbaxes, orientation='vertical')
# cb.set_label(label='[ppm]', size=12)
# cb.set_ticks(np.linspace(0, 0.6, 7))
# cb.ax.tick_params(labelsize=9)
# plt.savefig('/Users/ghkerr/phd/tracer/figs/'+'fig3_std.pdf', dpi=500)
# plt.show() 


# # # # # Find standard deviation vs. PW-EW composites in mid-latitudes
# where40 = np.where(lat==40.)[0][0]
# where60 = np.where(lat==60.)[0][0]
# # Overall variability
# TRAC_10_20_jja_std = np.nanstd(np.nanmean(TRAC_10_20_jja, axis=1), axis=0)
# TRAC_40_50_jja_std = np.nanstd(np.nanmean(TRAC_40_50_jja, axis=1), axis=0)
# TRAC_70_80_jja_std = np.nanstd(np.nanmean(TRAC_70_80_jja, axis=1), axis=0)
# # Jet-related variability
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_TRAC_10_20_jja, \
#     eqjet_TRAC_10_20_jja = globalo3_calculate.segregate_field_bylat(
#     np.nanmean(TRAC_10_20_jja, axis=1), lng, edj_jja, 
#     np.arange(0, len(TRAC_10_20_jja), 1))
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_TRAC_40_50_jja, \
#     eqjet_TRAC_40_50_jja = globalo3_calculate.segregate_field_bylat(
#     np.nanmean(TRAC_40_50_jja, axis=1), lng, edj_jja, 
#     np.arange(0, len(TRAC_40_50_jja), 1))
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_TRAC_70_80_jja, \
#     eqjet_TRAC_70_80_jja = globalo3_calculate.segregate_field_bylat(
#     np.nanmean(TRAC_70_80_jja, axis=1), lng, edj_jja, 
#     np.arange(0, len(TRAC_70_80_jja), 1))
# # Fraction of overall variability related to the jet
# frac_10_20 = np.abs(pwjet_TRAC_10_20_jja-eqjet_TRAC_10_20_jja)/TRAC_10_20_jja_std
# frac_40_50 = np.abs(pwjet_TRAC_40_50_jja-eqjet_TRAC_40_50_jja)/TRAC_40_50_jja_std
# frac_70_80 = np.abs(pwjet_TRAC_70_80_jja-eqjet_TRAC_70_80_jja)/TRAC_70_80_jja_std
# np.nanmean(frac_10_20[where40:where60+1])
# np.nanmean(frac_40_50[where40:where60+1])
# np.nanmean(frac_70_80[where40:where60+1])
# land = globalo3_calculate.find_grid_overland(lat[:-3], lng)
# np.nanmean((frac_40_50[:-3]*land)[where40:where60+1])