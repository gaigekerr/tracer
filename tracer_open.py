#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Open files related to GEOSChem Rn–Pb–Be tracer simulations 

Module opens GEOSChem (and Replay) model simulations and MERRA-2 assimilated 
data and parses the desired horizontal and vertical extent, time period, and
variables.

Revision History
    26022020 -- initial version created
    27022020 -- function 'open_merra2_inst3_3d_asm_Nv_specifieddomain' added
"""

def open_geoschem(years, months, simulation, collection, varname, latmin, 
    latmax, lngmin, lngmax, pmin, pmax, operation=None):
    """function opens daily mean GEOSChem output for tracer/specie of interest
    and extracts columned diagnostics for the latitudes, longitudes, and levels
    of interest. A sum, average, standard deviation is performed over the 
    column, if desired.
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    simulation : str
        GEOS-Chem simulation; should be the same as the run directory
        (e.g., merra2_2x25_RnPbBe_co50, merra2_2x25_tropchem)
    collection : str
        GEOS-Chem history collection (e.g., SpeciesConc, AerosolMass)
    varname : str
        Variable name (e.g., SpeciesConc_TRAC50_40_50, SpeciesConc_O3)
    latmin : float    
        Latitude (degrees north) of bottom edge of bounding box for focus 
        region. For this parameter and others defining the bounding box, 
        function finds the closest index to bounding box edges
    latmax : float
        Latitude (degrees north) of upper edge of bounding box for focus region
    lngmin : float
        Longitude (degrees east, 0-360) of left edge of bounding box for focus 
        region        
    lngmax : float
        Longitude (degrees east, 0-360) of right edge of bounding box for focus 
        region
    pmin : float
        The upper pressure level of interest, units of hPa
    pmax : float
        The lower pressure level of interest, units of hPa    
    operation : str, optional
        Operation (mean, sum, std) to apply to the tracer diagnostic column
        (the default is None, which implies no operation will be applied to
        the column axis).
        
    Returns
    -------
    var : numpy.ndarray
        GEOSChem tracer diagnostics for the region of interest, units of dry 
        mixing ratio (mol mol-1), [time, lat, lng] or [time, lev, lat, lng]
    lat : numpy.ndarray
        GEOSChem latitude coordinates, units of degrees north, [lat,]    
    lng : numpy.ndarray
        GEOSChem longitude coordinates, units of degrees east, [lat,]    
    lev : numpy.ndarray
        GEOSChem pressure level coordinates corresponding to the top edge of 
        the model layer, units of hPa, [lev,]
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading GEOSChem %s...' %varname)        
    import numpy as np
    from netCDF4 import Dataset
    import sys
    sys.path.append('/Users/ghkerr/phd/utils/')
    from geo_idx import geo_idx
    # Convert month abbrevations to integers
    months_int = []
    for m in months:
        months_int.append(str(time.strptime(m,'%b').tm_mon).zfill(2))
    # List will filled with monthly GEOSChem output for variable of interest
    var = []    
    # Loop through years, months of interest    
    for year in years:
        PATH_GEOSCHEM='/Users/ghkerr/phd/tracer/data/'+\
            '%s/%d/'%(simulation,year)
        for month in months_int:
            infile = Dataset(PATH_GEOSCHEM+'GEOSChem.%s.%d%s.nc4'
                %(collection,year,month),'r') 
            # On first iteration, extract dimensions and find indicies 
            # corresponding to (closest to) desired domain
            if (year==years[0]) and (month==months_int[0]):
                lat = infile.variables['lat'][:]
                lng = infile.variables['lon'][:] 
                lev = infile.variables['lev'][:]*1000. # Convert to hPa
                # Convert longitude from (-180-180) to (0-360)
                lng = lng%360          
                # Shift grid such that it spans (0-360) rather than (180-360, 
                # 0-180)
                lng = np.roll(lng,int(lng.shape[0]/2))
                latmin = geo_idx(latmin,lat)
                latmax = geo_idx(latmax,lat)
                lngmin = geo_idx(lngmin,lng)
                lngmax = geo_idx(lngmax,lng)
                pmax = (np.abs(lev-pmax)).argmin()
                pmin = (np.abs(lev-pmin)).argmin()
                lnglen = lng.shape[0]
                # Restrict coordinates over focus region 
                lat = lat[latmin:latmax+1]
                lng = lng[lngmin:lngmax+1] 
                # Although this function indexes the GEOSChem output using the 
                # 'lev' dimension, it will return the 'ilev' dimension instead.
                # ilev represents the hybrid level at interfaces, whereas lev 
                # represents the hybrid level at midpoints. MERRA-2 dimensions 
                # represent the model layer top edge and are therefore the 
                # most directly comparable to ilev. 
                # i.e., 
                # GEOSCHem ilev             GEOSChem lev
                # --- 940 hPa --
                #                           --- 947.5 hPa --
                # --- 955 hPa --
                #                           --- 962.5 hPa --
                # --- 970 hPa --
                #                           --- 977.5 hPa --
                # --- 985 hPa --
                #                           --- 992.5 hPa --
                # -- 1000 hPa --
                # So if you have the MERRA-2 model layer corresponding to 
                # a model top edge at 985 hPa, you'd want to find the  the 
                # value of lev closest to this (lev[0]) and then take the 
                # ilev[0+1] to find tthe top edge.
                lev = infile.variables['ilev'][:]
                lev = np.round(lev[pmax+1:pmin+1+1]*1000.)
            # Extract variable for the month
            var_month = infile.variables[varname][:]
            # Roll grid similar to longitude grid
            var_month = np.roll(var_month, int(var_month.shape[-1]/2), 
                axis=np.where(np.array(var_month.shape)==lnglen)[0][0])
            # Extract output for horizontal and verticle domain of interest;
            # for this to work, dimensions must be (time, lev, lat, lng)
            var_month = var_month[:, pmax:pmin+1, latmin:latmax+1, 
                lngmin:lngmax+1]
            if operation == 'mean':
                var_month = np.nanmean(var_month, axis=1)
                var.append(var_month)
            elif operation == 'sum':
                var_month = np.nansum(var_month, axis=1)            
                var.append(var_month)
            elif operation == 'std':
                var_month = np.nanstd(var_month, axis=1)
                var.append(var_month)
            elif operation == None:
                var.append(var_month)
    # Stack 
    var = np.vstack(var)
    print('%s for %d-%d loaded in %.2f seconds!' %(varname, years[0], 
        years[-1], (time.time()-start_time)))
    return var.data, lat.data, lng.data, lev.data

def open_merra2_inst3_3d_asm_Nv_specifieddomain(years, months, varname, lngmin, 
    latmax, lngmax, latmin, pmin, pmax, operation=None):
    """function opens daily mean MERRA-2 inst3_3d_asm_Nv fields (3d, 3-Hourly,
    Instantaneous, Model-Level, Assimilation, Assimilated Meteorological Fields 
    V5.12.4) for the specified months and years. The variable of interest 
    is extracted for the region and pressure level(s) of interest. If 
    specified, function can also compute the column-averaged value.

    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    months : list
        Three letter abbreviations (lowercase) for months in measuring period    
    varname : str
        Variable of interest reanalysis Options include T (air temperature), 
        q (specific humidity), U (eastward wind), or V (northward wind)
    lngmin : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    latmax : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    
    lngmax : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    latmin : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north
    pmin : float
        The upper pressure level of interest, units of hPa
    pmax : float
        The lower pressure level of interest, units of hPa    
    operation : str, optional
            Operation (mean, sum, std) to apply to the tracer diagnostic column
            (the default is None, which implies no operation will be applied to
            the column axis).
        
    Returns
    -------
    var : numpy.ndarray
        Model output for specified variable, units of ppbv, if operation = 
        None then shape is [time, lat, lng] if else the shape is [time, lev, 
        lat, lng]. Note the 0th index of the level dimension corresponds to 
        the pressure level closest to the surface (i.e., closest to pmax)
    lat : numpy.ndarray
        Model latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Model numpy.ndarray coordinates, units of degrees east, [lng,]    
    lev : numpy.ndarray
        Model pressure levels corresponding to the top edge of the layer, units 
        of hPa, [lev]
    """
    import time
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading %s from MERRA-2 simulation...' %varname)
    start_time = time.time()
    import numpy as np
    from netCDF4 import Dataset
    import sys
    sys.path.append('/Users/ghkerr/phd/utils/')
    from geo_idx import geo_idx
    # Define pressures; for whatever reason the value of the pressues didn't 
    # download with the MERRA-2 model-level data, so I manually copy and 
    # pasted them from here from pg. 10 of 
    # https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
    # Within this dataset, the pressure corresponding to the vertical level 
    # are nominal for a 1000 hPa surface pressure and refer to the top edge of 
    # the layer. Note that the bottom level layer has a nominal thickness of 15 hPa.
    presslev = {43 : 208.152, 44 : 244.875, 45 : 288.083, 46 : 337.500,
        47 : 375.000, 48 : 412.500, 49 : 450.000, 50 : 487.500, 51 : 525.000, 
        52 : 562.500, 53 : 600.000, 54 : 637.500, 55 : 675.000, 56 : 700.000, 
        57 : 725.000, 58 : 750.000, 59 : 775.000, 60 : 800.000, 61 : 820.000,
        62 : 835.000, 63 : 850.000, 64 : 865.000, 65 : 880.000, 66 : 895.000, 
        67 : 910.000, 68 : 925.000, 69 : 940.000, 70 : 955.000, 71 : 970.000, 
        72 : 985.000}
    # Convert month abbrevations to integers
    months_int = []
    for m in months:
        months_int.append(str(time.strptime(m,'%b').tm_mon).zfill(2))
    # List will filled with monthly GEOSChem output for variable of interest
    var = []    
    # Loop through years, months of interest    
    for year in years:
        PATH_MERRA='/Users/ghkerr/phd/meteorology/data/inst3_3d_asm_Nv/'+\
            '%d/'%(year)
        for month in months_int:
            infile = Dataset(PATH_MERRA+'MERRA2_300.inst3_3d_asm_Nv.%d%s.nc'
                %(year,month),'r')
            # On first iteration, extract dimensions and find indicies 
            # corresponding to (closest to) desired domain
            if (year==years[0]) and (month==months_int[0]):
                lat = infile.variables['lat'][:]
                lng = infile.variables['lon'][:] 
                # Convert longitude from (-180-180) to (0-360)
                lng = lng%360          
                # Shift grid such that it spans (0-360) rather than (180-360, 
                # 0-180)
                lng = np.roll(lng,int(lng.shape[0]/2))
                # Kludgey 
                if (lng.data[0] > 359.5) & (lng.data[0] <= 360.):
                    lng[0] = 0.
                latmin = geo_idx(latmin,lat)
                latmax = geo_idx(latmax,lat)
                lngmin = geo_idx(lngmin,lng)
                lngmax = geo_idx(lngmax,lng)
                lev = np.fromiter(presslev.values(), dtype=float)
                pmax = (np.abs(lev-pmax)).argmin()
                pmin = (np.abs(lev-pmin)).argmin()
                lnglen = lng.shape[0]
                # Restrict coordinates over focus region 
                lat = lat[latmin:latmax+1]
                lng = lng[lngmin:lngmax+1] 
                lev = np.flip(lev[pmin:pmax+1], axis=0)
            # Extract variable for the month
            var_month = infile.variables[varname][:]
            # Roll grid similar to longitude grid
            var_month = np.roll(var_month, int(var_month.shape[-1]/2), 
                axis=np.where(np.array(var_month.shape)==lnglen)[0][0])
            # Extract output for horizontal and verticle domain of interest;
            # for this to work, dimensions must be (time, lev, lat, lng)
            var_month = var_month[:, pmin:pmax+1, latmin:latmax+1, 
                lngmin:lngmax+1] 
            # Flip level dimension such that the lowest layer corresponds to 
            # the pressure level closest to the ground (closest to pmax)
            var_month = np.flip(var_month, axis=1)
            if operation == 'mean':
                var_month = np.nanmean(var_month, axis=1)
                var.append(var_month)
            elif operation == 'sum':
                var_month = np.nansum(var_month, axis=1)            
                var.append(var_month)
            elif operation == 'std':
                var_month = np.nanstd(var_month, axis=1)
                var.append(var_month)
            elif operation == None:
                var.append(var_month)
    # Stack 
    var = np.vstack(var)
    print('%s for %d-%d loaded in %.2f seconds!' %(varname, years[0], 
        years[-1], (time.time()-start_time)))    
    return var, lat, lng, lev