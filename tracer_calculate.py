#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate fields related to GEOSChem Rn–Pb–Be tracer simulations   

Module calculates relevant quantities or fields for GEOSChem tracer simulations
such as the position and distance from the eddy-driven jet

Revision History
    26022020 -- initial version created and code from "tracer_plot.py" migrated
                to this module
"""

def find_jetlatitude(edj, lat, lng): 
    """Determine the daily latitude of the eddy-driven jet locally (at each 
    longitudinal grid cell) based on the latitude of maximum zonal wind (U). 
    Note this is usually done with ~500 hPa zonal winds restricted to 20-70˚N.
    
    Parameters
    ----------
    edj : numpy.ndarray
        Zonal wind (U) in region, units of m s-1, [time, lat, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]

    Returns
    -------
    edj_lat : numpy.ndarray
        The latitude of the eddy-driven jet, identifed by the latitude of 
        maximum zonal wind (U) at a particular pressure/model level, units 
        of degrees north, [time, lng]
    """
    import time
    start_time = time.time()    
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'Determining eddy-driven jet latitude...')   
    import numpy as np
    # Array will be filled the daily latitude of the jet (defined as the 
    # location with maximum U wind at ~500 mb) in the focus region
    edj_lat = np.empty(shape=(len(edj), len(lng)))
    edj_lat[:] = np.nan
    # Find latitude of jet location 
    for day in np.arange(0, len(edj), 1):
        edj_day = edj[day]
        # Loop through longitude
        for i in np.arange(0, len(lng), 1):  
            # U wind at 500 hPa for longitude/day of interest 
            edj_transect = edj_day[:, i]
            transectmax = np.where(edj_transect==np.nanmax(edj_transect))[0][0]   
            # Find latitude of jet
            edj_lat[day, i] = lat[transectmax]
    del i
    print('Latitude determined in %.2f seconds!'%(time.time() - start_time))
    return edj_lat

def determine_jetdist(edj, lat, lng):
    """Determine the distance, in degrees, from the lataitude of the eddy-
    driven jet on each day for each model grid cell. Here a positive (negative)
    distance implies that the jet is north (south) of a particular grid cell.

    Parameters
    ----------
    edj : numpy.ndarray
        The latitude of the eddy-driven jet, identifed by the latitude of 
        maximum zonal wind (U) at a particular pressure/model level, units 
        of degrees north, [time, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]

    Returns
    -------
    edj_dist : numpy.ndarray
        Distance from the eddy-driven jet where positive distances imply 
        a southward eddy-driven jet, units of degrees, [time, lat, lng]
    """
    import numpy as np
    lat_grid = np.tile(lat, (lng.shape[0], 1)).T
    # Grid will be populated with the distance from the eddy-driven jet (in 
    # degrees) at each lat, lon coordinate
    edj_dist = np.empty(shape=(len(edj),lat_grid.shape[0],lat_grid.shape[1]))
    edj_dist[:] = np.nan
    # Loop through days in measuring period; and, at each day, subtract the 
    # latitude of the jet from the local latitude  
    for day in np.arange(0,len(edj),1):
        edj_dist[day] = edj[day]-lat_grid
    return edj_dist

def tracer_mass_weight(tracer, lev):
    """Account for mass of different levels in column calculation by
    vertically integrating over pressure levels. 

    Parameters
    ----------
    tracer : numpy.ndarray
        Tracer mixing ratio, units of mol mol-1, [time, lev, lat, lng] 
    lev : numpy.ndarray
        Model pressure levels, units of hPa, [lev]               

    Returns
    -------
    tracer_integrate : numpy.ndarray    
        Vertically-integrated pressure-weighted tracer, units of kg m-1 s-2, 
        [time, 1, lat, lng]
    """
    import numpy as np
    # Difference between pressure levels
    dp = np.abs(np.diff(lev))
    # Append 15 hPa to beginning of array; this is the difference between the 
    # surface and second layer (see diagram in tracer_plot.py)
    dp = np.hstack((np.array([15.]), dp))
    # Convert to Pa 
    dp = dp*100.
    tracer_integrate = []
    for i, ilev in enumerate(lev):
        tracer_integrate.append(tracer[:,i]*dp[i])
    tracer_integrate = np.stack(tracer_integrate)
    tracer_integrate = np.nansum(tracer_integrate, axis=0)
    tracer_integrate = np.reshape(tracer_integrate, 
        [tracer_integrate.shape[0], 1, tracer_integrate.shape[1], 
         tracer_integrate.shape[2]])
    return tracer_integrate