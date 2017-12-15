# ========================================================================
# ========================================================================
# Filnemae; ncl_helper.py
# Purpose: To create functions that mimic NCL fucntionalities as close
# as possible
# ========================================================================
# ========================================================================


# ==================== Importing libraries ================================

import numpy as np
from datetime import datetime, timedelta
import warnings

# ======================== NCL Helper Clones =================================
def cd_calendar(time_array, time_start, option):

    """cd_calendar() is python equivalent of NCL's function for converting netcdf4 date formats such as 'Days Since'
    time_axis: Time in units of "Hours since" or "Days since"
    time_start: The starting point specified in the netcdf file. Use ncdump -h to find the units. Pass all 5 components (H, M, D, HR, S)
    options:
    'h' if Hours since
    'd' if Days since
    """
    
    if(option not in ('h','d')):
        return ValueError("option: Not a valid input. Refer function help for valid inputs")
    
    # Compatibility with some later timedelta. Conver to float64
    if(not isinstance(time_array, np.ndarray)):
        raise ValueError("Timedeltas must be a numpy array")
        
    time_array = time_array.astype(np.float64)

    # import datetime and timedelta before running this code
    return([datetime(time_start[0], time_start[1], time_start[2], time_start[3], time_start[4]) + (timedelta(days = tm)
        if (option == 'd') else timedelta(hours = tm)) for tm in time_array])
    
    
def subset_time(data, time, time0, time1):
    
    """subset_time() is used to subset arrays along a time dimension (usually the first dimension)
    inputs:
        data: A np.ndarray which has to be subsetted along the time dimension
        time: A datetime list of values along the time axis of the data to be subsetted
        [time0, time1]: Start and end of the time axis
    outputs:
        sub_data: A subset of the original dataset along the time dimesion from time0 to time1"""
        
        
    if(not isinstance(data, np.ndarray)):
        raise ValueError("subset_time: Data not a numpy array")
        
    shp = data.shape
    if(shp[0] != len(time)):
        raise ValueError("subset_time: first axis might not be defined as time or improper time coordinate")
    
    if((not isinstance(time0, datetime)) or (not isinstance(time1, datetime))):
        raise ValueError("subset_time: datetime object expected for time0 and time1")
    
    if((time0 < min(time)) or (time1 > max(time))):
        warnings.warn("subset_time: One or both of the bounds exceeds the time axis")
    
    dt_filt = np.array([(tm <= time1) and (tm >= time0) for tm in time])
    sub_time = [time[i] for i in range(len(time)) if dt_filt[i]]
    sub_data = data[dt_filt]
    
    return sub_data, sub_time
    


# Find closest value    
def find_closest(arr, elm, option):

    """find_closest() tries to find the closest value to a given x and returns the index of that value
    arr: a 1D array,
    elm: a value for which closest value exists in array arr
    option: -1: smallest possible match
             0: Any match as long as it is least
             1: Largest possible match"""

    diff = np.array(map(lambda x: (x - elm), arr))

    msk = np.empty_like(diff, bool)


    if(option == 0):
        diff = abs(diff)
        msk[:] = False

    elif(option == -1):
        msk = diff > 0

    elif(option == 1):
        msk = diff < 0

    else:
        raise ValueError("Expected only 0, -1 or 1 for option")

    diff_ma = np.ma.masked_array(abs(diff), mask = msk)
    
    return np.argmin(diff_ma)


# Find closest date (Datetime objects)
def find_closest_date(arr, elm, option):

    """find_closest() tries to find the closest value to a given x and returns the index of that value
    arr: a 1D array,
    elm: a value for which closest value exists in array arr
    option: -1: smallest possible match
             0: Any match as long as it is least
             1: Largest possible match"""

    diff = np.array(map(lambda x: (x - elm).days, arr))

    msk = np.empty_like(diff, bool)

    if(option == 0):
        diff = abs(diff)
        msk[:] = False

    elif(option == -1):
        msk = diff > 0

    elif(option == 1):
        msk = diff < 0

    else:
        raise ValueError("Expected only 0, -1 or 1 for option")
 
    diff_ma = np.ma.masked_array(diff, mask = msk)

    return np.argmin(diff_ma)


# Weighted area averaging
def wgt_area_avg(data, lat_wgt, lon_wgt):

    """wgt_area_avg() performas weighted area averaging over a geographical area.
    data: data of which last 2 dimensions are lat and lon. Strictly needs to be a masked array
    lat_wgt: weights over latitude of area (usually cos(lat * pi/180))
    lon_wgt: weights over longitude of area (usually 1)
    
    Returns, Numpy array with 2 less dimensions. 
    Masking features are preserved for consistence"""

    # Get data shape
    shp = data.shape
    ndims = data.ndim
    
    if(isinstance(lat_wgt, float)):
        lat_wgt = [lat_wgt] * shp[ndims - 2]
    if(isinstance(lon_wgt, float)):
        lon_wgt = [lon_wgt] * shp[ndims - 1]

    lat_wgt = np.array(lat_wgt).reshape(len(lat_wgt), 1)
    lon_wgt = np.array(lon_wgt)


    # Make grid of lon_wgt, lat_wgt with lat and lon coordinates (last 2 axis of data)
    wy = np.broadcast_to(lon_wgt, data.shape[ndims - 2:ndims])
    wx = np.broadcast_to(lat_wgt, data.shape[ndims - 2:ndims])

    # Mask the array
    # Get 2D mask from the array
    ds = data[0]
    for el in shp[1:ndims-2]:
        ds = ds[0]
    if(isinstance(ds, np.ma.masked_array)):
        msk = ds.mask
    else:
        msk = False

    wy = np.ma.masked_array(wy, msk)
    wx = np.ma.masked_array(wx, msk)
    
    data_wgt = data * wy * wx
    sm_wgt = data_wgt.sum(axis = (ndims - 2, ndims - 1))
    sm_wgt = sm_wgt/np.sum(wy * wx)

    if(isinstance(ds, np.ma.masked_array)):
        return sm_wgt
    else:
        return np.array(sm_wgt)

# Weighted area sum for geospatial computations
def wgt_area_sum(data, lat_wgt, lon_wgt):

    """wgt_area_sum() performas weighted area addition over a geographical area.
    data: data of which last 2 dimensions are lat and lon. Strictly needs to be a masked array
    lat_wgt: weights over latitude of area (usually cos(lat * pi/180))
    lon_wgt: weights over longitude of area (usually 1)
    
    Returns, Numpy array with 2 less dimensions (Masked array. 
    Mask is False if no mask was supplied with the input data. 
    Else mask is derived from the input data)"""

    # Get data shape
    shp = data.shape
    ndims = data.ndim
    
    if(isinstance(lat_wgt, float)):
        lat_wgt = [lat_wgt] * shp[ndims - 2]
    if(isinstance(lon_wgt, float)):
        lon_wgt = [lon_wgt] * shp[ndims - 1]

    lat_wgt = np.array(lat_wgt).reshape(len(lat_wgt), 1)
    lon_wgt = np.array(lon_wgt)


    # Make grid of lon_wgt, lat_wgt with lat and lon coordinates (last 2 axis of data)
    wy = np.broadcast_to(lon_wgt, data.shape[ndims - 2:ndims])
    wx = np.broadcast_to(lat_wgt, data.shape[ndims - 2:ndims])

    # Mask the array
    # Get 2D mask from the array
    ds = data[0]
    for el in shp[1:ndims-2]:
        ds = ds[0]
    if(isinstance(ds, np.ma.masked_array)):
        msk = ds.mask
    else:
        msk = False

    wy = np.ma.masked_array(wy, msk)
    wx = np.ma.masked_array(wx, msk)
    
    data_wgt = data * wy * wx
    sm_wgt = data_wgt.sum(axis = (ndims - 2, ndims - 1))
    # sm_wgt = sm_wgt/np.sum(wy * wx)

    return sm_wgt


# Geospatial volume calculations
def wgt_volume(data, lat_wgt, lon_wgt, lat_res = 1., lon_res = 1.):

    """wgt_volume() performs volume calculation over a geographical area.
    data: data of which last 2 dimensions are lat and lon. Strictly needs to be a masked array
    data represents depth in m
    lat_wgt: weights over latitude of area (usually cos(lat * pi/180))
    lon_wgt: weights over longitude of area (usually 1)
    lat_res: [Optional] Latitude resolution. Defaults to 1.
    lon_res: [Optional] Longitude resolution. Defaults to 1.
    
    Returns, Numpy array with 2 less dimensions. Units are in m^3"""
    
    # 1deg lat/lon is assumed to be 111 km
    dist_deg = 111e3
    
    # Calculate the bounds of the rectangle at equator
    lat_dist = dist_deg * lat_res
    lon_dist = dist_deg * lon_res

    # Get data shape
    shp = data.shape
    ndims = data.ndim
    
    if(isinstance(lat_wgt, float)):
        lat_wgt = [lat_wgt] * shp[ndims - 2]
    if(isinstance(lon_wgt, float)):
        lon_wgt = [lon_wgt] * shp[ndims - 1]

    lat_wgt = np.array(lat_wgt).reshape(len(lat_wgt), 1)
    lon_wgt = np.array(lon_wgt)


    # Make grid of lon_wgt, lat_wgt with lat and lon coordinates (last 2 axis of data)
    wy = np.broadcast_to(lon_wgt, data.shape[ndims - 2:ndims])
    wx = np.broadcast_to(lat_wgt, data.shape[ndims - 2:ndims])

    # Mask the array
    # Get 2D mask from the array
    ds = data[0]
    for el in shp[1:ndims-2]:
        ds = ds[0]
    msk = ds.mask

    wy = np.ma.masked_array(wy, msk)
    wx = np.ma.masked_array(wx, msk)
    
    # Calculate in m^3
    data_wgt = data * wy * wx * lat_dist * lon_dist
    sm_wgt = data_wgt.sum(axis = (ndims - 2, ndims - 1))
    # sm_wgt = sm_wgt/np.sum(wy * wx)

    return sm_wgt


# Latitude averages of data
def lat_avg(data, lat_wgt):
    
    """Perform latitude average of data:
    Inputs: 
    data - n dimensional spatial data. The last 2 dimensions are assumed to lat and
    lon respectively
    lat_wgt - weights by latitudes"""
    
    lat_shape = lat_wgt.shape
    data_shape = data.shape
    
    # If one dimensional:
    if(len(lat_wgt) == 1):
        lat_wgt_re = np.broadcast_to(lat_wgt, data.shape)
        
    elif(len(lat_shape) > 1):
        raise ValueError ("lat_wgt must be 1 dimensional latitude weights")

    else:
        lat_2d = np.broadcast_to(lat_wgt.reshape(len(lat_wgt), 1), data_shape[-2:])
        lat_wgt_re = np.broadcast_to(lat_2d, data_shape)
        
    return (data * lat_wgt_re).mean(axis = -2)


def subset_latlon(data, lat, lon, lat_lim, lon_lim):
    """Latitude longitude subsetting for netCDF4 data similar to the coordinate
    subscripting avialable in NCL
    Inputs:
        data - geospatial data. Last 2 dimensions assumed to be lat and lon
        lat - latitude coordinate array (1D)
        lon - longitude coordinate array (1D)
        lat_lim - limits for latitude; a 2 element list
        lon_lim - limits for longitude; a 2 element list
    Outputs:
        subsetted data, lat and lon""" 
    
    from netCDF4 import Variable
    
    # check for last 2 dimensions    
    shp = data.shape

    # Check for data types
    # ~ (latnd and lonnd) is the same as checking for (~latnd and ~lonnd)
    if( not (isinstance(lat, np.ndarray) or isinstance(lon, np.ndarray))):
        raise ValueError("Latitude or Longitude coordinates not a numpy array")
    if(not (isinstance(data, np.ndarray) or isinstance(data, Variable))):
        raise ValueError("Data not a numpy array or masked array")
    
    if(len(lat.shape) > 1):
        raise ValueError("Latitude coordinates must be 1D")
    if(len(lon.shape) > 1):
        raise ValueError("Longitude coordinates must be 1D")
    if(shp[-2] != len(lat)):
        raise ValueError("Dimension mismatch along axis: Latitude")
    elif(shp[-1] != len(lon)):
        raise ValueError("Dimension mismatch along axis: Longitude")
    
    if((len(lat_lim) != 2) or (len(lon_lim) != 2)):
        raise ValueError("Latitude or longitude bounds are not 2 element list/array")
    
    # Lat and longitude filtering
    lat_filt = np.logical_and(lat <= lat_lim[1], lat >= lat_lim[0])
    lon_filt = np.logical_and(lon <= lon_lim[1], lon >= lon_lim[0])
    
    if(np.sum(lat_filt) == 0) :
        raise ValueError("Coordinate out of bounds: Latitude")
    if(np.sum(lon_filt) == 0) :
        raise ValueError("Coordinate out of bounds: Longitude")
    
    # First filtering along lat and then along lon using
    # np.take slower and more resource heavy
    # data_sub = np.take(data,np.where(lat_filt)[0],axis = -2)
    # data_sub = np.take(data_sub, np.where(lon_filt)[0], axis = -1)
    data_sub = data[...,np.where(lat_filt)[0],np.where(lon_filt)[0]]
    
    lat_sub = lat[lat_filt]
    lon_sub = lon[lon_filt]
    
    return data_sub, lat_sub, lon_sub

if(__name__ == "__main__"):
    print("Import Module and run!")