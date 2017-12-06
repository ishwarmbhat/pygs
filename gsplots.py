# =======================================================
# =======================================================
# Filename: gsplots.py
# Purpose: To simplify plotting using basemaps and
# matplotlib in python
# Created: Nov 22 2017
# =======================================================
# =======================================================

# Import required libraries
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp
import numpy as np
import nclcmaps as ncm
from netCDF4 import Dataset
import os

# Private method for getting E or W labels for longitudes in -180 to 180 range
def __get_EW__(num):
    
    """Get east or west for longitudes
    num: Numeric value of longitude"""

    if((num > 0) & (num < 180)):
        return "E"
    elif((num < 0) & (num > -180)):
        return "W"
    else:
        return ""


# Redefine basemap plots with shapefile from Natural Earth
class Basemap2(Basemap):

    """ Modify Basemap to use Natural Earth data instead of GSHHG data"""

    def drawcoastlines(self):
        import os.path
        libpath = os.path.dirname(os.path.abspath(__file__))
        shapefile = libpath + '/ne/coastlines/ne_%sm_coastline' % \
                    {'l':110, 'm':50, 'h':10}[self.resolution]
        self.readshapefile(shapefile, 'coastline', linewidth=1.)

# Plot contour maps using basemap package
def plot_contour_map(contour_data, lats, lons,
                     minlev, maxlev, levspace,
                     lat_lim = [-90,90], lon_lim = [0, 360],
                     drawls = False, cmap = "testcmap", ax = None, 
                     conf = None, ms = 0.1, rs = 3):
    
    
    """Plot a contour map on a basemap imported from mpl_toolkits
    contour_data - Data for contour map should be lat x lon
    lats - latitude values
    lons - longitude values
    minlev - minimum contour level
    maxlev - maximum contour level
    levspace - spacing between contours
    lat_lim [optional] - Limits for latitude [0, 360] by default. Must be a 2 value array-like
    lon_lim [optional] - Limits for longitude [-90,90] by default. Must be a 2 value array-like
    drawls [Optional] - Toggles Landsea masks. Set to False by default
    cmap[Optional] - ncl colormap to use. Default is testcmap
    ax [Optional] - which axis to draw upon
    conf [Optional] - Draw a scatter plot of 1s and 0s signifiying confidence 
    levels"""   
 
    if (len(lat_lim) != 2) or (len(lon_lim)!=2):
        raise ValueError("Only 2 values expected for axis limits")
    
    m = Basemap(projection = 'cyl', llcrnrlat = lat_lim[0], \
        urcrnrlat = lat_lim[1], \
        area_thresh = 10000,
        llcrnrlon = lon_lim[0], # urcrnrlon = 360,
        urcrnrlon = lon_lim[1],
        resolution = 'l', suppress_ticks = False, ax = ax)
        
    m.drawcoastlines()

    # If lsmask is required    
    if(drawls):
        m.drawlsmask()        
    # Data for the contour map    
    clevs = np.arange(minlev,maxlev+levspace,levspace)
    
    # Private method delta to subtract a small number to account for error in floating point rounding off
    _delta_ = levspace/1000 * 4
        
    x,y = np.meshgrid(lons, lats)
    cmap = ncm.cmap(cmap)
    
    if(np.ma.is_masked(contour_data)):
        contour_data[np.logical_and(~contour_data.mask, contour_data < min(clevs))] = min(clevs) + _delta_
        contour_data[np.logical_and(~contour_data.mask,contour_data > max(clevs))] = max(clevs) - _delta_
    else:
        contour_data[contour_data < min(clevs)] = min(clevs) + _delta_
        contour_data[contour_data > max(clevs)] = max(clevs) - _delta_

    # contour_data = contour_data[:,lon_pos]
    cs = m.contourf(x, y, contour_data, clevs, cmap = cmap, origin = 'lower')
    
    if(conf is not None):
        conf[~ conf.mask] = ms
        m.scatter(x[::rs,::rs],y[::rs,::rs],conf[::rs,::rs], 
                  marker = '.', color = 'k', edgecolor = None, lw = 0)
    
    return m,cs
    

def hov_diagram(ax,lon,y,contour_data,levels,col = "testcmap", 
                rng = None, step = 60, ew = False, norm = False,
                vmin = None, vmax = None):
    
    """Draw hovmoller diagrams (time vs lat)
    Inputs:
    ax = an axis on a plot to draw the plot on
    lon, y = x and y axes respectively. lon is usually longitude, y is time,
    contour_data = data to be plotted (lon x time)
    levels = levels for the contour of the hovmoller diagram
    [optional]
    col = ncl colormap. By default test cmap
    rng = minimum axis tick and maximum x axis tick (an array of 2 values)
    defaults to None and choose the min and max of the longitudes passed
    step = steps taken on longitde axis"""
    
    l1 = np.min(lon)
    l2 = np.max(lon)
    
    X,Y = np.meshgrid(lon,y)

#    contour_data[contour_data > np.max(levels)] = np.max(levels)
#    contour_data[contour_data < np.min(levels)] = np.min(levels)
    
    if(isinstance(col, str)):
        cmap = ncm.cmap(col)
        CS = ax.contourf(X,Y,contour_data,levels,cmap=cmap, origin='lower')
    else:
        CS = ax.contourf(X,Y,contour_data,levels,colors = col, origin='lower')
        

    ax.minorticks_on()
    ax.tick_params(axis = "both", which = "both", direction = "out")
    ax.tick_params(axis = "y", which = "minor", left = "off")
    ax.tick_params(axis = "both", which = "both", top = "off", right = "off")
    
    if (rng is None):
        tck = np.arange(l1, l2, step)
    else:
        if(len(rng) == 2):
            tck = np.arange(rng[0], rng[1], step)
        else:
            raise ValueError("Expected only a 2 value array for rng")
        
    if(ew):
        tck_lb = [str(np.abs(t)) + __get_EW__(t) for t in tck]
    else:
        tck_lb = [str(t) for t in tck]
        
    ax.set_xticks(tck)
    ax.set_xticklabels(tck_lb, fontsize = 7)

    return CS


# Geospatial masking function
def lsmask(data, lat, lon, mask, keep_mask = False):
    """Function to mask data using landsea_mask.nc 
    1x1 resolution file from NCEP
    Inputs: 
        data - data to be masked. n dimensional with lat and lon being last 2 dimension
        lat - latitudes associated with the data
        lon - longitudes associated with the data    
        mask - mask to be applied (can be array)- 
        0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf 
        eg. if mask = 0, ocean will be masked out
        Refer to (https://www.ncl.ucar.edu/Document/Functions/Shea_util/landsea_mask.shtml)
        keep_mask [Optional] - Whether old mask has to be retained on top of new landsea mask
        
    Outputs:
        data - masked array of same dimensionality as input
    """
    
    #  Preload the maskfile
    # NOTE: landsea mask has range between -90 to 90 for lat
    # and 0 to 360 for lon
    min_lat = np.min(lat)
    max_lat = np.max(lat)
    min_lon = np.min(lon)    
    max_lon = np.max(lon)
    
    if(min_lon < 0):
        raise ValueError("lon: Longitude not in range 0:360!")
    
    mloc = os.path.join(os.path.dirname(__file__), 'data/landsea.nc')
    
    f = Dataset(mloc, "r")
    lsm = f.variables["LSMASK"]
    lat_ls = f.variables["lat"][:]
    lon_ls = f.variables["lon"][:]
    
    lat_filt = np.logical_and(lat_ls <= max_lat, lat_ls >= min_lat)
    lon_filt = np.logical_and(lon_ls <= max_lon, lon_ls >= min_lon)
    
    # Sometimes there is an issue with simultaneously subsetting both 
    # lat and lon. CHECK
    # lsm = lsm[lat_filt,:]    
    lsm = lsm[lat_filt,lon_filt]
    
    lat_in = lat_ls[lat_filt]
    lon_in = lon_ls[lon_filt]
    
    lon_out, lat_out = np.meshgrid(lon, lat)
    
    
    f.close()
    
    # Interpolate using the basemap.interp function
    lsm_interp = interp(lsm, xin = lon_in, yin = lat_in, 
                                xout = lon_out, yout = lat_out, order = 3)
    
    # After interpolation, mask the dataset
    shp = data.shape
    
    if(isinstance(mask, int)):
        msk = ~(lsm_interp == mask)
    else:
        msk = np.empty(lsm.shape, np.bool)
        msk[:] = False
        for ms in mask:
            msk = np.logical_or(msk, lsm_interp == ms)
        msk = ~msk
        
    msk3d = np.broadcast_to(msk, shp)
    
    # Mask the input dataset
    if(keep_mask):
        if(isinstance(data, np.ma.masked_array)):
            msk3d = np.logical_or(msk3d, data.mask)
        else:
            raise ValueError("keep_mask: Previous Mask not found")
    return np.ma.masked_array(data, msk3d)


if(__name__ == "__main__"):
    print "Import module and run"