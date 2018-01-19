# =======================================================
# =======================================================
# Filename: gsplots.py
# Purpose: To simplify plotting using basemaps and
# matplotlib in python
# Created: Nov 22 2017
# =======================================================
# =======================================================

"""Simplifying plotting results from atmospheric analysis in python"""

# Import required libraries
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from mpl_toolkits.basemap import interp
import numpy as np
import nclcmaps as ncm
import matplotlib.pyplot as plt
import matplotlib.colors as mcl
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
    
def __get_NS__(num):
    
    """Get north or south for latitude
    num: Numeric value of longitude"""

    if((num > 0) & (num <= 90)):
        return "N"
    elif((num < 0) & (num >= -90)):
        return "S"
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
                     minlev, maxlev, levspace, add_cyclic = False,
                     lat_lim = [-90,90], lon_lim = [0, 360],
                     proj = "cyl", drawgrids = None, gridfont = 10, gridlat = [-90,90],
                     gridlon = [0,360],                
                     center_lon = 0, 
                     drawls = False, cmap = "testcmap", ax = None, extend = 'both',
                     conf = None, ms = 0.1, rs = 3, scatter = True):
    
    
    """Plot a contour map on a basemap imported from mpl_toolkits
    contour_data - Data for contour map should be lat x lon
    lats - latitude values
    lons - longitude values
    minlev - minimum contour level
    maxlev - maximum contour level
    levspace - spacing between contours
    extend [optional] - Takes 4 values 'min', 'max', 'both', 'neither'. Used for the contour map. 'both' by default
    add_cyclic [optional] - Add a cyclic point to the plot. Useful for full longitude ranges
    proj [optional] - Projection of the map.
    Takes only three possible values from the range of projections available in matplotlib ['hammer', 'robin', 'cyl'].
    Rectangle region bounds are ignored for both hammer and robin. Only lat_0 is considered. Default is 'cyl'
    drawgrids [optional] - Draw grids on the map at this spacing -  a 2 element tuple. Is only used for non-cylindrical maps. Default None
    gridfont [optional] - Fontsize of grid label
    lat_lim [optional] - Limits for latitude [0, 360] by default. Must be a 2 value array-like
    lon_lim [optional] - Limits for longitude [-90,90] by default. Must be a 2 value array-like
    lonshift [optional] - Shift the longitude to match lon_0. This argument takes the amount by which to shift
    drawls [Optional] - Toggles Landsea masks. Set to False by default
    cmap[Optional] - ncl colormap to use. Default is testcmap, or a matplotlib map or a list of colors
    ax [Optional] - which axis to draw upon
    conf [Optional] - Draw a scatter plot of 1s and 0s signifiying confidence 
    levels"""   
 
    if (len(lat_lim) != 2) or (len(lon_lim)!=2):
        raise ValueError("Only 2 values expected for axis limits")
        
    
    # Data for the contour map    
    clevs = np.arange(minlev,maxlev+levspace/2.,levspace)
    
    if(proj not in ('cyl', 'hammer', 'robin')):
        raise ValueError("plot_contour_map: proj not in one of ['robin', 'hammer', cyl'")
    if(proj == 'cyl'):
        m = Basemap(projection = 'cyl', llcrnrlat = lat_lim[0], \
                    urcrnrlat = lat_lim[1], \
                    area_thresh = 10000,
                    llcrnrlon = lon_lim[0], # urcrnrlon = 360,
                    urcrnrlon = lon_lim[1],
                    resolution = 'l', suppress_ticks = False, ax = ax)
    elif(proj == 'hammer' or proj == 'robin'):
        m = Basemap(projection = proj,
                    area_thresh = 10000,
                    lon_0 = center_lon,
                    resolution = 'l', suppress_ticks = True, ax = ax)
            
    m.drawcoastlines()
    
    nlon = len(lons)
    
    # If lsmask is required    
    if(drawls):
        m.drawlsmask()
    
    
    if(proj != 'cyl'):
        
        # If central longitude is not zero or the central value of the longitude array, the longitude grid 
        # also has to be shifted to accomodate for the new grid structure        
        if(lons[nlon/2] != center_lon): # which means we need to ensure central longitude is set to lat_0
            lon_start = center_lon - 180
            contour_data, lons = shiftgrid(lon_start,contour_data, lons)
            
        # Addcyclic is applied after centering the logitude since it can break the code because of the monotonicity issue
        if(add_cyclic):
            contour_data, lons = addcyclic(contour_data, lons)
        
        x,y = np.meshgrid(lons, lats)
        x,y = m(x,y) # Project meshgrid onto map
        if(drawgrids is not None):
            if(len(drawgrids) != 2):
                raise ValueError("plot_contour_map: Expected 2 values as tuple for drawgrids")

            m.drawparallels(np.arange(gridlat[0], gridlat[1]+1, drawgrids[0]), labels = [1,0,0,0], fontsize = gridfont)
            m.drawmeridians(np.arange(gridlon[0], gridlon[1]+1, drawgrids[1]), fontsize = gridfont)
    else:
        
        if(add_cyclic):
            contour_data, lons = addcyclic(contour_data, lons)
        x,y = np.meshgrid(lons, lats)
        
    
    # contour_data = contour_data[:,lon_pos]
    if(extend not in ['both', 'min', 'max','neither']):
        raise ValueError("plot_contour_map: extend only takes 3 values - ['both', 'min', 'max']")
        
    if(isinstance(cmap, str)):
        cmap = ncm.cmap(cmap)
    
    if(isinstance(cmap, str) or isinstance(cmap, mcl.ListedColormap) or isinstance(cmap, mcl.LinearSegmentedColormap)):
        cs = m.contourf(x, y, contour_data, clevs, cmap = cmap, extend = extend)
    elif(isinstance(cmap, np.ndarray)):
        cs = m.contourf(x, y, contour_data, clevs, colors = cmap, extend = extend)
    else:
        raise ValueError("Invalid Colormap provided")
    
    # Hack for map limits. Note: there needs to be a better way or simply migrate to cartopy
    if(proj != 'cyl'):
        
        # Map projections for latitude
        llcrnrlon, llcrnrlat = m(min(lons), lat_lim[0])
        urcrnrlon, urcrnrlat = m(max(lons), lat_lim[1])
        
        if(ax is None):
            ax = plt.gca()
        
        ax.set_ylim(llcrnrlat, urcrnrlat)
        
        # Drawing map boundary with black background. 
        # Relies on resolution being a little off. Is an unavoidable hack. Need to find something better
        m.drawmapboundary(fill_color = 'k')
    
    # Hatching the plot with significance levels
    if(conf is not None):
        conf[~ conf.mask] = ms
        if(scatter):
            m.scatter(x[::rs,::rs],y[::rs,::rs],conf[::rs,::rs], 
                      marker = '.', color = 'k', edgecolor = None, lw = 0)
        else:
            m.contourf(x, y, conf, n_levels = 2, 
                       hatches = [None, '\\\+///'], colors= 'none', extend = 'lower',
                       linecolor = 'grey')
            
    return m,cs
    

def hov_diagram(lon,y,contour_data,levels,col = "testcmap", ax = None,
                rng = None, step = 60, ew = False, norm = False, 
                extend = 'both', vmin = None, vmax = None):
    
    """Draw hovmoller diagrams (time vs lat)
    Inputs:
    lon, y = x and y axes respectively. lon is usually longitude, y is time,
    contour_data = data to be plotted (lon x time)
    levels = levels for the contour of the hovmoller diagram
    [optional]
    ax = an axis on a plot to draw the plot on
    extend = Which side of the colorbar to extend: Must be one of 'both', 'min', 'max', 'neither'
    ew = Add "E" or "W" tag to longitude axis
    col = ncl colormap. By default test cmap
    rng = minimum axis tick and maximum x axis tick (an array of 2 values)
    defaults to None and choose the min and max of the longitudes passed
    step = steps taken on longitde axis"""
    
    l1 = np.min(lon)
    l2 = np.max(lon)
    
    if(extend not in ['both', 'min', 'max', 'neither']):
        raise ValueError("hov_diagram: Unexpected value for extend")
    
    X,Y = np.meshgrid(lon,y)
    
    if(isinstance(col, str)):
        cmap = ncm.cmap(col)
        if(ax is not None):
            CS = ax.contourf(X,Y,contour_data,levels,cmap=cmap, extend = extend)
        else:
            CS = plt.contourf(X,Y,contour_data,levels,cmap=cmap, extend = extend)
    else:
        if(ax is not None):
            CS = ax.contourf(X,Y,contour_data,levels,colors = col, extend = extend)
        else:
            CS = plt.contourf(X,Y,contour_data,levels,colors = col, extend = extend)
            
    if(ax is None):
        ax = plt.gca()
        
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
        tck_lb = [str(np.int(np.abs(t))) + __get_EW__(t) for t in tck]
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


def fliplon(data, lon):
    """Toggle longitude coordinate between -180 to 180 and 0 to 360
    Formula used: lon1 = (lon+180)%360 - 180
        Check the following link
        http://vikas-ke-funde.blogspot.in/2010/06/convert-longitude-0-360-to-180-to-180.html
    Inputs:
        data - The data with longitude axis to be flipped. Rightmost longitude is assumed to be longitude
        lon - The longitude coordinates as a 1D numpy array
    Outputs:
        data_re - data with toggled coordinates
    """
    
    if(not isinstance(data, np.ndarray)):
        raise ValueError("fliplon: Array passed is not numpy ndarray")

    data_dim = data.shape
    
    if(not isinstance(lon, np.ndarray)):
        raise ValueError("fliplon: longitude coordinate is not a numpy array")
    
    if(len(lon.shape) > 1):
        raise ValueError("fliplon: logitude coordinate is not 1 Dimensional")
    
    nlon = len(lon)
    
    if(data_dim[-1] != nlon):
        raise ValueError("fliplon: Longitude coordinate size mismatch")
    
    # If negative value is present, data is already in the range -180 to 180
    if((np.min(lon) < 0) and (np.max(lon) < 180)):
        flag = 1
    else:
        flag = 0 # Otherwise data is assumed to be in 0-360 range
       
    if(not flag): # if not in -180 to 180 form
        lon1 = (lon + 180) % 360 - 180
    else:
        lon1 = lon % 360 # If it is in -180 to 180 form
    
    lon_sort = np.argsort(lon1)
    
    lon_re = lon1[lon_sort]
    data_re = data[...,lon_sort]
    
    return data_re, lon_re

if(__name__ == "__main__"):
    print "Import module and run"