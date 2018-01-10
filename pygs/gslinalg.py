# ======================================================
# ======================================================
# Filename: gslinalg.py
# Purpose: Miscellaneous functions. See docstrings for 
# more details
# ======================================================
# ======================================================

"""Some simple linear algebraic functions for geospatial computations"""


import numpy as np
import scipy.stats as st

# ND matrix correlations
def corr3d_oned(matnd, vec): 
     
    """Calculate correlation between a 3D matrix and a vector 
    In the nD matrix, the first dimension needs to be the dimension 
    on which the correlation is calculated 
    """ 
    
    if(len(vec.shape) != 1): 
        raise ValueError("vec needs to be one dimensional") 
 
    # Calculate mean and standard deviation of  
    ndm = len(matnd.shape)
    mat_std = np.std(matnd, axis=0)
    mat_mean = np.mean(matnd, axis=0)
    tlen = matnd.shape[0]
    matnd = (matnd - mat_mean)

    # Move time axis to end from first position
    matnd = np.rollaxis(matnd, 0, ndm)
 
    # Calculate standard deviation of vector 
    vec_len = len(vec) 
    if(tlen != vec_len): 
        raise ValueError("Time dimension of matrix and lenght of vector \
                d not match, Matrix has " + str(tlen) + " units, vector \
                has " +str(vec_len) + " units") 
    vec_std = np.std(vec) 
    vec = (vec - np.mean(vec))/vec_len 
 
    # Matrix correlation 
    mat_corr = np.sum(matnd * vec, axis = ndm-1)
    mat_corr = mat_corr/(mat_std * vec_std) 
    return mat_corr


# ND matrix correlations
def corr3d(matnd1, matnd2): 
     
    """Calculate correlation between a 3D matrix and another matrix of same
    dimension. 
    In the nD matrices, the first dimension needs to be the dimension 
    on which the correlation is calculated 
    """ 
 
    if(len(matnd1.shape) != len(matnd2.shape)):
        raise ValueError("Both matrices must be of same dimensionality")
 
    # Calculate mean and standard deviation of  
    ndm = len(matnd1.shape)
    mat_std1 = np.std(matnd1, axis=0)
    mat_mean1 = np.mean(matnd1, axis=0)
    tlen = matnd1.shape[0]
    matnd1 = (matnd1 - mat_mean1)

    # Move time axis to end from first position
    matnd1 = np.rollaxis(matnd1, 0, ndm)
 
    # Calculate standard deviation of vector 
    mat_std2 = np.std(matnd2, axis = 0) 
    matnd2 = (matnd2 - np.mean(matnd2, axis = 0))/tlen 
    matnd2 = np.rollaxis(matnd2, 0, ndm)
 
    # Matrix correlation 
    mat_corr = np.sum(matnd1 * matnd2, axis = ndm-1)
    mat_corr = mat_corr/(mat_std1 * mat_std2) 
    return mat_corr


def reg3d(array1, array2):
    """reg3d() performs linear regression of 2 arrays in 3 dimension
    Inputs:
    array1 - a 3D array of spatial data. First axis is time (Usually the regressand)
    array2 - a 1D array of data (Usually regressor)
    
    returns:
    slp_, intc_, r_, p_, sterr_ - 2D arrays of slope, intercept, 
                                Rsquared, p values and standard errors"""
    
    a1 = array1.shape
    a2 = array2.shape
    
    if(len(a1) != 3):
        raise ValueError("Array 1 needs to be 3 dimensional")
    if(len(a2) != 1):
        raise ValueError("Array 2 needs to be 1 dimensional")
    
    def regfunc(a, b):
        sl, inc, rsq, pv, sterr = st.linregress(b, a)
        return sl, inc, rsq, pv, sterr
        
    return np.apply_along_axis(regfunc, 0, array1, array2)
    

def pad_nd(data, k = 1, value = 0.):
    
    """Function to pad nD geospatial like data on all sides with value
    for the size defined by k
    Inputs:
        data: Assumed to be geospatial. Last 2 dimensions are lat x lon
        k [Optional]: The size of the padding on the edges. Defaults to 1 for a 3x3 kernel
        value [Optional]: The value of the padding. Defaults to 0
    """
    
    
    shp = data.shape
    if(not isinstance(data, np.ndarray)):
        raise ValueError("pad_nd: Data not a numpy array")
    elif(len(shp) < 2):
        raise ValueError("pad_nd: Data must have atlease lat and lon dimenions")
    
    lat_pad_shape = list(data.shape)
    lat_pad_shape[-2] = k
    lat_pad = np.zeros(lat_pad_shape)
    data_pad = np.concatenate((lat_pad, data, lat_pad), axis = -2)       
    
    # Padding the y-axis
    lon_pad_shape = list(data_pad.shape)
    lon_pad_shape[-1] = k
    lon_pad = np.zeros(lon_pad_shape)
    data_pad = np.concatenate((lon_pad, data_pad, lon_pad), axis = 2) 
    
    return data_pad


# Spatial kernel autocorrelation function
def spatial_autocorr_kernel(datax, k = 1, datay = None):
    
    """Function to perform local spatial autocorrelation using a kernel of 
    specified size
    Input:
        datax - Spatio-temporally varying data. First axis is assumed to be time. Last 2 axes are assumed to latitude and longitude
        k - k is a parameter for kernel size. The actual size of the kernel is 2k+1. Always > 1
        datay [Optional]: If a cross correlation is required. Can be used in the case of exploring multivariate relationships
    Output:
        Correlation matrix of size nlat x nlon - Each point is the mean correlation of one point and its neighbours defined by the kernel
        of size (2k+1) x (2k+1)"""
        
    shpx = datax.shape
    
    if len(shpx)!= 3:
        raise ValueError("spatial_autocorr_kernel: Input array needs to be 3D")
        
    ntime = shpx[0]
    nlat = shpx[1]
    nlon = shpx[2]
        
    if(datay is not None):
        shpy = datay.shape
        if len(shpy)!= 3:
            raise ValueError("spatial_autocorr_kernel: Target array datay must be 3D")
        ntimey, nlaty, nlony = shpy
        if(ntime != ntimey or nlon != nlony or nlat != nlaty):
            raise ValueError("spatial_autocorr_kernel: Arrays datax and datay mismatched")
        
        # Subtract mean
        datay = (datay - datay.mean(axis = 0))/datay.std(axis = 0)
    
    # Normalizing the data for calculating correlations
    datax = (datax - datax.mean(axis = 0))/datax.std(axis = 0)
    
    corr_mat = np.empty([nlat, nlon], np.float32)
    
    # add padding    
    data_pad = pad_nd(datax, k = k, value = 0.)
    
    # At each grid
    for i in range(k,nlat+k):
        for j in range(k,nlon+k):
            
            # Kernel movement for a single point
            acc = 0
            
            if(datay is not None):
                xij = datay[:,i-k, j-k]
            else:
                xij = data_pad[:,i,j]
            
            # Loop over the kernel
            for l in range(-k, k+1):
                for m in range(-k, k+1):            
                    xshift = data_pad[:,i+l,j+m]
                    corr = np.dot(xij,xshift)/ntime
                    acc = acc + corr
            corr_mat[i-k,j-k] = acc/((2*k+1)**2)
    
    return corr_mat

if(__name__ == "__main__"):
    print "Import module and run"