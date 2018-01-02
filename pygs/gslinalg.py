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
    

# Spatial kernel autocorrelation function
def spatial_autocorr_kernel(data, k = 1):
    
    """Function to perform local spatial autocorrelation using a kernel of 
    specified size
    Input:
        data - Spatio-temporally varying data. First axis is assumed to be time. Last 2 axes are assumed to latitude and longitude
        k - k is a parameter for kernel size. The actual size of the kernel is 2k+1. Always > 1
    Output:
        Correlation matrix of size nlat x nlon - Each point is the mean correlation of one point and its neighbours defined by the kernel
        of size (2k+1) x (2k+1)"""
        
    shp = data.shape
    
    if len(shp)!= 3:
        raise ValueError("spatial_autocorr_kernel: Input array needs to be 3D")
    ntime = shp[0]
    nlat = shp[1]
    nlon = shp[2]
    
    # Normalizing the data for calculating correlations
    data = (data - data.mean(axis = 0))/data.std(axis = 0)
    
    corr_mat = np.empty([nlat, nlon], np.float32)
    
    # Padding on the y-axis
    lat_pad = np.zeros([ntime, k, nlon])
    data_pad = np.concatenate((lat_pad, data, lat_pad), axis = 1)
    
    # Padding the y-axis
    nlat_new = nlat + 2 * k
    lon_pad = np.zeros([ntime, nlat_new, k])
    data_pad = np.concatenate((lon_pad, data_pad, lon_pad), axis = 2) 

    
    # At each grid
    for i in range(k,nlat+k):
        for j in range(k,nlon+k):
            
            # Kernel movement for a single point
            acc = 0
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