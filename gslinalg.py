# ======================================================
# ======================================================
# Filename: misc.py
# Purpose: Miscellaneous functions. See docstrings for 
# more details
# ======================================================
# ======================================================


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
    
# Calculating confidence intervals using t-test
def cnfint(xstd,nx,clev=0.95,tail=1):
    """ Calculates confidence interval using t-test.
        xstd: numpy array of standard deviation of the ensemble.
        nx: numpy array of number of members of the ensemble.
        clev: the confidence interval (default 0.95 => 95%).
        tail: single of double-sided tail. (default 1).
        return value: the width of the interval (symmetric in both up and down directions) is returned."""

    if tail==1:
        cl=clev
    else:
        cl=1. - (1-clev)/2.0

    sdbar = xstd/np.sqrt(nx);
    tinv_val = st.t.ppf(cl,nx-1);
    diff = tinv_val*sdbar;

    return diff


if(__name__ == "__main__"):
    print "Import module and run"