'''
monteCarloBeta.src.informationTheory.py
Basic information and probability theoretic functions.
Created on Aug 7, 2012
@author: eliotpbrenner
'''

import numpy as np

def KL_div_term(x,y):
    """
    One term in the KL divergence
    """
    if np.fabs(x) < 1e-30:
        return 0
    else:
        return x*np.log(x/y)
    
def ChernoffRadius(CentralProbability, N):
    """
    According to Chernoff bounds (version with Taylor estimate of KL-divergence), 
    if the Chernoff radius is greater than the returned t
    then the probability of empirical sequence of size N lying within 
    t of the expected value is at least CentralProbability
    """
    return np.sqrt(-np.log(1-CentralProbability)/N)



