'''
Created on Jul 29, 2012

@author: eliotpbrenner
'''

import probabilityDistributionPath as pdp
import csv
import logging
import numpy as np
import probabilityDistributionFactory as pdf
import probabilityDistributionPathFactory as pdpf

class ProbabilityCalculator(object):
    '''
    Creates object for repeated calculation of probabilities of emission from a fixed distribution:
    the point of this class is to avoid having to reconstruct p_eta_distribution each time.
    '''


    def __init__(self, eta, k,l):
        '''
        Constructor
        '''
        p_eta_distribution = pdf.probabilityDistributionFactory(k,l).get_p_eta(eta)  
        self.underlyingDistribution = p_eta_distribution
        self.k = k
        self.l = l
        self.m = k*l
        
   
   
    def emissionProbabilityFromP_eta_ofProductLikeTypeSizeN(self, p_gamma_distribution, N):
        """
        Estimate of probability of emission of p_gamma from p_eta (self): 
        as determined by Robbins' sharpening of Stirling's formula
        (derived on p. 39 of Csiszar-Korner's Information Theory: Coding Theorems 
        for Discrete Memoryless Systems)
        p_gamma_distribution is of an object type probabilityDistribution, 
        defined within this python module.
        """       
        m = self.m
        KL_DivergenceTermOfLog = -N*self.underlyingDistribution.KL_divergence_as_base(
            p_gamma_distribution.distribution)
        N_TermOfLog = -((m-1)/(2.0))*np.log(2.0*np.pi*N)
        sumLogParametersTermOfLog = -(0.5)*p_gamma_distribution.sumOfLogParameters()
        return np.exp(KL_DivergenceTermOfLog + N_TermOfLog + sumLogParametersTermOfLog)
        
        
  