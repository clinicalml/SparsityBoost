'''
monteCarloBeta.src.probabilityDistribution.py

Structure holding the m parameters of joint probability distribution
over cartesian product of two sets of cardinality k,l respectively.
Created on Aug 8, 2012
@author: eliotpbrenner
'''

import numpy as np
import informationTheory as it
tolerance = 1e-10
import probabilityDistributionFactory as pdf 
import scipy.special

def sumOfLogParametersOfDistributionAssociatedToType(aType,N):
    return np.sum(np.log(np.array(aType))) - len(aType)*np.log(N)

class ProbabilityDistribution(object):
    
    #initialize to uniform distribution
    def __init__(self, first_valset_cardinality=2, second_valset_cardinality=2):   #initialize parameters to specific distribution, uniform is default
        num_entries = first_valset_cardinality*second_valset_cardinality
        self.distribution = (1.0/num_entries)*np.ones((first_valset_cardinality, second_valset_cardinality))
    
    #true iff the parameters are equal within tolerance    
    def __eq__(self, other):
        return np.allclose(self.distribution, other.distribution, atol = tolerance)
    
    def min_parameter(self):
        return np.min(self.distribution)
    
    def max_parameter(self):
        return np.max(self.distribution)
    
    def k(self):
        return np.shape(self.distribution)[0]
    
    def l(self):
        return np.shape(self.distribution)[1]
    
    def m(self):
        return self.k()*self.l() 
    
    def first_variable_marginals(self):
        return np.sum(self.distribution, axis = 1)
    
    def second_variable_marginals(self):
        return np.sum(self.distribution, axis = 0)
    
    def flattenedMProjectionParameters(self):
        pA0, pA1 = self.first_variable_marginals()
        pB0, pB1 = self.second_variable_marginals()
        return np.array([pA0*pB0,pA0*pB1,pA1*pB0,pA1*pB1])
    
    def setParametersProductDistAssocWithMarginals(self, firstMarginalsMat, secondMarginalsMat):
        """
        Sets parameters of the probability distribution so that the distribution
        is the product distribution with specified marginals
        The marginals are given as matrices of 1 row and k, l entries respectively,
        and the output is k by l.
        """
        self.distribution = np.array(np.transpose(firstMarginalsMat)*secondMarginalsMat)
        
        
    def params_of_product_distribution_sharing_marginals(self):
        product_dist = np.zeros_like(self.distribution)
        row_marginals = self.first_variable_marginals()
        col_marginals= self.second_variable_marginals()
        for row in range(np.shape(product_dist)[0]):
            for column in range(np.shape(product_dist)[1]):
                product_dist[row][column] = row_marginals[row]*col_marginals[column]
        return product_dist
 
    
    def KL_divergence_as_base(self, p_prime_distribution):        
        """
        For another probability distribution P^{\prime} specified by the parameter list gives the KL Divergence
                                H(P^{\prime} | P))
        """
        p_distribution = self.distribution
        if np.shape(p_distribution) != np.shape(p_prime_distribution):
            raise Exception('p_distribution has shape' + str(np.shape(p_distribution)) + 
                            'while p_prime_distribution has shape' + str(np.shape(p_prime_distribution)))
        num_rows = np.shape(p_prime_distribution)[0]
        num_cols = np.shape(p_prime_distribution)[1]
        term_array = np.zeros_like(p_prime_distribution)
        for row in range(num_rows):
            for col in range(num_cols):
                term_array[row][col] = it.KL_div_term(p_prime_distribution[row][col], self.distribution[row][col])
        return np.sum(term_array)
 

    def emissionProbabilityFromP_eta_ofProductLikeTypeSizeN(self, p_gamma_distribution, N):
        """
        Estimate of probability of emission of p_gamma from p_eta (self): as determined by Robbins' sharpening of Stirling's formula
        (exercise ?? of Csiszar-Korner)
        """       
        m = self.m()
        KL_DivergenceTermOfLog = -N*self.KL_divergence_as_base(p_gamma_distribution.distribution)
        N_TermOfLog = -((m-1)/(2.0))*np.log(2.0*np.pi*N)
        sumLogParametersTermOfLog = -(0.5)*p_gamma_distribution.sumOfLogParameters()
        return np.exp(KL_DivergenceTermOfLog + N_TermOfLog + sumLogParametersTermOfLog)
    
    def sumOfLogParameters(self):
        """
         The sum of the logs of the parameters (needed to calculate the "g" function, for example)
         """
        return np.sum(self.logParametersFlattened())

    def logParametersFlattened(self):
        """
        the log of the parameters flattened (row major order)
        """
        return np.log(self.distribution).flatten()
    
    def KLDivergenceOfP_gammaFromDist(self, gamma, k=2, l=2):
        p_gamma  = pdf.probabilityDistributionFactory(k, l).get_p_eta(gamma)
        return self.KL_divergence_as_base(p_gamma.distribution)
    
    def exactEmissionProbability(self, aType,N):
        """Type is m-tuple of integers"""
        #if self.m() != 4 or len(aType) != 4:
        #    raise Exception("Not implemented.")
        return np.exp((scipy.special.gammaln(N+1) - sum(scipy.special.gammaln(np.add(aType,[1,1,1,1])))) + np.dot(self.logParametersFlattened(), aType))
    
    def RobbinsEstimatedEmissionProbability(self, aType, N):
        """
        Upper bound coming from Robbins' sharpening
        of Sterling's formula (for esimating the cardinality |T|)
        exp(-N(H(p_T || p)))*(2*pi*N)^(-(|X|-1)/2)*(product of p_T's parameters)^(-1/2)
        If any parameters of T (equivalently p_T) are 0, revert to exactEmissionProbability
        """
        if any(elt == 0 for elt in aType):
            return self.exactEmissionProbability(aType, N)
        #print (self.k(), -1)
        #print np.reshape(aType, (self.k(), -1))
        #print -N*self.KL_divergence_as_base((1.0/N)*np.reshape(aType, (self.k(),-1))) 
        #print (0.5)*(-self.m() + 1)*(np.log(2*np.pi*N))
        #print (-0.5)*sumOfLogParametersOfDistributionAssociatedToType(aType, N)
        return np.exp(-N*self.KL_divergence_as_base((1.0/N)*np.reshape(aType, (self.k(),-1))) 
                      +(0.5)*(-self.m() + 1)*(np.log(2*np.pi*N))
                      +(-0.5)*sumOfLogParametersOfDistributionAssociatedToType(aType, N)
            )
        
    
    def tau(self):
        flattened_parameters = self.distribution.flatten()
        return np.dot(flattened_parameters, (np.log(flattened_parameters) - np.log(self.flattenedMProjectionParameters())))
    