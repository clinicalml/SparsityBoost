'''
Created on Aug 9, 2012

@author: eliotpbrenner
'''

import probabilityDistribution as pd
import probabilityDistributionPath as pdp

class probabilityDistributionFactory(object):
 
    def __init__(self, k,l):
        self.k = k
        self.l = l
        
        
    def get_p_eta(self,eta):
        '''
        Makes a distribution p_eta with uniform marginals and test statistic eta.
        '''
        l = self.l
        k = self.k
        prob_dist = pd.ProbabilityDistribution(k,l)
        uniform_dist = pd.ProbabilityDistribution(k,l)
        prob_dist_path = pdp.probabilityDistributionPath(uniform_dist)
        prob_dist_parameters = prob_dist_path.distribution_at_specified_divergence_from_base_pos_t(eta)
        prob_dist.distribution = prob_dist_parameters
        return prob_dist
    
    def distributionWrappingParameters(self, parameters):
        wrappingDistribution = pd.ProbabilityDistribution() 
        wrappingDistribution.distribution = parameters
        return wrappingDistribution
