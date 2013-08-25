'''
Created on May 25, 2012

@author: eliotpbrenner

A path p(t) of distributions in the space of joint distributions
over  pair of variables, all of which share fixed marginals, as in the paper/thesis.
'''

max_iterations_for_binary_search = 40000
tolerance = 1E-10
toleranceForPGammaRecovery = 1E-15

import probabilityDistribution as pd
import numpy as np
import binarySearch as bS



class probabilityDistributionPath(object):
    '''
    Assumes that the path is of the form p_0 + t* Delta, where 
    p_0 is a product distribution and Delta is a 
    checkerboard pattern of +1 and -1
    '''

    def __init__(self, probability_distribution_on_path):
        self.base_probability_distribution = pd.ProbabilityDistribution()
        self.base_probability_distribution.distribution = probability_distribution_on_path.params_of_product_distribution_sharing_marginals()
        self.t_max = self.t_max()
        self.t_min = self.t_min()
        #for situations where it is convenient to mark a special p_eta as the base for calculating KL-divergnece:
        self.markedProbabilityDist = None
        self.binarySearchOnKLDivergence = bS.binarySearch(
        tolerance_exp=10,maxDepth=500,increasingFunction=self.KL_divergence_at_t)
        self.binarySearchOnMinusKLDivergence = bS.binarySearch(
        tolerance_exp=10,maxDepth=500,increasingFunction= self.minus_KL_divergence_at_t)

    def t_max(self):
        base_parameters = self.base_probability_distribution.distribution
        num_rows = np.shape(base_parameters)[0]
        num_cols = np.shape(base_parameters)[1]
        t_max_found = 1.0
        for row in range(num_rows - num_rows%2):
            for col in range(num_cols - num_cols%2):
                if (row+col)%2 == 0:
                    t_max_this_coord = 1 - base_parameters[row][col]
                else:
                    t_max_this_coord = base_parameters[row][col]
                if t_max_this_coord < t_max_found:
                    t_max_found = t_max_this_coord
        return t_max_found
    
    def t_min(self):
        base_parameters = self.base_probability_distribution.distribution
        num_rows = np.shape(base_parameters)[0]
        num_cols = np.shape(base_parameters)[1]
        t_max_found = 1.0
        for row in range(num_rows - num_rows%2):
            for col in range(num_cols - num_cols%2):
                if (row+col)%2 == 0:
                    t_max_this_coord = base_parameters[row][col]
                else:
                    t_max_this_coord = 1 - base_parameters[row][col]
                if t_max_this_coord < t_max_found:
                    t_max_found = t_max_this_coord
        return -t_max_found
    
    
    def distribution_at_t(self,t):
        """return the array for the distribution at t:
         Assumes that the path is of the form p_0 + t* Delta, where p_0 is a product distribution and Delta is a checkerboard pattern of +1 and -1
        """
        t_max = self.t_max
        t_min = self.t_min
        if t > t_max + tolerance or t < t_min - tolerance:
            raise ValueError('t= ' + str(t) + ' but t_max=' + str(t_max) + ' and t_min=' + str(t_min)) 
        base_parameters = self.base_probability_distribution.distribution
        num_rows = np.shape(base_parameters)[0]
        num_cols = np.shape(base_parameters)[1]
        distribution_at_t_parameters=np.copy(base_parameters)
        num_changing_rows = num_rows - num_rows%2
        num_changing_cols = num_cols - num_cols%2  #should probably store these somewhere so they don't have to be recomputed
        for row in range(num_changing_rows):
            for col in range(num_changing_cols):
                distribution_at_t_parameters[row][col] += (-1)**(row+col)*t
        return distribution_at_t_parameters
    
    
    def distribution_at_t_as_distribution(self,t):
        """return a distribution object version of the distribution_at_t"""
        k,l = np.shape(self.base_probability_distribution.distribution)
        new_distribution = pd.ProbabilityDistribution(k,l)
        new_distribution.distribution = self.distribution_at_t(t)
        return new_distribution
    
    def KL_divergence_at_t(self,t):
        base_distribution = self.base_probability_distribution
        distribution_at_t_for_this_path = self.distribution_at_t(t)
        return base_distribution.KL_divergence_as_base(distribution_at_t_for_this_path)
    
    def minus_KL_divergence_at_t(self,t):
        return -self.KL_divergence_at_t(t)
    
    def KL_divergence_at_max_t(self):
        return self.KL_divergence_at_t(self.t_max)
    
    def KL_divergence_at_min_t(self):
        return self.KL_divergence_at_t(self.t_min)
  
    def distribution_at_specified_divergence_from_base_pos_t(self, eta):   
        """This only looks in the direction of "positive t" from the base to find a distribution of the specified KL-divergence from the
        "base", where base means product distribution
        """
        upper_limit = self.t_max-tolerance
        lower_limit = 0
        KL_divergence_at_current_t = self.KL_divergence_at_t(upper_limit)
        if (eta > KL_divergence_at_current_t):
            raise ValueError("KL_divergence_at_t_max= " + str(KL_divergence_at_current_t) + " is less than eta= " + str(eta))
        return self.distribution_at_t(self.binarySearchOnKLDivergence.search(lower_limit,upper_limit,eta))
     
    def distribution_at_speicified_divergence_from_base_pos_t_as_distribution(self,eta):
        upper_limit = self.t_max-tolerance
        lower_limit = 0
        KL_divergence_at_current_t = self.KL_divergence_at_t(upper_limit)
        if (eta > KL_divergence_at_current_t):
            raise ValueError("KL_divergence_at_t_max= " + str(KL_divergence_at_current_t) + " is less than eta= " + str(eta))
        return self.distribution_at_t_as_distribution(self.binarySearchOnKLDivergence.search(lower_limit,upper_limit,eta))
    

    def distribution_at_specified_divergence_from_base_neg_t(self, eta):
        """
        This only looks in the direction of "negative t" from the base to find a distribution of the specified KL-divergence from the
        "base", where base means product distribution
        """
        upper_limit = -self.t_min-tolerance  #positive upper limit on the abs. val of t to search
        lower_limit = 0  #lower limit on abs val
        current_guess_t = upper_limit
        #print current_guess_t
        KL_divergence_at_current_t = self.KL_divergence_at_t(-upper_limit)
        if (eta > KL_divergence_at_current_t):
            raise ValueError("KL_divergence_at_t_max= " + str(KL_divergence_at_current_t) + " is less than eta= " + str(eta))
        return self.distribution_at_t(self.binarySearchOnKLDivergence.search(-upper_limit,lower_limit,eta))



    def t_at_specified_divergence_from_base_pos_t(self, eta):
        """This only looks in the direction of "positive t" from the base to find a distribution of the specified KL-divergence from the
        "base", where base means product distribution
        """
        upper_limit = self.t_max-tolerance
        lower_limit = 0
        KL_divergence_at_current_t = self.KL_divergence_at_t(upper_limit)
        if (eta > KL_divergence_at_current_t):
            raise ValueError("KL_divergence_at_t_max= " + str(KL_divergence_at_current_t) + " is less than eta= " + str(eta))
        return self.binarySearchOnKLDivergence.search(lower_limit,upper_limit,eta)        
 
 
    def t_at_specified_divergence_from_base_pos_t_orMax_t(self, eta):
        """
        This only looks in the direction of "positive t" from the base to find a distribution of the specified KL-divergence from the
        "base", where base means product distribution: if eta is greater than the KL_divergence at t_max, returns t_max-tolerance
        """
        upper_limit = self.t_max-tolerance
        lower_limit = 0
        current_guess_t = upper_limit
        #print current_guess_t
        KL_divergence_at_current_t = self.KL_divergence_at_t(current_guess_t)
        if (eta > KL_divergence_at_current_t):
            return current_guess_t #which equals self.t_max-tolerance in this case
        return self.binarySearchOnKLDivergence.search(lower_limit,upper_limit,eta)
    
    
    
    """
    This only looks in the direction of "negative t" from the base to find a distribution of the specified KL-divergence from the
    "base", where base means product distribution
    """
    def t_at_specified_divergence_from_base_neg_t(self, eta):
        upper_limit = -self.t_min-tolerance  #positive upper limit on the abs. val of t to search
        lower_limit = 0
        current_guess_t = upper_limit
        #print current_guess_t
        KL_divergence_at_current_t = self.KL_divergence_at_t(-current_guess_t)
        if (eta > KL_divergence_at_current_t):
            raise ValueError("KL_divergence_at_t_max= " + str(KL_divergence_at_current_t) + " is less than eta= " + str(eta))
        return self.binarySearchOnMinusKLDivergence.search(-upper_limit,lower_limit,-eta)
    
    def lengthOfSegmentofKLDivergenceLessThanSpecified(self,eta):
        """
        $\ell_{\gamma}$, in the notation of the paper
        """
        if self.KL_divergence_at_max_t() < eta:
            rightEndpoint = self.t_max
        else:
            rightEndpoint = self.t_at_specified_divergence_from_base_pos_t(eta) 
        if self.KL_divergence_at_min_t() < eta:
            leftEndpoint = self.t_min
        else:
            leftEndpoint = self.t_at_specified_divergence_from_base_neg_t(eta)
        return rightEndpoint - leftEndpoint
    
    def largestPos_t_atWhichKLDivergenceFromBaseIsLessThanEta(self, eta):
        """
        Similar to t_at_specified_divergence_from_base_pos_t(eta), but if the KL-divergence at t_max is
        less than eta, return t_max
        """
        if self.KL_divergence_at_max_t() < eta:
            rightEndpoint = self.t_max
        else:
            rightEndpoint = self.t_at_specified_divergence_from_base_pos_t(eta) 
        return rightEndpoint
    
    def smallestNeg_t_atWhichKLDivergenceFromBaseIsLessThanEta(self, eta):
        """
        Similar to t_at_specified_divergence_from_base_pos_t(eta), but if the KL-divergence at t_min is
        less than eta, return t_min
        """
        if self.KL_divergence_at_min_t() < eta:
            leftEndpoint = self.t_min
        else:
            leftEndpoint = self.t_at_specified_divergence_from_base_neg_t(eta)
        return leftEndpoint
    
    def markP_eta(self, eta):
        '''
        Set p_eta to be the marked distribution
        '''
        t_eta_plus = self.t_at_specified_divergence_from_base_pos_t(eta)
        self.markedProbabilityDist = self.distribution_at_t_as_distribution(
            t_eta_plus)
        
    def convertTauToKLDivergenceFromMarkedDistribution(self, tauValue):
        '''
        If we are given the gamma (=tau) used to form a p_gamma, this will 
        return the KL Divergence from p_eta (marked distribution)
        to p_gamma, i.e. KL(p_\tau^+ \| p_\eta^+)
        '''
        if not self.markedProbabilityDist:
            raise ValueError("No marked distribution.")
        t_tau_plus = self.t_at_specified_divergence_from_base_pos_t(tauValue) 
        distributionAt_t_tau_plus =  self.distribution_at_t_as_distribution(
            t_tau_plus) #p_\tau^+
        return self.markedProbabilityDist.KL_divergence_as_base(
            distributionAt_t_tau_plus.distribution) #KL(p_tau^+ \| p_eta^+)
    
    def tOfMarkedDistribution(self):
        if self.base_probability_distribution.k() != 2 or self.base_probability_distribution.l() != 2:
            raise NotImplementedError("k=%s, l=%s"%(self.base_probability_distribution.k(),self.base_probability_distribution.l()))
        if not self.markedProbabilityDist:
            raise ValueError("No marked distribution.")
        differenceOfParameters = self.markedProbabilityDist.distribution - \
            self.base_probability_distribution.distribution
        return differenceOfParameters[0,0]
        
     
    
    def t_at_specifiedDivergenceFromMarkedDistInDirectionOfBase(self, specifiedDivergence):
        
        """
        This only looks in the direction of negative t from the marked distribution
        to find a distribution of the specified KL-divergence from the
        "marked distribution" with the marked distribution as base.
        The method supposes that the marked distribution is in the 
        positive-t direction from the base:
        See the testProbabilityDistributionPath.
            test_t_atSpecified_KL_DivergenceFromMarkedDistribution
        for more documentation (through the test)
        """
        if not self.markedProbabilityDist:
            raise ValueError("No marked distribution.")
        upper_limit = self.tOfMarkedDistribution()
        lower_limit = 0
        #KL_divergence_at_upper_limit = KL(p^0 \| p^\eta)
        KL_divergence_at_upper_limit = self.markedProbabilityDist.\
            KL_divergence_as_base(  #KL(\cdot \| p^\eta)
        self.distribution_at_t_as_distribution(
            self.tOfMarkedDistribution() - upper_limit #0
            ).distribution) #p^0
        
        if (specifiedDivergence > KL_divergence_at_upper_limit):
            raise ValueError("Specified divergence " + str(specifiedDivergence) + 
            " is greater than KL_divergence of uniform dist based at p^\eta = " 
            + str(KL_divergence_at_upper_limit))
            


        increasingFunct = lambda t : self.markedProbabilityDist.\
            KL_divergence_as_base(  #KL(\cdot \| p^\eta)
        self.distribution_at_t_as_distribution(
            self.tOfMarkedDistribution() - t # t_eta^+ - t
            ).distribution) #p(t_\eta^+ - t)

        self.binarySearchOnKLDivFromMarkedDistribution = bS.binarySearch(
            tolerance_exp=15,maxDepth=500,increasingFunction=increasingFunct)
            
        s = self.binarySearchOnKLDivFromMarkedDistribution.search(
            lower_limit,upper_limit,specifiedDivergence)
        return self.tOfMarkedDistribution() - s
        
        
      
      
    def t_at_specifiedDivergenceFromMarkedDistAwayFromBase(self, specifiedDivergence):
        """
        This only looks in the direction of positive t from the marked distribution
        to find a distribution of the specified KL-divergence from the
        "marked distribution" with the marked distribution as base.
        The method supposes that the marked distribution is in the
        positive-t direction from the base:
        See the 
        testProbabilityDistributionPath.test_t_at_specifiedDivergenceFromMarkedDistAwayFromBase
        for more documentation (through the test)
        """
        if not self.markedProbabilityDist:
            raise ValueError("No marked distribution.")
        lower_limit = self.tOfMarkedDistribution()
        upper_limit = self.t_max
        
        #KL(p(t_max) || p^\eta)
        KL_divergence_at_upper_limit = self.markedProbabilityDist.\
            KL_divergence_as_base(  #KL(\cdot \| p^\eta)
        self.distribution_at_t_as_distribution(
            upper_limit #0
            ).distribution) #p(t_max)
        
        if (specifiedDivergence > KL_divergence_at_upper_limit):
            raise ValueError("Specified divergence " + str(specifiedDivergence) + 
            " is greater than KL_divergence of uniform dist based at p^\eta = " 
            + str(KL_divergence_at_upper_limit))
            

        #increasingFunct(t)=KL(p(t) ||  p^\eta)
        increasingFunct = lambda t : self.markedProbabilityDist.\
            KL_divergence_as_base(  #KL(\cdot \| p^\eta)
        self.distribution_at_t_as_distribution(t).distribution) #p(t)

        self.binarySearchOnKLDivFromMarkedDistribution = bS.binarySearch(
            tolerance_exp=15,maxDepth=500,increasingFunction=increasingFunct)
            
        return self.binarySearchOnKLDivFromMarkedDistribution.search(
            lower_limit,upper_limit,specifiedDivergence)
    
def main():
    pass
    
if __name__ == '__main__':
    main()    

    