'''
Created on Aug 8, 2012

@author: eliotpbrenner
'''


max_t_comparison_tolerance = 1e-10
import probabilityDistributionPathFactory as pdpf
import probabilityDistributionFactory as pdf
import probabilityCalculator as pc

class emissionProbabilityCalculator(object):
    '''
    For construction of objects giving the exact and estimated probability of emission of a probability distribution or type
    from a fixed distribution
                        p^\eta
    The probability calculator object is one of the arguments of the constructors.
    '''



    def __init__(self,eta,k,l, N):
        '''
        Constructor: probabilityCalculatorObject has an __init__ method of a probabilityCalculatorObject which takes (eta,k,l) as arguments
        '''
        self.probabilityCalculatorObject = pc.ProbabilityCalculator(eta, k, l)  #point of this is to avoid reconstructing p_ta repeatedly
        self.p_etaDistribution = self.probabilityCalculatorObject.underlyingDistribution
        self.N = N
        self.k = k
        self.l = l
        self.m = k*l
        self.gamma = None
        
    def setGamma(self, gamma):
        self.gamma = gamma 
        
    def setN(self, N):
        self.N = N
    
    def RobbinsEstimateOfEmissionProbability(self, firstMarginal, secondMarginal, t):
        """
        Evaluate the function estimating, from above, the probability of emission
        of a type of size N "close to" 
        the the probability distribution parameterized by the triple 
        (firstMarginal, secondMarginal, t)
        Before calling the actual evaluation of the probabilityCalculatorObject, 
        checks that the three parameters (firstMarginal, secondMarginal, t) 
        actually parameterize a valid probability distribution
        Only implemented for binary/binary (k=l=2 case) yet.
        """
        N = self.N
        #checking the marginals are within bounds or probability simplex: in case not or if the parameter t
        #is out of bounds we return 0
        if firstMarginal > 1 - max_t_comparison_tolerance or firstMarginal < 0 + max_t_comparison_tolerance:
            return 0
        if secondMarginal > 1 - max_t_comparison_tolerance or secondMarginal < 0 + max_t_comparison_tolerance:
            return 0
        pathBasedAtMarginals = pdpf.probabilityDistributionPathFactory([firstMarginal, secondMarginal], self.k,self.l).construct()
        max_t = pathBasedAtMarginals.t_max
        min_t = pathBasedAtMarginals.t_min
        if t > max_t - max_t_comparison_tolerance or t < min_t + max_t_comparison_tolerance:
            return 0
        p_gammaDistribution = pdf.probabilityDistributionFactory(self.k, self.l).distributionWrappingParameters(pathBasedAtMarginals.distribution_at_t(t))
        return self.probabilityCalculatorObject.emissionProbabilityFromP_eta_ofProductLikeTypeSizeN( p_gammaDistribution, N)
    
    
    def RobbinsEstimateOfEmissionProbabilityTimesCharFunctionOfTauMinusGamma(self, firstMarginal, secondMarginal, t):
        """
        Evaluate the function estimating, from above, the probability of emission of a type of size N "close to" 
        the the probability distribution parameterized by the triple (firstMarginal, secondMarginal, t)
        Before calling the actual evaluation of the probabilityCalculatorObject, checks that
        the three parameters (firstMarginal, secondMarginal, t) actually parameterize a valid probability distribution
        Only implemented for binary/binary (k=l=2 case) yet
        """
        N = self.N
        #checking the marginals are within bounds or probability simplex: in case not or if the parameter t
        #is out of bounds we return 0
        if firstMarginal > 1 - max_t_comparison_tolerance or firstMarginal < 0 + max_t_comparison_tolerance:
            return 0
        if secondMarginal > 1 - max_t_comparison_tolerance or secondMarginal < 0 + max_t_comparison_tolerance:
            return 0
        pathBasedAtMarginals = pdpf.probabilityDistributionPathFactory([firstMarginal, secondMarginal], self.k,self.l).construct()
        max_t = pathBasedAtMarginals.t_max
        min_t = pathBasedAtMarginals.t_min
        if t > max_t - max_t_comparison_tolerance or t < min_t + max_t_comparison_tolerance:
            return 0
        KLDivergenceAt_t = pathBasedAtMarginals.KL_divergence_at_t(t)
        if KLDivergenceAt_t > self.gamma:
            return 0
        else:
            p_gammaDistribution = pdf.probabilityDistributionFactory(self.k, self.l).distributionWrappingParameters(pathBasedAtMarginals.distribution_at_t(t))
            return self.probabilityCalculatorObject.emissionProbabilityFromP_eta_ofProductLikeTypeSizeN( p_gammaDistribution, N)  
        
        """
        KLDivergenceAt_t = pathBasedAtMarginals.KL_divergence_at_t(t)
        if KLDivergenceAt_t > self.gamma:
            return 0
        else:
             return 1
        """