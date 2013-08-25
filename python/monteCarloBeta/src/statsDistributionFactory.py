'''
Created on Aug 12, 2012


Contains methods of creating objects of type
<class 'scipy.stats.distributions.rv_frozen'>
from persistent and transient data
    
    

@author: eliot
'''
import probabilityDistributionPathFactory as pdpf
from scipy import stats
import functionAlgorithms as fa
import emissionProbabilityCalculator as epc



import logging

# create logger with 'mc_application'
logger = logging.getLogger('mc_application')
logger.setLevel(logging.ERROR)

# create file handler which logs even debug messages
import os
homeDirectory = os.getenv("HOME", "/Users/eliotpbrenner")
fh = logging.FileHandler(homeDirectory + '/Documents/sontag/logs/mcLog.log')
fh.setLevel(logging.ERROR)

# create console handler with a higher log level
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(fh)

#helper functions
def t_gammaPlusMinus_l_gamma(marginalPair, gamma, k,l):
    """
    Input: k = |Val(A)|, l=|Val(B)|, marginalPair = [firstMarginal, secondMarginal]
    Return list of t_gamma_plus, t_gamma_minus and the "length" l_gamma := t_gamma_plus - t_gamma_minus """
    probDistPath = pdpf.probabilityDistributionPathFactory(marginalPair, k, l).construct()
    t_gamma_plus = probDistPath.largestPos_t_atWhichKLDivergenceFromBaseIsLessThanEta(gamma)
    t_gamma_minus = probDistPath.smallestNeg_t_atWhichKLDivergenceFromBaseIsLessThanEta(gamma)
    return t_gamma_plus, t_gamma_minus, t_gamma_plus - t_gamma_minus #length of relevant segment
    
def constructIntegrandFunctionObject(marginalPair, eta, k,l,N):
    """
    Input: marginalPair, e.g. [0.5,0.3]
           eta, k,l, N
    Output:
    integrand function which object maps t to Estimated emission probability of p(baseMarginals)(t) from p^\eta
    and which has some additional methods/properties, facilitating binary search, associated with it.
    """
    #functionFromTwoMarginalsAndParameterToIntegrand(x,y,t) = Estimated emission probability of p(x,y,t) from p^\eta
    functionFromTwoMarginalsAndParameterToIntegrand = epc.emissionProbabilityCalculator(eta, k, l, N).RobbinsEstimateOfEmissionProbability
    
    
    t=0.001  #used for logger
    logger.info("For marginals %s, Robbins function takes value %s at t=%s"%(
            marginalPair, functionFromTwoMarginalsAndParameterToIntegrand(marginalPair[0], marginalPair[1], t), t))
    functionFromParameterToIntegrandObject = fa.functionAlgorithms(functionFromTwoMarginalsAndParameterToIntegrand)
    
    functionFromParameterToIntegrandObject.setFixedArgumentList(marginalPair)
    logger.info("Fixed argument list set to %s"%(str(functionFromParameterToIntegrandObject.fixedArgumentList)))
    #functionFromParameterToIntegrand(t) = Estimated emission probability of p(baseMarginals)(t) from p^\eta
    functionFromParameterToIntegrand = functionFromParameterToIntegrandObject.functionOfOneVariable
    logger.info("As func. of one variable, takes value %s at t=%s"%(functionFromParameterToIntegrand(t), t))
    
    return fa.functionAlgorithms(functionFromParameterToIntegrand)
    
def calculateScaleRatio(t_gamma_plus, computedScale, scriptL_gamma ):
    """
    Input parameters:
    t_gamma_plus: t_parameter at boundary of A_{0}^\gamma
    computedScale:  a point between t_gamma_minus and t_gamma_plus:
    scriptL_gamma = t_gamma_plus - t_gamma_minus
    Output:
    scaleRatio:  (t_gamma_plus - computedScale)/scriptL_gamma
    """
    #for scaleInputtoStatsNorm we compute actual width of gaussian
    scaleInputToStatsNorm = t_gamma_plus - computedScale  
    logger.info("Normal dist. constructed with loc=%s and scale=%s"%(
             t_gamma_plus, scaleInputToStatsNorm))
    #for scaleRatio we take the ratio of scaleInputToStatsNorm and divide by scriptL_gamma 
    scaleRatio = scaleInputToStatsNorm/scriptL_gamma
    logger.info("Normal dist. constructed with loc=%s and scaleRatio=%s"%(
                t_gamma_plus, scaleRatio))
    return scaleRatio

class statsDistributionFactory(object):
    '''
    Contains methods of creating objects of type
    <class 'scipy.stats.distributions.rv_frozen'>
    from persistent and transient data
    '''
    
    def __init__(self,gamma):
        '''
        Store persistent data
        '''
        self.gamma = gamma
        self.eta = None
        self.N = None
        self.k = 2
        self.l = 2
        self.theProportion = stats.norm.pdf(1)/stats.norm.pdf(0)
        self.scaleRatio = None   #a persistent scale_ratio that is used for all marginals...if 0.3, 
        #say then the scale is always 0.3*length of relevant segment of path
        
      
    def set_eta(self, eta):
        self.eta = eta
        
    def set_N(self, N):
        self.N = N
    
    def set_gamma(self, gamma):
        self.gamma = gamma
    
    def evaluateMethodAtMinOfEtaAndGamma(self, marginals, aMethod):
        """
        Evaluate aMethod with gamma set to the min of gamma and eta
        If gamma changes in the process, restore gamma to its original value
        """
        if self.gamma <= self.eta:
            return aMethod(marginals)
        else:                #Gamma is larger than eta
            oldGamma = self.gamma   #store real value of gamma
            self.set_gamma(self.eta)     #temporarily replace gamma with the smaller value eta
            distribution = aMethod(marginals)
            self.set_gamma(oldGamma)   #restore real value of gamma 
            return distribution
        
    
    def CalculateAndSetScaleRatio(self, baseMarginals):
        '''
        baseMarginals: list of the form (say) [0.5,0.5] or [0.1,0.5]
        Calculate and set the scaleRatio based on the baseMarginals: presupposes that gamma, eta, and N have already been set
        '''
        #constants
        k,l=self.k, self.l  #currently hardcoding for binary random variables
        gamma, eta, N = self.gamma, self.eta, self.N #for readability
        firstMarginal, secondMarginal = baseMarginals
    
        logger = logging.getLogger('mc_application')
        logger.setLevel(logging.ERROR)
        logger.info('Computing the scale ratio with gamma=%s, eta=%s, N=%s, baseMarginals=%s'%(gamma,eta,N,str(baseMarginals)))
        
        #STEP 1: set up function objects: integrandFunctionObject
        integrandFunctionObject = constructIntegrandFunctionObject([firstMarginal, secondMarginal], eta, k,l,N)
        
        #STEP 2: compute scriptL_gamma: length of segment of probability distribution path based at p(baseMarginals) inside the set A_0^\gamma
        t_gamma_plus, t_gamma_minus, scriptL_gamma  = t_gammaPlusMinus_l_gamma([firstMarginal, secondMarginal], gamma, k,l)
        
        #STEP 3: Compute scale
        logger.info("Searching for where the function is proportion %s of the max between %s and %s"%(
                self.theProportion, t_gamma_minus, t_gamma_plus))
        computedScale = integrandFunctionObject.searchArgWhereIncreasingFunctionTakesProportionOfMaxVal(
                self.theProportion, t_gamma_minus, t_gamma_plus)
        logger.info("For marginals (%s,%s), N=%s, computed scale is %s"%(firstMarginal, secondMarginal, N, computedScale))
        
        #STEP 4: transform computedScale to scaleRatio
        self.scaleRatio =  calculateScaleRatio(t_gamma_plus, computedScale, scriptL_gamma )
    
    def CalculateAndSetScaleRatioRobustToGammaExceedingEta(self, baseMarginals):
        '''
        Calculate and set the scaleRatio based on the baseMarginals: presupposes that gamma, eta, and N have already been set
        Uses t_eta_plus in place of t_gamma if gamma exceeds eta
        '''
        self.evaluateMethodAtMinOfEtaAndGamma(baseMarginals, self.CalculateAndSetScaleRatio)


    
    def GaussianCenteredAttGammaPlusfromMarginals(self, marginals):
        """
        Inputs a (2-element) list of marginals, outputs a <class 'scipy.stats.distributions.rv_frozen'>
        """
        #constants
        k,l = self.k, self.l
        
        probabilityDistPathBasedAtMarginals = pdpf.probabilityDistributionPathFactory(marginals, k, l).construct()
        logger.info("For marginals=%s, gamma=%s, about to compute t_gamma_plus"%(str(marginals), self.gamma))
        t_gamma_plus = probabilityDistPathBasedAtMarginals.t_at_specified_divergence_from_base_pos_t_orMax_t(self.gamma)
        logger.info('For gamma=%s, marginals =%s, returning normal dist of loc and scale = %s'%(self.gamma, str(marginals), t_gamma_plus))
        return stats.norm(loc=t_gamma_plus, scale=t_gamma_plus)
    
    def GaussianCenteredAttGammaPlusOrtEtafromMarginals(self, marginals):
        """
        Inputs a (2-element) list of marginals, outputs a <class 'scipy.stats.distributions.rv_frozen'>
        Differs from GaussianCenteredAttGammaPlusfromMarginals because it uses eta instead of gamma if gamma>eta,
        i.e. the min(gamma,eta) instead of gamma
        """
        return self.evaluateMethodAtMinOfEtaAndGamma(marginals, self.GaussianCenteredAttGammaPlusfromMarginals)

    def GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio(self, marginals):
        """
        This procedure does NOT set the scale ratio.  It assumes the scale ratio has already been set.
        Inputs a (2-element) list of marginals, outputs a <class 'scipy.stats.distributions.rv_frozen'>
        Currently, k=l=2 (binary-binary) hardcoded into this method.  Assumed that gamma, eta, and N have been set
        """
        if self.scaleRatio is None:
            raise TypeError("The scale ratio of the statsDistributionFactory has to be set.")
            
        #constants
        k,l, gamma, eta, N = self.k, self.l, self.gamma, self.eta, self.N
        firstMarginal, secondMarginal = marginals
        
        logger = logging.getLogger('mc_application')
        logger.setLevel(logging.ERROR)
        logger.debug('Recovering the scale from the scale ratio with gamma=%s, eta=%s, N=%s, marginals=%s'%(gamma,eta,N,str(marginals)))
        
        t_gamma_plus, t_gamma_minus, scriptL_gamma = t_gammaPlusMinus_l_gamma(marginals, gamma, k,l)
        scaleInputToStatsNorm = scriptL_gamma*self.scaleRatio
        logger.debug('With scaleRatio =%s and scriptL_gamma=%s, scaleInputToStatsNorm=%s'%(self.scaleRatio, scriptL_gamma, scaleInputToStatsNorm))
        logger.debug("For marginals (%s,%s), N=%s, normal dist. constructed with loc=%s and scale=%s"%(firstMarginal, secondMarginal, N, t_gamma_plus, scaleInputToStatsNorm))
        return stats.norm(loc=t_gamma_plus, scale=scaleInputToStatsNorm)
        
    def GaussianCenteredAtMinOftGammaPlusOrEtafromMarginalsScaledFromScaleRatio(self, marginals):
        """
        This procedure does NOT set the scale ratio.  It assumes the scale ratio has already been set.
        Inputs a (2-element) list of marginals, outputs a <class 'scipy.stats.distributions.rv_frozen'>
        Currently, k=l=2 (binary-binary) hardcoded into this method.  Assumed that gamma, eta, and N have been set
        """
        return self.evaluateMethodAtMinOfEtaAndGamma(marginals, self.GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio)


    def GaussianCenteredAttGammaPlusfromMarginalsScaledDependingOn_N(self, marginals):
        """
        Inputs a (2-element) list of marginals, outputs a <class 'scipy.stats.distributions.rv_frozen'>
        Currently, k=l=2 (binary-binary) hardcoded into this method.  Assumed that gamma, eta, and N have been set
        """
        self.CalculateAndSetScaleRatioRobustToGammaExceedingEta(marginals)
        return self.GaussianCenteredAtMinOftGammaPlusOrEtafromMarginalsScaledFromScaleRatio(marginals)
        
        
      