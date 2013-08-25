'''
Created on Aug 13, 2012

@author: eliot
'''
import unittest
from scipy import stats
import sys
sys.path.insert(0, '../src')

import unittest
import numpy as np
import functionAlgorithms as fa
import emissionProbabilityCalculator as epc
import probabilityDistributionPathFactory as pdpf
import probabilityDistributionPath as pdp
import informationTheory as it

import logging
logger = logging.getLogger('myapp')
import os
homeDirectory = os.getenv("HOME", "/Users/eliotpbrenner")
hdlr = logging.FileHandler(homeDirectory + '/logs/monteCarloBeta.log')
formatter = logging.Formatter("%(asctime)s [%(funcName)s: %(filename)s,%(lineno)d] %(message)s")
#formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.CRITICAL)


class testFunctionAlgorithms(unittest.TestCase):


    def testBinarySearch(self):
        newFunctionAlgorithmObject = fa.functionAlgorithms(np.sin)
        np.testing.assert_almost_equal(
          newFunctionAlgorithmObject.searchArgWhereIncreasingFunctionTakesVal(0,
          -np.pi/2.0, np.pi/2.0),0.0)
        np.testing.assert_almost_equal(
          newFunctionAlgorithmObject.searchArgWhereIncreasingFunctionTakesVal(-1,
          -np.pi/2.0, np.pi/2.0), -1.57079632094)
        np.testing.assert_almost_equal(
          newFunctionAlgorithmObject.searchArgWhereIncreasingFunctionTakesVal(1,
          -np.pi/2.0, np.pi/2.0),  1.57079632094)
        
    def testSearchArgWhereIncreasingFunctionTakesProportionOfMaxVal(self):
        theProportion = stats.norm.pdf(1)/stats.norm.pdf(0)
        N_list = list(10*np.array(range(1,10)))
        N_list.extend((100*np.array(range(1,10))))
        k,l=2,2
        eta=0.01
        gamma = 0.001
        logger.debug("Set eta=%s, gamma=%s"%(eta,gamma))
        #firstMarginalsDist = stats.uniform(loc = .4,scale = .2)
        #secondMarginalsDist = stats.uniform(loc = .4,scale = .2)
        for N in N_list:
            for iteration in range(1):
                #firstMarginal, secondMarginal = [firstMarginalsDist.rvs(), secondMarginalsDist.rvs()]
                firstMarginal, secondMarginal = 1.0/l, 1.0/k #0.5, 0.5 for binary-binary
                logger.debug("Randomly chosen marginals: (%s,%s)"%(firstMarginal, secondMarginal))
        
                
                functionFromTwoMarginalsAndParameterToIntegrand = epc.emissionProbabilityCalculator(eta, k, l, N).RobbinsEstimateOfEmissionProbability
                
                
                logger.debug("For marginals (%s,%s), Robbins function takes value %s at t=%s"%(firstMarginal, secondMarginal, functionFromTwoMarginalsAndParameterToIntegrand(firstMarginal, secondMarginal, 0.001), 0.001))
                functionFromParameterToIntegrandObject = fa.functionAlgorithms(functionFromTwoMarginalsAndParameterToIntegrand)
                
                functionFromParameterToIntegrandObject.setFixedArgumentList([firstMarginal, secondMarginal])
                logger.debug("Fixed argument list set to %s"%(str(functionFromParameterToIntegrandObject.fixedArgumentList)))
                functionFromParameterToIntegrand = functionFromParameterToIntegrandObject.functionOfOneVariable
                logger.debug("As func. of one variable, takes value %s at t=%s"%(functionFromParameterToIntegrand(0.001), 0.001))
                
                probDistPath = pdpf.probabilityDistributionPathFactory([firstMarginal, secondMarginal], k, l).construct()
                t_gamma_plus = probDistPath.largestPos_t_atWhichKLDivergenceFromBaseIsLessThanEta(gamma)
                t_gamma_minus = probDistPath.smallestNeg_t_atWhichKLDivergenceFromBaseIsLessThanEta(gamma)
                logger.debug("Robbins function takes value %s at t_gamma_plus=%s"%(functionFromParameterToIntegrandObject.theFunction(firstMarginal, secondMarginal, t_gamma_plus), t_gamma_plus))
                logger.debug("Robbins function takes value %s at t_gamma_minus=%s"%(functionFromParameterToIntegrandObject.theFunction(firstMarginal, secondMarginal, t_gamma_minus), t_gamma_minus))
                
                integrandFunctionObject = fa.functionAlgorithms(functionFromParameterToIntegrand)
                logger.info("Searching for where the function is proportion %s of the max between %s and %s"%(theProportion, t_gamma_minus, t_gamma_plus))
                computedScale = integrandFunctionObject.searchArgWhereIncreasingFunctionTakesProportionOfMaxVal(theProportion, t_gamma_minus, t_gamma_plus)
                logger.info("For marginals (%s,%s), N=%s, computed scale is %s"%(firstMarginal, secondMarginal, N, computedScale))
                

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()