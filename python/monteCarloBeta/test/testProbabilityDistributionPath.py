'''
Created on Oct 22, 2012

@author: eliotpbrenner
'''


import sys
sys.path.insert(0, '../montecarlobeta')
from nose import tools

import unittest
import probabilityDistributionPathFactory as pdpf
import probabilityDistribution as pd
import probabilityDistributionPath as pdp
import probabilityDistributionFactory as pdf

class Test(unittest.TestCase):


    def setUp(self):
        self.factory = pdpf.probabilityDistributionPathFactory([0.1,0.9], 2, 2)
        self.path = self.factory.construct()
        self.factoryUniform = pdpf.probabilityDistributionPathFactory([0.5,0.5], 2, 2)
        self.pathUniform = self.factoryUniform.construct()
        self.pathUniform.markP_eta(0.01)
        self.distributionFactory = pdf.probabilityDistributionFactory(2,2)
        
    def tearDown(self):
        del(self.factory)
        del(self.factoryUniform)

    def testMinAndMax_t(self):
        self.failUnlessAlmostEqual(self.path.t_max, 0.0099999999999999985)
        self.failUnlessAlmostEqual(self.path.t_min, -0.089999999999999983)
    
    def testKL_DivergenceCalcs(self):
        self.failUnlessAlmostEqual(self.path.KL_divergence_at_t(0.005), 0.0018225000889761831)
        self.failUnlessAlmostEqual(self.path.KL_divergence_at_max_t(), 0.011134087132719479)
        self.failUnlessAlmostEqual(self.path.KL_divergence_at_min_t(), 0.32508297339144726)
    
    def testFinding_t_from_KL_Divergence(self):
        self.failUnlessAlmostEqual(
        self.path.t_at_specified_divergence_from_base_pos_t_orMax_t(0.0018225000889761831), 0.005)
        self.failUnlessAlmostEqual(
        self.path.t_at_specified_divergence_from_base_pos_t_orMax_t(0.0018225000889761831+.1), 0.0099999999999999985)
        self.failUnlessAlmostEqual(
        self.path.smallestNeg_t_atWhichKLDivergenceFromBaseIsLessThanEta(0.005), -0.0099156144151186665 )
        
        
    @tools.raises(ValueError)
    def testRaisingExceptionForEtaTooLarge(self):
         self.pathUniform.t_at_specified_divergence_from_base_neg_t(1.0)
    
    @tools.raises(ValueError)
    def testRaisingExceptionForEtaTooLarge(self):
         self.pathUniform.distribution_at_t(1.0)
    
    def test_t_at_specified_divergence_from_base_pos_t_orMax_t(self):
         t=self.pathUniform.t_at_specified_divergence_from_base_pos_t_orMax_t(1.0)
         self.failUnlessAlmostEqual(t,0.25,7)         
         
    #@unittest.skip("demonstrating skipping")    
    def testLengthOfSegment(self):
        self.failUnlessAlmostEqual(self.path.lengthOfSegmentofKLDivergenceLessThanSpecified(0.005), 0.017617308260382991 )
    
    def testMarkedDistributions(self):
        self.failUnlessAlmostEqual(self.pathUniform.tOfMarkedDistribution(), 0.035296285045184561)
        self.failUnlessAlmostEqual(self.pathUniform.convertTauToKLDivergenceFromMarkedDistribution(1e-3), 0.0047119020864278272)


    @tools.raises(NotImplementedError)
    def testRaiseNotImplementedError(self):
        k,l = 3,3
        prob_dist = pd.ProbabilityDistribution(k,l)
        probDistPath = pdp.probabilityDistributionPath(prob_dist)
        probDistPath.tOfMarkedDistribution()
    
    def test_t_atSpecified_KL_DivergenceFromMarkedDistribution(self):
        tauSpecified = 1e-3
        self.failUnlessAlmostEqual(self.pathUniform.t_at_specifiedDivergenceFromMarkedDistInDirectionOfBase(tauSpecified), 0.024206270787431289)
        recoveredDistribution = self.pathUniform.distribution_at_t_as_distribution(self.pathUniform.t_at_specifiedDivergenceFromMarkedDistInDirectionOfBase(tauSpecified))
        markedDistribution = self.pathUniform.markedProbabilityDist
        self.failUnlessAlmostEqual(markedDistribution.KL_divergence_as_base(recoveredDistribution.distribution), tauSpecified)

    def test_t_at_specifiedDivergenceFromMarkedDistAwayFromBase(self):
        tauSpecified = 1e-3
        self.failUnlessAlmostEqual(self.pathUniform.t_at_specifiedDivergenceFromMarkedDistAwayFromBase(tauSpecified), 0.046339227999196209)
        recoveredDistribution = self.pathUniform.distribution_at_t_as_distribution(self.pathUniform.t_at_specifiedDivergenceFromMarkedDistAwayFromBase(tauSpecified))
        markedDistribution = self.pathUniform.markedProbabilityDist
        self.failUnlessAlmostEqual(markedDistribution.KL_divergence_as_base(recoveredDistribution.distribution), tauSpecified)

    @tools.raises(ValueError)
    def test_exceptional_t_at_specifiedDivergenceFromMarkedDistAwayFromBase(self):
        tauSpecified = 1.0
        self.failUnlessAlmostEqual(
            self.pathUniform.t_at_specifiedDivergenceFromMarkedDistAwayFromBase(
                tauSpecified), 0.046339227999196209)
    
    def testTau(self):
        eta0 = 0.001
        eta1 = 0.69
        pAtKLDivergence0 = self.distributionFactory.distributionWrappingParameters(
            self.pathUniform.distribution_at_specified_divergence_from_base_pos_t(eta0))
        pAtKLDivergence1 = self.distributionFactory.distributionWrappingParameters(
            self.pathUniform.distribution_at_specified_divergence_from_base_pos_t(eta1))
        self.failUnlessAlmostEqual(pAtKLDivergence0.tau(), eta0)
        self.failUnlessAlmostEqual(pAtKLDivergence1.tau(), eta1)
    
        
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()