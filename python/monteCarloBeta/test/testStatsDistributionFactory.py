'''
Created on Mar 18, 2013

@author: eliotpbrenner
'''
import sys
sys.path.insert(0, '../src')
import unittest
import statsDistributionFactory as sdf

class Test(unittest.TestCase):


    def setUp(self):
        gamma, eta, N = 0.01, 0.1, 30
        self.statsDistFact = sdf.statsDistributionFactory(gamma)
        self.statsDistFact.set_eta(eta)
        self.statsDistFact.set_N(N)
        

    def tearDown(self):
        pass


    def testCalculateAndSetScaleRatio(self):
        uniformMarginals = [0.5,0.5]
        self.statsDistFact.CalculateAndSetScaleRatio(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.16348715431721675)
        self.statsDistFact.set_N(10*self.statsDistFact.N)
        self.statsDistFact.CalculateAndSetScaleRatio(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.017720649896901121)
    
    def testCalculateAndSetScaleRatioRobustToGammaExceedingEta(self):
        uniformMarginals = [0.5,0.5]
        self.statsDistFact.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.16348715431721675)
        self.statsDistFact.set_N(10*self.statsDistFact.N)
        self.statsDistFact.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.017720649896901121)
        self.statsDistFact.gamma = 0.099
        self.statsDistFact.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.05448988305890981)
        self.statsDistFact.gamma = 0.11
        self.statsDistFact.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.056455307910827257)
        self.statsDistFact.gamma = 0.2
        self.statsDistFact.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)
        self.failUnlessAlmostEqual(self.statsDistFact.scaleRatio, 0.056455307910827257)
        self.failUnlessEqual(self.statsDistFact.gamma, 0.2)  #check that the procedure didn't alter gamma
        
    def testGaussianCenteredAttGammaPlusOrtEtafromMarginals(self):
        uniformMarginals = [0.5,0.5]
        distribution = self.statsDistFact.GaussianCenteredAttGammaPlusfromMarginals(uniformMarginals)
        self.failUnlessAlmostEqual(distribution.mean(), 0.0352962850452)
        self.failUnlessAlmostEqual(distribution.std(), 0.0352962850452)
        distribution = self.statsDistFact.GaussianCenteredAttGammaPlusOrtEtafromMarginals(uniformMarginals)
        self.failUnlessAlmostEqual(distribution.mean(), 0.0352962850452)
        self.failUnlessAlmostEqual(distribution.std(), 0.0352962850452)
        self.statsDistFact.gamma = 0.2
        distribution = self.statsDistFact.GaussianCenteredAttGammaPlusOrtEtafromMarginals(uniformMarginals)
        self.failUnlessAlmostEqual(distribution.mean(), 0.10989731312990136)
        self.failUnlessAlmostEqual(distribution.std(), 0.10989731312990136)
        self.failUnlessEqual(self.statsDistFact.gamma,0.2)
        self.statsDistFact.gamma = 0.1
        distribution = self.statsDistFact.GaussianCenteredAttGammaPlusOrtEtafromMarginals(uniformMarginals)
        self.failUnlessAlmostEqual(distribution.mean(), 0.10989731312990136)
        self.failUnlessAlmostEqual(distribution.std(), 0.10989731312990136)
    
    def testGaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio(self):
        uniformMarginals = [0.5,0.5]
        self.statsDistFact.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)
        distribution = self.statsDistFact.GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio(uniformMarginals)
        self.failUnlessAlmostEqual(distribution.mean(), 0.0352962850452)
        self.failUnlessAlmostEqual(distribution.std(), 0.0115409784)
        
    def testGaussianCenteredAttGammaPlusfromMarginalsScaledDependingOn_N(self):
        uniformMarginals = [0.5,0.5]
        distribution = self.statsDistFact.GaussianCenteredAttGammaPlusfromMarginalsScaledDependingOn_N(uniformMarginals)
        self.failUnlessAlmostEqual(distribution.mean(), 0.0352962850452)
        self.failUnlessAlmostEqual(distribution.std(), 0.0115409784)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()