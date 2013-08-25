'''
Created on Nov 1, 2012

@author: eliotpbrenner
'''
import unittest

import sys
sys.path.insert(0, '../src')
import emissionProbabilityCalculator as epc

class Test(unittest.TestCase):


    def setUp(self):
        k,l=2,2
        N=100
        gamma=1e-3
        eta=1e-3
        probabilityCalculatorAssociatedToEta = epc.emissionProbabilityCalculator(eta, k, l, N)
        probabilityCalculatorAssociatedToEta.setGamma(gamma)
        self.functionProbabilityOfEmissionByP_eta_Robbins = probabilityCalculatorAssociatedToEta.RobbinsEstimateOfEmissionProbabilityTimesCharFunctionOfTauMinusGamma


    def tearDown(self):
        del(self.functionProbabilityOfEmissionByP_eta_Robbins)


    def testName(self):
        self.failUnlessAlmostEqual(self.functionProbabilityOfEmissionByP_eta_Robbins(0.5,0.5,1e-1),0)
        self.failUnlessAlmostEqual(self.functionProbabilityOfEmissionByP_eta_Robbins(0.5,0.5,1e-2),0.0010163942209161392)
        self.failUnlessAlmostEqual(self.functionProbabilityOfEmissionByP_eta_Robbins(0.5,0.5,1e-3),0.00093502678060530602)
        self.failUnlessAlmostEqual(self.functionProbabilityOfEmissionByP_eta_Robbins(0.5,0.5,1e-4),0.00092080067221769139)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()