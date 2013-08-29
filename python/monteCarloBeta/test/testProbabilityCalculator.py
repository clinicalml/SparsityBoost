'''
Created on Oct 31, 2012

@author: eliotpbrenner
'''
import sys
sys.path.insert(0, '../montecarlobeta')
import unittest
import probabilityCalculator as pc
import probabilityDistributionFactory as pdf

class Test(unittest.TestCase):


    def setUp(self):
        k,l=2,2
        eta=0.01
        gamma=1e-3
        self.probabilityCalculatorObject = pc.ProbabilityCalculator(eta, k, l)
        self.p_gamma_distribution = pdf.probabilityDistributionFactory(2,2).get_p_eta(gamma)
        self.p_eta_distribution = pdf.probabilityDistributionFactory(2,2).get_p_eta(eta)

    def tearDown(self):
        del(self.probabilityCalculatorObject)


    def testProbabilityCalculator(self):
        N=100
        self.failUnlessAlmostEqual(
            self.probabilityCalculatorObject.emissionProbabilityFromP_eta_ofProductLikeTypeSizeN(
            self.p_gamma_distribution, N),0.00063545388776551162)
        N=1000
        self.failUnlessAlmostEqual(
            self.probabilityCalculatorObject.emissionProbabilityFromP_eta_ofProductLikeTypeSizeN(
            self.p_gamma_distribution, N),2.89311883867e-07)
        N=100
        self.failUnlessAlmostEqual(
            self.probabilityCalculatorObject.emissionProbabilityFromP_eta_ofProductLikeTypeSizeN(
            self.p_eta_distribution, N),0.0010365601817726436)
        

if __name__ == "__main__":
    unittest.main()