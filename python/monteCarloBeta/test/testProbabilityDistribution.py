'''
Created on Oct 21, 2012

@author: eliotpbrenner
'''


import sys
sys.path.insert(0, '../src')

import unittest
import probabilityDistribution as pd
import probabilityDistributionFactory as pdf
import numpy as np

class TestProbabilityDistribution(unittest.TestCase):


    def testFirstAndSecondMarginals(self):
        new_distribution = pd.ProbabilityDistribution(3,2)
        new_distribution.distribution = np.array([[0.1,0.15],[0.05,0.20],[.2,.3]])
        np.testing.assert_almost_equal(new_distribution.first_variable_marginals(),np.array([ 0.25 , 0.25,  0.5 ]))      
        np.testing.assert_almost_equal(new_distribution.second_variable_marginals(), np.array([ 0.35,  0.65 ]))
        
    def testProductDistributionFormation(self):
        product_distribution = pd.ProbabilityDistribution(3,2)
        product_distribution.distribution = np.array([[ 0.0875 , 0.1625],[ 0.0875,  0.1625] ,[ 0.175 ,  0.325 ]])
        product_distribution_reconstruct = pd.ProbabilityDistribution(3,2) 
        product_distribution_reconstruct.setParametersProductDistAssocWithMarginals(np.matrix([ 0.25 , 0.25,  0.5 ]), np.matrix([ 0.35,  0.65 ]))
        self.assertEqual(product_distribution,product_distribution_reconstruct) #uses __eq__ defined in the probabilityDistribution class

    
    def testProductDistributionSharingMarginals(self):
        new_distribution = pd.ProbabilityDistribution(3,2)
        new_distribution.distribution = np.array([[0.1,0.15],[0.05,0.20],[.2,.3]])
        product_distribution = pd.ProbabilityDistribution(3,2)
        product_distribution.distribution = np.array([[ 0.0875 , 0.1625],[ 0.0875,  0.1625] ,[ 0.175 ,  0.325 ]])
        product_distribution_reconstruct = pd.ProbabilityDistribution(3,2) 
        product_distribution_reconstruct.distribution = new_distribution.params_of_product_distribution_sharing_marginals()
        self.assertEqual(product_distribution, product_distribution_reconstruct)
    
    def testKL_divergenceAsBase(self):
        new_distribution2 = pd.ProbabilityDistribution(3,2)
        new_distribution2.distribution = np.array([[0.1,0.15],[0.05,0.20],[.2,.3]])
        another_distribution_params = np.array([ [.25,0],  [0,.25],  [.25,.25] ])
        self.assertEqual(new_distribution2.KL_divergence_as_base(another_distribution_params),
                               .25*np.log(.25/.1) + .25*np.log(.25/.20) + .25*np.log(.25/.2) + .25*np.log(.25/.3))
        
    def testEmissionProbabilityFromP_eta_ofProductLikeTypeSizeN(self):
        new_distribution = pd.ProbabilityDistribution(3,2)
        new_distribution.distribution = np.array([[0.1,0.15],[0.05,0.20],[.2,.3]])
        product_distribution = pd.ProbabilityDistribution(3,2)
        product_distribution.distribution = np.array([[ 0.0875 , 0.1625],[ 0.0875,  0.1625] ,[ 0.175 ,  0.325 ]])
        np.testing.assert_almost_equal(new_distribution.emissionProbabilityFromP_eta_ofProductLikeTypeSizeN(product_distribution,100), 4.37171106784e-06)
    
    def testKLDivergenceOfP_gammaFromDist(self):
        new_distribution = pd.ProbabilityDistribution(2,2)
        new_distribution.distribution = np.array([[0.25,0.25],[0.25,0.25]])
        np.testing.assert_almost_equal(new_distribution.KLDivergenceOfP_gammaFromDist(0.001), 0.000999999947487)
        
    def testProbabilityDistributionFactory(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        np.testing.assert_almost_equal(p_eta.distribution, np.array([[ 0.35989731,  0.14010269],[ 0.14010269,  0.35989731]]) ) 
    
    def testExactEmissionProbability(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        np.testing.assert_almost_equal(p_eta.exactEmissionProbability([2,1,1,1],5),0.054900966738391482)
        np.testing.assert_almost_equal(p_eta.exactEmissionProbability([1,2,1,1],5),0.021372132192157535)

    def testRobbinsEstimatedEmissionProbability(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        aType = [200,100,100,100]
        N = 500
        np.testing.assert_almost_equal(p_eta.exactEmissionProbability(aType, N), 2.64474060719e-19)
        np.testing.assert_almost_equal(p_eta.RobbinsEstimatedEmissionProbability(aType, N), 2.65202363049e-19)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testProbabilityDistribution']
    unittest.main()