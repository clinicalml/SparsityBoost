'''
Created on Aug 10, 2012

@author: eliotpbrenner
'''

import sys
sys.path.insert(0, '../src')

import unittest
import numpy as np
import informationTheory as it
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

class Test(unittest.TestCase):
        
    def test_KL_div_term(self):
        np.testing.assert_almost_equal(it.KL_div_term(0, 1),0,decimal=10)
        np.testing.assert_almost_equal(it.KL_div_term(.3, .7), -0.254189358116,decimal=10)
        

    def testChernoffRadius(self):
        for N in 100*np.array(range(1,11)):
            for CentralProbability in 0.1*0.68*np.array(range(1,11)):
                t = it.ChernoffRadius(CentralProbability, N)
                if N == 100 and CentralProbability == 0.68:
                    np.testing.assert_almost_equal(t,0.106744287116,decimal=10)
                if N == 100 and CentralProbability == 0.476:
                    np.testing.assert_almost_equal(t, 0.0803905214973,decimal=10)
                if N == 900 and CentralProbability == 0.68:
                    np.testing.assert_almost_equal(t,0.03558142903863639,decimal=10)
                logging.info("""With N=%s and CentralProbability=%s, the Chernoff radius is %s"""%(N,CentralProbability, t))

     
if __name__ == "__main__":
    unittest.main()