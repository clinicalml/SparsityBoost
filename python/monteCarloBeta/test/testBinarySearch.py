# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 14:54:45 2013

@author: eliotpbrenner
"""

'''
Created on Aug 10, 2012

@author: eliotpbrenner
'''

import sys
from nose import tools
sys.path.insert(0, '../src')
import math
import unittest
import binarySearch as bS
import logging
logging.basicConfig(format='%(levelname)s:%(filename)s:%(message)s', level=logging.INFO)

logging.basicConfig(
    filename = fileName,
    format = "%(levelname) -10s %(asctime)s %(module)s:%(lineno)s %(funcName)s %(message)s",
    level = logging.DEBUG
)

class testBinarySearch(unittest.TestCase):
        
    def testSqrtBinarySearch(self):        
        newBinarySearchObject = bS.binarySearch(tolerance_exp=7,maxDepth = 10000, increasingFunction=math.sqrt)
        x_hat = newBinarySearchObject.search(1,9,2,0)
        self.assertAlmostEqual(x_hat,4,7)
    
    @tools.raises(ValueError)
    def testCatchingNoFunct(self):        
        newBinarySearchObject = bS.binarySearch(tolerance_exp=7, maxDepth=1000)
        newBinarySearchObject.search(1,9,2,0)

    @tools.raises(ValueError)
    def testNonInterval(self):        
        newBinarySearchObject = bS.binarySearch(tolerance_exp=7, maxDepth=1000, increasingFunction=math.sqrt)
        newBinarySearchObject.search(9,1,2,0)
        
       
    def testMaxDepth(self):        
        newBinarySearchObject = bS.binarySearch(tolerance_exp=7, maxDepth=500, increasingFunction=math.sqrt)
        self.failUnlessAlmostEqual(newBinarySearchObject.search(0,1,2,0),1,200)
    
    @tools.raises(ValueError)
    def testNonIncreasingFunction(self):
        def f(posReal):
            return -math.sqrt(posReal)            
        newBinarySearchObject = bS.binarySearch(tolerance_exp=7, maxDepth=500, increasingFunction=f)
        newBinarySearchObject.search(1,9,2,0)
        
 
if __name__ == "__main__":
    unittest.main()