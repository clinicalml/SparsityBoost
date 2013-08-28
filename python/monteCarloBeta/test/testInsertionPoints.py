#!/system/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
Created on June 10, 2013

@author: eliotpbrenner
'''

import sys
sys.path.insert(0, '../src')

import unittest
import insertionPoints as iP
from pprint import pprint

class Test(unittest.TestCase):

    def setUp(self):
		self.entrySet = [100,200,300,400,1000,2000,10000]
		self.increasingList = [1,100,100,200,200,200,300,1000,1000,1000,2000,20000]

    def tearDown(self):
        self.entrySet = None
        self.increasingList = None


    def testAccountForType(self):
		result = iP.firstIndexHash(self.entrySet,self.increasingList)
		self.assertEqual(result,{100: 1, 200: 3, 300: 6, 400: 7, 1000: 7, 2000: 10, 10000: 11})

    def test_fine_le_withBounds(self):
		a=[1,2,3,4,1,2,3,4,1,2,3,4]
		x=2.6
		self.failUnlessEqual(iP.find_le_withBounds(a, x,4,8),(2,6))

if __name__ == "__main__":
    unittest.main()
