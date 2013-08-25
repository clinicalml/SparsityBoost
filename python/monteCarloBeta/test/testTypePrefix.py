'''
Created on Mar 21, 2013

@author: eliotpbrenner
'''
import unittest
import sys
sys.path.insert(0, '../src')
import typePrefix as tp


class TestTypePrefix(unittest.TestCase):


    def testHasChildren(self):
        rootPrefix = tp.typePrefix(3)  #root prefix is [] and has children [0],[1],[2],[3]
        self.failUnless(rootPrefix.hasChildren())
        fullPrefix = tp.typePrefix(3,[3,0,0,0])
        self.failIf(fullPrefix.hasChildren())
    
    def testProcess(self):
        rootPrefix = tp.typePrefix(3)
        totalNumberOfTypesLen4Size3 = rootPrefix.numDescendants()
        globalList = []
        rootPrefix.DFSProcess(globalList.append)
        ExpectedGlobalList = [[0, 0, 0, 3], [0, 0, 1, 2], [0, 0, 2, 1], [0, 0, 3, 0], 
                              [0, 1, 0, 2], [0, 1, 1, 1], [0, 1, 2, 0], [0, 2, 0, 1], 
                              [0, 2, 1, 0], [0, 3, 0, 0], [1, 0, 0, 2], [1, 0, 1, 1], 
                              [1, 0, 2, 0], [1, 1, 0, 1], [1, 1, 1, 0], [1, 2, 0, 0], 
                              [2, 0, 0, 1], [2, 0, 1, 0], [2, 1, 0, 0], [3, 0, 0, 0]]
        self.failUnlessEqual(globalList, ExpectedGlobalList)
        self.failUnlessEqual(len(globalList), totalNumberOfTypesLen4Size3) #20
        globalList = []
        fullPrefix = tp.typePrefix(3,[3,0,0,0])
        fullPrefix.DFSProcess(globalList.append)
        self.failUnlessEqual(globalList,[[3,0,0,0]])
        self.failUnless(len(globalList),fullPrefix.numDescendants()) #1
        
    def testChildrenListLastEntry_k_Mod_m(self):
        rootPrefix = tp.typePrefix(100)
        childrenWhichAre2Mod10 = rootPrefix.childrenListLastEntry_k_Mod_m(2, 10)
        prefixes = [prefix.data for prefix in childrenWhichAre2Mod10]
        self.failUnlessEqual(prefixes, [[2], [12], [22], [32], [42], [52], [62], [72], [82], [92]] )
        
if __name__ == "__main__":
    unittest.main()