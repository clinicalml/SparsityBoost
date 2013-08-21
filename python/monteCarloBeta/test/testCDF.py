'''
Created on Feb 22, 2013

@author: eliotpbrenner
'''

import sys
sys.path.insert(0, '../src')
import unittest
import CDF as CDF
import probabilityDistributionFactory as pdf
from pprint import pprint
from joblib import Parallel, delayed
import datetime
import pickle
import pandas
import src.typePrefix as tp

def timeStamped(fname, fmt='%Y-%m-%d-%H-%M-%S_{fname}'):
    return datetime.datetime.now().strftime(fmt).format(fname=fname)

class Test(unittest.TestCase):


    def setUp(self):
        self.CDF = CDF.CDF()
        self.CDF2 = CDF.CDF()

    def tearDown(self):
        pass


    def testAccountForType(self):
        """
        When the (partial) CDF accounts for only one, type
          T=[1,1,1,4]
        there is only one discontinuity point, which is locted at tau(T),
        and is a jump in the CDF from 0 to the emission probability of T
        from the reference distribution
        """
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(7)
        self.CDF.accountForType([1,1,1,4])
        self.failUnlessAlmostEqual(CDF.tauOfType([1,1,1,4], 7),0.0427972344694)
        self.failUnlessAlmostEqual(self.CDF.Dictionary, {0.042797234469424295: 0.024888873765445504})
    
    def _testAccountForAllTypesWithTwoElementPrefix(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(5)
        prefix = [2,1]
        remainingMass = 2
        self.CDF.accountForAllTypesWithTwoElementPrefix(prefix, remainingMass)
        pprint(self.CDF.Dictionary)
        self.failUnlessAlmostEqual(self.CDF.Dictionary, {0.013844293808390619: 0.054900966738391482,
                                                         0.11849392256130019: 0.065587032834470232,
                                                         0.2911031660323688: 0.13610213450308134} )
               
    def _testAccountForAllTypes(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(30)
        self.CDF.setn(4)
        self.CDF.accountForAllTypes()
        pprint(self.CDF.Dictionary)
        print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        pprint(self.CDF.AscendingDiscontinuityList)
        print "**********************"
        self.failUnlessEqual(len(self.CDF.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF.assignCumulativeProbability(0.1), 0.470627961298)
        
        self.CDF2.referenceDistribution = p_eta
        self.CDF2.setN(30)
        self.CDF2.setn(4)
        #rootPrefix = tp.typePrefix(30)
        self.CDF2.accountForAllTypes()
        self.failUnlessEqual(len(self.CDF2.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF2.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF2.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF2.assignCumulativeProbability(0.1), 0.470627961298)
        
        """
        self.CDF2.referenceDistribution = p_eta
        self.CDF2.setN(30)
        self.CDF2.accountForTypesByModuloClass10()
        pprint(self.CDF2.Dictionary)
        print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        pprint(self.CDF2.AscendingDiscontinuityList)
        print "**********************"
        self.failUnlessEqual(len(self.CDF2.Dictionary), 2009)
        self.failUnlessAlmostEqual(self.CDF2.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF2.assignCumulativeProbability(0.1), 0.470627961298)
        """
    def _testMerge(self):
        self.CDF.accountForEvent(1, 0.3)
        self.CDF.accountForEvent(2, 0.05)
        self.CDF2.accountForEvent(1.5, 0.05)
        self.CDF2.accountForEvent(2.5, 0.6)
        pprint(self.CDF.Dictionary)
        pprint(self.CDF2.Dictionary)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList),2)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 0.35)
        self.CDF.mergeInto(self.CDF2)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList),4)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        
    def _testaccountForListOfTypes(self):
        self.CDF.setReferenceDistribution_p_eta(0.1)
        self.CDF.setN(30)
        self.CDF.setn(4)
        rootPrefix = tp.typePrefix(self.CDF.N, [], self.CDF.n)
        modulus = 10
        rootsForParallelJobs = [rootPrefix.childrenListLastEntry_k_Mod_m(k, modulus) for k in range(modulus)]    
        #print rootsForParallelJobs[0]
        #print self.CDF.makeNewTypeAccountingForPrefix(rootsForParallelJobs[0][0])
        ListOfCDFLists = [self.CDF.accountForListOfTypes(rootList) for rootList in rootsForParallelJobs]
        ListOfCDFs = [CDF for sublist in ListOfCDFLists for CDF in sublist]
        self.CDF.mergeListInto(ListOfCDFs)
        self.failUnlessEqual(len(self.CDF.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF.assignCumulativeProbability(0.1), 0.470627961298)
        
    def _testNewParallel(self):
        self.CDF.setReferenceDistribution_p_eta(0.1)
        self.CDF.setN(30)
        self.CDF.setn(4)
        self.CDF.accountForAllTypesParallelized(10)
        self.failUnlessEqual(len(self.CDF.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF.assignCumulativeProbability(0.1), 0.470627961298)
    
    
       
    def _testParallel(self):
        import os
        pathFromHomeToPickleFile = '/Documents/sontag/timing_experiments/'
        pathToFile = os.getenv('HOME') + pathFromHomeToPickleFile
        #check if path exists and is writable
        if not os.access(pathToFile, os.W_OK):
            raise IOError("Path " + pathFromHomeToPickleFile + "is not writable.")
        eta = 0.01
        N_List = [10, 20]#  120, 130, 140, 150, 160, 170, 180, 190]
        for N in N_List:
            parameterList = [[N,eta,k] for k in range(10)]
            res = Parallel(n_jobs=-1, verbose=50)(delayed(CDF.returnCDFAccountingForTypesOfModuloClassk)(*listForAnInstance) 
                                           for listForAnInstance in parameterList)
            mergedCDF = CDF.CDF()
            for resultCDF in res:
                mergedCDF.mergeInto(resultCDF)
            
            fileName = timeStamped('_' + str(N) + '_mergedCDF.p')
            pickle.dump(mergedCDF.Dictionary, open(pathToFile + fileName, 'wb'))
            if N == 10: self.failUnlessEqual(len(mergedCDF.Dictionary),96)
            if N == 20: self.failUnlessEqual(len(mergedCDF.Dictionary),583)
            pprint(mergedCDF.Dictionary)
    
    def _testFindExactFromCDFFile(self):
        import os
        pathFromHomeToPickleFile = '/Documents/sontag/timing_experiments/mergedCDFs'
        pathToFile = os.getenv('HOME') + pathFromHomeToPickleFile
        fileName= pathToFile + 'N'  +  '_mergedCDF.p'
       
        #get the CDF stored in the file
        theDict = pickle.load(open( fileName, 'rb'))
        dictKeys = theDict.keys()
        dictKeys.sort()   #put in ascending order
        exactCDF = CDF.CDF()
        exactCDF.Dictionary = theDict
        exactCDF.AscendingDiscontinuityList = dictKeys #already sorted in ascending order
        
        
        #load the dataFrame with the estimated value
        approxDF = pandas.read_csv(pathToFile + '2013-02-25-03-48-47_MC_results_N200.csv')
        f = lambda x : exactCDF.assignCumulativeProbability(x, debug = True)
        approxDF['exactBeta'] = approxDF['gamma'].map(f)
        fileName = timeStamped('MC_exact_comparison_N200.csv')
        approxDF.to_csv(pathToFile + fileName)
        pprint(approxDF)
    
    def _testAccountForAllTypesRobbins(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(30)
        self.CDF.setn(4)
        self.CDF.accountForAllTypesRobbins()
        pprint(self.CDF.Dictionary)
    
if __name__ == "__main__":
    unittest.main()