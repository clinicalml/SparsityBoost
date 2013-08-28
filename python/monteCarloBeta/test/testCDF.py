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
import typePrefix as tp
import probabilityDistributionPathFactory as pdpf 

def timeStamped(fname, fmt='%Y-%m-%d-%H-%M-%S_{fname}'):
    return datetime.datetime.now().strftime(fmt).format(fname=fname)



class Test(unittest.TestCase):


    def setUp(self):
        self.CDF = CDF.CDF()
        self.CDF2 = CDF.CDF()

    def tearDown(self):
        pass


    @unittest.skip("showing class skipping")
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


    @unittest.skip("showing class skipping")       
    def testAccountForAllTypesWithTwoElementPrefix(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(5)
        prefix = tp.typePrefix(5,data=[2,1],n=4)
        self.CDF.accountForTypesWithPrefix(prefix)
        self.failUnlessAlmostEqual(self.CDF.Dictionary, {0.013844293808390619: 0.054900966738391482,
                                                         0.11849392256130019: 0.065587032834470232,
                                                         0.2911031660323688: 0.13610213450308134} )

    @unittest.skip("showing class skipping")
    def testAccountForAllTypes(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(30)
        self.CDF.setn(4)
        self.CDF.accountForAllTypes()
        self.failUnlessEqual(len(self.CDF.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF.assignCumulativeProbability(0.1), 0.470627961298)
        
        self.CDF2.referenceDistribution = p_eta
        self.CDF2.setN(30)
        self.CDF2.setn(4)
        self.CDF2.accountForAllTypes()
        self.failUnlessEqual(len(self.CDF2.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF2.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF2.Dictionary[self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(self.CDF2.assignCumulativeProbability(0.1), 0.470627961298)
 
    @unittest.skip("showing class skipping")       
    def testMerge(self):
        self.CDF.accountForEvent(1, 0.3)
        self.CDF.accountForEvent(2, 0.05)
        self.CDF2.accountForEvent(1.5, 0.05)
        self.CDF2.accountForEvent(2.5, 0.6)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList),2)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[
          self.CDF.AscendingDiscontinuityList[-1]], 0.35)
        self.failUnlessAlmostEqual(self.CDF2.Dictionary[
          self.CDF2.AscendingDiscontinuityList[-1]], 0.65)        
        self.CDF.mergeInto(self.CDF2)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList),4)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[
          self.CDF.AscendingDiscontinuityList[-1]], 1.0)

    @unittest.skip("showing class skipping")    
    def testNewParallel(self):
        self.CDF.setReferenceDistribution_p_eta(0.1)
        self.CDF.setN(30)
        self.CDF.setn(4)
        self.CDF.accountForAllTypesParallelized(10)
        self.failUnlessEqual(len(self.CDF.Dictionary), 2009)
        self.failUnlessEqual(len(self.CDF.AscendingDiscontinuityList), 2009)
        self.failUnlessAlmostEqual(self.CDF.Dictionary[
          self.CDF.AscendingDiscontinuityList[-1]], 1.0)
        self.failUnlessAlmostEqual(
          self.CDF.assignCumulativeProbability(0.1), 0.470627961298)
    
    @unittest.skip("showing class skipping")
    def testAccountForAllTypesRobbins(self):
        p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(0.1)
        self.CDF.referenceDistribution = p_eta
        self.CDF.setN(30)
        self.CDF.setn(4)
        self.CDF.accountForAllTypesRobbins()
        self.failUnlessAlmostEqual(
          self.CDF.assignCumulativeProbability(0.1), 0.49556072704210913)

    def testParallelNonUniformMarginals(self):
        import itertools 
        import numpy as np
        eta=0.01  #for reference distribution
        displacements = np.arange(0.05,0.55,0.05)
        for displacedMarginals in itertools.product(displacements,displacements):
            probabilityDistributionPath = pdpf.probabilityDistributionPathFactory(
              displacedMarginals,2,2).construct()
            p_eta_sub1 = \
              probabilityDistributionPath.\
              distribution_at_speicified_divergence_from_base_pos_t_as_distribution(eta)
            CDF1 = CDF.CDF()
            CDF1.referenceDistribution = p_eta_sub1
            CDF1.setN(30)
            CDF1.setn(4)
            CDF1.accountForAllTypes()
            print "*************************************"
            print displacedMarginals
            print len(CDF1.Dictionary)
            print CDF1.assignCumulativeProbability(eta)
            
             
    
    @unittest.skip("showing class skipping")    
    def testFindExactFromCDFFile(self):
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
    

    
if __name__ == "__main__":
    unittest.main()