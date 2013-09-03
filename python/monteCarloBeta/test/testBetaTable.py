#!/system/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python

'''
Created on Aug 19, 2012

@author: eliotpbrenner
'''

import sys

sys.path.insert(0, '../montcarlobeta')

import unittest
import numpy as np
import betaTable as bt
import probabilityDistributionPathFactory as pdpf

pathToPickleFiles = '/Users/eliotpbrenner/Documents/sontag/betaTables'


class testBetaTable(unittest.TestCase):
    @unittest.skip("demonstrating skipping")
    def testNList(self):
        eta = 0.01
        newBetaTable = bt.betaTable(eta)
        NMaxList = [100, 200, 500, 1000, 10000, 110000]
        stepSizeList = [5, 10, 50, 100, 1000, 10000]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        self.assertEquals(len(newBetaTable.NList), 59)
        self.assertEquals(newBetaTable.NList[:15], [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])
        self.assertEquals(newBetaTable.NList[16:31],
                          [85, 90, 95, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250])
        self.assertEquals(newBetaTable.NList[-15:],
                          [5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
                           100000])

    @unittest.skip("demonstrating skipping")
    def testNormalized_KLDivergenceList(self):
        eta = 0.01
        newBetaTable = bt.betaTable(eta)
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 10, 2)
        np.testing.assert_almost_equal(newBetaTable.normalizedKL_divergenceList,
                                       np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                                                 0.8, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96,
                                                 0.97, 0.98, 0.99]))
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 2, 5)
        np.testing.assert_almost_equal(newBetaTable.normalizedKL_divergenceList,
                                       [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.9875, 0.99375])

    @unittest.skip("demonstrating skipping")
    def testGenerate_N_ToMinimumGammaDict(self):
        eta = 0.01
        k, l = 2, 2
        newBetaTable = bt.betaTable(eta)
        NMaxList = [100, 200, 500, 1000, 10000, 110000]
        stepSizeList = [5, 10, 50, 100, 1000, 10000]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        newBetaTable.generate_N_toMinimumGammaDict()
        uniformMarginalsBasedPath = pdpf.probabilityDistributionPathFactory([1.0 / k, 1.0 / l],
                                                                            k, l).construct()

        for index in [1, 39, 50]:
            N = newBetaTable.NList[index]
            Nnext = newBetaTable.NList[index + 1]
            gammaN = newBetaTable.NToMinimumGammaDict[N]
            scriptL_gammaN = uniformMarginalsBasedPath.lengthOfSegmentofKLDivergenceLessThanSpecified(gammaN)
            recoveredN = np.int(np.ceil(1.0 / (scriptL_gammaN)))
            print N, recoveredN, Nnext
            self.failUnless(recoveredN > N)
            self.failUnless(recoveredN < Nnext)

    @unittest.skip("demonstrating skipping")
    def testGenerate_N_ToMaxKLDiv(self):
        eta = 0.01
        newBetaTable = bt.betaTable(eta)
        NMaxList = [100, 200, 500, 1000, 10000, 110000]
        stepSizeList = [5, 10, 50, 100, 1000, 10000]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        newBetaTable.generate_N_toMinimumGammaDict()
        newBetaTable.convertMinimumGammaToMaxKL_Divergence(eta)
        self.failUnlessAlmostEqual(newBetaTable.NToMaxKL_DivergenceDict[5], 0.0)
        self.failUnlessAlmostEqual(newBetaTable.NToMaxKL_DivergenceDict[200], 0.00875832083579)
        self.failUnlessAlmostEqual(newBetaTable.NToMaxKL_DivergenceDict[500], 0.00952414885724)
        self.failUnlessAlmostEqual(newBetaTable.NToMaxKL_DivergenceDict[9000], 0.010036955125)


    @unittest.skip("demonstrating skipping")
    def testgenerateDictFromNToKLDiv(self):
        eta = 0.1
        newBetaTable = bt.betaTable(eta)
        newBetaTable.unPickleIt(pathToPickleFiles)
        NMaxList = [100, 200, 500, 1000, 10000, 110000]
        stepSizeList = [5, 10, 50, 100, 1000, 10000]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        newBetaTable.generate_N_toMinimumGammaDict()
        newBetaTable.convertMinimumGammaToMaxKL_Divergence(eta)
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 2, 4)
        newBetaTable.generateDictFromNTo_KL_divergenceListAndGammaList()
        np.testing.assert_almost_equal(newBetaTable.NToKL_divergenceList[5], [])
        np.testing.assert_almost_equal(newBetaTable.NToKL_divergenceList[200], np.array(
            [0.01073601, 0.02147203, 0.03220804, 0.04294405, 0.05368006, 0.06441608, 0.07515209, 0.0858881, 0.09662412,
             0.10199212]))
        np.testing.assert_almost_equal(newBetaTable.NToKL_divergenceList[900], np.array(
            [0.01073601, 0.02147203, 0.03220804, 0.04294405, 0.05368006, 0.06441608, 0.07515209, 0.0858881, 0.09662412,
             0.10199212, 0.10467613, 0.10601813]))
        np.testing.assert_almost_equal(newBetaTable.NToKL_divergenceList[9000], np.array(
            [0.01073601, 0.02147203, 0.03220804, 0.04294405, 0.05368006, 0.06441608, 0.07515209, 0.0858881, 0.09662412,
             0.10199212, 0.10467613, 0.10601813]))
        np.testing.assert_almost_equal(newBetaTable.NToGammaList[5], np.array([]))
        np.testing.assert_almost_equal(newBetaTable.NToGammaList[200], np.array(
            [4.72911163e-02, 3.10520743e-02, 2.08606503e-02, 1.38177009e-02, 8.79790876e-03, 5.22301296e-03,
             2.74885554e-03, 1.15079180e-03, 2.72508398e-04, 6.64500683e-05]))
        np.testing.assert_almost_equal(newBetaTable.NToGammaList[900],
                                       np.array([4.72911163e-02, 3.10520743e-02, 2.08606503e-02,
                                                 1.38177009e-02, 8.79790876e-03, 5.22301296e-03,
                                                 2.74885554e-03, 1.15079180e-03, 2.72508398e-04,
                                                 6.64500683e-05, 1.64295047e-05, 4.09326296e-06]))
        np.testing.assert_almost_equal(newBetaTable.NToGammaList[9000],
                                       np.array([4.72911163e-02, 3.10520743e-02, 2.08606503e-02,
                                                 1.38177009e-02, 8.79790876e-03, 5.22301296e-03,
                                                 2.74885554e-03, 1.15079180e-03, 2.72508398e-04,
                                                 6.64500683e-05, 1.64295047e-05, 4.09326296e-06]))
        #TODO: replace with assertions
        print len(newBetaTable.NToKL_divergenceList[200])
        print len(newBetaTable.NToGammaList[200])

    @unittest.skip("demonstrating skipping")
    def testGenerateDictFromNTo_KL_divergenceListAndGammaListIncludingGammaGreaterThanEta(self):
        eta = 0.1
        newBetaTable = bt.betaTable(eta)
        NMaxList = [100, 200, 500, 1000, 10000, 110000]
        stepSizeList = [5, 10, 50, 100, 1000, 10000]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        newBetaTable.generate_N_toMinimumGammaDict()
        newBetaTable.convertMinimumGammaToMaxKL_Divergence(eta)
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 2, 4)

        newBetaTable.generateDictFromNTo_KL_divergenceListAndGammaListIncludingGammaGreaterThanEta()
        np.testing.assert_almost_equal(newBetaTable.NToKL_divergenceList[200][-10:],
                                       np.array([0.03287861, 0.06575721, 0.09863582, 0.13151442, 0.16439303,
                                                 0.19727164, 0.23015024, 0.26302885, 0.29590745, 0.32878606]))
        np.testing.assert_almost_equal(newBetaTable.NToGammaList[200][-10:],
                                       np.array([0.23624684, 0.30816165, 0.36908873, 0.42411776, 0.47520695,
                                                 0.52333623, 0.56905008, 0.61264038, 0.65417864, 0.69314369]))

    @unittest.skip("demonstrating skipping")
    def testProduceParameterListAndNumComputations(self):
        """


        """
        eta = 0.1
        newBetaTable = bt.betaTable(eta)
        NMaxList = [100, 200, 500, 1000, 10000, 110000]
        stepSizeList = [5, 10, 50, 100, 1000, 10000]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        newBetaTable.generate_N_toMinimumGammaDict()
        newBetaTable.convertMinimumGammaToMaxKL_Divergence(eta)
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 2, 4)
        newBetaTable.generateDictFromNTo_KL_divergenceListAndGammaListIncludingGammaGreaterThanEta()
        parameterList = newBetaTable.produceParameterList()
        print len(parameterList)
        np.testing.assert_almost_equal(parameterList[::50],
                                       [[0.1, 5, 0.10000001584978346], [0.1, 20, 0.36908872832865341],
                                        [0.1, 35, 0.10000001584978346], [0.1, 50, 0.020860650332333031],
                                        [0.1, 60, 0.52333622567850091], [0.1, 75, 0.0011507917956953204],
                                        [0.1, 90, 0.047291116251685306], [0.1, 100, 0.23624684232748688],
                                        [0.1, 130, 0.047291116251685306], [0.1, 150, 0.23624684232748688],
                                        [0.1, 180, 0.047291116251685306], [0.1, 200, 0.00027250839811708907],
                                        [0.1, 300, 0.52333622567850091], [0.1, 450, 0.031052074303808991],
                                        [0.1, 600, 0.0011507917956953204], [0.1, 800, 0.10000001584978346],
                                        [0.1, 1000, 0.42411775903703708], [0.1, 3000, 0.61264038021325939],
                                        [0.1, 6000, 0.031052074303808991], [0.1, 8000, 0.0052230129584135421],
                                        [0.1, 10000, 6.6450068260958534e-05], [0.1, 30000, 0.23624684232748688],
                                        [0.1, 50000, 0.47520695397905599], [0.1, 70000, 0.65417863953094035],
                                        [0.1, 100000, 0.020860650332333031]])
        self.assertEquals(len(parameterList), 1221)
        self.assertEquals(len(parameterList), newBetaTable.numberOfComputations())
        commandsString = newBetaTable.produceCommandLineScript()
        self.failUnlessEqual(commandsString[0:45], "python betaCommandLine.py 0.1 5 0.10000001585")
        self.failUnlessEqual(commandsString[-53:], "\npython betaCommandLine.py 0.1 100000 0.693143691468\n")


    @unittest.skip("demonstrating skipping")
    def testComputeBeta(self):
        eta = 0.1
        nbt = bt.betaTable(eta)
        computationResult = bt.ComputeBeta(0.1, 1000, 0.01, 2, 2, nbt.maxIterations, nbt.frequencyOfRecordingResult,
                                           nbt.frequencyOfApplicationOfStoppingCriterion, nbt.accuracy_percent,
                                           nbt.probability_accuracy_achieved)[1:]
        print "computationResult :", computationResult
        np.testing.assert_almost_equal(computationResult, (0.1, 1000, 0.01, 2.8744098845630804e-24))

    @unittest.skip("demonstrating skipping")
    def testComputeBetasSerialize(self):
        eta = 0.1
        newBetaTable = bt.betaTable(eta)
        NMaxList = [300]
        stepSizeList = [100]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        newBetaTable.generate_N_toMinimumGammaDict()
        newBetaTable.convertMinimumGammaToMaxKL_Divergence(eta)
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 1, 1)
        newBetaTable.generateDictFromNTo_KL_divergenceListAndGammaListIncludingGammaGreaterThanEta()
        newBetaTable.ComputeBetasSerialize('testBetas.p')
        recoveredBetaTb = bt.betaTable(eta)
        recoveredBetaTb.unSerialize('testBetas.p')
        vect1 = recoveredBetaTb.betasByAscendingNAscendingGamma[-20]
        np.testing.assert_almost_equal(np.array((vect1[1], vect1[2], vect1[3])), np.array(
            (0.1, 200L, 0.0002725083981170891))) #, 2.3840021343733423e-10))47.49394100000001
        vect2 = recoveredBetaTb.betasByAscendingNAscendingGamma[-21]
        np.testing.assert_almost_equal(np.array((vect2[1], vect2[2], vect2[3])), np.array(
            ( 0.1, 100L, 0.693143691467917)))  #, 1.0444680167646372))12.575713000000007
        self.assertEquals(recoveredBetaTb.nextToLargestNLastIndex, -21)
        self.assertEquals(recoveredBetaTb.largestNFirstIndex, -21)
        self.assertEquals(recoveredBetaTb.nextToLargestNFirstIndex, -40)
        self.failUnlessEqual(recoveredBetaTb.NFirstIndexHash, {200: 20, 100: 0})

    @unittest.skip("demonstrating skipping")
    def testInterpolateAndReturnBeta(self):
        etaList = [0.001] #, 0.04, 0.02, 0.01, 0.005, 0.001]
        prefixList = [
            'etaOneThousandth'] #, 'etaFourHundredths', 'etaTwoHundredths', 'etaOneHundredth', 'etaFiveThousandths']
        for eta, prefix in zip(etaList, prefixList):
            recoveredBetaTb = bt.betaTable(eta)
            print prefix + 'AdaptiveScalingLargeIncludingLargeGamma.p'
            recoveredBetaTb.unSerialize(prefix + 'AdaptiveScalingLargeIncludingLargeGamma.p')
            recoveredBetaTb.unSerializeDumpToCSV(prefix + 'AdaptiveScalingLargeIncludingLargeGamma.p', prefix)
            #res1 = recoveredBetaTb.interpolateAndReturnLogBetaWithoutChecking(150,0.01)
            #res2 = recoveredBetaTb.interpolateAndReturnLogBetaWithoutChecking(50,0.01)
            #res3 = recoveredBetaTb.interpolateAndReturnLogBetaWithoutChecking(250,0.01)
            #res4 = recoveredBetaTb.interpolateAndReturnLogBetaWithoutChecking(150,0.12)
            #res5 = recoveredBetaTb.interpolateAndReturnLogBetaWithoutChecking(150,0.0001)
            #res6 = recoveredBetaTb.interpolateAndReturnLogBetaWithoutChecking(150,0.000001)
            #res7 = recoveredBetaTb.interpolateAndReturnLogBeta(150,0.01)
            #np.testing.assert_approx_equal(res1,-9.65588259336, significant=1)
            #np.testing.assert_approx_equal(res2,0, significant=1)
            #np.testing.assert_approx_equal(res3,-14.9776264458,  significant=1)
            #np.testing.assert_approx_equal(res4,-0.636875219165, significant=1)
            #np.testing.assert_approx_equal(res5,-16.8125513352, significant=2)

    def testInterpolateAndReturnBetaWithChecking(self):
        eta = 0.01 #0.1
        newBetaTable = bt.betaTable(eta)
        NMaxList = [150, 170]
        stepSizeList = [150, 10]
        newBetaTable.generate_N_List(NMaxList, stepSizeList)
        print newBetaTable.NList
        newBetaTable.generate_N_toMinimumGammaDict()
        newBetaTable.convertMinimumGammaToMaxKL_Divergence(eta)
        newBetaTable.generateNormalizedKL_divergenceList(0.1, 1, 1)
        newBetaTable.generateDictFromNTo_KL_divergenceListAndGammaListIncludingGammaGreaterThanEta()
        #newBetaTable.ComputeBetasSerialize('testBetas.p')
        recoveredBetaTb = bt.betaTable(eta)
        recoveredBetaTb.unSerialize('testBetas.p')
        res1 = recoveredBetaTb.interpolateAndReturnLogBeta(155, 0.001, 0)
        res2 = recoveredBetaTb.interpolateAndReturnLogBeta(50, 0.001, -1)
        res3 = recoveredBetaTb.interpolateAndReturnLogBeta(250, 0.001, -10)
        res4 = recoveredBetaTb.interpolateAndReturnLogBeta(150, 0.12, 0)
        res5 = recoveredBetaTb.interpolateAndReturnLogBeta(150, 0.0001, -20)
        res6 = recoveredBetaTb.interpolateAndReturnLogBeta(150, 0.000001, 0)
        res6another = recoveredBetaTb.interpolateAndReturnLogBeta(150, 0.000000001, 0)
        res7 = recoveredBetaTb.interpolateAndReturnLogBeta(150, 0.01, -11)
        np.testing.assert_approx_equal(res1, -20, significant=1)
        np.testing.assert_approx_equal(res2, -1, significant=1)
        np.testing.assert_approx_equal(res3, -14.9776264458, significant=1)
        np.testing.assert_approx_equal(res4, -0.636875219165, significant=1)
        np.testing.assert_approx_equal(res5, -20, significant=2)
        np.testing.assert_approx_equal(res6, -17.0776918599, significant=2)
        np.testing.assert_approx_equal(res6another, -17.0776918599, significant=2)
        np.testing.assert_approx_equal(res7, -11, significant=2)


if __name__ == "__main__":
    unittest.main()
