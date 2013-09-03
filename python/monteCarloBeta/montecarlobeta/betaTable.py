'''
Created on Aug 19, 2012

@author: eliotpbrenner

See Section 4.2 of the Thesis: "Building the Table"
'''
import sys
sys.path.insert(0,'../src')

import pickle
import csv
import math
import numpy as np
import time
import probabilityDistributionPathFactory as pdpf
import probabilityDistributionFactory as pdf
import emissionProbabilityCalculator as epc
import mcIntegrationWithStoppingCriteria as iwsc
import informationTheory as it
import statsDistributionFactory as sdf
from scipy.interpolate import griddata
from scipy import stats
from joblib import Parallel, delayed
import insertionPoints as iP

class Mytimer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start


tolerance = 1e-5  #scale certain KL-divergences so that they are "feasible" by 1-tolerance


class betaTable(object):
    '''
    A "table" object for storing the results of beta results for a large number of N's and gammas done for a single eta
    '''


    def __init__(self, eta, dataFilename = None, k=2,l=2):
        '''
        Constructor: sets nothing apart from eta
                     Also loads dataFilename if provided.
        '''
        self.eta = eta      #set before unpickling
        self.NList = None   #set from pickled results list
        self.NFirstIndexHash = None
        self.normalizedKL_divergenceList = None  #used in construction from lists
        self.NToMinimumGammaDict = None #used in construction from lists
        self.NToMaxKL_DivergenceDict = None #used in construction from lists
        self.NToKL_divergenceList = None #computed from pickled results list
        self.NToGammaList = None  #set from pickled results list
        self.NToBetaList = None  #set from pickled results list
        self.NumberOfGammasGreaterThanEta = 10
        self.largestN = None
        self.nextToLargestN = None
        self.largestNFirstIndex = None
        self.nextToLargestNFirstIndex = None
        self.nextToLargestNLastIndex = None
        self.k = k
        self.l = l
        self.uniformMarginals = [1.0/k, 1.0/l]
        self.CentralProbability = 0.68  #probability that a Gaussian distribution lies within 1 st. deviation of mean
        self.maxIterations = 10000
        self.frequencyOfRecordingResult = 100
        self.frequencyOfApplicationOfStoppingCriterion = 100
        self.accuracy_percent = 10
        self.probability_accuracy_achieved = 0.95
        if dataFilename is not None:
            print dataFilename
            self.unSerialize(dataFilename)
        if self.k != 2 or self.l !=2:
            raise NotImplementedError("Beta table class implemented only for k=l=2 so far")


    def generate_N_List(self, intervalMaxList, stepSizeList): #tested
        '''
        :Parameters:
        intervalMaxList : list of ints length n
        stepSizeList list of ints of length n
        
        :Effect:
        Fills the NList data member with a list of Ns.    
        Takes steps of size stepSizeList[i] between
        intervalMaxList[i-1] (=0 for i=-1) and intervalMaxList[i].  
        
		Fills the NFirstIndexHash data member with a dictionary
		from N values to their first occurrence in the NList 

        :Example:
        See testBetaTable.testNList
        '''
        intervalMaxList.insert(0,0)  #tack on 0 to the beginning of the intervalMaxList
        self.NList = []
        for index in range(len(stepSizeList)):
            self.NList.extend(range(intervalMaxList[index], intervalMaxList[index+1], stepSizeList[index]))
            self.NList = [entry for entry in self.NList if entry > 0]  #eliminate any zeros because they are invalid anyway


    def generateDictFromNToGammaList(self):
        """
        :Effect:
        Convert NToKL_divergenceList to NToGammaList
        """
        if not self.NToKL_divergenceList:
            raise ValueError("Beta table has no NToKL_divergenceList.")
        self.NToGammaList = {}
        probabilityDistPathBasedAtUniform = pdpf.probabilityDistributionPathFactory([1.0/self.k, 1.0/self.l], self.k, self.l).construct()
        probabilityDistPathBasedAtUniform.markP_eta(self.eta)
        for N in self.NList:
            GammaList = [ probabilityDistPathBasedAtUniform.KL_divergence_at_t(
                probabilityDistPathBasedAtUniform.t_at_specifiedDivergenceFromMarkedDistInDirectionOfBase(KLDivergence))
                          for KLDivergence in self.NToKL_divergenceList[N]]
            self.NToGammaList[N] = np.array(GammaList)

    def generateDictFromNTo_KL_divergenceListAndGammaList(self):
        '''
        :Effect:
        Fills NToKL_divergenceList and NToGammaList data members based on the 
        normalizedKL_divergenceList with dictionaries having 
        as keys the integer elements of the NList and as values np.arrays
        '''
        self.generateDictFromNTo_KL_divergenceList()
        self.generateDictFromNToGammaList()
        return



    def generateDictFromNTo_KL_divergenceListAndGammaListIncludingGammaGreaterThanEta(self): #tested
        "Also gives numGammasGreaterThanEta Gammas"

        if self.normalizedKL_divergenceList is None or len(self.normalizedKL_divergenceList) == 0:
            raise Exception("Beta table has no normalized KL_divergnece list")
        if not self.NumberOfGammasGreaterThanEta:
            raise Exception("Beta table has no variable for number of gammas greater than eta")
        self.generateDictFromNTo_KL_divergenceList()
        self.generateDictFromNToGammaList()
        p_eta = pdf.probabilityDistributionFactory(self.k, self.l).get_p_eta(self.eta)

        uniformMarginals = [1.0/self.k,1.0/self.l]
        probabilityDistPathBasedAtUniformMarginals = pdpf.probabilityDistributionPathFactory(uniformMarginals, self.k, self.l).construct()
        t_max = probabilityDistPathBasedAtUniformMarginals.t_max
        distributionAt_t_max_OneUniformBasedPath = probabilityDistPathBasedAtUniformMarginals.distribution_at_t(t_max)
        KLDivergenceFromP_etaToDistributionAtTMaxOnPath = p_eta.KL_divergence_as_base(distributionAt_t_max_OneUniformBasedPath)

        probabilityDistPathBasedAtUniform = pdpf.probabilityDistributionPathFactory([1.0/self.k, 1.0/self.l], self.k, self.l).construct()
        probabilityDistPathBasedAtUniform.markP_eta(self.eta)

        numLgGam = int(self.NumberOfGammasGreaterThanEta)
        rawKLDivergenceListForGammaGreaterThanEta = KLDivergenceFromP_etaToDistributionAtTMaxOnPath* ( (1.0-tolerance)/numLgGam )*np.array(range(numLgGam+1) )


        for N in self.NList:
            self.NToKL_divergenceList[N].extend(rawKLDivergenceListForGammaGreaterThanEta)
            GammaListForGammaGreaterThanEta = [ probabilityDistPathBasedAtUniform.KL_divergence_at_t(
                probabilityDistPathBasedAtUniform.t_at_specifiedDivergenceFromMarkedDistAwayFromBase(KLDivergence)) for KLDivergence in rawKLDivergenceListForGammaGreaterThanEta]
            self.NToGammaList[N] = np.append(self.NToGammaList[N],np.array(GammaListForGammaGreaterThanEta))


    def generateDictFromNTo_KL_divergenceList(self):
        """
        :Effect:
        Convert normalizedKL_divergenceList to NToKL_divergenceList
        by multiplying by KLDistributionFromP_etaToUniform.
        """

        p_eta = pdf.probabilityDistributionFactory(self.k, self.l).get_p_eta(self.eta)
        uniformDistribution = pdf.probabilityDistributionFactory(self.k, self.l).get_p_eta(0.0)
        KLDistributionFromP_etaToUniform = p_eta.KL_divergence_as_base(uniformDistribution.distribution)
        rawKLDivergenceList = KLDistributionFromP_etaToUniform*np.array(self.normalizedKL_divergenceList)
        self.NToKL_divergenceList = {}
        for N in self.NList:
            self.NToKL_divergenceList[N] = [KLDivergence for KLDivergence in rawKLDivergenceList if KLDivergence < self.NToMaxKL_DivergenceDict[N]]


    def generateNormalizedKL_divergenceList(self, stepsize, levelRatio, numLevels):
        '''
        :Parameters:
        stepsize : float, typically < 1.0 (say 0.1)
        levelRatio: integer > 1, typically 2,3,...,10
        numLevels: integer >=1, typically around 4 or 5
        
        :Effect:
        Fill the normalizedKL_divergenceList datamember with
        list of numbers from 0 to 1, with density increasing, 
        i.e. tick size decreasing, as ticks move from 0 to 1, in general:
          *stepsize is the initial stepsize between 0 and 1/levelRatio (level 1)
          *levelRatio is the ratio between the inter-tick distance at level i and level i+1
          *numLevels is the number of times the inter-tick distance decreases by a factor of 1/levelRatio 
          
        :Example:
        See testBetaTable.testNoramalized_KLDivergenceList
        '''

        self.normalizedKL_divergenceList = np.array([])
        interval = levelRatio*stepsize
        print "interval %s"%(interval)
        for exponentIncrement in range(numLevels-1):
            self.normalizedKL_divergenceList = np.append(self.normalizedKL_divergenceList,levelRatio**(-numLevels+exponentIncrement)*
                                                                                          (np.arange(interval,levelRatio*interval,interval)))
        exponentIncrement = numLevels-1
        self.normalizedKL_divergenceList = np.append(self.normalizedKL_divergenceList,levelRatio**(-numLevels+exponentIncrement)*
                                                                                      (np.arange(interval,levelRatio,interval)))
        self.normalizedKL_divergenceList = 1.0 - self.normalizedKL_divergenceList
        self.normalizedKL_divergenceList = self.normalizedKL_divergenceList[::-1]


    def generate_N_toMinimumGammaDict(self):
        '''
        :Effect:
        Implements (145) in Section 4.2 of Chapter 3
        
        Generates a dictionary mapping each value of N to the minimum gamma
        and assigned to the NToMinimumGammaDict datamember
        
        :Example:
        testGenerate_N_ToMinimumGammaDict
        '''
        if not self.NList:
            raise ValueError("No NList for betaTable.")
        logGammas = (-1.0/10.0)*np.array(range(1,1000))
        Gammas = np.exp(logGammas)
        uniformMarginalsBasedPath = pdpf.probabilityDistributionPathFactory([1.0/self.k, 1.0/self.l],
                                                                            self.k,self.l).construct()
        Ns = [np.int(np.ceil(1.0/(uniformMarginalsBasedPath.lengthOfSegmentofKLDivergenceLessThanSpecified(gamma))))  for gamma in Gammas]
        NToMinimumGammaDict = {}
        for N in self.NList:
            index = 0
            while Ns[index] <= N:
                index += 1
            NToMinimumGammaDict[N] = Gammas[index]
        self.NToMinimumGammaDict = NToMinimumGammaDict


    def unSerialize(self, pickleFile):     #TODO: refactor to name "deSerialize" 
        #Load in results stored in serialized pickle format
        #auxiliary function
        def N_KLDiv_Pts(arrayOfResultsForThisEta,p_eta):
            N_KLDiv_Pts = [tuple([alist['N'], alist['gamma'], alist['beta']]) for alist in arrayOfResultsForThisEta]
            N_KLDiv_Pts = np.array( N_KLDiv_Pts, dtype=[('N', '<i8'), ('KL_Div','<f8'), ('beta', '<f8')])
            #import ipdb; ipdb.set_trace() 
            for row in N_KLDiv_Pts:
                row['KL_Div'] = p_eta.KLDivergenceOfP_gammaFromDist(row['KL_Div'])
            N_KLDiv_Pts.sort(order=['N','KL_Div'])
            return N_KLDiv_Pts

        #main body:
        if not self.eta:
            raise Exception("Eta must be defined prior to unserializing results of beta computation")
        probabilityDistFactoryUniformMarginals = pdf.probabilityDistributionFactory(self.k, self.l)
        p_eta = probabilityDistFactoryUniformMarginals.get_p_eta(self.eta)
        listOfResults= pickle.load(open(pickleFile,'rb'))
        print pickleFile
        listOfResultsAsTuples = [tuple(alist) for alist in listOfResults]
        arrayOfResults = np.array(listOfResultsAsTuples, dtype=[('time', '<f8'),('eta', '<f8'), ('N', '<i8'), ('gamma', '<f8'), ('beta', '<f8')])
        arrayOfResultsForThisEtaGammaLessThanEta = arrayOfResults[ arrayOfResults['eta'] == self.eta] #and arrayOfResults['gamma'] < self.eta]
        arrayOfResultsForThisEtaGammaLessThanEta = arrayOfResultsForThisEtaGammaLessThanEta[arrayOfResultsForThisEtaGammaLessThanEta['gamma'] < self.eta]
        arrayOfResultsForThisEtaGammaAtLeastEta =  arrayOfResults[ arrayOfResults['eta'] == self.eta]
        arrayOfResultsForThisEtaGammaAtLeastEta = arrayOfResultsForThisEtaGammaAtLeastEta[ arrayOfResultsForThisEtaGammaAtLeastEta ['gamma'] >= self.eta ]

        self.betasByAscendingNAscendingGamma = np.array(listOfResultsAsTuples, dtype=[('time', '<f8'),('eta', '<f8'), ('N', '<i8'), ('gamma', '<f8'), ('beta', '<f8')])
        self.betasByAscendingNAscendingGamma.sort(order = ['N', 'gamma'])

        self.N_KLDivPtsForInterpolationGammaLessThanEta =  N_KLDiv_Pts(arrayOfResultsForThisEtaGammaLessThanEta,p_eta)
        self.N_KLDivPtsForInterpolationGammaAtLeastEta = N_KLDiv_Pts(arrayOfResultsForThisEtaGammaAtLeastEta,p_eta)
        #import ipdb; ipdb.set_trace()     
        N_KLDivPtsForInterpolation_nd_GammaLessThanEta = self.N_KLDivPtsForInterpolationGammaLessThanEta[['N', 'KL_Div']].view(np.ndarray).reshape(len(self.N_KLDivPtsForInterpolationGammaLessThanEta), -1)
        self.points_GammaLessThanEta = np.array([list(arow[0].view(np.ndarray)) for arow in N_KLDivPtsForInterpolation_nd_GammaLessThanEta])
        N_KLDivPtsForInterpolation_nd_GammaAtLeastEta = self.N_KLDivPtsForInterpolationGammaAtLeastEta[['N', 'KL_Div']].view(np.ndarray).reshape(len(self.N_KLDivPtsForInterpolationGammaAtLeastEta), -1)
        self.points_GammaAtLeastEta = np.array([list(arow[0].view(np.ndarray)) for arow in N_KLDivPtsForInterpolation_nd_GammaAtLeastEta])
        self.valuesForInterpolationGammaLessthanEta = [np.log(row['beta']) for row in self.N_KLDivPtsForInterpolationGammaLessThanEta]

        self.valuesForInterpolationGammaAtLeastEta = [np.log(row['beta']) for row in self.N_KLDivPtsForInterpolationGammaAtLeastEta]

        self.largestN = self.betasByAscendingNAscendingGamma['N'][-1]
        #  if N is greater than the last (largest) N in the N_KLDivPtsForInterpolation['N'], compute value of function
        #  on this largest N and the given gamma, then second largest N and given gamma, and extrapolate from those two.
        self.nextToLargestN, self.nextToLargestNLastIndex = self.nextLargestNAndLastIndex(self.largestN)

        self.largestNFirstIndex = self.nextToLargestNLastIndex
        thirdLargestN,thirdLargestNLastIndex = self.nextLargestNAndLastIndex(self.nextToLargestN, self.nextToLargestNLastIndex)
        self.nextToLargestNFirstIndex = thirdLargestNLastIndex
        self.NList = sorted(set(self.betasByAscendingNAscendingGamma['N']))

        self.NFirstIndexHash = iP.firstIndexHash(set(self.NList),list(self.betasByAscendingNAscendingGamma['N']))

    def unSerializeDumpToCSV(self,pickleFileName, baseFileName):
    # :Input:
    # Path to the pickle file with info stored in python-friendly format
    # Basefilename, including path, of the csv's output by this method
    #
    # :Output:
    # --Four files:
    #      *baseFilenameN_KLDivPtsGammaLTEta.csv
    #      *baseFilenameN_KLDivPtsGammaGTEta.csv
    #      *baseFilenameLgBetaValsGammaLTEta.csv
    #      *baseFilenameLgBetaValsGammaGTEta.csv
    #  Each of format given in DumpToCSV

    #get the data from the pickle file into memory
        self.unSerialize(pickleFileName)

        #setup file names
        fileSuffixes = ['N_KLDivPtsGammaLTEta.csv', 'N_KLDivPtsGammaGTEta.csv', 'LgBetaValsGammaLTEta.csv', 'LgBetaValsGammaGTEta.csv' ]
        N_KL_LT_fname, N_KL_GT_fname, LgBeta_LT_fname, LgBeta_GT_fname = [baseFileName + suffix for suffix in fileSuffixes]
        #change betas to their logs
        #import ipdb; ipdb.set_trace()
        self.N_KLDivPtsForInterpolationGammaLessThanEta['beta'] = np.log(self.N_KLDivPtsForInterpolationGammaLessThanEta['beta'])
        self.N_KLDivPtsForInterpolationGammaAtLeastEta['beta'] = np.log(self.N_KLDivPtsForInterpolationGammaAtLeastEta['beta'])
        self.N_KLDivPtsForInterpolationGammaLessThanEta.dtype.names = 'N', 'KL_Div', 'LgBeta'
        self.N_KLDivPtsForInterpolationGammaAtLeastEta.dtype.names = 'N', 'KL_Div', 'LgBeta'
        DumpToCSV(N_KL_LT_fname, self.eta, self.NList, 'KL_Div', self.N_KLDivPtsForInterpolationGammaLessThanEta,)
        DumpToCSV(N_KL_GT_fname, self.eta, self.NList, 'KL_Div', self.N_KLDivPtsForInterpolationGammaAtLeastEta)
        DumpToCSV(LgBeta_LT_fname, self.eta, self.NList, 'LgBeta', self.N_KLDivPtsForInterpolationGammaLessThanEta)
        DumpToCSV(LgBeta_GT_fname, self.eta, self.NList, 'LgBeta', self.N_KLDivPtsForInterpolationGammaAtLeastEta)





    def convertMinimumGammaToMaxKL_Divergence(self, eta):
        """
        :Parameters:
        eta : float, used to construct p^\eta in the below
        
        :Explanation:
        Max KL divergence is KL(p^\gamma_minimum \| p^\eta)
        unless gamma_minimum is > eta, in which case it's zero
        
        self must have NToMinimumGammaDict
        """
        def MaxKL_Divergence(eta, p_eta, probabilityDistFactory, gamma_minimum):
            if gamma_minimum > eta:
                return 0
            else:
                return p_eta.KL_divergence_as_base(probabilityDistFactory.get_p_eta(gamma_minimum).distribution)

        if not self.NToMinimumGammaDict:
            raise ValueError("No minimum gamma dictionary")
        self.NToMaxKL_DivergenceDict = {}
        probabilityDistFactory = pdf.probabilityDistributionFactory(self.k,self.l)
        p_eta = probabilityDistFactory.get_p_eta(eta)
        for N in self.NToMinimumGammaDict.keys():
            self.NToMaxKL_DivergenceDict[N] = MaxKL_Divergence(
                eta,p_eta,probabilityDistFactory,self.NToMinimumGammaDict[N])


    def numberOfComputations(self):
        """
        :Returns:
        Number of computations of betas which have to be done
        """
        theNumber = 0
        for aList in self.NToKL_divergenceList.values():
            theNumber += len(aList)
        return theNumber

    def produceCommandLineScript(self):
        """
        :Output:
	  a list of command line python commands with the appropriate
          arguments useful for multicore processing, e.g.
				python betaCommandLine.py 0.01 100 0.001
        """
        if not self.NToGammaList:
            raise Exception("No dictionary from N to list of gammas in betaTable object")
        if not self.eta:
            raise Exception ("No eta value in betaTable object")
        commandsString = ""
        parameterList = self.produceParameterList()
        for eta, N, gamma in parameterList:
            commandsString += "python betaCommandLine.py %s %s %s\n"%(eta,N,gamma)
        return commandsString

    def produceParameterList(self):
        """
        :Output:
        list
        [...,[eta,N,gamma],...]
        """
        if not self.NToGammaList:
            raise ValueError("No dictionary from N to list of gammas in betaTable object")
        if not self.eta:
            raise ValueError("No eta value in betaTable object")
        parameterList = []
        for N in sorted(self.NToGammaList.keys()):
            for gamma in self.NToGammaList[N]:
                listForThisInstance = [self.eta, N, gamma]
                parameterList.append(listForThisInstance)
        return parameterList


    def ComputeBetasSerialize(self,pickleFileName, n_jobs = -1, verbose = 50):
        """
        :Input:
          pickleFileName: serialized version goes here
          n_jobs: number of jobs for Parallel processing: default (-1) is to use all available cores
          verbose: for reporting progress of Parallel.  50 by default.
        """
        parameterList = self.produceParameterList()
        res = ComputeBetasParallel(parameterList,self.k,self.l,self.maxIterations,self.frequencyOfRecordingResult,self.frequencyOfApplicationOfStoppingCriterion, self.accuracy_percent, self.probability_accuracy_achieved)
        pickle.dump(res,open(pickleFileName, 'wb'))

    def nextLargestNAndLastIndex(self,largerN, currIndexFromEnd = -1):
        nextLargestN  = largerN
        while nextLargestN ==  largerN and math.fabs(currIndexFromEnd) < len(self.betasByAscendingNAscendingGamma['N']):
            currIndexFromEnd -= 1
            nextLargestN =  self.betasByAscendingNAscendingGamma['N'][currIndexFromEnd]
        return nextLargestN, currIndexFromEnd


    def extrapolatedValue(self, N, N_KLDivPtsForInterpolation, largestN, gamma,currSmallest =0):
        valueCalculatedWithLargestN = self.interpolateAndReturnLogBetaWithoutChecking( largestN, gamma)
        nextToLargestN =   largestN  #TODO: optimize this by storing it once!!
        currIndexFromEnd = -1
        while nextToLargestN ==  largestN:
            currIndexFromEnd -= 1
            nextToLargestN =  N_KLDivPtsForInterpolation['N'][currIndexFromEnd]
        valueCalculatedWithNextToLargestN = self.interpolateAndReturnLogBetaWithoutChecking(nextToLargestN, gamma)
        slope = (valueCalculatedWithLargestN - valueCalculatedWithNextToLargestN)/(1.0*( largestN - nextToLargestN))
        if (gamma > self.eta and slope < 0) or (gamma < self.eta and slope > 0):  #correct for spurious values: the slope shouldn't be negative if gamma > eta or positive if gamma < 0
            slope = 0
            print 'not stopping early. Returning with latest N interpolation:',  min(0,valueCalculatedWithLargestN + slope*(N-largestN)), currSmallest, gamma, N
        result = valueCalculatedWithLargestN + slope*(N-largestN) #TODO: write wrapper function that mins this with 0
        return result

    def interpolateAndReturnLogBeta(self, N, gamma, currSmallest=0):
    #TODO: redesign API so that only this method and a few other similar ones are exposed.
        if not self.eta:
            raise ValueError("Eta must be defined prior to interpolating and returning log beta")

            #check against currSmallest
        #test case: N=150, gamma=0.01
        leftNVal = iP.find_ge(self.NList,N)[0]  #200
        if leftNVal is self.NList[-1]:
            rightNVal = leftNVal
        else:
            rightNVal = iP.find_gt(self.NList, leftNVal)[0] #300
        rightNIndex = self.NFirstIndexHash[rightNVal]  #20
        leftNIndex = self.NFirstIndexHash[leftNVal]   #40
        clampedGamma = max(gamma, self.betasByAscendingNAscendingGamma['gamma'][leftNIndex])
        clampedGamma = min(clampedGamma, self.betasByAscendingNAscendingGamma['gamma'][rightNIndex])
        candidateSmallerLogBeta = np.log(self.betasByAscendingNAscendingGamma['beta'][iP.find_le_withBounds(
            self.betasByAscendingNAscendingGamma['gamma'], clampedGamma, leftNIndex, rightNIndex)[1]])
        if candidateSmallerLogBeta > currSmallest:
            return currSmallest
        return min(currSmallest, self.interpolateAndReturnLogBetaWithoutChecking(N, gamma))  #TODO: write wrapper function that mins this with 0

    def N_KLDivPtsAndBetaValuesForInterpolation(self,gamma):
        #checks that necessary values are defined as data members
        #If so, returns the N_KLDivPts and Values for interpolation,
        #depending on whether gamma > eta or not.
        if not self.eta:
            raise Exception("Eta must be defined prior to interpolating and returning log beta")

        if gamma >= self.eta:
            N_KLDivPtsForInterpolation = self.N_KLDivPtsForInterpolationGammaAtLeastEta
            valuesForInterpolation = self.valuesForInterpolationGammaAtLeastEta
            points = self.points_GammaAtLeastEta
        else:  #gamma < self.eta
            N_KLDivPtsForInterpolation = self.N_KLDivPtsForInterpolationGammaLessThanEta
            valuesForInterpolation = self.valuesForInterpolationGammaLessthanEta
            points = self.points_GammaLessThanEta
        if not np.shape(N_KLDivPtsForInterpolation)[0]>0:
            raise ValueError("N_KLDivPtsForInterpolation must be defined prior to interpolating and returning log beta")
        if not np.shape(valuesForInterpolation)[0]>0:
            raise ValueError("valuesForInterpolation must be defined prior to interpolating and returning log beta")
        return N_KLDivPtsForInterpolation, valuesForInterpolation, points

    def KL_DivFromP_etaOfP_gamma(self, gamma,k,l, largestKL_DivergenceForThisN=10000):
        probabilityDistFactoryUniformMarginals = pdf.probabilityDistributionFactory(k, l)
        p_eta = probabilityDistFactoryUniformMarginals.get_p_eta(self.eta)
        KL_DivergenceFromP_etaOfP_gamma = p_eta.KLDivergenceOfP_gammaFromDist(gamma)
        #largestKL_DivergenceForThisN =  10000 #TODO: eliminate
        return min(KL_DivergenceFromP_etaOfP_gamma, largestKL_DivergenceForThisN)

    def interpolateAndReturnLogBetaWithoutChecking(self, N, gamma):

        def thresholdWithGamma_0(points,N,KL_DivergenceFromP_etaOfP_gamma):
            LastIndexHoldingKLDivergencesForN = np.searchsorted(points[:,0], N+0.5)-1
            if KL_DivergenceFromP_etaOfP_gamma > points[LastIndexHoldingKLDivergencesForN][1]:
                KL_DivergenceFromP_etaOfP_gamma = points[LastIndexHoldingKLDivergencesForN][1]
            return KL_DivergenceFromP_etaOfP_gamma

            #Pick out the right part of the grid based on whether gamma > eta or gamma < eta
        N_KLDivPtsForInterpolation, valuesForInterpolation, points = self.N_KLDivPtsAndBetaValuesForInterpolation(gamma)

        #######################
        # EXTRAPOLATION CASES #
        ########################
        #  if given N is LESS THAN THE FIRST (smallest) N in the 
        #N_KLDivPtsForInterpolation['N'] value, return 0.
        if N < N_KLDivPtsForInterpolation['N'][0]:
            return 0
            #  if N is GREATER THAN THE LAST (largest) N in the N_KLDivPtsForInterpolation['N'],
            #compute value of function
            #  on this largest N and the given gamma, then second largest N and given gamma,
            #and extrapolate from those two.
        if N > self.largestN:
            return self.extrapolatedValue(N, N_KLDivPtsForInterpolation, self.largestN, gamma)

        ######################
        # INTERPOLATION CASE #
        ######################
        KL_DivergenceFromP_etaOfP_gamma = self.KL_DivFromP_etaOfP_gamma(gamma, self.k, self.l)
        #Threshhold with gamma_zero
        KL_DivergenceFromP_etaOfP_gamma = thresholdWithGamma_0(points, N,KL_DivergenceFromP_etaOfP_gamma)
        #points: size 27 by 2,
        """
        [[  1.50000000e+02   1.00662110e-03]
 [  1.50000000e+02   2.01324231e-03]
 [  1.50000000e+02   3.01986347e-03]
 x[  1.50000000e+02   4.02648443e-03]
 x[  1.50000000e+02   5.03310554e-03]
 [  1.50000000e+02   6.03972653e-03]
 [  1.50000000e+02   7.04634791e-03]
 [  1.50000000e+02   8.05296826e-03]
 [  1.60000000e+02   1.00662110e-03]
 [  1.60000000e+02   2.01324231e-03]
 [  1.60000000e+02   3.01986347e-03]
 x[  1.60000000e+02   4.02648443e-03]
 x[  1.60000000e+02   5.03310554e-03]
 [  1.60000000e+02   6.03972653e-03]
 [  1.60000000e+02   7.04634791e-03]
 [  1.60000000e+02   8.05296826e-03]]

         values for interpolation [-1.2144560173716061, -1.4511231883100346, -1.7313625289331374, x-1.9615352001829527, x-2.3478184438998579, -2.6377119764078074, -2.9469896408941301, -3.4276382731799142, -1.1472162478207573, -1.5210917187085564, -1.7544707257217145, x-2.0798175100315799, x-2.3384896837338087, -2.6977560657153621, -3.0536992279369453, -3.5608678191779584]
         N_KLDivPtsForInterpolation: N, KL, beta
        [(150, 0.0010066210995419156, 0.29687146329623676)
         (150, 0.002013242312216161, 0.23430696939425433)
         (150, 0.0030198634712750097, 0.17704301932011676)
         x(150, 0.004026484433649273, 0.1406423409529187)
         x(150, 0.005033105540203621, 0.09557744249810314)
         (150, 0.006039726534445542, 0.07152473275721835)
         (150, 0.007046347905582378, 0.05249750465505174)
         (150, 0.008052968262266202, 0.032463520273886486)
         (160, 0.0010066210995419156, 0.31751943566794116)
         (160, 0.002013242312216161, 0.2184732453818809)
         (160, 0.0030198634712750097, 0.17299878187022594)
         x(160, 0.004026484433649273, 0.12495301278943756)
         x(160, 0.005033105540203621, 0.09647323334905836)
         (160, 0.006039726534445542, 0.06735648681850963)
         (160, 0.007046347905582378, 0.04718405657234157)
         (160, 0.008052968262266202, 0.02841415566238527)]

         KL_DivergenceFromP_etaOfP_gamma 0.00471190208643
         N 155
        """
        return griddata(points, valuesForInterpolation, (N, KL_DivergenceFromP_etaOfP_gamma), method='linear') #min(0, )
    #End of Class Definition


    def ComputeBetasDumpToCSV(self,CSVFileName, n_jobs = -1, verbose = 50):
        """
        TODO: test
        :Input:
          CSVFileName: dump result into CSV
          n_jobs: number of jobs for Parallel processing: default (-1) is to use all available cores
          verbose: for reporting progress of Parallel.  50 by default.
        """
        parameterList = self.produceParameterList()
        res = ComputeBetasParallel(parameterList,self.k,self.l,self.maxIterations,self.frequencyOfRecordingResult,self.frequencyOfApplicationOfStoppingCriterion, self.accuracy_percent, self.probability_accuracy_achieved)
        with open(CSVFileName, 'wb') as f:
            csvWriter = csv.writer(f, delimiter=',')
            for row in res:
                #row=[runningTime, eta, N, gamma, beta]
                runningTime, eta, N, gamma, beta = row
                csvWriter.writerow([N,gamma,beta])

def functionFromMarginalsTo_t_Distribution(eta,N,gamma,k,l):
    statsDistFactory = sdf.statsDistributionFactory(gamma)
    statsDistFactory.set_eta(eta)
    statsDistFactory.set_N(N)
    uniformMarginals = [1.0/k,1.0/l]
    statsDistFactory.CalculateAndSetScaleRatioRobustToGammaExceedingEta(uniformMarginals)

    return statsDistFactory.GaussianCenteredAtMinOftGammaPlusOrEtafromMarginalsScaledFromScaleRatio

def functionProbabilityOfEmissionByP_eta_Robbins(eta,N,gamma,k,l):
    probabilityCalculatorAssociatedToEta = epc.emissionProbabilityCalculator(eta, k, l, N)
    probabilityCalculatorAssociatedToEta.setGamma(gamma)
    return probabilityCalculatorAssociatedToEta.RobbinsEstimateOfEmissionProbabilityTimesCharFunctionOfTauMinusGamma

def ComputeBeta(eta,N,gamma,k,l,maxIterations, frequencyOfRecordingResult,frequencyOfApplicationOfStoppingCriterion, accuracy_percent, probability_accuracy_achieved):
    the_functionFromMarginalsTo_t_Distribution = functionFromMarginalsTo_t_Distribution(eta,N,gamma,k,l)
    the_functionProbabilityOfEmissionByP_eta_Robbins = functionProbabilityOfEmissionByP_eta_Robbins(eta,N,gamma,k,l)
    Gaussianscale = it.ChernoffRadius(0.68, N) #TODO: eliminate this magic number, central probability
    normalDistLocUniformScaleDepOnChernoff = stats.norm(loc=0.5,scale=Gaussianscale)


    with Mytimer() as t:
        iteration_values, rho_hat_values = iwsc.IntegrateWithStoppingCriterion(k,l,
                                                                               the_functionProbabilityOfEmissionByP_eta_Robbins,
                                                                               normalDistLocUniformScaleDepOnChernoff,
                                                                               the_functionFromMarginalsTo_t_Distribution,
                                                                               maxIterations, frequencyOfRecordingResult,
                                                                               frequencyOfApplicationOfStoppingCriterion,
                                                                               accuracy_percent, probability_accuracy_achieved)
    runningTime = t.interval
    beta = rho_hat_values[-1]*N**3  #TODO: fix magic number (3=k*l-1)
    return [runningTime, eta, N, gamma, beta]



def ComputeBetasParallel(parameterList,k,l,maxIterations, frequencyOfRecordingResult,frequencyOfApplicationOfStoppingCriterion, accuracy_percent, probability_accuracy_achieved):
    for alist in parameterList:
        alist.extend([k,l,maxIterations, frequencyOfRecordingResult,frequencyOfApplicationOfStoppingCriterion, accuracy_percent, probability_accuracy_achieved])
    res = Parallel(n_jobs=-1, verbose=50)(delayed(ComputeBeta)(*alist) for alist in parameterList)
    return res

def DumpToCSV(filename, eta, NList, fieldName, recArray):
    # :Output Format:
    # 0.01 #eta
    # 100,200,300 #NList
    # 100,[comma-separated values for N=100] #values are values in the recArray[fieldName] for records with recArray['N'] == 100
    # ...
    # 300,[comma-separated values for N=300]
    with open(filename,'wb') as f:
        csvWriter = csv.writer(f)
        csvWriter.writerow([eta])
        csvWriter.writerow(NList)
        for N in NList:
            csvWriter.writerow([N] + list(recArray[recArray['N'] == N][fieldName]))
