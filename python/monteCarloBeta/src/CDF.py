'''
Created on Feb 17, 2013

@author: eliotpbrenner

Holds a CDF (cumulative probability distribution), and has methods for
augmenting the CDF stored so far to take account of the calculated
probability of a new event.
Data members (key):
 --AscendingDiscontinuityList: the locations of discontinuities of the CDF.
    In current implementation, only grows as events are accounted for.
 --Dictionary: the cumulative probability at each discontinuity point.

Data members (auxiliary):
 --Reference Distribution: when E is the event of emission of a sequence belongii
ng to type class T,
     the the probability of emission of T from p, where p is the reference distrr
ibution.
 --N: fixed size of the type classes T to be considered.
 --n: fixed length of the type classes T to be considered.
 --probDistributionFactory: builds reference distributions.

  The currently implemented methods for generating CDF's all relate
  to the Mutual Information (MI/tau) random variable of two binary (Bernoulli) rr
andom variables.
'''


import bisect
import csv
import probabilityDistributionFactory as pdf
import numpy as np
import typePrefix as tp
import pandas
from joblib import Parallel, delayed

def accountForListOfTypePrefixes(listOfTypePrefixes, N, n, referenceDistribution):
        """
        Return a single CDF accounting for the entire list of Type Prefixes
        """
        #Make the new CDF
        newCDF = CDF()
        newCDF.setN(N)
        newCDF.setn(n)
        newCDF.setReferenceDistributionFromDistribution(referenceDistribution)
        for aTypePrefix in listOfTypePrefixes:
            newCDF.accountForTypesWithPrefix(aTypePrefix)
        return newCDF

def accountForListOfTypePrefixesRobbins(listOfTypePrefixes, N, n, referenceDistribution):
        """
        Return a single CDF accounting for the entire list of Type Prefixes
        """
        #Make the new CDF
        newCDF = CDF()
        newCDF.setN(N)
        newCDF.setn(n)
        newCDF.setReferenceDistributionFromDistribution(referenceDistribution)
        for aTypePrefix in listOfTypePrefixes:
            newCDF.accountForTypesWithPrefixRobbins(aTypePrefix)
        return newCDF
        
def tauOfType(aType,N):
    """
    For aType, obtain the probability distribution p associated to T by dividing
    entries of aType by N.
    Compute and return the mutual information statistic MI(p).
    """
    firstMarginalUnnormalized = aType[0] + aType[1]
    secondMarginalUnnormalized = aType[0] + aType[2]
    if (firstMarginalUnnormalized == 0 or firstMarginalUnnormalized == N) or (secondMarginalUnnormalized == 0) or (secondMarginalUnnormalized == N):
        return 0
    pA0 = firstMarginalUnnormalized/(1.0*N)
    pB0 = secondMarginalUnnormalized/(1.0*N)
    pA1 = 1.0 - pA0
    pB1 = 1.0 - pB0
    MProjectionParameters = [pA0*pB0, pA0*pB1, pA1*pB0, pA1*pB1]
    result = 0
    for i in range(4):
        if aType[i] > 0:
            typeParameter = aType[i]/(1.0*N)
            result += typeParameter*np.log(typeParameter/MProjectionParameters[i])
    return result
        

def returnCDFAccountingForTypesOfModuloClassk(N,eta,k):
    resultCDF = CDF()
    p_eta = pdf.probabilityDistributionFactory(2,2).get_p_eta(eta)
    resultCDF.referenceDistribution = p_eta
    resultCDF.setN(N)
    resultCDF.accountForTypesForWhichFirstEntryIs_k_Mod_10(k)
    return resultCDF

class CDF(object):

    def __init__(self):
        '''
        Parameterless Constructor: Start with empty Discontinuity list and Probability list
        '''
        self.AscendingDiscontinuityList = []
        self.Dictionary = dict([])  #keys are elements of AscendingDiscontinuityList, values are the probabilities
        self.referenceDistribution = None
        self.N = None
        self.n = None
        self.probDistributionFactory = pdf.probabilityDistributionFactory(2,2)
        
    def setReferenceDistribution(self, eta):
        """
        Set the reference distribution, for calculating emission probabilities of types,
        equal to p^\eta, the distribution of uniform marginals with MI(p^\eta)=\eta.
        """
        p_eta = self.probDistributionFactory.get_p_eta(eta)
        self.referenceDistribution = p_eta
    
    def setReferenceDistributionFromDistribution(self, distribution):
        """
        If the distribution is already formed
        """
        self.referenceDistribution = distribution
        
    
    def accountForEvent(self,value, probability):
        """
        Account for the probability by adding probability to all the existing probabilities
        in the dictionary self.Dictionary for discontinuity jumps greater than or equal to the value.
        N.B. that because of numerical issues, this method may not be suitable for use in
        the circumstance where many small probabilities add to significantly change the cumulative probability,
        """
        insertionOfDisconinuity = True
        positionAtWhichValueWouldBeInserted = bisect.bisect(
             self.AscendingDiscontinuityList, value) 
        if (self.AscendingDiscontinuityList) and\
            self.AscendingDiscontinuityList[positionAtWhichValueWouldBeInserted-1] == value: #insertion not needed
            firstPositionForAppendingProbability = positionAtWhichValueWouldBeInserted - 1 
            insertionOfDisconinuity = False
        else:     #insertion is needed
            bisect.insort(self.AscendingDiscontinuityList, value)
            firstPositionForAppendingProbability = positionAtWhichValueWouldBeInserted
        if insertionOfDisconinuity:
            if firstPositionForAppendingProbability > 0:
                valueToLeftOfNewDiscontinuity = self.AscendingDiscontinuityList[
                    firstPositionForAppendingProbability-1]
                #print self.Dictionary
                self.Dictionary[value] = self.Dictionary[valueToLeftOfNewDiscontinuity]
            else:
                self.Dictionary[value] = 0
        for discontinuityListIndex in range(firstPositionForAppendingProbability,
                len(self.AscendingDiscontinuityList)):
            valueToBeIncremented = self.AscendingDiscontinuityList[discontinuityListIndex]
            self.Dictionary[valueToBeIncremented] += probability
    
    def mergeInto(self, another):
        """
        Account for the probability mass currently accounted for by another CDF object.
        """
        if not self.AscendingDiscontinuityList:  #nothing to merge into: put another's data into self
            self.AscendingDiscontinuityList = another.AscendingDiscontinuityList
            self.Dictionary = another.Dictionary
            return
        if not another.AscendingDiscontinuityList:  #nothing to merge: self stays the same
            return
        #another has at least one discontinuity to merge into self
        lowElementIndex = len(another.AscendingDiscontinuityList) - 1
        highElementIndex = lowElementIndex + 1
        while lowElementIndex >=0: #len(another.AscendingDiscontinuityList):
            lowValue = another.AscendingDiscontinuityList[lowElementIndex]
            lowProbability = another.Dictionary[lowValue]
            positionAtWhichValueWouldBeInserted = bisect.bisect(self.AscendingDiscontinuityList, lowValue) 
            if self.AscendingDiscontinuityList[positionAtWhichValueWouldBeInserted - 1] == lowValue:
                firstPositionForAppendingProbability = positionAtWhichValueWouldBeInserted - 1
                insertionOfDiscontinuity = False
            else:
                bisect.insort(self.AscendingDiscontinuityList, lowValue)
                firstPositionForAppendingProbability = positionAtWhichValueWouldBeInserted
                insertionOfDiscontinuity = True
            if insertionOfDiscontinuity:
                if firstPositionForAppendingProbability > 0:
                    valueToLeftOfNewDiscontinuity = self.AscendingDiscontinuityList[firstPositionForAppendingProbability-1]
                    #print self.Dictionary
                    self.Dictionary[lowValue] = self.Dictionary[valueToLeftOfNewDiscontinuity]
                else:
                    self.Dictionary[lowValue] = 0
            if highElementIndex < len(another.AscendingDiscontinuityList):
                highValue = another.AscendingDiscontinuityList[highElementIndex]
                #increment only up to the high value
                upperIndex = bisect.bisect_left(self.AscendingDiscontinuityList, highValue)
                #print "Upper index = %s"%(upperIndex)
            else:
                upperIndex = len(self.AscendingDiscontinuityList)
                #increment all elements to the right
                #print "Upper Index length of self.AscendingDiscontinuityList = %s"%(upperIndex)
            for discontinuityListIndex in range(firstPositionForAppendingProbability,upperIndex):
                valueToBeIncremented = self.AscendingDiscontinuityList[discontinuityListIndex]
                self.Dictionary[valueToBeIncremented] += lowProbability
            lowElementIndex -= 1
            highElementIndex -= 1
        return
    
    def mergeListInto(self,listOfOthers):
        for CDF in listOfOthers:
            self.mergeInto(CDF)
        

    def reconstructFromFile(self, filepath):
        """Read from a csv file where lines are of the form
        discontinuity_value_i, cumulativeProbability_i
        """
        f = open(filepath, 'rt')
        try:
            reader = csv.reader(f)
            reader.next()
            for row in reader:
                value = float(row[0])
                cumulativeProbability = float(row[1])
                bisect.insort(self.AscendingDiscontinuityList,value)
                self.Dictionary[value] = cumulativeProbability
        finally:
            f.close()

    #only works if the AscendingDiscontinuityList and Dictionary data members are filled
    def assignCumulativeProbability(self,intermediateValue, debug = False):
        """ 
        Assignment is made on basis of the stored values of discontinuities
        and dictionary of probabilities at those discontinuities.
        """
        positionAtWhichValueWouldBeInserted = bisect.bisect(self.AscendingDiscontinuityList, intermediateValue) 
        if positionAtWhichValueWouldBeInserted > 0:
            if debug == True:
                print "Intermediate val %s would be inserted between %s and %s"%(intermediateValue,self.AscendingDiscontinuityList[positionAtWhichValueWouldBeInserted-1],
                                                                                 self.AscendingDiscontinuityList[positionAtWhichValueWouldBeInserted] )
                print "and is assigned value %s"%(self.Dictionary[self.AscendingDiscontinuityList[positionAtWhichValueWouldBeInserted-1]])
            return self.Dictionary[self.AscendingDiscontinuityList[positionAtWhichValueWouldBeInserted-1]]
        else:
            return self.Dictionary[self.AscendingDiscontinuityList[0]]
    
    def assignCumulativeProbabilityToLinspace(self,linSpaceMin, linSpaceStep, numSteps):
        linspace = [linSpaceMin+i*linSpaceStep for i in range(numSteps)]
        return [self.assignCumulativeProbability(intermediateValue) for intermediateValue in linspace]
    
    def setReferenceDistribution_p_eta(self,eta):
        k,l = 2,2  #hardcoding binary variables for now
        self.referenceDistribution = pdf.probabilityDistributionFactory(k,l).get_p_eta(eta)
    
    def setN(self, N):
        self.N = N
    
    def setn(self,n):
        self.n = n
    
    
    def accountForTypeRobbins(self, aType):
        """
        Modify the CDF of the Mutual Information (MI) random variable
        to account for the probability of emission of aType using Robbins' approximation
        """
        if not self.N:
            raise ValueError("N not set, so we cannot account for type.")
        N=self.N
        self.accountForEvent(tauOfType(aType, N), 
            self.referenceDistribution.RobbinsEstimatedEmissionProbability(aType, N))
    
    def accountForType(self, aType):
        """
        Modify the CDF of the Mutual Information (MI) random variable
        to account for the probability of emission of aType
        """
        if not self.N:
            raise ValueError("N not set, so we cannot account for type.")
        N=self.N
        self.accountForEvent(tauOfType(aType, N), 
            self.referenceDistribution.exactEmissionProbability(aType, N))
    
    
    def accountForTypesWithPrefix(self, aTypePrefix):
        """
        For an object of class typePrefix with data (02) (when N=3, n=4),
        processes the type (0201) and (0210)
        """    
        aTypePrefix.DFSProcess(self.accountForType)
    
    def accountForTypesWithPrefixRobbins(self, aTypePrefix):
        """
        For an object of class typePrefix with data (02) (when N=3, n=4),
        processes the type (0201) and (0210): uses Robbins' approximation
        as the method for processing a leaf node (true type)
        """    
        aTypePrefix.DFSProcess(self.accountForTypeRobbins)

    def accountForAllTypes(self):
        """
        Generate all the type classes of size N and account for all types in the CDF:
        currently assumes the random variable whose CDF we are generating is the MI r.v.
        """
        if not self.N or not self.n:
            raise ValueError("Cannot generate types unless both n and N are set.")
        rootPrefix = tp.typePrefix(self.N, [], self.n)
        self.accountForTypesWithPrefix(rootPrefix)

    def accountForAllTypesRobbins(self):
        """
        Generate all the type classes of size N and account for all types in the CDF:
        currently assumes the random variable whose CDF we are generating is the MI r.v.
        """
        if not self.N or not self.n:
            raise ValueError("Cannot generate types unless both n and N are set.")
        rootPrefix = tp.typePrefix(self.N, [], self.n)
        self.accountForTypesWithPrefixRobbins(rootPrefix)
     
    def accountForAllTypesParallelized(self, modulus):
        """
        Break the job into account for all types which have the prefix (0), (1),...(modulus - 1) % modulus
        """
        if not self.N or not self.n:
            raise ValueError("Cannot generate types unless both n and N are set.")
        rootPrefix = tp.typePrefix(self.N, [], self.n)
        rootsForParallelJobs = [rootPrefix.childrenListLastEntry_k_Mod_m(k, modulus) for k in range(modulus)]
        ListOfCDFs = Parallel(n_jobs=-1, verbose=50)(delayed(accountForListOfTypePrefixes)(
                                        aRootList,self.N, self.n, self.referenceDistribution) 
                                        for aRootList in rootsForParallelJobs)
        self.mergeListInto(ListOfCDFs)
        
    def accountForAllTypesParallelizedRobbins(self, modulus):
        """
        Break the job into account for all types which have the prefix (0), (1),...(modulus - 1) % modulus
        """
        if not self.N or not self.n:
            raise ValueError("Cannot generate types unless both n and N are set.")
        rootPrefix = tp.typePrefix(self.N, [], self.n)
        rootsForParallelJobs = [rootPrefix.childrenListLastEntry_k_Mod_m(k, modulus) for k in range(modulus)]
        ListOfCDFs = Parallel(n_jobs=-1, verbose=50)(delayed(accountForListOfTypePrefixesRobbins)(
                                        aRootList,self.N, self.n, self.referenceDistribution) 
                                        for aRootList in rootsForParallelJobs)
        self.mergeListInto(ListOfCDFs)
        
    def toDataFrame(self):
        """
        Export a calculated CDF to a dataframe: convenient for storing in 
        and recovering from .csv format
        Columns:
        gamma,beta
        where gamma,beta are the gamma,beta values at all the discontinuities
        """
        df = pandas.DataFrame(np.array(self.Dictionary.items()), columns=['gamma','beta'])
        return df        
    
        
    
    