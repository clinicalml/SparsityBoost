'''
Created on Aug 12, 2012

@author: eliot
'''
import numpy as np

default_intervalTolerance = 1e-15
default_valueTolerance = 1e-200



import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
import os
homeDirectory = os.getenv("HOME", "/Users/eliotpbrenner")


class functionAlgorithms(object):
    '''
    persistent class with a function: perform various algorithms on it such as binary search
    for argument with given value
    '''


    def __init__(self,theFunction):
        '''
        initialize the function.  Interval tolerance and value tolerance are set to 1e-10 by default
        '''
        self.theFunction = theFunction
        self.intervalTolerance = default_intervalTolerance
        self.valueTolerance = default_valueTolerance
        self.fixedArgumentList = None
    
    def setFixedArgumentList(self, fixedArgumentList):
        self.fixedArgumentList = fixedArgumentList
    
    def searchArgWhereIncreasingFunctionTakesVal(self, theValue, lowerBound, upperBound):
        '''
        Assuming theIncreasingFunction is increasing, finds the argument where
        the function takes the value (or returns the upper or lower bound of the interval if)
        the function value is outside the interval and above/below it, resp.
        '''
        #logging.debug("Searching for theValue %s between %s and %s"%(theValue, lowerBound, upperBound))
        if np.abs(upperBound - lowerBound) < self.intervalTolerance:
            return lowerBound
        if theValue < self.theFunction(lowerBound) + self.valueTolerance:
            return lowerBound
        if theValue > self.theFunction(upperBound) - self.valueTolerance:
            return upperBound
        midpoint = (lowerBound + upperBound)/2.0
        midpointValue = self.theFunction(midpoint)
        if np.abs(midpointValue - theValue) < self.valueTolerance:
            return midpoint
        if midpointValue < theValue:
            #logging.debug("With midpointvalue=%s < theValue=%s, resetting lower bound to midpoint=%s"%(midpointValue, theValue, midpoint))
            return self.searchArgWhereIncreasingFunctionTakesVal(theValue, midpoint, upperBound)
        else:
            #logging.debug("With midpointvalue=%s > theValue=%s, resetting upper bound to midpoint=%s"%(midpointValue, theValue, midpoint))
            return self.searchArgWhereIncreasingFunctionTakesVal(theValue, lowerBound, midpoint)
         
    def searchArgWhereIncreasingFunctionTakesProportionOfMaxVal(self, theProportion, lowerBound, upperBound):
        '''
        [Binary] search within [lowerBound, upperBound] for the argument where the function takes the value
                f(found argument) = theProportion*(f(upperBound))
        '''
        logging.debug("Searching for value of theProportion*self.theFunction(upperBound)=%s between lowerBound=%s and upperBound=%s"%(theProportion*self.theFunction(upperBound), 
                                                             lowerBound, upperBound))
        return self.searchArgWhereIncreasingFunctionTakesVal(theProportion*self.theFunction(upperBound), 
                                                             lowerBound, upperBound)
    
    def functionOfOneVariable(self,varyingArgument):
        '''
        Assuming the function is a function of n+1 variables and fixedArgumentList
        is a list of n fixed arguments, this is a function of 1 argument where only the final
        argument varies
        '''
        allArguments = [argument for argument in self.fixedArgumentList]
        allArguments.append(varyingArgument)
        #logging.debug("All arguments=%s"%(str(allArguments)))
        return self.theFunction(*allArguments)
