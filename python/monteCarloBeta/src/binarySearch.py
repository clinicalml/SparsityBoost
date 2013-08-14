# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 12:45:30 2013

@author: eliotpbrenner

Utility class for performing binary search

"""
import math
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


class binarySearch(object):
    
    def __init__(self,tolerance_exp=7,maxDepth=500, increasingFunction=None):
        self.increasingFunction = increasingFunction
        self.tolerance = 10**(-tolerance_exp)
        self.maxDepth = maxDepth
        
    def search(self,x0,x1,y_hat,depth=0):
        logging.debug(
        "Performing binary search with x0=%s, x1=%s, y_hat=%s, depth=%s,"%(
        x0,x1,y_hat,depth))
        if self.increasingFunction is None:
            raise ValueError("Cannot perform binary search without an increasing function.")
        if depth > self.maxDepth:
            logging.warn("Maxdepth of %s in binary search reached with x0=%s, x1=%s, y_hat=%s"%(depth,x0,x1,y_hat))
            return x0
        if x0 > x1:
            raise ValueError("To define interval for search cannot have x0=%s > x1=%s"%(x0,x1))
        y0 = self.increasingFunction(x0)
        y1 = self.increasingFunction(x1)
        if y0 > y1:
            raise ValueError("Since x0=%s, x1=%s, y0=%s, y1=%s, the function is not increasing"%(x0,x1,y0,y1))
        x_mid = (x0+x1)/2.0
        y_mid = self.increasingFunction(x_mid)
        rho = y_hat - y_mid
        if math.fabs(rho) < self.tolerance:
            logging.debug("y_mid=%s was within tolerance=%s of y_hat=%s at search depth %s, so returning x_mid=%s"%(y_mid,self.tolerance,y_hat,depth,x_mid))
            return x_mid
        if rho > 0:
            return self.search(x_mid,x1,y_hat,depth+1)
        else:
            return self.search(x0,x_mid,y_hat,depth+1)
            
        

