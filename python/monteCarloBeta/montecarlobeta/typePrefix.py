'''
Created on Mar 21, 2013

@author: eliotpbrenner
'''

import scipy.special
import numpy as np

class typePrefix(object):
    '''
    A prefix of a type of length n and size N
    Data members:
      N         size
      n         length
      data      [] for root; in general list of length <=n and sum <=N
      processed boolean used in traversals such as DFS
    '''


    def __init__(self, N, data = [], n = 4):
        '''
        Constructor
        '''
        self.N = N  #size of the type of which this is a prefix
        self.n = n  #length of the type of which this is a prefix
        self.data = data  #default is [], root: otherwise a list of length between and n 
        self.processed = False

    def hasChildren(self):
        '''
        Boolean saying whether other prefixes start with this prefix:
          the condition is just that the length of the data is STRICTLY less than n
        '''
        return len(self.data) < self.n
    
    def childrenList(self):
        '''
        If prefix is [x_1,...x_r] then [[x_1,...,x_r,i]] where i ranges from 0 to N-x_1-...x_r.
        '''
        if not self.hasChildren():
            raise ValueError("Prefix is complete; does not have children.")
        remainingMass = self.N - sum(self.data)
        if len(self.data) < self.n - 1:
            return [ typePrefix(self.N,self.data + [i], self.n ) for i in range(remainingMass+1) ]  #TODO: change to yield
        else:
            return [ typePrefix(self.N, self.data + [remainingMass],self.n) ]
    
    def childrenListLastEntry_k_Mod_m(self,k,m):
        '''
        If prefix is [x_1,...x_r] then [[x_1,...,x_r,i]] where i%m = k and i ranges from 0 to N-x_1-...x_r.
        '''
        if not self.hasChildren():
            raise ValueError("Prefix is complete; does not have children.")
        remainingMass = self.N - sum(self.data)
        if len(self.data) < self.n - 1:
            return [ typePrefix(self.N,self.data + [i], self.n ) for i in range(remainingMass+1) if i % m == k]  #TODO: change to yield
        else:
            if remainingMass % m == k:
                return [ typePrefix(self.N, self.data + [remainingMass],self.n) ]
            else:
                return []

    def childrenGenerator(self):
        '''
        If prefix is [x_1,...x_r] then [[x_1,...,x_r,i]] where i ranges from 0 to N-x_1-...x_r.
        '''
        if not self.hasChildren():
            raise ValueError("Prefix is complete; does not have children.")
        remainingMass = self.N - sum(self.data)
        if len(self.data) < self.n - 1:
            return (typePrefix(self.N,self.data + [i], self.n ) for i in range(remainingMass+1) )  #TODO: change to yield
        else:  
            return ( typePrefix(self.N, self.data + [remainingMass],self.n) for i in range(1))
        
    def DFSPrint(self):
        '''
        DFSProcess where the process is simply to print the data of the prefix:
          primarily useful for debugging.
        '''
        def printMethod(anObject):
            print anObject
        self.DFSProcess(printMethod)
        
    def DFSProcess(self, processMethod):
        '''
        Depth first processing.  Example:
          globalList = []
          rootPrefix = tp.typePrefix(3)
          rootPrefix.DFSProcess(globalList.append)
        Afterward, globalList has all types of length 4 and size 3
        '''
        if self.hasChildren():
            for child in self.childrenGenerator():
                if not child.processed:
                    child.DFSProcess(processMethod)
            self.processed = True
        else:
            processMethod(self.data)
            self.processed = True
            
         
    def numDescendants(self):
        '''
        Useful for performing fast tests/checks/timing estimates
          the number of children and children of children...
        '''
        remainingMass = self.N - sum(self.data)
        remainingLength = self.n - len(self.data)
        #special cases
        if not self.hasChildren():
            return 0
        if remainingLength < 1:
            raise ValueError("Recorded a prefix with children but no remainingLength!")
        return int(scipy.special.binom(remainingMass + remainingLength - 1,remainingLength - 1))
    
    def reshapeAsDistribution(self,k):
        return (1.0/self.N)*np.reshape(self.data, (k,-1))
    
            