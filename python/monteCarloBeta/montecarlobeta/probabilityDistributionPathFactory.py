'''
Created on Aug 8, 2012

@author: eliotpbrenner

TODO: explain purpose
'''

import probabilityDistributionPath as pdp
import probabilityDistribution as pd
import numpy as np


class probabilityDistributionPathFactory(object):
    '''
    ProbabilityDistributionPath Factory class
    '''

    def __init__(self, feasibleMarginals, k,l):
        self.feasibleMarginals = feasibleMarginals
        self.k = k
        self.l = l
        
    def construct(self):
        firstMarginalsAsMatrix = np.matrix([self.feasibleMarginals[0], 1.0 - self.feasibleMarginals[0]])
        secondMarginalsAsMatrix = np.matrix([self.feasibleMarginals[1], 1.0 - self.feasibleMarginals[1]])
        #construct the path through the product distribution with those marginals
        productDistributionWithFeasibleMarginals = pd.ProbabilityDistribution(self.k,self.l)
        productDistributionWithFeasibleMarginals.setParametersProductDistAssocWithMarginals(firstMarginalsAsMatrix,
                                                                                             secondMarginalsAsMatrix)
        return pdp.probabilityDistributionPath(productDistributionWithFeasibleMarginals)

def main():
    import itertools    
    displacements = np.arange(0.05,0.45,0.05)
    for displacedMarginals in itertools.product(displacements,displacements):
        probabilityDistribution = probabilityDistributionPathFactory(
          displacedMarginals,2,2).construct()
        
        
        
        
        
if __name__ == "__main__":
    main()     