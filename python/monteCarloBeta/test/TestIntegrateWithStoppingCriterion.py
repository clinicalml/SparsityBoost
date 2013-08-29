'''
Created on Aug 16, 2012

@author: eliotpbrenner
'''


import sys
sys.path.insert(0, '../src')
import unittest
from scipy import stats
import emissionProbabilityCalculator as epc
import mcIntegrationWithStoppingCriteria as iwsc
import informationTheory as it
import Timer as Timer
from numpy import genfromtxt
import pickle
import csv

import statsDistributionFactory as sdf

#logging
import logging


import os
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

CentralProbability = 0.68  #probability that a Gaussian distribution lies within 1 st. deviation of mean
benchMarkInterval_N = 100 #the kth line in the benchmark file contains results for N = (k+1)*benchMarkInterval_N


class IntegrateWithStoppingCriterionTest(unittest.TestCase):


    def _testIntegratetoFlatFile(self):        
        #load in benchmark data
        benchmarks = genfromtxt('cdf_for_valsets_of_variables_of_size2and2', delimiter=',', skip_header=4)
        benchmarks_without_stepsize = benchmarks[:,1:]
        lineInBenchMarkFileList = [0, 9]
        
        #prepare the csv file
        f=open('/home/eliot/Documents/sontag/benchmarkcomparisons/N@100_N@1000_eta@onehundredth.csv', 'wt')
        writer = csv.writer(f)
        for lineInBenchMarkFile in lineInBenchMarkFileList:
            benchmark_N = (lineInBenchMarkFile+1)*benchMarkInterval_N
            benchmarksForThisN = benchmarks_without_stepsize[lineInBenchMarkFile]
            stepSizeInBenchmarks = benchmarks[0,0]
            numStepsInBenchmarks = len(benchmarks[0])
            logging.info("Read in %s lines of benchmarks with stepsize=%s and %s steps"%(len(benchmarks_without_stepsize), stepSizeInBenchmarks, numStepsInBenchmarks))
            
            #constants common to all of the runs
            k, l = 2,2
            eta = 0.01
            
            N = benchmark_N
            maxIterations = 10000
            dictionaryFromNameOfMethodToResults = {}
            
            #build the inputs to IntegrateWithStoppingCriterion
            #the integrand, namely probability-of-emission calculator from p-eta, which is common to ALL methods
            probabilityCalculatorAssociatedToEta = epc.emissionProbabilityCalculator(eta, k, l, N)
            
            for gammaNumber in range(numStepsInBenchmarks-1): 
                gamma = (gammaNumber+1)*stepSizeInBenchmarks
                benchMarkBeta = benchmarksForThisN[gammaNumber]
                probabilityCalculatorAssociatedToEta.setGamma(gamma)
                functionProbabilityOfEmissionByP_eta_Robbins = probabilityCalculatorAssociatedToEta.RobbinsEstimateOfEmissionProbabilityTimesCharFunctionOfTauMinusGamma
                
                #the probability distribution used to choose marginals distribution, also common to ALL methods (except for completely-uniform, if implemented)
                Gaussianscale = it.ChernoffRadius(CentralProbability, N)
                normalDistLocUniformScaleDepOnChernoff = stats.norm(loc=0.5,scale=Gaussianscale)
                
                
                #e.g. 5, .95 means "seeking 95 percent certainty answer is within 5 percent of truth , typically also common to all runs
                accuracy_percent = 10
                probability_accuracy_achieved = 0.95
                
                #the method for picking the parameter-distribution based on the marginals, specific to this run
                #
                statsDistFactory = sdf.statsDistributionFactory(gamma)
                dictionaryFromNameOfMethodToMethod = {}
                dictionaryFromNameOfMethodToMethod['GaussianScaleAdaptiveTo_l_gamma']=statsDistFactory.GaussianCenteredAttGammaPlusfromMarginals
                #dictionaryFromNameOfMethodToMethod['GaussianScaleAdaptiveTo_N']=statsDistFactory.GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio
                
               
                #functionFromMarginalsTo_t_Distribution = statsDistFactory.GaussianCenteredAttGammaPlusfromMarginals
                
                uniformMarginals = [1.0/k, 1.0/l]
                
                for nameOfMethod in dictionaryFromNameOfMethodToMethod.keys():
                    if dictionaryFromNameOfMethodToResults.has_key( nameOfMethod):
                        resultsDictionary = dictionaryFromNameOfMethodToResults[nameOfMethod]
                    else:
                        resultsDictionary = {}
                    statsDistFactory.set_eta(eta)
                    statsDistFactory.set_N(N)
                    statsDistFactory.CalculateAndSetScaleRatio(uniformMarginals)
                    #functionFromMarginalsTo_t_Distribution = statsDistFactory.GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio
                    functionFromMarginalsTo_t_Distribution = dictionaryFromNameOfMethodToMethod[nameOfMethod]
                     
                    #start the calculation proper
                    with Timer.Timer() as t:
                        iteration_values, rho_hat_values = iwsc.IntegrateWithStoppingCriterion(k,l,
                                               functionProbabilityOfEmissionByP_eta_Robbins, 
                                               normalDistLocUniformScaleDepOnChernoff, functionFromMarginalsTo_t_Distribution,  #!!!MAIN THING TO CHANGE!!!
                                               maxIterations, 100, 100, #frequencyOfRecodingreslt, frequencyOfApplicationOfStoppingCriterion 
                                               accuracy_percent, probability_accuracy_achieved)
                        
                    logging.info("Time elapsed is %.03f sec." % t.interval)
                    logging.info("Iteration_values for the method %s for eta=%s, gamma=%s, N=%s, accpercent=%s, prob_achieved=%s"%(nameOfMethod, eta,gamma,N,accuracy_percent, probability_accuracy_achieved))
                    logging.info(iteration_values)
                    logging.info("rho_hat_values for the method %s for eta=%s, gamma=%s, N=%s, accpercent=%s, prob_achieved=%s"%(nameOfMethod, eta,gamma,N,accuracy_percent, probability_accuracy_achieved))
                    logging.info( rho_hat_values)
                    betaResult = rho_hat_values[-1]*N**3
                    if resultsDictionary.has_key('time'):
                        resultsDictionary['time'].append(t.interval)
                    else:
                        resultsDictionary['time'] = [t.interval]
                        
                    if resultsDictionary.has_key('CalculatedResult'):
                        resultsDictionary['CalculatedResult'].append(betaResult)
                    else:
                        resultsDictionary['CalculatedResult'] = [betaResult]
                    
                    if resultsDictionary.has_key('BenchmarkBeta'):
                        resultsDictionary['BenchmarkBeta'].append(benchMarkBeta)
                    else:
                        resultsDictionary['BenchmarkBeta'] = [benchMarkBeta]
                        
                        
                    dictionaryFromNameOfMethodToResults[nameOfMethod] = resultsDictionary
                    fileNameForPickleDump = "/home/eliot/Documents/sontag/benchmarkcomparisons/testRun" + str(N)
                    pickle.dump( dictionaryFromNameOfMethodToResults, open( fileNameForPickleDump, "wb" ) )
            
            #write to the CSV file
            N_header = 'N=%s'%(N)
            writer.writerow( [N_header])
            for resultsDictionary in dictionaryFromNameOfMethodToResults.values():
                for key in resultsDictionary.keys():
                    print key
                    print resultsDictionary[key]
                    rowList = resultsDictionary[key]
                    rowList.insert(0,key)
                    print rowList
                    writer.writerow(rowList)
        f.close()



    def testIntegratetoPickle(self):      
        #load in benchmark data
        benchmarks = genfromtxt('cdf_for_valsets_of_variables_of_size2and2', delimiter=',', skip_header=5)
        benchmarks_without_stepsize = benchmarks[:,1:]
        lineInBenchMarkFileList = [4, 9, 14, 18]
        for lineInBenchMarkFile in lineInBenchMarkFileList:
            benchmark_N = (lineInBenchMarkFile+1)*benchMarkInterval_N
            benchmarksForThisN = benchmarks_without_stepsize[lineInBenchMarkFile]
            stepSizeInBenchmarks = benchmarks[0,0]
            numStepsInBenchmarks = len(benchmarks[0])
            logging.info("Read in %s lines of benchmarks with stepsize=%s and %s steps"%(len(benchmarks_without_stepsize), stepSizeInBenchmarks, numStepsInBenchmarks))
            
            #constants common to all of the runs
            k, l = 2,2
            eta = 0.01
            
            N = benchmark_N
            maxIterations = 10000
            dictionaryFromNameOfMethodToResults = {}
            
            #build the inputs to IntegrateWithStoppingCriterion
            #the integrand, namely probability-of-emission calculator from p-eta, which is common to ALL methods
            probabilityCalculatorAssociatedToEta = epc.emissionProbabilityCalculator(eta, k, l, N)
            
            for gammaNumber in range(numStepsInBenchmarks-1): 
                gamma = (gammaNumber+1)*stepSizeInBenchmarks
                benchMarkBeta = benchmarksForThisN[gammaNumber]
                probabilityCalculatorAssociatedToEta.setGamma(gamma)
                functionProbabilityOfEmissionByP_eta_Robbins = probabilityCalculatorAssociatedToEta.RobbinsEstimateOfEmissionProbabilityTimesCharFunctionOfTauMinusGamma
                
                #the probability distribution used to choose marginals distribution, also common to ALL methods 
                #(except for completely-uniform, if implemented)
                Gaussianscale = it.ChernoffRadius(CentralProbability, N)
                normalDistLocUniformScaleDepOnChernoff = stats.norm(loc=0.5,scale=Gaussianscale)
                
                
                #e.g. 5, .95 means "seeking 95 percent certainty answer is within 5 percent of truth , typically also common to all runs
                accuracy_percent = 10
                probability_accuracy_achieved = 0.95
                
                #the method for picking the parameter-distribution based on the marginals, specific to this run
                #
                statsDistFactory = sdf.statsDistributionFactory(gamma)
                dictionaryFromNameOfMethodToMethod = {}
                dictionaryFromNameOfMethodToMethod['GaussianScaleAdaptiveTo_l_gamma']=statsDistFactory.GaussianCenteredAttGammaPlusfromMarginals
                dictionaryFromNameOfMethodToMethod['GaussianScaleAdaptiveTo_N']=statsDistFactory.GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio
                
               
                #functionFromMarginalsTo_t_Distribution = statsDistFactory.GaussianCenteredAttGammaPlusfromMarginals
                
                uniformMarginals = [1.0/k, 1.0/l]
                
                for nameOfMethod in dictionaryFromNameOfMethodToMethod.keys():
                    if dictionaryFromNameOfMethodToResults.has_key( nameOfMethod):
                        resultsDictionary = dictionaryFromNameOfMethodToResults[nameOfMethod]
                    else:
                        resultsDictionary = {}
                    statsDistFactory.set_eta(eta)
                    statsDistFactory.set_N(N)
                    statsDistFactory.CalculateAndSetScaleRatio(uniformMarginals)
                    #functionFromMarginalsTo_t_Distribution = statsDistFactory.GaussianCenteredAttGammaPlusfromMarginalsScaledFromScaleRatio
                    functionFromMarginalsTo_t_Distribution = dictionaryFromNameOfMethodToMethod[nameOfMethod]
                     
                    #start the calculation proper
                    with Timer.Timer() as t:
                        iteration_values, rho_hat_values = iwsc.IntegrateWithStoppingCriterion(k,l,
                                               functionProbabilityOfEmissionByP_eta_Robbins, 
                                               normalDistLocUniformScaleDepOnChernoff, functionFromMarginalsTo_t_Distribution,  #!!!MAIN THING TO CHANGE!!!
                                               maxIterations, 100, 100, #frequencyOfRecodingreslt, frequencyOfApplicationOfStoppingCriterion 
                                               accuracy_percent, probability_accuracy_achieved)
                        
                    logging.info("Time elapsed is %.03f sec." % t.interval)
                    logging.info("Iteration_values for the method %s for eta=%s, gamma=%s, N=%s, accpercent=%s, prob_achieved=%s"%(nameOfMethod, eta,gamma,N,accuracy_percent, probability_accuracy_achieved))
                    logging.info(iteration_values)
                    logging.info("rho_hat_values for the method %s for eta=%s, gamma=%s, N=%s, accpercent=%s, prob_achieved=%s"%(nameOfMethod, eta,gamma,N,accuracy_percent, probability_accuracy_achieved))
                    logging.info( rho_hat_values)
                    betaResult = rho_hat_values[-1]*N**3
                    if resultsDictionary.has_key('time'):
                        resultsDictionary['time'].append(t.interval)
                    else:
                        resultsDictionary['time'] = [t.interval]
                    if resultsDictionary.has_key('ratioOfResultToBenchmark'):
                        resultsDictionary['ratioOfResultToBenchmark'].append(betaResult/benchMarkBeta)
                    else:
                        resultsDictionary['ratioOfResultToBenchmark'] = [betaResult/benchMarkBeta]
                    dictionaryFromNameOfMethodToResults[nameOfMethod] = resultsDictionary
                    fileNameForPickleDump = "testRun" + str(N)
                    pickle.dump( dictionaryFromNameOfMethodToResults, open( fileNameForPickleDump, "wb" ) )

if __name__ == "__main__":
    unittest.main()