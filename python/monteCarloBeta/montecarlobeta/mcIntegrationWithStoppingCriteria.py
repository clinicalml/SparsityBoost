'''
Created on Aug 15, 2012

@author: eliotpbrenner

Functions which is mean to perform Monte Carlo integration on sets defined by
        \tau(p) < \gamma
of integrands, with stopping criterion based on the formula (4.6) of Bucklew's book
'''

import logging
# create logger
logging.basicConfig(format='%(levelname)s:%(filename)s: %(message)s', level=logging.INFO)


from scipy import stats
def IntegrateWithStoppingCriterion(k,l,
                                   integrand, 
                                   marginalsSamplingDistribution, functionFromMarginalsToProbabilityDistributionOnFibers,
                                   maxIterations, frequencyOfRecordingEstimate, frequencyOfApplicationOfStoppingCriterion,
                                   desiredPrecisionPercent, desiredConfidence):
    """Returns:
    iteration_values, rho_hat_values
    Parameters:
     k,l: |valset|'s of marginals
     integrand: typically a probability over the characteristic set defined by gamma, e.g. probabilityCalculatorAssociatedToEta.RobbinsEstimateOfEmissionProbabilityTimesCharFunctionOfTauMinusGamma
     marginalsSamplingDistribution: of data type of type <class 'scipy.stats.distributions.rv_frozen'>, with an rvs and pdf method, e.g. Gaussian with radius determined by Chernoff bounds
     functionFromMarginalsToProbabilityDistributionOnFibers: e.g. output a Gaussian in the parameter t located at t_gamma_plus and scaled according to decay of integrand
     maxIterations: in no case will we use more than this many trials
     desiredPrecisionPercent: the stopping criterion tries to avoid an error of greater than this percent in the result, x in Bucklew's notation
     desiredConfidence: the stopping criterion is that the precision is achieved with at least this probability, y in Bucklew's notation
     frequencyOfApplicationOfStoppingCriterion: test every this-number-of-trials to see if the stopping criterion is actually satisfied. 
    """
    logger = logging.getLogger('mc_application')

    logger.debug('Entering a call to IntegrateWithStoppingCriterion')
    
    #Calculate t_y from y=desiredConfidence
    t_y = stats.norm.isf((1-desiredConfidence)/2.0)
    #constant part of k*
    k_starConstantFactor = (t_y*100/desiredPrecisionPercent)**2
    
    marginalsRandomDraw = marginalsSamplingDistribution.rvs
    marginalsPDF = marginalsSamplingDistribution.pdf
    
    
    
    # store every N intermediate integral approximations in an # array I and record the corresponding k value
    rho_hat_values = []  
    F_hat_values = []
    iteration_values = []
    iteration = 1
    stoppingCriterion = False
    s=0
    F_s = 0  #cumulative sum to form the quantity F on the bottom of p. 71 of Bucklew book
    while iteration < maxIterations+1 and stoppingCriterion == False:
        #randomly draw the marginals
        marginal1 = marginalsRandomDraw()
        marginal2 = marginalsRandomDraw()
        if marginal1 <= 0 or marginal1 >= 1 or marginal2 <=0 or marginal2 >= 1:
            s_summand = 0
            F_s_summand = 0
        else:
            randomDrawDistribution_t_Parameter = functionFromMarginalsToProbabilityDistributionOnFibers([marginal1, marginal2])
            t_ParamterRandomDraw = randomDrawDistribution_t_Parameter.rvs
            t_ParameterPDF =  randomDrawDistribution_t_Parameter.pdf
            #randomly draw the t-parameter
            t_parameter = t_ParamterRandomDraw()
            #x3 = randomDrawFunction3()
            #print x1,x2,x3
            s_summand = integrand(marginal1,marginal2,t_parameter)/(marginalsPDF(marginal1)*marginalsPDF(marginal2)*t_ParameterPDF(t_parameter))
            F_s_summand = s_summand**2
        s += s_summand
        F_s += F_s_summand
        if iteration % frequencyOfRecordingEstimate == 0:
            I = (1.0/iteration)*s
            rho_hat_values.append(I)
            iteration_values.append(iteration)
        if iteration % frequencyOfApplicationOfStoppingCriterion == 0:
            F = (1.0/iteration)*F_s
            k_starOfIterations = k_starConstantFactor*(F/(I*I) - 1.0 )  #(4.6) in Bucklew, top on p. 72
            if iteration >= k_starOfIterations: 
                stoppingCriterion = True
                logger.warning('Stopping integration with iteration=%s exceeding k_startOfIterations=%s'%(iteration, k_starOfIterations))
            else:
                logger.info('iteration=%s less than k_startOfIterations=%s'%(iteration, k_starOfIterations))
        iteration += 1
    return iteration_values, rho_hat_values

