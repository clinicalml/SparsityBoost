SparsityBoost
=============



A repository holding the code implementing the methods of 


"SparsityBoost: A New Scoring Function for Learning Bayesian
Network Structure"
by Eliot Brenner and David Sontag


Presented at Uncertainty In Artificial Intelligence 2013, Friday July 12, Bellevue, Washington, USA

Published in <i>Uncertainty in Artificial Intelligence: Proceedings of the Twenty-Eighth Conference (2013)</i>.


http://cs.nyu.edu/~dsontag/papers/BrennerSontag_uai13.pdf

Basic Usage:
==============
The library exists in two independent parts, in different languages and with different purposes: the parts under c/, in C/C++, are for learning networks, and the parts under python/monteCarloBeta
are for computing the values of beta needed for the interpolation.  We recommend starting by trying to learn the structure of a network, by assuming that eta (a lower bound for the edge strength) takes one of the following values:

0.04, 0.02, 0.01, 0.005, 0.001

In this case, the user only needs to run the C/C++ parts of the library.  A typical invocation of the main method, ../c/bscore, looks like this:

../c/bscore -bic -nodes 37 -maxpa 4 -data ../data/synthetic_examples/experiments/0/alarm4400.dat -mod_out ../results/experiments/0/sparsity_eps005/model4400.mod -edge_scores .005
 gobnilp1.3/bin/gobnilp -g../results/experiments/0/sparsity_eps005/gobnilp4400.set ../results/experiments/0/sparsity_eps005/model4400.mod > ../results/experiments/0/sparsity_eps005/output4400.txt


Dependencies:
=============

run_experiments.pl has been tested on OSX 10.8 and Ubuntu Linux.

Gobnilp1.3 and its dependency Scip

In order to run perl/run_experiments.pl you have to have 

gobnilp1.3/bin/gobnilp

under the perl directory.u

In ther perl directory itself, you have to have ../perl/gobnilp.set

In run_experiments.pl, comment in/out the sed commands which will work on your OS.  These lines are labelled "OSX" or "Linux". 

