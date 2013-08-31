SparsityBoost
=============



A repository holding the code implementing the methods of 


"SparsityBoost: A New Scoring Function for Learning Bayesian
Network Structure"
by Eliot Brenner and David Sontag


Presented at Uncertainty In Artificial Intelligence 2013, Friday July 12, Bellevue, Washington, USA

Published in <i>Uncertainty in Artificial Intelligence: Proceedings of the Twenty-Eighth Conference (2013)</i>.


http://cs.nyu.edu/~dsontag/papers/BrennerSontag_uai13.pdf

Installation:
==============

To install the c/ part of the library (the only part needed for learning networks) first type:

make

from the c/ directory.  This should create the bscore executable file in your c/ directory.

Basic Usage:
==============
The library exists in two independent parts, in different languages and with different purposes: the parts under c/, in C/C++, are for learning networks, and the parts under python/monteCarloBeta
are for computing the values of the CDF of the Mutual Information ("beta") needed for the interpolation.  

We recommend starting by trying to learn the structure of a network, by assuming that eta (a lower bound for the edge strength) takes one of the following values:

0.04, 0.02, 0.01, 0.005, 0.001

In this case, the user only needs to run the C/C++ parts of the library.  A typical invocation of the main method, ../c/bscore, looks like this:

../c/bscore -bic -nodes 37 -maxpa 4 -data ../data/synthetic_examples/experiments/0/alarm1000.dat -mod_out ../results/experiments/0/sparsity_eps005/model1000.mod -edge_scores .005

The .dat file consists of raw data, in the format of a space-separated file.
Each row holds one observation which is a line of 0s and 1s, one column for each variable in the system.
The .mod file is in GOBNILP Score File Format, described on p. 11 of the GOBNILP manual:

http://www.cs.york.ac.uk/aig/sw/gobnilp/manual.pdf

In particular, the first line is the number of variables in the network, and what follows are the local parts of our score corresponding to a pair of edges 
(on a line by themsevles) and possible parent sets (at the beginning of the following lines with the scores). 

In order to use this mod file to actually learn the nework, you can use GOBNILP or any other system that learns network structure from local scores.  For example, using gobnilp1.3:

gobnilp1.3/bin/gobnilp ../results/experiments/0/sparsity_eps005/model4400.mod > ../results/experiments/0/sparsity_eps005/output4400.txt

The script perl/run_experiments.pl reproduces the experiments in the UAI paper. 

Reproducing UAI 2013 Results:
==============================
In order to run perl/run_experiments.pl you have to have 

gobnilp1.3/bin/gobnilp

under the perl directory.

In run_experiments.pl, comment in/out the sed commands which will work on your OS.  These lines are labelled "OSX" or "Linux". 
In the perl directory itself, you have to have the GOBNILP settings file gobnilp.set, which comes with this repository.  The sed commands will copy this GOBNILP settings file to the appropriate locations, such as../results/experiments/i/sparsity_eps005/, etc. 

The data for the experiment can be downloaded from the following location:

http://cs.nyu.edu/~dsontag/data/alarm_data_uai13.tar.gz

The archive file contains a structure like this:

data/synthetic_examples/experiments/9/alarm6600.dat

so that if you untar it directly in the toplevel ("SparsityBoost") folder of the project, the paths will be what the run_experiments.pl script expects. 

Computing the CDF of the Mutual Information ("beta")
====================================================
The reason for doing this is if you wish to learn a network assuming a minimal edge strength ("epsilon") other than the four values cited above.  Use the class betaTable.py in MonteCarlo/montecarlo/.  Instructions to follow.


Dependencies:
=============
Required for the second step of learning the network and for reproducing our results from the UAI paper:

Gobnilp1.3 and its dependency Scip

General installation and usage instructions for gobnilp1.3 are given in the GOBNILP manual linked above.

run_experiments.pl has been tested on OSX 10.8 and Ubuntu Linux.


License:
========
SparsityBoost is free software; you can redistribute it
and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation

SparsityBoost is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with SparsityBoost; see the file gpl-3.0.txt.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

####################################################
