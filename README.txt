Seismic-2d
==========
This is a simplification of a model originally used for global-scale
seismology.

NET-VISA: Network Processing Vertically Integrated Seismic
Analysis. Nimar S. Arora, Stuart Russell, Erik Sudderth. Bulletin of the
Seismological Society of America (BSSA) April 2013, vol. 103 no. 2A
pp709-729.

Files
=====

description.odt -- a description of the model
generate.py   -- generates the physics of a 2-D world and some episodes
util.py       -- some geophysical utility functions
solve.py      -- a sample solver that learns the physics and solves the episodes
evaluate.py   -- evaluates a solution versus a reference
mwmatching.py -- utility script for max-weight max cardinality matching
training.data -- 10K episodes for training
test.data     -- another 10K episodes for testing
test.blind    -- the test data with the event-to-detection mapping omitted 
test.solution -- the sample solution on the test data

Overview
========

The model is completely described in ```description.odt``` and this should be
translated in the Probabilistic Programming Language of your
choosing. The unlabeled data in ```test.blind``` (and optionally the labeled
data in ```training.data```) comprises the observations to the model. The
query of interest is the seismic bulletin for each of the observed
episodes.

Once all the bulletins have been produced, the script ```evaluate.py```
can be used to produce the reports on the accuracy versus the reference
script ```test.data```. One can also compare the results versus the
baseline in ```test.solution```.

The files ```generate.py``` and ```solve.py``` have only been provided
for convenience they shouldn't normally be used. However, if you want to
check the performance of your model on more than just the provided data
you may generate more as needed. The sample solver is based loosely on
the published greedy algorithm, and may be used as a competitive
baseline.

Authors
=======
Nimar S. Arora, Bayesian Logic Inc., nimar.arora@gmail.com
Stuart Russell, Deptt. of Computer Science, Berkeley.
