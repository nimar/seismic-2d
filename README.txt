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

description.tex -- a description of the model
generate.py   -- generates the physics of a 2-D world and some episodes
util.py       -- some geophysical utility functions
solve.py      -- a sample solver that learns the physics and solves the episodes
pysolve.py    -- a python-based solver
csolve.c      -- a C-based solver, identical to pysolve.py but much faster
evaluate.py   -- evaluates a solution versus a reference
mwmatching.py -- utility script for max-weight max cardinality matching

short_data/ or large_data/
                  -- data directory where the following files
                     need to be downloaded
 physics.data     -- physics for the training and test episodes
 training.data    -- 100 (or 10K) episodes for training
 test.data        -- another 100 (or 10K) episodes for testing
 test.blind       -- the test data with the event-to-detection mapping omitted
 test.solution    -- a sample solution on the test data
 
The data files needed for the PPAML CP4 challenge problem are located here:
http://ppaml.kitware.com/midas/item/4543 

Overview
========

The model is completely described in ```description.tex``` and this
should be translated in the Probabilistic Programming Language of your
choosing. The unlabeled test data in ```test.blind``` and the labeled
training data in ```training.data``` comprises the observations to the
model. The query of interest is the seismic bulletin for each of the
observed episodes in the test data.

Once all the bulletins have been produced, the script ```evaluate.py```
can be used to produce the reports on the accuracy versus the reference
script ```test.data```. One can also compare the results versus the
baseline in ```test.solution```.

The files ```generate.py``` and ```solve.py``` have only been provided
for convenience. These files shouldn't normally be used. However, if you
want to check the performance of your model on more than just the
provided data you may generate more as needed. The sample solver is
based loosely on the published greedy algorithm, and may be used as a
very simple baseline.

The sample solver first tries to run the C-based solver and if that is
not found it attempts the python version. In order to build the C-based
solver the following command must be executed:

    python setup.py build_ext --inplace

This will generate a file csolve.so (or csolve.pyd on
Windows). Compiling Python C extensions on Windows usually involves
additional steps. We strongly recommend using a 32-bit Python even on
64-bit Windows, as the installation is a lot easier.

- Compiling on 32-bit Python using MinGW:

    o.  Install MinGW

    o.  Compile using the --compiler flag:

        python setup.py build_ext --inplace --compiler=mingw32
    
- Compiling on 64-bit Python:

    o. Install Visual Studio 10.0 Express:
        http://www.microsoft.com/visualstudio/eng/downloads#d-2010-express

    o. Install Microsoft Windows SDK for Windows 7 and .NET Framework 4 
        http://www.microsoft.com/en-us/download/details.aspx?id=8279
 
    o. Create the following file:
    
  C:\Program Files (x86)\Microsoft Visual Studio 10.0\vc\bin\amd64\vcvars64.bat

       with the single-line content:
       
  CALL "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /x64

    o. Compile as usual

        python setup.py build_ext --inplace

Authors
=======
Nimar S. Arora, Bayesian Logic Inc., nimar.arora@gmail.com
Stuart Russell, Deptt. of Computer Science, Berkeley.
