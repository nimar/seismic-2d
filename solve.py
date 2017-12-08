#!/usr/bin/python
"""
Solves the seismic-2d episodes using a very simple algorithm.

./solve.py data/training.data data/test.blind data/test.solution

This script first learns an approximation of the physics, and then uses it to
associate arrivals to candidate events.

Author: Nimar Arora (feel free to modify, use, or redistribute)
"""
from __future__ import print_function, division

import sys

from scipy.stats import laplace
from scipy.stats import gamma, norm, invgamma, poisson, uniform, expon,\
     bernoulli, laplace, cauchy
from numpy import array, sqrt, log, exp, pi, arcsin, degrees, linspace, seterr,\
     inf
seterr(all='raise')

from util import *
from generate import sample_physics

try:
  from csolve import *
except ImportError:
  print("Warning: csolve module not found. Please build it with:\n"
        "  python setup.py build_ext --inplace\n"
        "Reverting to python-based solver module, which is very slow!!")
  from pysolve import *

def main():
  if len(sys.argv) != 4:
    print("Usage: solve.py training.data test.blind test.solution",
          file = sys.stderr)
    sys.exit(1)
    
  train_file_name, blind_file_name, solution_file_name = sys.argv[1:]
  
  traindata = read_episodes(train_file_name)
  
  testdata = read_episodes(blind_file_name)

  physics = learn_physics(traindata)
  
  # We will write out each episode as we solve it, so that it can
  # be concurrently analyzed.
  fp = open(solution_file_name, "w")
  
  print("Episodes:", file=fp, end="\n\n")
  fp.flush()
  
  for episode in testdata:
    
    events, assocs = solve_episode(physics, STATIONS, episode.detections)
    
    write_single_episode(Episode(events, [], assocs), fp)
    
    fp.flush()
    
  fp.close()

def learn_physics(traindata):
  # Sample a physics object to get the values for the constants. All the other
  # attributes of this object will be learned.
  physics = sample_physics()

  # learn the event rate
  avgevents = sum(len(epi.events) for epi in traindata) / len(traindata)
  
  lambda_e = avgevents / (4 * pi * physics.R**2 * physics.T)

  # learn the detection probability coefficients

  physics = physics._replace(lambda_e = lambda_e)

  for field in physics._fields:
    print("{0} = {1}".format(field, getattr(physics, field)))

  import sys
  sys.exit(1)

  return physics

if __name__ == "__main__":
  try:
    main()
  except SystemExit:
    raise
  except:
    import pdb, traceback, sys
    traceback.print_exc(file=sys.stdout)
    pdb.post_mortem(sys.exc_traceback)
    raise
