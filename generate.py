#!/usr/bin/python
"""
Generate the physics of a hypothetical 2-D spherical world, and then generate
seismic events and detections.

The data for this problem was generated using the following command:

./generate.py 10000 physics.data training.data test.data test.blind

Author: Nimar Arora (feel free to modify, use, or resdistribute)
"""
from __future__ import print_function, division

## IMPORTS

import sys

from scipy.stats import gamma, norm, invgamma, poisson, uniform, expon,\
     bernoulli, laplace
from numpy import array, sqrt, log, exp, pi, arcsin, degrees
from collections import namedtuple

from util import *

def main():
  
  # first we process the command line arguments
  
  if len(sys.argv) != 6:
    print("Usage: generate.py <num episodes> <physics-file> <training-file>"
          " <test-file> <blind-test-file>", file=sys.stderr)
    sys.exit(1)

  numepisodes, phys_file_name, train_file_name, test_file_name, blind_file_name\
               = sys.argv[1:]
  numepisodes = int(numepisodes)

  # next we generate the physics of the world
  
  physics = sample_physics()

  # and then we generate the test and training data based on the physics
  
  training = sample_episodes(numepisodes, physics)
  
  test = sample_episodes(numepisodes, physics)

  # the blind test data has no events
  test_blind = strip_events(test)

  
  # finally we write out the physics and all the data files
  
  write_namedtuple(physics, phys_file_name)

  write_episodes(training, train_file_name)

  write_episodes(test, test_file_name)

  write_episodes(test_blind, blind_file_name)
  
  
def sample_physics():
  """
  Return a physics object with all the model parameters sampled from their
  appropriate hyperpriors.
  """

  T = 3600                              # an episode is 3600 seconds

  R = 6371                             # the radius of the earth in km
  
  # event rate
  lambda_e = gamma.rvs(6.0, loc=0, scale = 1/(4 * pi * R**2 * T))
  
  # event magnitude
  lambda_m = log(10)
  
  # station-specific attributes
  
  mu_d0, mu_d1, mu_d2 = [], [], []
  mu_t, theta_t, mu_z, theta_z, mu_s, theta_s = [], [], [], [], [], []
  mu_a0, mu_a1, mu_a2, sigma_a = [], [], [], []
  lambda_f, mu_f, sigma_f = [], [], []
  
  for sta in STATIONS:
    
    # per-station detection probability
    (d0, d1, d2) = mvar_norm_sample(array([-10.4, 3.26, -0.0499]),
                                    array([[13.43,-2.36, -.0122],
                                           [-2.36, 0.452, .000112],
                                           [-.0122, .000112, .000125]]))
    mu_d0.append(d0)
    mu_d1.append(d1)
    mu_d2.append(d2)

    # time
    mu_t.append(0.0)
    theta_t.append(invgamma.rvs(120, 0.0, 118))
    
    # azimuth
    mu_z.append(0.0)
    theta_z.append(invgamma.rvs(5.2, 0.0, 44.0))

    # slowness
    mu_s.append(0.0)
    theta_s.append(invgamma.rvs(6.7, 0.0, 7.5))

    # amplitude
    (a0, a1, a2) = mvar_norm_sample(array([-7.3, 2.03, -.00196]),
                                    array([[1.23, -.227, -.000175],
                                           [-.227, .0461, .0000245],
                                           [-.000175, .0000245,  .000000302]]))
    mu_a0.append(a0)
    mu_a1.append(a1)
    mu_a2.append(a2)
    
    sigma_a.append(sqrt(invgamma.rvs(21.1, 0, 12.6)))
    
    # false arrivals
    lambda_f.append(gamma.rvs(2.1, 0, 0.0013))
    mu_f.append(norm.rvs(-0.6, 0.6))
    sigma_f.append(sqrt(invgamma.rvs(10.34, 0, 9.49)))
    
  return Physics(T, R, lambda_e, lambda_m, mu_d0, mu_d1, mu_d2,
                 mu_t, theta_t, mu_z, theta_z, mu_s, theta_s,
                 mu_a0, mu_a1, mu_a2, sigma_a,
                 lambda_f, mu_f, sigma_f)

def sample_episodes(numepisodes, physics):

  episodes = []
  
  for epinum in xrange(numepisodes):
    
    events = []
    detections = []
    assocs = []

    # first generate all the events

    numevents = poisson.rvs( physics.lambda_e * 4 * pi * physics.R ** 2
                             * physics.T )
    
    for evnum in xrange(numevents):
      
      evlon = uniform.rvs(-180, 360)
      evlat = degrees(arcsin(uniform.rvs(-1, 2)))
      evmag = expon.rvs(physics.lambda_m, 2.0)
      evtime = uniform.rvs(0, physics.T)

      event = Event(evlon, evlat, evmag, evtime)
      
      events.append(event)

      truedets = []

      #print ("event mag %f" % event.mag)
      
      # for each event generate its set of true detections
      for stanum, station in enumerate(STATIONS):

        dist = compute_distance((station.lon, station.lat),
                                (event.lon, event.lat))
        sta_to_ev_az = compute_azimuth((station.lon, station.lat),
                                       (event.lon, event.lat))

        detprob = logistic(physics.mu_d0[stanum]
                           + physics.mu_d1[stanum] * event.mag
                           + physics.mu_d2[stanum] * dist)

        #print ("stanum %d dist %f detprob %f" % (stanum, dist, detprob))
        
        # is detected ?
        if bernoulli.rvs(detprob):
          dettime = laplace.rvs(event.time + compute_travel_time(dist)
                                + physics.mu_t[stanum],
                                physics.theta_t[stanum])

          degdiff = laplace.rvs(physics.mu_z[stanum], physics.theta_z[stanum])
          detaz = (sta_to_ev_az + degdiff + 360) % 360

          detslow = laplace.rvs(compute_slowness(dist) + physics.mu_s[stanum],
                                physics.theta_s[stanum])

          detamp = exp(norm.rvs(physics.mu_a0[stanum]
                                + physics.mu_a1[stanum] * event.mag
                                + physics.mu_a2[stanum] * dist,
                                physics.sigma_a[stanum]))

          truedets.append(len(detections))
          detections.append(Detection(stanum, dettime, detaz, detslow, detamp))

      assocs.append(truedets)


    # now generate the false detections
    for stanum in xrange(len(STATIONS)):
      numfalse = poisson.rvs(physics.lambda_f[stanum] * physics.T)
      
      for dnum in xrange(numfalse):
        dettime = uniform.rvs(0, physics.T)
        detaz = uniform.rvs(0, 360)
        detslow = uniform.rvs(compute_slowness(180),
                              compute_slowness(0) - compute_slowness(180))
        detamp = exp(norm.rvs(physics.mu_f[stanum], physics.sigma_f[stanum]))
        
        detections.append(Detection(stanum, dettime, detaz, detslow, detamp))
    
    episodes.append(Episode(events, detections, assocs))
    
  return episodes

def strip_events(episodes):
  """
  Returns episodes with identical detections but no events
  """
  
  return [Episode([], episode.detections, []) for episode in episodes]


if __name__ == "__main__":
  main()
  
