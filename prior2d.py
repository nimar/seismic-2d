# Learn the hyper-priors for the seismic-2d model using real IDC data.

from __future__ import division

import pickle
import numpy as np
from scipy.stats import norm, expon, gamma, invgamma, laplace, cauchy, mode

from idcnetvisa.database.dataset import *
from idcnetvisa import netvisa
from idcnetvisa.utils.LogisticModel import LogisticModel
from idcnetvisa.utils import geog
from idcnetvisa.priors import SecDetPrior

DATAFILE="/home/nimar/netvisa/sandbox/models/model-1314835200-1325376000/Data.pic"

def iasp_ttime(delta): return -.023 * delta ** 2 + 10.7 * delta + 5
def iasp_slow(delta): return -.046 * delta  + 10.7

T = 3600

def main():
  data = pickle.load(file(DATAFILE))
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  earthmodel = netvisa.EarthModel(sitenames, sites, phasenames, phasetimedef,
                                  phaseprop, start_time,
                                  ttime_prefix, ddrange_file,
                                  qfvc_file, hydro_dir, infra_dir)

  print "data file unpickled"

  learn_ttime_approx(earthmodel)
  
  topstations = learn_stations(data, earthmodel)

  learn_evrate(data)
  
  learn_detprob(data, earthmodel, topstations)

  learn_arraz(data, earthmodel, topstations)
  
  learn_arrtime(data, earthmodel, topstations)

  learn_arrslow(data, earthmodel, topstations)
  
  learn_arramp(data, earthmodel, topstations)

  false_dets = SecDetPrior.compute_false_detections(detections, leb_evlist)
  
  learn_falserate(data, earthmodel, false_dets)
  
  learn_falseamp(data, earthmodel, false_dets)

def learn_ttime_approx(earthmodel):
  print "Learning Travel Time Approximation"
  
  X, Y = [], []
  
  for dist in range(180):
    try:
      ttime = earthmodel.TravelTime(0, 0, dist)
    except:
      continue
    
    X.append(dist)
    Y.append(ttime)
  
  X = np.array(X)
  Y = np.array(Y)
  
  p = np.polyfit(X, Y, 2)
  print p

def learn_stations(data, earthmodel):

  print "Top Stations:"
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()
  truecnt = [0 for s in xrange(numsites)]

  phaseid = 0
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10 \
           or event[EV_MB_COL] < 2:
      continue

    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:
        
        siteid = int(detections[detnum, DET_SITE_COL])
        
        truecnt[siteid] += 1
        
  for siteid in xrange(numsites):
    if truecnt[siteid] > 1000:
      print "siteid", siteid, sitenames[siteid]
      print sites[siteid]
      
  return [siteid for siteid in xrange(numsites) if truecnt[siteid] > 1000]

def learn_falseamp(data, earthmodel, false_dets):
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data
  
  numsites = earthmodel.NumSites()
  truecnt = [0 for s in xrange(numsites)]
  site_raw = [[] for s in xrange(numsites)]
  
  for detnum in false_dets:
    siteid = int(detections[detnum, DET_SITE_COL])
    amp = detections[detnum, DET_AMP_COL]
    if amp > 0:
      site_raw[siteid].append(np.log(amp))
      
  phaseid = 0
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10 \
           or event[EV_MB_COL] < 2:
      continue
    
    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:

        siteid = int(detections[detnum, DET_SITE_COL])

        truecnt[siteid] += 1

  print "False Amp"

  all_loc, all_scale = [], []
  
  for siteid, raw in enumerate(site_raw):
    if truecnt[siteid] > 1000:
      print "siteid", siteid,

      loc, scale = cauchy.fit(raw)

      print "siteid", siteid, "Cauchy", loc, scale

      all_loc.append(loc)
      all_scale.append(scale)

  print "Gaussian of scale:", norm.fit(all_loc)
  print "InvGamma of scale:", invgamma_fit(all_scale)

def learn_falserate(data, earthmodel, false_dets):

  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()

  falsecnt = [0 for s in xrange(numsites)]
  truecnt = [0 for s in xrange(numsites)]
  
  for detnum in false_dets:
    siteid = int(detections[detnum, DET_SITE_COL])
    falsecnt[siteid] += 1

  phaseid = 0
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10 \
           or event[EV_MB_COL] < 2:
      continue

    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:

        siteid = int(detections[detnum, DET_SITE_COL])

        truecnt[siteid] += 1

  print "False Rate"
  
  all_false = []
  for siteid in xrange(numsites):
    if truecnt[siteid] > 1000:
      rate = falsecnt[siteid] / (end_time - start_time)
      all_false.append(rate)
      print "siteid", siteid,
      print "false rate", rate

  print "overall false rate gamma:", gamma.fit(all_false)
  print "overall false rate gamma (mom):", gamma_fit(all_false)

def learn_arramp(data, earthmodel, topstations):
  print "Arrival Amplitude:"
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()

  phaseid = 0

  site_raw = [[] for site in xrange(numsites)]
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10 \
           or event[EV_MB_COL] < 2:
      continue

    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:

        siteid = int(detections[detnum, DET_SITE_COL])

        if earthmodel.SiteMedium(siteid) != MEDIUM_SEISMIC:
          continue

        (arrtime, arraz, dist, oop_angle, arr_henergy) \
                  = earthmodel.ArrivalParams(event[EV_LON_COL],
                                             event[EV_LAT_COL],
                                             event[EV_DEPTH_COL],
                                             event[EV_TIME_COL],
                                             event[EV_HYDRO_ENERGY_COL],
                                             phaseid, siteid)
        
        if arrtime > 0 and detections[detnum, DET_AMP_COL] > 0:

          logamp = np.log(detections[detnum, DET_AMP_COL])
          evmag = event[EV_MB_COL]
          ttime = arrtime - event[EV_TIME_COL]
          
          site_raw[siteid].append((logamp, evmag, ttime))

  all_coeff, all_std = [], []
  
  for siteid in topstations:

    raw = np.array(site_raw[siteid])

    # solve ax = b by least squares
    b = np.array(raw[:,0])
    a = np.c_[np.ones(len(b)), raw[:,1:3]]
    x = np.linalg.lstsq(a, b)[0]
    res = b - np.dot(a, x)
    res_std = np.std(res)

    coeff = np.array(x)
    all_coeff.append(coeff)
    all_std.append(res_std)

    print "siteid:", siteid,
    print "coeff:", coeff,
    print "std:", res_std

  mu, Sigma = mvar_norm_fit(all_coeff)
  print "coeff mean:", mu
  print "coeff covariance:", Sigma
  print "residual variance: Inverse-Gamma",invgamma_fit([s**2 for s in all_std])

def learn_arrslow(data, earthmodel, topstations):
  print "Arrival Slowness:"
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()

  phaseid = 0

  site_raw = [[] for site in xrange(numsites)]
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10:
      continue

    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:

        siteid = int(detections[detnum, DET_SITE_COL])

        if earthmodel.SiteMedium(siteid) != MEDIUM_SEISMIC:
          continue

        (arrtime, arraz, dist, oop_angle, arr_henergy) \
                  = earthmodel.ArrivalParams(event[EV_LON_COL],
                                             event[EV_LAT_COL],
                                             event[EV_DEPTH_COL],
                                             event[EV_TIME_COL],
                                             event[EV_HYDRO_ENERGY_COL],
                                             phaseid, siteid)
        
        if arrtime > 0:
          
          arrslow = earthmodel.ArrivalSlowness(event[EV_LON_COL],
                                               event[EV_LAT_COL],
                                               event[EV_DEPTH_COL],
                                               phaseid, siteid)
          
          res = arrslow - detections[detnum, DET_SLO_COL]
          
          site_raw[siteid].append(res)

  locs, scales = [], []
  for siteid in topstations:
    raw = site_raw[siteid]
    
    loc, scale = laplace.fit(raw)
    locs.append(loc)
    scales.append(scale)
    print siteid, loc, scale
      
  print "Inverse-Gamma:", invgamma_fit(scales)
      

def learn_arraz(data, earthmodel, topstations):
  print "Arrival Azimuth"
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()

  phaseid = 0

  site_raw = [[] for site in xrange(numsites)]
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10:
      continue

    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:

        siteid = int(detections[detnum, DET_SITE_COL])

        if earthmodel.SiteMedium(siteid) != MEDIUM_SEISMIC:
          continue

        (arrtime, arraz, dist, oop_angle, arr_henergy) \
                  = earthmodel.ArrivalParams(event[EV_LON_COL],
                                             event[EV_LAT_COL],
                                             event[EV_DEPTH_COL],
                                             event[EV_TIME_COL],
                                             event[EV_HYDRO_ENERGY_COL],
                                             phaseid, siteid)

        if arrtime > 0:
          res = geog.degdiff(arraz, detections[detnum, DET_AZI_COL])

          site_raw[siteid].append(res)

  locs, scales = [], []
  for siteid in topstations:
    
    raw = site_raw[siteid]
    
    loc, scale = laplace.fit(raw)
    locs.append(loc)
    scales.append(scale)
    print siteid, loc, scale
      
  print "Inverse-Gamma:", invgamma_fit(scales)
      

def learn_arrtime(data, earthmodel, topstations):
  print "Arrival Time:"
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()

  phaseid = 0

  site_raw = [[] for site in xrange(numsites)]
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10:
      continue

    for ph, detnum in leb_evlist[evnum]:
      
      if ph==phaseid:

        siteid = int(detections[detnum, DET_SITE_COL])

        if earthmodel.SiteMedium(siteid) != MEDIUM_SEISMIC:
          continue

        (arrtime, arraz, dist, oop_angle, arr_henergy) \
                  = earthmodel.ArrivalParams(event[EV_LON_COL],
                                             event[EV_LAT_COL],
                                             event[EV_DEPTH_COL],
                                             event[EV_TIME_COL],
                                             event[EV_HYDRO_ENERGY_COL],
                                             phaseid, siteid)

        if arrtime > 0:
          res = detections[detnum, DET_TIME_COL] - arrtime

          site_raw[siteid].append(res)

  locs, scales = [], []
  for siteid in topstations:
    raw = site_raw[siteid]
    
    loc, scale = laplace.fit(raw)
    locs.append(loc)
    scales.append(scale)
      
    print siteid, loc, scale

  print "Inverse-Gamma:", invgamma_fit(scales)
      
def learn_evrate(data):
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  NUMINT = (end_time - start_time) // T

  for minmag in [2, 2.5, 3, 3.5, 4]:
    print "minimum magnitude", minmag, ":"
    
    # count the number of events in each interval
    intcounts = dict((i,0) for i in xrange(NUMINT))

    mags = []
    
    for evnum in xrange(len(leb_events)):
      
      evmag = leb_events[evnum, EV_MB_COL]

      if evmag < minmag:
        continue

      mags.append(evmag)
      
      evtime = leb_events[evnum, EV_TIME_COL]
  
      if evtime < end_time - T:
      
        intcounts[ (evtime - start_time) // T ] += 1
  
    counts = intcounts.values()
  
    print "  mean event rate", np.mean(counts)
    print "  gamma loc=0 dist", gamma_fit(counts)
    print "  mean magnitude", np.mean(mags)

def learn_detprob(data, earthmodel, topstations):
  print "Detection Probability:"
  
  (start_time, end_time, detections, leb_events, leb_evlist,
   site_up, sites, phasenames, phasetimedef, phaseprop,
   sitenames, ttime_prefix, ddrange_file, qfvc_file,
   hydro_dir, infra_dir) = data

  numsites = earthmodel.NumSites()

  phaseid = 0

  site_raw = [[] for site in xrange(numsites)]
  
  for evnum, event in enumerate(leb_events):

    if event[EV_MEDIUM_COL] != MEDIUM_SEISMIC or event[EV_DEPTH_COL] > 10:
      continue

    det_site = set(detections[detnum, DET_SITE_COL]\
                   for ph, detnum in leb_evlist[evnum] if ph==phaseid)


    for siteid in range(numsites):
      if earthmodel.SiteMedium(siteid) != MEDIUM_SEISMIC:
        continue
      
      (arrtime, arraz, dist, oop_angle, arr_henergy) \
                = earthmodel.ArrivalParams(event[EV_LON_COL], event[EV_LAT_COL],
                                           event[EV_DEPTH_COL],
                                           event[EV_TIME_COL],
                                           event[EV_HYDRO_ENERGY_COL],
                                           phaseid, siteid)
      
      oop_angle = abs(oop_angle)
      
      # check if the site is in the shadow zone of this phase (or the
      # phase is blocked)
      if arrtime < 0:
        continue
      
      # check if the site was up at the expected arrival time
      if arrtime < start_time or arrtime >= end_time \
          or not site_up[siteid, (arrtime - start_time) // UPTIME_QUANT]:
        continue
      
      isdet = int(siteid in det_site)
      
      site_raw[siteid].append((isdet, event[EV_MB_COL], dist))
    

  all_coeff = []
  
  for siteid in topstations:

    raw = site_raw[siteid]
    
    numevents = sum(x[0] for x in raw)
    if numevents >= 10:
      print siteid

      model = LogisticModel("isdet", ["mb", "dist"],
                            [[x[1] for x in raw], [x[2] for x in raw]],
                            [x[0] for x in raw])

      coeff = np.array([model.coeffs[-1], model.coeffs[0], model.coeffs[1]])
      
      print "siteid %d coeffs:" % siteid, coeff

      all_coeff.append(coeff)

  mu, Sigma = mvar_norm_fit(all_coeff)
  print "coeff mean:", mu
  print "coeff covariance:", Sigma

def invgamma_fit(data):
  """
  match mean and var assuming location=0
  """
  mu = np.mean(data)
  var = np.var(data)

  alpha = (1/var) * mu ** 2 + 2

  beta = mu * (alpha - 1)

  return alpha, beta

def gamma_fit(data):
  """
  match mean and var assuming location=0
  """
  mu = np.mean(data)
  var = np.var(data)

  scale = var / mu
  shape = mu / scale

  return shape, scale

def mvar_norm_fit(data):
  """
  Fits a multivariate normal using row vectors as data.

  Returns a mean vector (also a row vector) and a covariance matrix.
  """
  assert(len(data) > 1)
  
  mu = sum(data) / len(data)

  # data should be a row vectors
  assert(len(mu.shape)==1)

  # ML estimator of covariance

  Sigma = sum( np.outer(x - mu, x-mu) for x in data ) / len(data)

  return mu, Sigma

def mvar_norm_rvs(mu, Sigma):
  """
  Generates a random variable from a multivariate normal distribution with
  given mean vector and covariance matrix.
  
  Returns a vector of the same shape as the mean vector
  
  see http://en.wikipedia.org/wiki/Multivariate_normal_distribution  
  """
  # mu better be a row vector
  assert(len(mu.shape)==1)

  # the square root (roughly speaking) of Sigma
  # i.e. root * root.T = Sigma
  root = np.linalg.cholesky(Sigma)

  z = np.array([norm.rvs() for _ in mu])
  
  return mu + np.dot(root, z)

if __name__ == "__main__":
  main()
