'''
Created on 2010-06-17

@author: adam
'''

from qfault.noise import Bound
from qfault.util import concurrency
from qfault.util.iteration import PartitionIterator
import functools
import gmpy
import logging
import operator


__all__ = ['counts_as_poly', 
           'summed_counts_as_poly',
           'pr_at_least_k_failures']

logger = logging.getLogger('counting.probability')

def counts_as_poly(counts, locations, noise, bound=Bound.UpperBound):
    '''
    Convert weighted error likelihood counts into a polynomial that
    gives the (possibly un-normalized) failure probability as a function
    of the original error strength.
    '''
    coeffs = [sum(countsK.values()) for countsK in counts]
    return summed_counts_as_poly(coeffs, locations, noise, bound)

def summed_counts_as_poly(summed_counts, locations, noise, bound=Bound.UpperBound):
    '''
    Returns a polynomial for the given (summed) weighted likelihood counts.
    
    :param summed_counts: A list of weighted-counts, indexed by fault-order.
    :param locations: The locations to which the counts correspond.
    :param noise: The noise model.
    :param bound: (optional) The type of bound, either UpperBound or LowerBound.
    '''
    prefactor = _likelihood_prefactor(locations, noise, bound)
    likelyhood = noise.likelyhood(bound)
    coeff_sum = sum(summed_counts[k] * (likelyhood ** k) 
                    for k in range(len(summed_counts)))
    
    return prefactor * coeff_sum

def pr_k_of_n_ub(k, n, pr):
    '''
    Returns an upper bound on the probability that some event (e.g. failure)
    occurs for at least k out of n objects.  Assumes each object is independent.
    
    :param k: The number of objects for which the event occurs
    :param n: The total number of objects
    :param pr: The probability that the event occurs for a single object.
    
    >>> pr_k_of_n_ub(2, 4, 0.25)
    0.375
    '''
    B = long(gmpy.comb(n, k))
    p = B * pow(pr, k)
    return p

def pr_failure_partition(totals, failure_probs, partition):
    '''
    Returns an upper bound on the probability that the given
    failure partition occurs.  The partition is defined by
    a list integers.  Each integer specifies the number of
    failures of a category of locations.
    
    :param totals: The total number of locations for each category
    :param failure_probs: The failure probability of a location, by category.
    :param partition: A list of integers, one for each category, that specifies
                      the number of failures in each category.
                      
      >>> pr_failure_partition([4, 5], [0.1, 0.2], [1, 2])
      0.16000000000000003
    '''
    return reduce(operator.mul, [pr_k_of_n_ub(*args)
                                 for args in zip(partition, totals, failure_probs)])

def _location_totals_by_type(locations):
    '''
    Returns a tuple of length two.  The first
    element is a list of representative locations of
    each type.  The second element is a corresponding
    tally of all of the locations of each type.
    '''
    
    # Sort locations by type
    locs_by_type = {}
    for loc in locations:
        ltype = loc['type']
        locs_by_type[ltype] = locs_by_type.get(ltype, []) + [loc]
        
    # Take a representative of each location type.  Assume
    # that locations of the same type have the same failure
    # statistics.
    loc_types = [locs_of_type[0] for locs_of_type in locs_by_type.values()]
    loc_totals = map(len, locs_by_type.values())
    
    return loc_types, loc_totals
    

def pr_at_least_k_failures(kMin, locations, noiseModel, kMax=None):
    '''
    Returns Pr[kMin <= k <= kMax], an upper bound on the probability that  between kMin and kMax
    failures (of any kind) occur at the given locations.
    If kMax=None, then returns Pr[kMin <= k] an upper bound on the probability that at least kMin
    failures occur.
    
    Arguments
    ---------
    kMin        -- The minimum number of failures.
    loc_totals    -- (LocationCount) The total number of each type of location.
    noiseModel    -- The noise model.
    kMax        -- (optional) The maximum number of failures.
    
    >>> import qfault.circuit.location as location
    >>> from qfault.noise import NoiseModelXZSympy
    >>> from qfault.qec.error import Pauli
    >>> locations = location.Locations([location.prep(Pauli.Z, 'test', 0),
    ...                                 location.prep(Pauli.X, 'test', 1),
    ...                                 location.cnot('test', 1, 'test', 0),
    ...                                 location.meas(Pauli.Z, 'test', 1)])
    >>> noise_model = NoiseModelXZSympy()
    >>> pr = pr_at_least_k_failures(1, locations, noise_model, kMax=2)
    >>> str(pr)
    '(-15*x + 1)*(-4*x + 1)**3*(228*x**2/(-15*x + 1)**2 + 27*x/(-15*x + 1))'
    >>> pr(0.01)
    0.262610462117647
    >>> pr = pr_at_least_k_failures(1, locations, noise_model)
    >>> str(pr)
    '(-15*x + 1)*(-4*x + 1)**3*(960*x**4/(-15*x + 1)**4 + 784*x**3/(-15*x + 1)**3 + 228*x**2/(-15*x + 1)**2 + 27*x/(-15*x + 1))'
    >>> pr(0.01)
    0.263584338015876
    '''
    
    #===============================================================================
    # The calculation performed here is of the form:
    #
    # A_\vec{n} sum_\vec{k} L^k \prod_i B_i Pr[fail_i]
    #
    # A_n            -- A prefactor based on the number of each type of location.
    #                   (e.g. (1-12g)^10 (1-8g)^5 (1-4g)^8)
    # L                -- The likelyhood term from the noise model (e.g. g/(1-12g))
    # B_i            -- The binomial coefficient \choose{n_i}{k_i}
    # Pr[fail_i]    -- The probability of failure for location type i. (e.g. 8g)
    #
    # The sum is over all kMin <= k <= kMax.  The product is over all (six) location
    # types. Each product term is computed in parallel.
    #===============================================================================

    boundType = Bound.UpperBound
    
    if None == kMax:
        # 10 is arbitrary, but seems to work well.
        kMax = kMin + 10
        bound = True
    else:
        bound = False
        
    kMax = min(kMax, len(locations)) + 1
    
    loc_types, loc_totals = _location_totals_by_type(locations)
    
    logger.debug('Computing Pr[{0} <= k < {1}] for {2}, {3}'.format(kMin, 
                                                                    kMax, 
                                                                    loc_totals, 
                                                                    noiseModel))
    
    loc_type_weight = lambda loc_type: sum(noiseModel.getWeight(loc_type, e, boundType) 
                                           for e in noiseModel.errorList(loc_type))
    loc_weights = map(loc_type_weight, loc_types)

    # Check for identical weights.  These can be grouped together which
    # will significantly reduce the number of partition iterations.
    weight_set = set(loc_weights)
    new_loc_weights = []
    new_loc_totals = []
    for w in weight_set:
        n = sum(loc_totals[i] for i in range(len(loc_weights)) if loc_weights[i] == w)
        new_loc_weights.append(w)
        new_loc_totals.append(n)
        
    loc_totals = new_loc_totals
    loc_weights = new_loc_weights
                    
    likelyhood = noiseModel.likelyhood(boundType)
    
    # For each k, compute the failure probability for each partition
    # of k faults.  There may be many partitions, so we compute
    # them in parallel.
    pr = 0
    for k in range(kMin, kMax):
        partitioner = PartitionIterator(k, len(loc_totals), loc_totals)
        map_func = functools.partial(pr_failure_partition, loc_totals, loc_weights)
        weight = concurrency.mapreduce_concurrent(map_func, sum, partitioner)
        pr += weight * (likelyhood ** k)

    prefactor = _likelihood_prefactor(locations, noiseModel, bound)
    logger.debug('A=%s, pr=%s', prefactor, pr)
    pr *= prefactor
        
    if bound and (len(locations) >= kMax):
        # Not all possible failure configurations were computed.  Bound
        # the probability by ignoring the prefactor and using probabilities (instead of
        # likelyhoods)
        
        iterator = PartitionIterator(kMax, len(loc_totals), loc_totals) 
        prFailList = [noiseModel.prFail(l, boundType) for l in loc_types]
        map_func = functools.partial(pr_failure_partition, loc_totals, prFailList)
        cap = concurrency.mapreduce_concurrent(map_func, sum, iterator)
        logger.debug('adding bounding cap: %s', cap)
        pr += cap
    
    return pr

#def pr_bad(kGood, locations, noiseModel, kMax=None):
#    # Count up all of the locations
#    
#    prBad = pr_at_least_k_failures(kGood+1, locations, noiseModel, kMax)
#    if int == type(prBad):
#        prBad = SymPolyWrapper(sympoly1d([prBad]))
#    
#    return prBad

def _likelihood_prefactor(locations, noise_model, bound):
    '''
    Returns the prefactor for calculating probabilities based
    on likelihoods.
    Should look something like (1-12g)^nCnot * (1-8g)^nRest * ...
    '''
    loc_types, loc_totals = _location_totals_by_type(locations)
    prefactor = reduce(operator.mul, 
                       (noise_model.prIdeal(l, bound) ** total 
                        for l, total in zip(loc_types, loc_totals)), 
                       1)
    return prefactor
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    