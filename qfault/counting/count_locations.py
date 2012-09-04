'''

Functions counting error propagation within a given set of (physical) locations. 

@author: adam
'''

from qfault.qec.error import Pauli
from qfault.util import listutils, concurrency
import logging
import itertools
from copy import copy
from qfault import noise
import functools

__all__ = ['count_errors_of_order_k', 'count_location_set', 'propagate_location_errors']
LOGGER = logging.getLogger('Counting')


def propagate_location_errors(locations):
    '''
    Propagates all bit errors for all locations.  Returns a list of resulting
    errors, one for each location.  Each propagated error is indexed by the
    Pauli error (i.e., X, Y, Z) from which it originated. Locations 
    containing two qubits (i.e., CNOT) are indexed by two-qubit Pauli group 
    elements.
    
    >>> import qfault.circuit.location as location
    >>> cnot = location.cnot('test', 0, 'test', 1)
    >>> measX = location.meas(Pauli.X, 'test', 1)
    >>> propagate_location_errors(location.Locations([cnot, measX]))
    [{IZ: {'test': IZ}, ZX: {'test': ZI}, YX: {'test': YI}, ZY: {'test': ZZ}, YY: {'test': YZ}, XX: {'test': XI}, XY: {'test': XZ}, XZ: {'test': XZ}, ZI: {'test': ZI}, XI: {'test': XI}, YZ: {'test': YZ}, YI: {'test': YI}, IX: {'test': II}, ZZ: {'test': ZZ}, IY: {'test': IZ}}, {Z: {'test': IZ}}]
    '''
    # First extract a list of all the block IDs used in the computation
    blocklengths = locations.blocklengths()
    #if len(blocklengths) != 14:    # used to have only 10 blocklengths, before adding data snapshots
    #    print "Error!  Too few blocklengths!"
    
    # now for each location, consider the types of errors that can occur there
    #     for a 1-qubit location, it is sufficient to consider X and Z errors separately
    #     for a 2-qubit location, it is sufficient to consider XI, IX, ZI, IZ
    # for each type of error, propagate it forward, and store a dictionary in the location
    def initialize_errors(blocklengths):
        """This returns a properly formatted errors array, initialized to 0."""
        return {block_id: Pauli.I ** length
                for block_id, length in blocklengths.iteritems()}
    
    def propagate_through_remaining(errors, locations):
        for loc in locations:
            propagate_errors_through_loc(errors, loc)
            
    def error_parts(error):
        xpart_list = ['X'+str(i+1) for i in error.partial(Pauli.X).support()]
        zpart_list = ['Z'+str(i+1) for i in error.partial(Pauli.Z).support()]
        return xpart_list + zpart_list
    
    def error_product(loc_errors, err):
        eparts = [loc_errors[e_part] for e_part in error_parts(err)]
        return reduce(lambda epart_a, epart_b: {block: epart_a[block] * epart_b[block]
                                                for block in epart_a.keys()}, 
                      eparts)
        
    propagated = []    
    
    for k in range(len(locations)):
        loc = locations[k]
        remaining_locs = locations[k+1:]
        loc_errors = {}
        
        loc_errors['X1'] = initialize_errors(blocklengths)
        if loc['type'] not in ('prepX', 'measX'):
            loc_errors['X1'][loc['block1']][loc['bit1']] = Pauli.X
            propagate_through_remaining(loc_errors['X1'], remaining_locs)
        
        loc_errors['Z1'] = initialize_errors(blocklengths)
        if loc['type'] not in ('prepZ', 'measZ'):
            loc_errors['Z1'][loc['block1']][loc['bit1']] = Pauli.Z
            propagate_through_remaining(loc_errors['Z1'], remaining_locs)
        
        if 'block2' in loc:
            loc_errors['X2'] = initialize_errors(blocklengths)
            loc_errors['X2'][loc['block2']][loc['bit2']] = Pauli.X
            propagate_through_remaining(loc_errors['X2'], remaining_locs)
            
            loc_errors['Z2'] = initialize_errors(blocklengths)
            loc_errors['Z2'][loc['block2']][loc['bit2']] = Pauli.Z
            propagate_through_remaining(loc_errors['Z2'], remaining_locs)
                
        # Now compute the error for each possible Pauli error
        # that can occur at this location.  Each error is a product
        # of the individual X/Z parts computed above.
        loc_errors = {error: error_product(loc_errors, error) 
                      for error in noise.errorListXZ[loc['type']]}
            
        propagated.append(loc_errors)
        
    return propagated
    
def propagate_errors_through_loc(errors, loc):
    '''
    Propagates the given errors though the location.
    :param errors: A dictionary of Pauli errors indexed block name.
    :param loc: The primitive circuit location.
    
    >>> import qfault.circuit.location as location
    >>> from qfault.qec.error import PauliError
    >>> errors = {'test': PauliError.fromstring('YX')}
    >>> cnot = location.cnot('test', 0, 'test', 1)
    >>> propagate_errors_through_loc(errors, cnot)
    {'test': YI}
    >>> measZ = location.meas(Pauli.Z, 'test', 0)
    >>> propagate_errors_through_loc(errors, measZ)
    {'test': XI}
    >>> prepZ = location.prep(Pauli.Z, 'test', 0)
    >>> propagate_errors_through_loc(errors, prepZ)
    {'test': II}
    >>> errors = {'test1': Pauli.X, 'test2': PauliError.fromstring('IZ')}
    >>> cnot = location.cnot('test1', 0, 'test2', 1)
    >>> propagate_errors_through_loc(errors, cnot)
    {'test1': Y, 'test2': IY}
    '''
    if loc['type'] in ('prepX', 'prepZ'):
        # Preparation resets any incoming errors
        errors[loc['block1']][loc['bit1']] = Pauli.I
    elif loc['type'] == 'measZ':
        # Z-basis measurement erases Z errors
        errors[loc['block1']][loc['bit1']] = \
            errors[loc['block1']][loc['bit1']].partial(Pauli.X)
    elif loc['type'] == 'measX':
        # X-basis measurement erases X errors
        errors[loc['block1']][loc['bit1']] = \
            errors[loc['block1']][loc['bit1']].partial(Pauli.Z)
    elif loc['type'] == 'cnot':
        ctrl_block, targ_block = loc['block1'], loc['block2']
        ctrl_bit, targ_bit = loc['bit1'], loc['bit2']
        # First X errors
        errors[targ_block][targ_bit] *= errors[ctrl_block][ctrl_bit].partial(Pauli.X)
        #Now Z errors
        errors[ctrl_block][ctrl_bit] *= errors[targ_block][targ_bit].partial(Pauli.Z)
        
    return errors

def count_errors_of_order_k(k, 
                            locations, 
                            noise_model, 
                            block_order=None,
                            block_error_maps=None):
    '''
    Counts the the errors that occur in the given set of locations with order 'k'
    according to the given noise model.  It returns a dictionary indexed by error.
    Each error is a tuple that represents the error on each logical block. Block
    errors are sorted according to the 'block_order', or according to the order
    of locations.blocknames() if no order is given. 
    The optional block_error_maps parameter may be used to map errors on each block
    to other kinds of keys (e.g. error syndromes).
    
    :param k: The error order (i.e., the number of faults)
    :param locations: The set of physical locations
    :param noise_model: The noise model
    :param block_order: (optional) An ordered list of block names.
    :param block_error_maps: A list of maps, one for each block.
    
    >>> import qfault.circuit.location as location
    >>> import qfault.noise as noise
    >>> cnot = location.cnot('test1', 0, 'test2', 0)
    >>> measX = location.meas(Pauli.X, 'test2', 0)
    >>> locations = location.Locations([cnot, measX])
    >>> noise_model = noise.CountingNoiseModelXZ()
    >>> count_errors_of_order_k(1, locations, noise_model)
    {(Z, Z): 2, (I, I): 1, (I, Z): 3, (X, Z): 2, (X, I): 2, (Y, Z): 2, (Y, I): 2, (Z, I): 2}
    >>> count_errors_of_order_k(1, locations, noise_model, block_order=('test2','test1'))
    {(Z, Z): 2, (I, I): 1, (I, Z): 2, (I, Y): 2, (Z, I): 3, (Z, Y): 2, (I, X): 2, (Z, X): 2}
    '''
    if None == block_order:
        block_order = locations.blocknames()
    
    if None == block_error_maps:
        block_error_maps = [concurrency.Noop()] * len(block_order)
    
    if 0 == k:
        # Only possible error with 0 failures is the trivial one
        lengths = locations.blocklengths()
        error = tuple(Pauli.I ** lengths[name] for name in block_order)
        return {error: 1}

    propagated_errors = propagate_location_errors(locations)
    location_index_sets = tuple(itertools.combinations(range(len(locations)), k))
    counts = concurrency.mapreduce_concurrent(functools.partial(_count_func,
                                                                locations,
                                                                propagated_errors,
                                                                noise_model,
                                                                block_order,
                                                                block_error_maps),
                                              _merge_counts,
                                              location_index_sets)
                
    return counts

def count_location_set(propagated_errors,
                       error_weights, 
                       block_order,
                       block_error_maps):
    '''
    Counts all error configurations for the given locations according 
    to the specified weights.
    :param propagated_errors: The set of propagated errors for each
                              location.
    :param error_weights: A list errors and weights for each location.  
                          The list should be of the form:
                          [[(error1, weight1), (error2, weight2), ...], 
                           [(error1, weight1), ...], ...].
                          It must have the same length as 
                          'propagated_errors'.
    :param key_generator: The error key generator.
    :param blocklengths: A dictionary of block lengths, indexed by 
                         block name.
    :rtype dict:  A dictionary of counts, indexed by error key.
    '''

    counts = {}
        
    for error_config in itertools.product(*error_weights):
        errors = [ec[0] for ec in error_config]
        weights = [ec[1] for ec in error_config]
            
        total_weight = listutils.mul(weights)
        block_errors = (listutils.mul(propagated_errors[i][err][name] 
                                            for i, err in enumerate(errors))
                        for name in block_order)
                
        error_key = tuple(block_error_maps[i](error) 
                          for i, error in enumerate(block_errors))
#        print 'blockErrors=', blockErrors, 'key=', errorKey, 'weight=', total_weight
        counts[error_key] = counts.get(error_key, 0) + total_weight

    return counts

#def count_blocks_by_syndrome(locations, blocks, noise, kMax):
#    counterUtils.propagateAllErrors(locations)
#    
#    keyGenerator = MultiBlockSyndromeKeyGenerator(blocks)
#    counts = [count_errors_of_order_k(k, locations, noise, keyGenerator)
#              for k in range(kMax+1)]
#    
#    return counts

#def mapCounts(counts, keymap):
#
#    newCounts = []
#    for countsK in counts:
#        newCountsK = {}
#        for key,count in countsK.iteritems():
#            mappedKey = keymap(key)
#            newCountsK[mappedKey] = newCountsK.get(mappedKey, 0) + count
#            
#        newCounts.append(newCountsK)
#        
#    return newCounts
#
#def maxCount(*countss):
#    '''
#    Converts the list of counts into a single set of counts that
#    maximizes all other counts.  More precisely, the returned counts
#    counts_max have the property:
#    counts_max[k][key] >= countss[i][k][key]
#    for every i, k and key.
#    '''
#    lengths = [len(counts) for counts in countss]
#    countsMax = [{} for _ in range(max(lengths))]
#    for counts in countss:
#        for k in range(len(counts)):
#            for key, val in counts[k].iteritems():
#                countsMax[k][key] = max(countsMax[k].get(key, val), val)
#                
#    return countsMax


def _count_func(locations,
                propagated_errors,
                noise_model,
                block_order,
                block_error_maps,
                indices):
    locs = [locations[i] for i in indices]
    prop_errs = [propagated_errors[i] for i in indices]
    
    error_weights = [[(e, noise_model.getWeight(l, e)) 
                      for e in noise_model.errorList(l)] 
                      for l in locs]

    return count_location_set(prop_errs, 
                              error_weights, 
                              block_order, 
                              block_error_maps)
    
def _merge_counts(counts):
    master = copy(counts[0])
    for count in counts[1:]:
        for key, val in count.iteritems():
            master[key] = master.get(key, 0) + val
            
    return master

if __name__ == '__main__':
    concurrency.initialize_concurrency(1)
    import doctest
    doctest.testmod()