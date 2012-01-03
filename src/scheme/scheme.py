'''
Created on 2011-11-14

@author: adam
'''
from counting import probability, levelOne
from qec.error import Pauli, PauliError
from settings.noise import NoiseModelXSympy, NoiseModelZSympy, NoiseModelXZSympy
import logging
from counting.key import MultiBlockSyndromeKeyGenerator
from util import iteration, listutils
import itertools
from counting.component.adapter import InputAdapter

logger = logging.getLogger('scheme')


class Scheme(object):
    '''
    classdocs
    '''
    
    defaultNoiseModels = {Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy(),
              }


    def __init__(self):
        '''
        Constructor
        '''
        
    def count(self):
        raise NotImplementedError
                        
                        
    def countExRec(self, exRec, noiseModels, pauli=Pauli.Y, prInput=1):
        logicalProbabilities = {}
    
        result = exRec.count(noiseModels, pauli)

        prBad = exRec.prBad(noiseModels[pauli], pauli)
        prAccept = exRec.prAccept(noiseModels)
        
        locTotals = exRec.locations(pauli).getTotals()
        noise = noiseModels[pauli]
        
        print 'result.counts=', result.counts
        prAccepts = {p: prAccept(p/15) for p in (0, 0.0001, 0.001, 0.01)}
        prBads = {p: prBad(p/15) for p in (0, 0.0001, 0.001, 0.01)}
        prInputs = {p: prInput(p/15) for p in (0, 0.0001, 0.001, 0.01)}
        print 'Pr[accept]=', prAccepts
        print 'Pr[bad]=', prBads
        print 'Pr[input]=', prInputs
        pthresh = levelOne.pseudoThreshold(result.counts, locTotals, prBad, prAccept, noise, prInput)
            
        p = pthresh
        gamma = pthresh/15

        print 'pseudothreshold=', pthresh
        print 'Pr[bad]({0})'.format(p), prBad(gamma)
        prAccept = exRec.prAccept(noiseModels)
        print 'Pr[accept]({0})'.format(p), prAccept(gamma)
        
        # Construct a polynomial for each logical error
        # exRec counts are indexed by [k][error]
        # rearrange so that the index is [error][k]
        logicalCounts = {}
        for k, counts in enumerate(result.counts):
            for e, count in counts.iteritems():
                logicalCounts[e] = logicalCounts.get(e, [{e:0}]*k) + [{e:count}]
                
        for e, counts in logicalCounts.iteritems():
            pr = probability.countsToPoly(counts, locTotals, noise)
            pr = prInput * (pr / prAccept) + prBad
            logicalProbabilities[e] = pr
            
        for e, pr in logicalProbabilities.iteritems():
            print('Pr[E={0}] <= {1}'.format(e, pr(gamma)))
                
        return logicalProbabilities
                
    def getInputs(self, component):
        inBlocks = component.inBlocks()
        code = inBlocks[0].getCode()
        if not all(code == block.getCode() for block in inBlocks):
            raise Exception('All input blocks must use the same code')
        
        logger.info('Computing inputs for %s', component)
        
        eRange = range(1 << code.n)
        errorList = [PauliError(xbits=eX, zbits=eZ) for eX in eRange for eZ in eRange]
        keyGen = MultiBlockSyndromeKeyGenerator(inBlocks)
        
        inputs = {}
        for errors in itertools.combinations_with_replacement(errorList, len(inBlocks)):
            errors = {inBlocks[k].name: errors[k] for k in range(len(errors))}
            syndromes = tuple(code.getSyndrome(e) for e in errors.values())
            #pr = listutils.mul(self.prInputSyndrome(code.getSyndrome(e)) for e in errors.values())
            inputs[syndromes] = keyGen.getKey(errors)
            
        return inputs


        
    def getExRecs(self):
        raise NotImplementedError
    
    def prInputSyndrome(self, s):
        raise NotImplementedError