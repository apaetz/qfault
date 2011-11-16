'''
Created on 2011-11-14

@author: adam
'''
from counting import probability, levelOne
from qec.error import Pauli
from settings.noise import NoiseModelXSympy, NoiseModelZSympy, NoiseModelXZSympy
import logging

logger = logging.getLogger('scheme')


class Scheme(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        
    def count(self):
        noiseModels = {Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy(),
              }
        
        countKey = Pauli.Y
        exrecLibrary = self.getExRecs()
        for name, exRecs in exrecLibrary.iteritems():
            
            
            for i, exRec in enumerate(exRecs):
            
                logger.info('Counting %s exRec %d of %d', name, i+1, len(exRecs))
                
                result = exRec.count(noiseModels, countKey)
            
                prBad = exRec.prBad(noiseModels[countKey], countKey)
                prAccept = exRec.prAccept(noiseModels)
                
                locTotals = exRec.locations(countKey).getTotals()
                noise = noiseModels[countKey]
                
                print 'result.counts=', result.counts
                pthresh = levelOne.pseudoThreshold(result.counts, locTotals, prBad, prAccept, noise)
                    
                p = pthresh
                gamma = pthresh/15
    
                print 'pseudothreshold=', pthresh
                print 'Pr[bad]({0})'.format(p), prBad(gamma)
                prAccept = exRec.prAccept(noiseModels)
                print 'Pr[accept]({0})'.format(p), prAccept(gamma)

        
    def getExRecs(self):
        raise NotImplementedError