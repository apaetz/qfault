'''
Created on 2011-11-14

@author: adam
'''
from counting import probability, levelOne
from qec.error import Pauli
from settings.noise import NoiseModelXSympy, NoiseModelZSympy, NoiseModelXZSympy

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
        
        pauli = Pauli.Y
        exrecLibrary = self.getExRecs()
        for name, exRecs in exrecLibrary.iteritems():
            for exRec in exRecs:
                result = exRec.count(noiseModels, pauli)
            
                prBad = exRec.prBad(noiseModels[pauli], pauli)
                prAccept = exRec.prAccept(noiseModels)
                
                locTotals = exRec.locations(pauli).getTotals()
                noise = noiseModels[pauli]
                
                print 'result.counts=', result.counts
                pthresh = 1
                for pauli in (Pauli.X, Pauli.Z, Pauli.Y):
                    counts = [{pauli:c.get(pauli, 0)} for c in result.counts]
                    countPoly = probability.countsToPoly(counts, locTotals, noise)
                    pr = countPoly / prAccept + prBad
                    #print 'Pr[', pauli, ']({0})='.format(p), pr(gamma)
                    print counts
                    pthreshNew = levelOne.pseudoThreshold(counts, exRec.locations().getTotals(), prBad, prAccept, noise)
                    pthresh = min(pthresh, pthreshNew)
                    
                p = pthresh
                gamma = pthresh/15
    
                print 'Pr[bad]({0})'.format(p), prBad(gamma)
                prAccept = exRec.prAccept(noiseModels)
                print 'Pr[accept]({0})'.format(p), prAccept(gamma)

        
    def getExRecs(self):
        raise NotImplementedError