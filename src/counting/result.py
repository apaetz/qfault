'''
Created on 2011-10-22

@author: adam
'''    

class CountResult(object):
    ''' 
    '''
    
    def __init__(self, counts, keyMeta, blocks, rejectedCounts=None, name=None):
                    
        self.blocks = blocks
        self.counts = counts
        self.rejected = rejectedCounts
        #self.prAccept = prAccept
        self.keyMeta = keyMeta
        
    def __get__(self, etype):
        return self.counts[etype]
    
#    def counts(self):
#        return self.counts
#    
#    def rejectedCounts(self):
#        return self.rejected
#    
#    def prAccept(self):
#        return self.prAccept
#    
#    def blocks(self):
#        return self.blocks
#    
#    def keyMeta(self):
#        return self.keyMeta