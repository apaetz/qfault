'''
Data structures for quantum circuit locations.

The primitive gate set implemented here is as follows,
{|0>, |+>, I (rest), CNOT, Z-basis measurement, X-basis measurement}.

The intended use of these circuits is for Pauli-error propagation and
fault-tolerant threshold analysis.  Pauli gates are not included in 
the gate set because they can be tracked classically through the 
Pauli frame.  Hadamard gates (or S or T, etc.) could be added, but 
have thus far been unnecessary due to availability of |+>.

Circuit and gate representation
--------------------------------
A circuit is simply a sequence of primitive gate locations.  The (qu)bits of
the circuit can be partitioned into logical blocks. Each gate input is assigned
a block and a bit position within the block.

Each primitive gate is represented by a Python dictionary with the following
keys: type, block1, bit1. Fanin-2 gates (i.e. CNOT) also contain keys: block2
bit2.

@author: Adam Paetznick
Some parts are adapted from code by Ben Reichardt.
'''


class Locations(object):
    '''
    Container for a set of circuit locations.  The primary motivation for this
    class is to allow a name to be attached to the list of locations.
    '''
    
    def __init__(self, sequence, name=''):
        self.list = list(sequence)
        self.name = name
        
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return str(self)

    def __getitem__(self, index):
        return self.list[index]
    
    def __len__(self):
        return len(self.list)
    
    def __add__(self, other):
        return Locations(self.list + other.list, str(self) + ':' + str(other))
                        
    def blocknames(self):
        '''
        Returns a list of all of the block names.
        '''
        return self.blocklengths().keys()
        
    def blocklengths(self):
        '''
        Returns a dictionary of block lengths,
        indexed by block name.
        '''
        blocks = {}
        for loc in self:
            blockname = loc['block1']
            bit = loc['bit1']
            blocks[blockname] = max(blocks.get(blockname, 0), bit+1)
            if 'block2' in loc: 
                blockname = loc['block2']
                bit = loc['bit2']
                blocks[blockname] = max(blocks.get(blockname, 0), bit+1)
        return blocks


def rest(block, bit):
    return {'type': 'rest', 'block1': block, 'bit1': bit}

def prep(pauli, block, bit):
    return {'type': 'prep'+str(pauli), 'block1': block, 'bit1': bit}

def cnot(block1, bit1, block2, bit2):
    return {'type': 'cnot', 'block1': block1, 'bit1': bit1, 'block2': block2, 'bit2': bit2}

def meas(pauli, block, bit):
    return {'type': 'meas'+str(pauli), 'block1': block, 'bit1': bit}

def supported_types():
    '''
    Returns a list of all of the supported location types.
    '''
    return ('cnot', 'rest', 'prepX', 'prepZ', 'measX', 'measZ')