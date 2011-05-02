'''
Created on 2010-05-30

@author: adam
'''
from cache import DataManager
import os

if __name__ == '__main__':
    import sys
    key = sys.argv[1]
    
    dm = DataManager(os.path.curdir + os.path.sep)
    
    obj = dm.load(key)
    
    print obj