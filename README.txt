This collection of Python modules implements the malignant set counting procedure described in TBD.
Its primary purpose is to compute a lower bound on the tolerable depolarizing noise threshold for fault-tolerant circuits based on the Golay code. 

1. Software requirements 
------------------------
Python modules were written and tested for Python version 2.7.  Some or all modules may be compatible with previous versions, but this has not been tested.
The following third-party packages are also required.  Tested versions are listed.
	- Sympy 0.7.1.rc1 http://sympy.org
	- matplotlib 1.0.0 http://matplotlib.sourceforge.net
	- gmpy 1.11 http://gmpy.sourceforge.net
	- Cython 0.12.1 http://cython.org
	- numpy 1.3.0 http://numpy.scipy.org
	- gcc 3.4.3 http://gcc.gnu.org (for build purposes only)
	
	
TODO:  The build and run instructions are out of date.

2. Build instructions
---------------------
Functions for the most time consuming part of the counting procedure are written in C and bound to a Python interface using Cython.
To compile the C source and build Python bindings run,
  cd src/component/cython
  sh buildCython
  
3. How to run the code
----------------------
The top-level Python module is located at src/GolayCounting.py.  It is recommended run this module by executing the wrapper script 'countGolay'.
Usage of countGolay is as follows:
  countGolay [nslots] [prep] [countRests]

  Options:
    nslots       -- The number of threads (slots/CPUs/cores) to use
    prep         -- The ancilla preparation, either 'rand' or 'overlap'
    countRests   -- 'true' to count rest locations, 'false' to ignore rest locations 
    
Running time on 31 cores is approximately four days.  Checkpoint results are saved as compressed pickle files in the 'data/' directory and allow execution to be be terminated and then restarted at a later time.