#
# Cython wrapper for Golay exRec counting.  This wrapper exposes (to Python) the C implementation
# of countExRecForConfigs().
# See the corresponding setup.py distutils file.
# To compile, run
# python setup.py build_ext; mv build/lib.linux-i686-2.6/CountGolayExRecCython.so .
#

cdef extern from "exrec.c":
	object _countExRecForConfigs "countExRecForConfigs"  (object sLECLists,
						 object sCLists,
						 object countsLEC,
						 object countsC,
						 object lcountsTEC,
						 object configs,
						 object kMax)


def countExRecForConfigs_c(sLECLists, sCLists, countsLEC, countsC, lcountsTEC, configs):
	'''
	Counts the IX, XI, and XX malignant sets for the given failure configurations and
	for each of the possible combinations of trailing error corrections.  The
	underlying implementation is written in C which makes this function much faster than
	the equivalent Python implementation.
	
	Returns a list of malignant counts indexed as [ec][error][k], where:
	ec: value 0-3 indicating which trailing ECs are present (in binary)
	error: value 0-2 indicating the logical error (IX=0, XI=1, XX=2)
	k: the number of failing locations.
	
	Also returns the total number of sets that were counted.
	
	Arguments:
	----------
	sLECLists  -- A table indexed by [k][i]=syndrome, containing the syndromes for which countsLEC[k][syndrome] != 0
	sCLists    -- A table indexed by [k][i]=syndrome, containing the syndromes for which countsC[k][syndrome] != 0
	countsLEC  -- A table indexed by [k][syndrome] containing weighted counts for each syndrome at the output
	              of the leading error correction for fault order k.
	countsC    -- A table indexed by [k][syndrome] containing weighted counts for each syndrome at the output
                  of the transversal CNOT for fault order k.
	lcountsTEC -- A table indexed by [k][syndrome][logical] containing weighted counts for each input syndrome
	              to the trailing error correction.  logical=TRUE is the count for the event that ideal decoding
	              of the output of the TEC yeilds a logical X-error.
	configs    -- The configurations (kLECa, kLECb, kCnot, kTECa, kTECb) to count
	'''
	kMax = max([sum(c) for c in configs])
	totalSets = sum([len(sLECLists[c[0]]) * len(sLECLists[c[1]]) * len(sCLists[c[2]]) for c in configs])
	return _countExRecForConfigs(sLECLists, sCLists, countsLEC, countsC, lcountsTEC, configs, kMax), totalSets
