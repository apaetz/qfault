#
# Some basic routines for preparing encoded stabilizer states.
#
# Adam Paetznick (9/12/2011)
from qfault.circuit import location
from qfault.qec.error import PauliError, Pauli


def permuteCircuit(circuit, rounds):
	return [circuit[r] for r in rounds]

def ancillaZPrep(schedule, roundPermutation=None, name='0', reverseQubitOrder=True):
	"""
	Returns the circuit for preparing encoded |0> based on the given CNOT schedule (given as a list of rounds).
	By default, the qubit ordering of the circuit is reversed.  For example, if there are n qubits, then 
	n <-> 0, (n-1) <-> 1, etc.  This is usually desired because the location counting code uses a
	descending (little endian) bit ordering.  That is, an error on qubit k is represented by a 1 on bit
	k.  However, it is more natural to express operators, stabilizer generators in particular, in ascending
	(big endian) bit ordering.  For example XII is an X on qubit 0, not on qubit 2.
	To disable qubit reversal, set reverseQubitOrder=False.
	"""
	if None != roundPermutation:
		cnots = permuteCircuit(schedule, roundPermutation)
	else:
		cnots = schedule
	
	locations = []
	def addLoc(loc):
		locations.append(loc)

	#flatcnots = reduce(lambda x,y: x+y, cnots)
	flatcnots = [cnot for round in cnots for cnot in round]
	allQubits = set(qubit for cnot in flatcnots for qubit in cnot)
	
	numQubits = max(allQubits)
	bitLookup = [i for i in range(numQubits+1)]
	if reverseQubitOrder:
		bitLookup = list(reversed(bitLookup))
	
	# Qubits that are used first as a control are prepared as |+>.
	# Qubits taht are used first as a target are prepared as |0>.
	controls = [c[0] for c in flatcnots]
	targets = [c[1] for c in flatcnots]
	prepAsZ  = set()
	for q in range(len(allQubits)):
		if (q not in controls) or (q in targets and controls.index(q) > targets.index(q)):
			prepAsZ.add(q)

	prepAsX = allQubits - prepAsZ

	# the initial qubit preparations
	for i in prepAsX:
		addLoc(location.prep(Pauli.X, name, bitLookup[i]))
	for i in prepAsZ:
		addLoc(location.prep(Pauli.Z, name, bitLookup[i]))

	touched = set()
	
	for rnd in xrange(len(cnots)):			# simulate each round sequentially
		# Rests will only be added for qubits which have already been touched.
		resting = set(touched)
		
		for src, tgt in cnots[rnd]:
			cnotQubits = set([src, tgt])
			addLoc(location.cnot(name, bitLookup[src], name, bitLookup[tgt]))
			resting -= cnotQubits
			touched.update(cnotQubits)
			
		# now add the rest locations
		for r in resting:
			addLoc(location.rest(name, bitLookup[r]))
			
	#propagateAllErrors(locations)
	return locations