#
# Some basic routines for preparing encoded stabilizer states.
#
# Adam Paetznick (9/12/2011)
from util.counterUtils import locXprep, locZprep, loccnot, locrest, \
	propagateAllErrors


def permuteCircuit(circuit, rounds):
	return [circuit[r] for r in rounds]

def ancillaZPrep(schedule, roundPermutation=None, name='0'):
	"""
	Returns the circuit for preparing encoded |0> based on the given CNOT schedule (given as a list of rounds).
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
	
	# A qubit that is a target of any CNOT in the prep circuit should
	# be prepared as |0>.  All other qubits should be prepared as |+>
	prepAsZ  = set([c[1] for c in flatcnots])
	prepAsX = allQubits - prepAsZ

	# the initial qubit preparations
	for i in prepAsX:
		addLoc(locXprep(name, i))
	for i in prepAsZ:
		addLoc(locZprep(name, i))

	touched = set()
	
	for rnd in xrange(len(cnots)):			# simulate each round sequentially
		resting = set(touched)
		for src, tgt in cnots[rnd]:
			cnotQubits = set([src, tgt])
			addLoc(loccnot(name, src, name, tgt))
			resting -= cnotQubits
			touched.update(cnotQubits)
			
		# now add the rest locations, provided we are in at least the second round
		for r in resting:
			addLoc(locrest(name, r))
			
	propagateAllErrors(locations)
	return locations