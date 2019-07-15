from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quil import DefGate
import numpy as np

from pyquil.api import WavefunctionSimulator
wf_sim = WavefunctionSimulator()

# Create core operation for diffusion step
def create_diffusion_op(n):
	a = 2.0 * np.full((2**n, 2**n), 1/(2**n)) - np.identity(2**n)
	diff_def = DefGate("DIFF", a)
	DIFF = diff_def.get_constructor()
	return diff_def, DIFF

def grovers(N: int, begin: str, end: str, rounds = -1):
	ancilla_qbits = 0
	debug_qbits = 0
	total_qbits = N + ancilla_qbits + debug_qbits

	qc = get_qc(str(total_qbits) + 'q-qvm')
	p = Program()

	### Make Gate Def'ns ###
	diff_def, diffusion_op = create_diffusion_op(N)
	
	### Make program def'ns ###
	oracle = make_oracle(begin, end, N)
	
	### Add gates def'n ###
	p += diff_def

	# initalize state of non-ancilary qubits
	for i in range(N):
		p += H(i)

	# default number of rounds
	if rounds == -1:
		rounds = int(np.pi/4 * np.sqrt(2**N))
	print("Rounds: " + str(rounds))

	for r in range(2): #debug, rounds -> 1
		p += oracle #Oracle isn't gate, contains subcircuits
		p += diffusion_op(*range(N)) 

		# # simulator debug
		# wavefunction = wf_sim.wavefunction(p)
		# print(wavefunction)

	# Real run
	results = qc.run_and_measure(p, trials = 30)
	for i in range(total_qbits):
		print("qubit " + str(i) +  ": " + str(results[i]))

# Construct diag operator that returns -1 if all 
# gates are 1
def create_Nbit_CZ(N: int):
	my_mat = np.identity(2 ** N)
	my_mat[2 ** N - 1, 2 ** N - 1] = -1
	Nbit_CZ_def = DefGate("Nbit_CZ", my_mat)
	return Nbit_CZ_def, Nbit_CZ_def.get_constructor()

# NZ: I propose Oracle should be a circuit. If it's a numpy array, we
# can't take advantage of the abstraction of building circuits using 
# pre-existing pyquil gates
def make_oracle(begin: str, end: str, N: int):

	oracle = Program()

	# Compute key XOR begin
	for i in range(N):
		if begin[i] == "0":
			oracle += I(int(i))

		elif begin[i] == "1":
			oracle += X(int(i))

	# Check prev result == end
	for i in range(N):
		if end[i] == "0":
			oracle += X(int(i))

		elif end[i] == "1":
			oracle += I(int(i))

	# If we have all one's, i.e, we have a match, invert phase
	Nbit_CZ_def, Nbit_CZ = create_Nbit_CZ(N)
	oracle += Nbit_CZ_def

	# Gate arguments must be non-empty list, even if the gates 
	# takes all qubits as input. See https://stackoverflow.com/questions/3941517/converting-list-to-args-when-calling-function
	oracle += Nbit_CZ(*(range(N)))

	# Undo our previous modifications to input
	for i in range(N):
		if end[i] == "0":
			oracle += X(int(i))

		elif end[i] == "1":
			oracle += I(int(i))

	for i in range(N):
		if begin[i] == "0":
			oracle += I(int(i))

		elif begin[i] == "1":
			oracle += X(int(i))

	return oracle

grovers(4, "0011", "0110")
