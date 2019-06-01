from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quil import DefGate
import numpy as np

from pyquil.api import WavefunctionSimulator
wf_sim = WavefunctionSimulator()

# Builds Grover's algorithm program
def run_grovers(N: int, begin: str, end: str, rounds = -1):
	
	qc = get_qc(str(N) + 'q-qvm')

	p = Program()

	# Initialize Non-Ancillary QBits
	for i in range(N):
		p += H(i)

	# Main Loop Setup
	if rounds == -1: # Default number of rounds
		rounds = int(np.pi/4 * np.sqrt(2**N))
	print("Rounds: " + str(rounds))

	diff_def, diffusion_op = create_diffusion_op(N)
	p += diff_def
	oracle = make_oracle(begin, end, N)

	# Main Loop
	for r in range(rounds): #debug, rounds -> 1
		p += oracle #Oracle isn't gate, contains subcircuits
		p += diffusion_op(*range(N)) 

	# Execute and Output
	results = qc.run_and_measure(p, trials = 30)
	for i in range(N):
		print("qubit " + str(i) +  ": " + str(results[i]))

# Builds oracle
def make_oracle(begin: str, end: str, N: int):

	oracle = Program()

	### Gate defs 
	Nbit_CZ_def, Nbit_CZ = create_Nbit_CZ(N)
	oracle += Nbit_CZ_def

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

# Create core operation for diffusion step
def create_diffusion_op(n):
	a = 2.0 * np.full((2**n, 2**n), 1/(2**n)) - np.identity(2**n)
	diff_def = DefGate("DIFF", a)
	DIFF = diff_def.get_constructor()
	return diff_def, DIFF

# Construct diag operator that returns -1 if all gates are 1
def create_Nbit_CZ(N: int):
	my_mat = np.identity(2 ** N)
	my_mat[2 ** N - 1, 2 ** N - 1] = -1
	Nbit_CZ_def = DefGate("Nbit_CZ", my_mat)
	return Nbit_CZ_def, Nbit_CZ_def.get_constructor()

run_grovers(4, "0011", "0110")

