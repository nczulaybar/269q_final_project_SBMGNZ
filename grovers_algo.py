from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quil import DefGate
import numpy as np

# Create core operation for diffusion step
def create_diffusion_op(n):
	a = 2.0 * np.full((2**n, 2**n), 1/(2**n)) - np.eye(2**n)
	diff_def = DefGate("DIFF", a)
	DIFF = diff_def.get_constructor()
	return diff_def, DIFF

def grovers(num_qubits, oracle, rounds = -1):
	qc = get_qc(str(num_qubits) + 'q-qvm')
	p = Program()

	diff_def, diffusion_op = create_diffusion_op(num_qubits)
	p += diff_def

	oracle_def = DefGate("ORACLE", oracle)
	ORACLE = oracle_def.get_constructor()
	p += oracle_def

	diff_args = range(num_qubits)

	# initalize state of non-ancilary qubits
	for i in range(num_qubits):
		p += H(i)

	# default number of rounds
	if rounds == -1:
		rounds = int(np.pi/4 * np.sqrt(2**num_qubits))
	print("Rounds: " + str(rounds))

	for r in range(rounds):
		# appy oracle
		p += ORACLE(*diff_args)
		p += diffusion_op(*diff_args)
		

	results = qc.run_and_measure(p, trials = 10)
	for i in range(num_qubits):
		print("qubit " + str(i) +  ": " + str(results[i]))

"""
Notes:
	1) oracle should be a numpy array
	2) make sure not to run for more than the suggested number of rounds, 
	   because the accuracy falls off (I'm not quite sure why)

"""
SEARCHED_STRING = "1101"
N = len(SEARCHED_STRING)
oracle = np.zeros(shape=(2 ** N, 2 ** N))
for b in range(2 ** N):
    if np.binary_repr(b, N) == SEARCHED_STRING:
        oracle[b, b] = -1
    else:
        oracle[b, b] = 1

grovers(N, oracle)
