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

def grovers(num_qubits: int, oracle, rounds = -1):
	qc = get_qc(str(num_qubits) + 'q-qvm')
	p = Program()

	diff_def, diffusion_op = create_diffusion_op(num_qubits)
	p += diff_def

	# oracle_def = DefGate("ORACLE", oracle)
	# ORACLE = oracle_def.get_constructor()
	# p += oracle_def
	p += oracle

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
		p += oracle #ORACLE(*diff_args)
		p += diffusion_op(*diff_args)
		

	results = qc.run_and_measure(p, trials = 10)
	for i in range(num_qubits):
		print("qubit " + str(i) +  ": " + str(results[i]))

# """
# Notes:
# 	1) oracle should be a numpy array
# 	2) make sure not to run for more than the suggested number of rounds, 
# 	   because the accuracy falls off (I'm not quite sure why)

# """
# SEARCHED_STRING = "1101"
# N = len(SEARCHED_STRING)
# oracle = np.zeros(shape=(2 ** N, 2 ** N))
# for b in range(2 ** N):
#     if np.binary_repr(b, N) == SEARCHED_STRING:
#         oracle[b, b] = -1
#     else:
#         oracle[b, b] = 1

# Construct diag operator that returns -1 if all 
# gates are 1
def create_Nbit_CZ(N: int):
	my_mat = np.identity(2 ** N)
	my_mat[2 ** N - 1, 2 ** N - 1] = -1
	#my_mat[0,0] = -1
	Nbit_CZ_def = DefGate("Nbit_CZ", my_mat)
	return Nbit_CZ_def, Nbit_CZ_def.get_constructor()

# NZ: I propose Oracle should be a circuit. If it's a numpy array, we
# can't take advantage of the abstraction of building circuits using 
# pre-existing pyquil gates
def make_oracle(begin: str, end: str, N: int):

	oracle = Program()

	# Compute key XOR begin
	for i in range(N):
		if begin[i] == 0:
			oracle += I(int(i))

		elif begin[i] == 1:
			oracle += X(int(i))

	# Check prev result == end
	for i in range(N):
		if end[i] == 0:
			oracle += X(int(i))

		elif end[i] == 1:
			oracle += I(int(i))

	# If we have all one's, i.e, we have a match, invert phase
	Nbit_CZ_def, Nbit_CZ = create_Nbit_CZ(N)
	oracle += Nbit_CZ_def

	# Gate arguments must be non-empty list, even if the gates 
	# takes all qubits as input. See https://stackoverflow.com/questions/3941517/converting-list-to-args-when-calling-function
	oracle += Nbit_CZ(*(range(N)))

	# Undo our previous modifications to input
	for i in range(N):
		if end[i] == 0:
			oracle += X(int(i))

		elif end[i] == 1:
			oracle += I(int(i))

	for i in range(N):
		if begin[i] == 0:
			oracle += I(int(i))

		elif begin[i] == 1:
			oracle += X(int(i))

	return oracle


oracle = make_oracle("11", "01", 2)

grovers(2, oracle)



