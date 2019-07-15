from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quil import DefGate
import numpy as np

#I have added a line

# Nbit CZ, but then stores result to N + 1 bit
# Needs to be unitary, so if we map 11..0 to 1, need to also map 11..1 to 0
def create_Np1_bit_CZ(N: int):
	my_mat = np.identity(2 ** (N + 1))

	# middle = ((2 ** N) // 2) - 1
	last = 2 ** (N + 1) - 1

	# my_mat[middle, middle] = 0 #Middle row, middle column 11..0 
	# my_mat[middle, last] = 1 # Bottom row, middle column 11..0 -> 11..1

	my_mat[last, last] = 0
	my_mat[last - 1, last - 1] = 0
	my_mat[last - 1, last] = 1
	my_mat[last, last - 1] = 1

	Np1_bit_CZ = DefGate("Nbit_CZ", my_mat)
	return Np1_bit_CZ, Np1_bit_CZ.get_constructor()

#Takes a superposition of lots of bits in, checks that the oracle chooses 11..1.
def debug_oracle_selector(N: int):
	total_qb = N + 1
	qc = get_qc(str(total_qb) + 'q-qvm')
	p = Program()

	for i in range(N):
		p += H(i)

    #Ancilla to tell us if we have 1111 state
	Np1_bit_CZ_def, Np1_bit_CZ = create_Np1_bit_CZ(N)
	p += Np1_bit_CZ_def

	p += Np1_bit_CZ(*range(total_qb))

	results = qc.run_and_measure(p, trials = 30)
	
	for i in range(total_qb):
		print("qubit " + str(i) +  ": " + str(results[i]))


def debug_oracle_calc
