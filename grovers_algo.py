from pyquil import Program, get_qc, list_quantum_computers
from pyquil.noise import add_decoherence_noise
from pyquil.gates import *
from pyquil.quil import DefGate
import numpy as np
from pyquil.api import QPU, QuantumComputer
import argparse
from pyquil.api import WavefunctionSimulator
wf_sim = WavefunctionSimulator()

parser = argparse.ArgumentParser()
parser.add_argument("-m", type=str, default="1011", help="message")
parser.add_argument("-c", type=str, default="1100", help="ciphertext")
parser.add_argument("-k", type=str, default="1011", help="correct key")
parser.add_argument("-offsets", nargs='*', default=[1,0], help="list of row offsets")
parser.add_argument("-noise_level_s", type=int, default=1, help="signals highest amount of noise (larger value=less noise)")
parser.add_argument("-noise_level_e", type=int, default=10000, help="signals lowest amount of noise")
parser.add_argument("-noise_level_incr", type=int, default=100, help="signals noise level incremention (larger=more detailed stats)")
parser.add_argument("-trials", type=int, default=50)
parser.add_argument("-qvms", nargs='*', default=['4q-qvm'], help="a list of qvm architectures to test")

args = parser.parse_args()

# Builds Grover's algorithm program
def run_grovers(qc:str, N: int, begin: str, end: str, rounds = -1, offsets = [1, 0], trials = 20, nlevel = 0):
	print("Using qvm: " + qc)
	
	ret = {}
	for i in range(N):
		ret[i] = [] 

	qvm = get_qc(qc, as_qvm=True)

	p = Program()
	ro = p.declare('ro', 'BIT', N)

	# Initialize Non-Ancillary QBits
	for i in range(N):
		p += H(i)

	# Main Loop Setup
	if rounds == -1: # Default number of rounds
		rounds = int(np.pi/4 * np.sqrt(2**N))
	print("Running for " + str(rounds) + " rounds")

	diff_def, diffusion_op = create_diffusion_op(N)
	p += diff_def
	oracle = make_shift_rows_oracle(begin, end, N, offsets)

	# Main Loop
	for r in range(rounds): 
		p += oracle 
		p += diffusion_op(*range(N)) 

	
	# compile
	p = qvm.compile(p)
	p = Program(p.program)
	p = add_decoherence_noise(p, nlevel*3e-5, nlevel*3e-5, 5e-8, 1.5e-7, 1)
	p.pop() # Remoes HALT instruction
	for i in range(N):
		p += MEASURE(i, ro[i])
	for n in range(trials):	
		results = qvm.run(p)[0]
		for j in range(N):
			ret[j].append(results[j])
	
	return ret

def make_shift_rows_oracle(begin: str, end: str, N: int, offsets: list):
	oracle = Program()

	### Gate defs 
	Nbit_CZ_def, Nbit_CZ = create_Nbit_CZ(N)
	oracle += Nbit_CZ_def

	rows = int(np.sqrt(N))
	cols = int(N/rows)

	rounds = 1 # Oracle is only designed to work when rounds=1; future work could
			   # make it general to any number of rounds

	for r in range(rounds):
		for i in range(rows):
			oracle += shift_row_p1(begin[i*cols: (i+1)*cols], end[i*cols: (i+1)*cols],
							cols, ((r+1)*offsets[i])%cols, i*cols)
	
	oracle += 	Nbit_CZ(*(range(N)))
	
	for r in range(rounds):
		for i in range(rows):
			oracle += shift_row_p2(begin[i*cols: (i+1)*cols], end[i*cols: (i+1)*cols],
						cols,((r+1)*offsets[i])%cols, i*cols)

	return oracle

# Builds oracle
def shift_row_p1(begin: str, end: str, N: int, offset: int, start: int):

	oracle = Program()

	# Compute key XOR begin
	for i in range(N):
		if begin[i] == "0":
			oracle += I(start + int((i+offset)%N))

		elif begin[i] == "1":
			oracle += X(start + int((i+offset)%N))

	# Check prev result == end
	for i in range(N):
		if end[i] == "0":
			oracle += X(start + int(i))

		elif end[i] == "1":
			oracle += I(start + int(i))

	return oracle


def shift_row_p2(begin: str, end: str, N: int, offset: int, start: int):
	oracle = Program()

	# Undo our previous modifications to input
	for i in range(N):
		if end[i] == "0":
			oracle += X(start + int(i))

		elif end[i] == "1":
			oracle += I(start + int(i))

	for i in range(N):
		if begin[i] == "0":
			oracle += I(start + int((i+offset)%N))

		elif begin[i] == "1":
			oracle += X(start + int((i+offset)%N))

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

assert(len(args.m) == len(args.c))
assert(float(len(args.offsets)) == float(np.sqrt(len(args.m))))
assert(args.noise_level_s < args.noise_level_e)
assert(args.noise_level_s > 0)

stats = []
for i in np.arange(args.noise_level_s, args.noise_level_e, args.noise_level_incr):
	for qc in args.qvms:
		print("Noise level: " + str(i))
		results = run_grovers(qc, len(args.m), args.m, args.c, offsets = args.offsets, trials= args.trials, nlevel=i)
		
		true = [int(n) for n in list(args.k)]
		stat = [float(sum(results[i]))/args.trials if true[i] else 
				 float(args.trials-sum(results[i]))/args.trials for i in range(len(true))]

		acc=0.0
		for i in range(args.trials):
			for j in range(len(true)):
				if(results[j][i] != true[j]):
					acc-=1
					break
			acc+=1

		stat += [acc/args.trials]
		stats.append(qc + "," + str(stat))

		for i in range(len(args.m)):
			print("Qubit " + str(i) +  ": " + str(results[i]))

	with open('stats.txt', 'w') as f:
		for s in stats:
			f.write(s+"\n")