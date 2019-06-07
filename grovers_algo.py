from pyquil import Program, get_qc, list_quantum_computers
from pyquil.noise import add_decoherence_noise
from pyquil.gates import *
from pyquil.quil import DefGate
import numpy as np
from pyquil.api import QPU, QuantumComputer


from pyquil.api import WavefunctionSimulator
wf_sim = WavefunctionSimulator()

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
	print("Rounds: " + str(rounds))

	diff_def, diffusion_op = create_diffusion_op(N)
	p += diff_def
	oracle = make_shift_rows_oracle(begin, end, N, offsets)
	#oracle = make_oracle(begin, end, N)

	# Main Loop
	for r in range(rounds): #debug, rounds -> 1
		p += oracle #Oracle isn't gate, contains subcircuits
		p += diffusion_op(*range(N)) 

	
	# compile
	p = qvm.compile(p)
	p = Program(p.program)
	p = add_decoherence_noise(p, nlevel*3e-5, nlevel*3e-5, 5e-8, 1.5e-7, 1)
	p.pop()
	for i in range(N):
		p += MEASURE(i, ro[i])
	for n in range(50):	
		results = qvm.run(p)[0]
		#results = qvm.run_and_measure(Program(p.program), trials = trials)
		for j in range(N):
			ret[j].append(results[j])
	
	return ret

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

def make_shift_rows_oracle(begin: str, end: str, N: int, offsets: list):
	oracle = Program()

	### Gate defs 
	Nbit_CZ_def, Nbit_CZ = create_Nbit_CZ(N)
	oracle += Nbit_CZ_def

	rows = int(np.sqrt(N))
	cols = int(N/rows)

	rounds = 1

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

# Builds oracle
def shift_row_offset(begin: str, end: str, N: int, offset: int, start: int):

	oracle = Program()

	### Gate defs 
	Nbit_CZ_def, Nbit_CZ = create_Nbit_CZ(N)
	oracle += Nbit_CZ_def

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

	# If we have all one's, i.e, we have a match, invert phase
	oracle += Nbit_CZ(*(range(start, start + N)))

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

""" 
Requires 4 arguments:
	0) number of qubits in system
	1) begin bistring (plaintext)
	2) end bitstring (ciphertext)
	4) list of offsets per row

Note that the number of qubits should be a square and that
the begin and end bitstrings will be interpreted as square 
arrays with sqrt(n) rows and sqrt(n) columns.

""" 

print(list_quantum_computers())

trials = 15
stats = []
for i in np.arange(10, 20):
	for qc in ['4q-qvm']:
		print("nlevel: " + str(i))
		results = run_grovers(qc, 4, "1011", "1100", offsets = [1, 0], trials=trials, nlevel=i)
		print(results)
		true = [1, 0, 1, 1]
		stat = [float(sum(results[i]))/trials if true[i] else 
				 float(trials-sum(results[i]))/trials for i in range(len(true))]

# 	acc=0.0
# 	for i in range(trials):
# 		for j in range(len(true)):
# 			if(results[j][i] != true[j]):
# 				acc-=1
# 				break
# 		acc+=1

# 	stat += [acc/trials]
# 	stats.append(qc + "," + str(stat))

# 	for i in range(4):
# 		print("qubit " + str(i) +  ": " + str(results[i]))

# with open('stats.txt', 'w') as f:
# 	for s in stats:
# 		f.write(s+"\n")



	#	print("worked for " + qc)
	#	break
	#except:
	#	continue
