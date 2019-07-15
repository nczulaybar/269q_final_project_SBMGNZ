# 269Q Final Project 
Sawyer Birnbaum
Maithreyi Gopalakrishnan
Nicolo Zulaybar

Instructions:

1. First, install PyQuil 2.0 (see http://docs.rigetti.com/en/stable/start.html for instructions)
2. Run "qvm -S" and "quilc -S" to setup QVM and Quilc Compiler instances
3. To test Grover's algorithm, run "python grovers_algo.py." The script takes a number of arguments.

* -m: this is the plaintext message
* -c: this is the encrypted ciphertext; the oracle will be built base on the message and ciphertext
* -k: this is the correct key used to encrypt the message
* -offsets: an array of per row offsets used in the ShiftRow step
* -noise_level_s, noise_level_e, noise_level_incr: the Noise Level varies from noise_level_s to noise_level_e in increments of noise_level_incr. A higher Noise Level means less noise.
* -trials: the number of trials to run per qvm/noise level
* -qvms: names of QVM's that will be tested
