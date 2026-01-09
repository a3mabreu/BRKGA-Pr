# BRKGA-Pr

C++ implementation of a Biased Random-Key Genetic Algorithm (BRKGA) for matrix profile reduction using CSR [1].


## Using the code

After cloning this repository, navigate to the project directory to test, compile, and run the BRKGA-Pr heuristic.


### Automated unit tests

To run the automated unit tests:

`make test`

### Cleaning

To clean the build files:

`make clean`

### Compiling

To compile the program:

`make`

### Running

To run the program with original parameters for 120 seconds:

`./bin/main --filename input/usps_norm_5NN.mtx --max_time 120 --pop 20 --elite 8 --mutants 4 --prob 0.7`

### Command-Line Arguments

Available options:

--filename <path>: path to the input graph file in Matrix Market (.mtx) format.

--max_time <integer>: the maximum execution time allowed for the heuristic, specified in seconds.

--pop <integer>: the size of the population.

--elite <integer>: the size of the elite set.

--mutants <integer>: the number of new individuals introduced in each generation.

--prob <float>: the probability of inheriting the key from the elite parent.

## Comments

- Tested with the C++23 standard and the g++ 14.2.0 compiler in Linux.


## References

[1] W. Tinney, J. Walker, Direct solutions of sparse network equations by optimally ordered triangular factorization, Proceedings of the IEEE 55 (11) (1967) 1801–1809. doi:10.1109/PROC.1967.6011