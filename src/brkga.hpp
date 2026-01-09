#ifndef BRKGA_H
#define BRKGA_H

#include "misc.hpp"
#include "random_keys.hpp"


/**** BRKGA-Pr
@param P: Population size
@param E: Elite set size
@param R: Number of mutant individuals (new solutions)
@param PROB: Probability of inheriting the key from the elite parent
@param INIT: Constructive method for the initial population
*/
void brkga(CSR& csr, const usize P, const usize E, const usize R, const realT PROB, const usize INIT);

void initPopulation(CSR& csr, const usize INIT_V, std::vector<SolutionRK>& population, const usize N);

#endif