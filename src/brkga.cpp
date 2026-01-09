#include "brkga.hpp"

/**** BRKGA-Pr
@param P: Population size
@param E: Elite set size
@param R: Number of mutant individuals (new solutions)
@param PROB: Probability of inheriting the key from the elite parent
@param INIT: Constructive method for the initial population
*/
void brkga(CSR& csr, const usize P, const usize E, const usize R, const realT PROB, const usize INIT) {
    std::vector<SolutionRK> population(P), next_population(P);
    csr.random_keys.resize(csr.m);
    csr.tmp_rk.resize(csr.m);
    csr.indexed_rk.resize(csr.m);
    std::uniform_int_distribution<usize> distributionElite(0, E - 1);
    std::uniform_int_distribution<usize> distributionPop(0, P - 1);

    /// Initial population
    initPopulation(csr, INIT, population, P);

    /// Alocating and next_population (Elite ones are moved)
    for (usize i = E; i < P; ++i) {
        next_population[i].random_keys.resize(csr.m);
        next_population[i].labels.resize(csr.m);
        std::iota(next_population[i].labels.begin(), next_population[i].labels.end(), 0);
    }

    const auto start = std::chrono::steady_clock::now();
    /// Main loop
    /// Time limit enforced at the bottom
    while (true) {
        /// Sort individuals based on OF
        std::sort(population.begin(), population.end());

        /// MUTANTS (Random solutions)
        for (usize i = E; i < (E + R); ++i) {
            if (INIT == 0) {
                csr.constructiveNSloanMGPS(realZeroOneInclusive());
            } else {
                csr.msWConstrutiveSM({});
            }

            encoder(csr);
            csr.evaluateProfile();

            next_population[i].profile = csr.profile;
            std::swap(next_population[i].labels, csr.labels);
            std::swap(next_population[i].random_keys, csr.random_keys);
        }

        /// MATING
        for (usize i = (E + R); i < P; ++i) {
            // Chose biased mates
            const usize parent1 = distributionElite(getMT());
            const usize parent2 = distributionPop(getMT());

            for (usize k = 0; k < csr.m; ++k) {
                /// Parametrized uniform crossover
                if (realRK() < PROB) {
                    csr.random_keys[k] = population[parent1].random_keys[k];
                } else {
                    csr.random_keys[k] = population[parent2].random_keys[k];
                }
            }

            decoder(csr);
            csr.evaluateProfile();

            next_population[i].profile = csr.profile;
            std::swap(next_population[i].random_keys, csr.random_keys);
            std::swap(next_population[i].labels, csr.labels);
        }

        /// ELITISM 
        std::move(population.begin(), population.begin() + E, next_population.begin());

        /// Check the time limit
        if (std::chrono::steady_clock::now() - start > csr.max_time) {
            sort(next_population.begin(), next_population.end());
            const usize current_profile = static_cast<usize>(next_population[0].profile);
            if (current_profile < csr.best_profile)
                csr.best_profile = current_profile;
            csr.profile = csr.best_profile;
            return;
        }

        /// EVOLVE (swap pointers)
        std::swap(population, next_population);
    }
}

// Initial population
void initPopulation(CSR& csr, const usize INIT, std::vector<SolutionRK>& population, const usize N) {
    encoder(csr);
    population[0].profile = csr.profile;
    population[0].labels = csr.labels;
    population[0].random_keys = csr.random_keys;

    /// Sloan-MGPS
    csr.sloanMGPS();
    encoder(csr);
    csr.evaluateProfile();
    population[1].profile = csr.profile;
    population[1].labels = csr.labels;
    population[1].random_keys = csr.random_keys;
    csr.sloanMGPS();
    encoder(csr);
    csr.evaluateProfile();
    population[2].profile = csr.profile;
    population[2].labels = csr.labels;
    population[2].random_keys = csr.random_keys;
    csr.sloanMGPS();
    encoder(csr);
    csr.evaluateProfile();
    population[3].profile = csr.profile;
    population[3].labels = csr.labels;
    population[3].random_keys = csr.random_keys;

    /// ML1W-SM (SloanMGPS label the coarsest graph)
    const usize algo_base = 0;
    csr.msW({}, algo_base);
    encoder(csr);
    csr.evaluateProfile();
    population[4].profile = csr.profile;
    population[4].labels = csr.labels;
    population[4].random_keys = csr.random_keys;
    csr.msW({}, algo_base);
    encoder(csr);
    csr.evaluateProfile();
    population[5].profile = csr.profile;
    population[5].labels = csr.labels;
    population[5].random_keys = csr.random_keys;
    csr.msW({}, algo_base);
    encoder(csr);
    csr.evaluateProfile();
    population[6].profile = csr.profile;
    population[6].labels = csr.labels;
    population[6].random_keys = csr.random_keys;

    /// MPG
    csr.mpg();
    encoder(csr);
    csr.evaluateProfile();
    population[7].profile = csr.profile;
    population[7].labels = csr.labels;
    population[7].random_keys = csr.random_keys;

    for (usize i = 8; i < N; ++i) {
        if (INIT == 0) {
            csr.constructiveNSloanMGPS(realZeroOneInclusive());
        } else {
            csr.msWConstrutiveSM({});
        }

        encoder(csr);
        csr.evaluateProfile();

        population[i].profile = csr.profile;
        population[i].random_keys = csr.random_keys;
        population[i].labels = csr.labels;
    }
}