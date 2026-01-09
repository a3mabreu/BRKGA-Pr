#include <getopt.h>
#include <csignal>
#include "csr.hpp"
#include "brkga.hpp"

void parseArguments(int argc, char *argv[], std::map<std::string, ParamValue> &params);

void brkgaEx(CSR& csr, std::map<std::string, ParamValue> &params) {
    if (get<int>(params["pop"]) < 10) {
        std::cerr << "\nYou need to use pop > 10;\n";
        exit(EXIT_FAILURE);
    }

    if (!get<int>(params["irace"]))
        std::cout << "\tBRKGA-Pr... \n"; 

    auto start = std::chrono::steady_clock::now();

    brkga(csr, get<int>(params["pop"]), get<int>(params["elite"]), get<int>(params["mutants"]), get<realT>(params["prob"]), get<int>(params["init"]));

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    if (!get<int>(params["irace"])) {
        std::cout << "Profile: " << csr.best_profile;
        std::cout << "\ttotal time: " << duration.count() << "s" << std::endl;
    } else {
        std::cout << csr.best_profile;
    }
}

void printParams(const std::map<std::string, ParamValue>& params) {
    for (const auto& [key, value] : params) {
        std::cout << key << ": ";
        std::visit([](auto&& arg) {
            std::cout << arg;
        }, value);
        std::cout << "\n";
    }
}

// Global pointer for handleSigterm
CSR* g_csr_pointer = nullptr;
// The signal handler function for SIGTERM
void handleSigterm(int /*signum*/) {
    // Write the message to standard output (file descriptor 1)
    if (g_csr_pointer != nullptr && g_csr_pointer->best_profile != 0) {
        std::cout << g_csr_pointer->best_profile << std::endl;
    } else {
        std::cout << std::numeric_limits<unsigned long>::max() << std::endl;
    }

    _exit(143); 
}

int main(int argc, char *argv[]) {
    // Register the signal handler for SIGTERM (signal 15)
    signal(SIGTERM, handleSigterm);
    // Initialize parameters with default values
    std::map<std::string, ParamValue> params = {
        {"irace", 0},
        {"filename", std::string("input/usps_norm_5NN.mtx")},
        {"init", 1},
        {"max_time", int64_t(10)},
        {"alpha", 0.0f},
        {"pop", 20},
        {"elite", 8},
        {"mutants", 4},
        {"prob", 0.75f},
    };
    parseArguments(argc, argv, params);

    /// Print parameters
    // if (!get<int>(params["irace"])) {
    //     printParams(params);
    // }
    
    // Create the CSR if it's not an experiment or MSH
    const auto filename = get<std::string>(params["filename"]);
    CSR csr(filename, true);    
    csr.max_time = std::chrono::seconds{get<int64_t>(params["max_time"])};
    // Set the global pointer to handle SIGTERM
    g_csr_pointer = &csr;

    csr.evaluateProfile();
    if (!get<int>(params["irace"]))
        std::cout << "\nInitial Profile: " << csr.profile << '\n';

    brkgaEx(csr, params);

    // Always verify feasibility and OF
    if (!csr.isFeasible()) {
        std::cerr << "Solution is not feasible." << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

void parseArguments(int argc, char *argv[], std::map<std::string, ParamValue> &params) {
    struct option long_options[] = {
        {"algo", required_argument, nullptr, 0},
        {"alpha", required_argument, nullptr, 0},
        {"alpha_sa", required_argument, nullptr, 0},
        {"bl", required_argument, nullptr, 0},
        {"cross", required_argument, nullptr, 0},
        {"crot_s", required_argument, nullptr, 0},
        {"delta", required_argument, nullptr, 0},
        {"delta_s", required_argument, nullptr, 0},
        {"elite", required_argument, nullptr, 0},
        {"exp_out", required_argument, nullptr, 0},
        {"filename", required_argument, nullptr, 0},
        {"hamming_t", required_argument, nullptr, 0},
        {"init", required_argument, nullptr, 0},
        {"irace", required_argument, nullptr, 0},
        {"k_step", required_argument, nullptr, 0},
        {"k_max", required_argument, nullptr, 0},
        {"l_0", required_argument, nullptr, 0},
        {"limit", required_argument, nullptr, 0},
        {"l_size", required_argument, nullptr, 0},
        {"level_d", required_argument, nullptr, 0},
        {"max_per", required_argument, nullptr, 0},
        {"max_it", required_argument, nullptr, 0},
        {"max_time", required_argument, nullptr, 0},
        {"min_zeros", required_argument, nullptr, 0},
        {"mi", required_argument, nullptr, 0},
        {"mp", required_argument, nullptr, 0},
        {"mutants", required_argument, nullptr, 0},
        {"n_pass", required_argument, nullptr, 0},
        {"per", required_argument, nullptr, 0},
        {"per_it", required_argument, nullptr, 0},
        {"pop", required_argument, nullptr, 0},
        {"pool_s", required_argument, nullptr, 0},
        {"prob", required_argument, nullptr, 0},
        {"prob_mut", required_argument, nullptr, 0},
        {"pr_prop", required_argument, nullptr, 0},
        {"prob_per", required_argument, nullptr, 0},
        {"prob_nex", required_argument, nullptr, 0},
        {"prob_rex", required_argument, nullptr, 0},
        {"psi", required_argument, nullptr, 0},
        {"repair", required_argument, nullptr, 0},
        {"r_max", required_argument, nullptr, 0},
        {"t_0", required_argument, nullptr, 0},
        {"t_f", required_argument, nullptr, 0},
        {nullptr, 0, nullptr, 0} // Terminating entry
    };

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
        if (opt == 0) {
            std::string option_name = long_options[option_index].name;
            if (option_name == "algo") {
                params["algo"] = std::stoi(optarg);
            } else if (option_name == "alpha") {
                params["alpha"] = std::stof(optarg);
            } else if (option_name == "bl") {
                params["bl"] = std::stoi(optarg);
            } else if (option_name == "alpha_sa") {
                params["alpha_sa"] = std::stof(optarg);
            } else if (option_name == "cross") {
                params["cross"] = std::stoi(optarg);
            } else if (option_name == "crot_s") {
                params["crot_s"] = std::stoi(optarg);
            } else if (option_name == "delta") {
                params["delta"] = std::stoi(optarg);
            } else if (option_name == "delta_s") {
                params["delta_s"] = std::stoi(optarg);
            } else if (option_name == "elite") {
                params["elite"] = std::stoi(optarg);
            } else if (option_name == "exp_out") {
                params["exp_out"] = optarg;
            } else if (option_name == "filename") {
                params["filename"] = optarg;
            }  else if (option_name == "hamming_t") {
                params["hamming_t"] = std::stof(optarg);
            } else if (option_name == "init") {
                params["init"] = std::stoi(optarg);
            } else if (option_name == "irace") {
                params["irace"] = std::stoi(optarg);
            } else if (option_name == "k_step") {
                params["k_step"] = std::stoi(optarg);
            } else if (option_name == "k_max") {
                params["k_max"] = std::stoi(optarg);
            } else if (option_name == "l_0") {
                params["l_0"] = std::stoi(optarg);
            } else if (option_name == "limit") {
                params["limit"] = std::stoi(optarg);
            } else if (option_name == "level_d") {
                params["level_d"] = std::stoi(optarg);
            } else if (option_name == "l_size") {
                params["l_size"] = std::stoi(optarg);
            } else if (option_name == "max_per") {
                params["max_per"] = std::stoi(optarg);
            } else if (option_name == "max_it") {
                params["max_it"] = std::stoi(optarg);
            } else if (option_name == "max_time") {
                params["max_time"] = std::stoll(optarg);
            } else if (option_name == "min_zeros") {
                params["min_zeros"] = std::stoi(optarg);
            } else if (option_name == "mutants") {
                params["mutants"] = std::stoi(optarg);
            } else if (option_name == "n_pass") {
                params["n_pass"] = std::stoi(optarg);
            } else if (option_name == "per") {
                params["per"] = std::stoi(optarg);
            } else if (option_name == "per_it") {
                params["per_it"] = std::stoi(optarg);
            } else if (option_name == "pop") {
                params["pop"] = std::stoi(optarg);
            } else if (option_name == "pool_s") {
                params["pool_s"] = std::stoi(optarg);
            } else if (option_name == "pr_prop") {
                params["pr_prop"] = std::stof(optarg);
            } else if (option_name == "prob_per") {
                params["prob_per"] = std::stof(optarg);
            } else if (option_name == "prob") {
                params["prob"] = std::stof(optarg);
            } else if (option_name == "prob_mut") {
                params["prob_mut"] = std::stof(optarg);
            } else if (option_name == "prob_nex") {
                params["prob_nex"] = std::stof(optarg);
            } else if (option_name == "prob_rex") {
                params["prob_rex"] = std::stof(optarg);
            } else if (option_name == "psi") {
                params["psi"] = std::stof(optarg);
            } else if (option_name == "repair") {
                params["repair"] = std::stoi(optarg);
            } else if (option_name == "r_max") {
                params["r_max"] = std::stoi(optarg);
            } else if (option_name == "t_0") {
                params["t_0"] = std::stof(optarg);
            } else if (option_name == "t_f") {
                params["t_f"] = std::stof(optarg);
            } else if (option_name == "mi") {
                params["mi"] = std::stoi(optarg);
            } else if (option_name == "mp") {
                params["mp"] = std::stof(optarg);
            }
        }
    }
}