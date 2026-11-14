#include "ASEPSimulation.hpp"
#include <ctime>

int main(int argc, char* argv[]) {
    // Default simulation parameters
    int N = 100;
    double alpha = 0.5;
    double p = 0.75;
    int trials = 10000;
    
    // Command-line argument passing
    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (i + 1 < argc) {
            if (arg == "--N" || arg == "-N") {
                N = std::stoi(argv[i + 1]);
            } else if (arg == "--alpha" || arg == "-a") {
                alpha = std::stod(argv[i + 1]);
            } else if (arg == "--p") {
                p = std::stod(argv[i + 1]);
            } else if (arg == "--trials" || arg == "-t") {
                trials = std::stoi(argv[i + 1]);
            } 
        }
    }
    // Simulation setup
    ASEPSimulation sim(N, alpha, p);   
    
    // Display simulation parameters
    std::cout << "=== ASEP simulation ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  N = " << N << std::endl;
    std::cout << "  k = " << std::ceil(N * alpha) << std::endl;
    std::cout << "  alpha = " << alpha << std::endl;
    std::cout << "  p = " << p << std::endl;
    std::cout << "  trials = " << trials << std::endl;
    std::cout << std::endl;
    
    // Execute simulations
    sim.generate_ht_data(trials);
    sim.compute_upper_bound_prob();
    sim.compute_stationary_gaps(trials);
    for (int l_0 = 0; l_0 < 10; ++l_0) {
        sim.generate_AN_data(trials,l_0);
    }
    sim.compute_lower_bound_prob();

    return 0;
}