#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <iomanip>
#include <numeric>
#include <sstream>

// This class implements simulations for ASEP on a one-dimensional
// lattice with N sites and k particles, providing methods to compute both
// upper and lower bounds for the cutoff phenomenon analysis.


class ASEPSimulation {
private:
    // Simulation parameters
    int N;          // System size
    double alpha;   // Particle density (k = ceil(N*alpha))
    int k;          // Number of particles
    double p;       // Probability moving to the right
    double q;       // Probability moving to the left
    double gamma;   // Asymmetry parameter (p - q)
    

    // Random number generators
    std::mt19937 gen;                               // Primary random generator
    std::uniform_int_distribution<> bond_dist;      // Bond selection [0,N-2]
    std::uniform_real_distribution<> prob_dist;     // Jump probability
    std::exponential_distribution<> time_dist;      // Exponential time increments 


public:
    // Construct ASEP with specified parameters
    ASEPSimulation(int N_val, double alpha_val, double p_val);

    // Initialize particle configuration
    // Parameter: If is_min is true, create a minimal configuration; if false, maximal
    // Return: Vector representing particle positions (1=occupied, 0=empty)
    std::vector<int> init_config(bool is_min) const;

    // Time scaling function g(k,c)
    // Parameter: c (Scaling parameter)
    // Return: Scaled time value
    double g(double c);

    // Hitting time scaling function h(x)
    // Parameter: x (Hitting time value to scale)
    // Return: Scaled hitting time value
    double h(double h);

    // -----------------------------------------------------------------------------------------------

    // Core simulation methods

    // Execute single ASEP time step
    // Parameters: config (Current particle configuration, modified in-place), 
    //             L (Position of leftmost particle, modified in-place),
    //             local_gen (Random number generator for this step)
    // Return: Time increment for this step
    double simulate_step(std::vector<int>& config, int& L, std::mt19937& local_gen);

    // Simulate hitting time from minimal to maximal configuration
    // Parameters: local_gen (Random number generator),
    //             max_time (Maximum simulation time, default: 10000)
    // Return: Hitting time (\infty if max_time exceeded)
    double simulate_single_ht(std::mt19937& local_gen, double max_time = 10000.0);

    // -----------------------------------------------------------------------------------------------
    
    // Upper bound: hitting time analysis

    // Generate hitting time data for upper bound analysis
    // Parameter: trials (Number of independent simulations to perform)
    void generate_ht_data(int trials);

    // Compute upper bound probabilities P(Z > c) from hitting time data
    void compute_upper_bound_prob();

    // -----------------------------------------------------------------------------------------------

    // Lower bound: event A_N analysis

    // Compute stationary distribution gaps for parameter selection
    // Parameter: trials (Number of independent simulations)
    void compute_stationary_gaps(int trials);

    // Generate event A_N(l_0) data for lower bound analysis
    // Parameters: trials (Number of independent simulations),
    //             l_0 (Gap parameter for event A_N(l_0))
    void generate_AN_data (int trials, int l_0);

    // Compute lower bound probabilities from A_N(l_0) event data
    void compute_lower_bound_prob();
};

