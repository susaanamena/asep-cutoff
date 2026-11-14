#include "ASEPSimulation.hpp"

// Constructor
ASEPSimulation::ASEPSimulation(int N_val, double alpha_val, double p_val) 
    : N(N_val), 
      alpha(alpha_val), 
      k(std::ceil(N_val * alpha_val)),
      p(p_val), 
      q(1-p_val), 
      gamma(p - q),
      gen(std::random_device{}()),
      bond_dist(0, N_val - 2),
      prob_dist(0.0, 1.0),
      time_dist(N_val)
      
{

}

// Initialize configuration
std::vector<int> ASEPSimulation::init_config(bool is_min) const {
    std::vector<int> config(N, 0);
    if (is_min) {
        // Minimal configuration: particles at leftmost positions
        for (int i = 0; i < k; ++i) {
            config[i] = 1;
        }
    } else {
        // Maximal configuration: particles at rightmost positions
        for (int i = N - k; i < N; ++i) {
            config[i] = 1;
        }
    }
    return config;
}

// Scaling functions
double ASEPSimulation::g(double c){
    return (pow(sqrt(k)+sqrt(N-k),2)+c*pow(N,1/3.))/gamma;
}

double ASEPSimulation::h(double x){  
    return (x-(pow(sqrt(k)+sqrt(N-k),2))/gamma)/(pow(N,1/3.)/gamma);
}

// -----------------------------------------------------------------------------------------------

// Core simulation methods

double ASEPSimulation::simulate_step(std::vector<int>& config, int& L, std::mt19937& local_gen) {
    // Select random bond and generate exponential time increment
    int bond = bond_dist(local_gen);
    double dt = time_dist(local_gen); 
    
    // Check bond configuration and attempt particle jump
    if (config[bond] == 0 && config[bond+1] == 1) {
        // Particle attempts to move left (probability q)
        if (prob_dist(local_gen) < q) {
            std::swap(config[bond], config[bond+1]);
            if (L == bond+1) {
                L--; // Update leftmost particle position
            }
        }
    } else if (config[bond] == 1 && config[bond+1] == 0) {
        // Particle attempts to move right (probability p)
        if (prob_dist(local_gen) < p) {
            std::swap(config[bond], config[bond+1]);
            if (L == bond) {
                L++; // Update leftmost particle position
            }
        }
    }
    return dt;
}

double ASEPSimulation::simulate_single_ht(std::mt19937& local_gen, double max_time) {
    // Initialize minimal configuration 
    std::vector<int> config = init_config(true);
    int L = 0;  // Position of leftmost particle
    double t = 0.0;  // Accumulated time
    
    // Simulate until reaching maximal configuration or timeout
    while (L != N - k && t < max_time) {
        auto dt=simulate_step(config, L, local_gen);
        t += dt;
    }
    return (t < max_time) ? t : std::numeric_limits<double>::infinity();
}

// -----------------------------------------------------------------------------------------------

// Upper bound: hitting time analysis

void ASEPSimulation::generate_ht_data(int trials) {
    // Open output file for hitting times
    std::string filename = "hittingtimes" + std::to_string(N) + ".dat";
    FILE* f = fopen(filename.c_str(), "a");
    if (!f) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Execute specified number of trials
    for (int i = 0; i < trials; ++i) {
        // Generate seed for each trial
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + i;
        std::mt19937 local_gen(seed);

        double hitting_time = simulate_single_ht(local_gen);
        fprintf(f, "%f\n", hitting_time);
    }
    fclose(f);
}

void ASEPSimulation::compute_upper_bound_prob(){
    std::string filename = "hittingtimes"+ std::to_string(N)+".dat"; 
    std::string outfilename = "ht_prob" + std::to_string(N)+ ".dat";
    
    std::ifstream infile(filename);
    std::ofstream outfile(outfilename);
     
    // Read and scale hitting times with h(x)
    std::vector<double> values;
    double num;
    while (infile >> num) {
        values.push_back(h(num)); 
    }
    infile.close();

    int total = values.size();

    // Calculate P(Z > c) for c in [-11.9, 4] with step 0.1
    for (double c = -11.9; c <= 4.0 + 1e-9; c += 0.1) {
        int counter = 0;
        for (double v : values) {
            if (v > c) counter++;
        }
        double prob = static_cast<double>(counter) / total;
        outfile << std::setprecision(6) << prob << std::endl;
    }
    outfile.close();
}

// -----------------------------------------------------------------------------------------------

// Lower bound: event A_N analysis

void ASEPSimulation::compute_stationary_gaps(int trials){
    std::string filename = "stationary" + std::to_string(N) + ".dat";
    FILE* f = fopen(filename.c_str(), "a");
    if (!f) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    // Simulate for twice the cutoff time
    double simulation_time = 2*ASEPSimulation::g(0);

    for (int i=0; i < trials; ++i){
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + i;
        std::mt19937 local_gen(seed);
        
        // Start from maximal configuration
        std::vector<int> config = init_config(false);
        int L=N-k; 
        double t=0; 

        // Simulate until equilibrium
        while (t<simulation_time){
            auto dt=simulate_step(config, L, local_gen);
            t+=dt;
        }
    fprintf(f, "%d\n", N-k-L); 
    }
    fclose(f);

}


//Lower bound probability calculations
void ASEPSimulation::generate_AN_data (int trials, int l_0){
    std::string filename = "prob" + std::to_string(N) + "_L"+ std::to_string(l_0) + ".dat";
    FILE* f = fopen(filename.c_str(), "a");
    if (!f) {
        std::cerr << "Error: Cannot open " << filename  << std::endl;
        return;
    }

    double g_k_0 = (pow(sqrt(k)+sqrt(N-k),2))/gamma;
    double scaling_fact = pow(N,1/3.)/gamma;

    for (int i=0; i < trials; ++i){
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + i;
        std::mt19937 local_gen(seed);

        // Start from minimal configuration
        std::vector<int>config = init_config(true);
        int L=0; 
        double t=0.0; 

        // Evaluate A_N(l_0) for each c value
        for (double c=-11.9;c<4.0;c+=0.1){
            double target_time = g_k_0 + c*scaling_fact; 
            
            // Simulate until target time and check event A_N(l_0)
            bool A_N = (l_0 + 1 < N - k - L); 
            while (t < target_time) {
                A_N = (l_0 + 1 < N - k - L); 
                double dt = simulate_step(config, L, local_gen);
                t += dt;
                }
            fprintf(f, "%d", A_N ? 1 : 0);      
        }
        fprintf(f, "\n");
    }
    fclose(f);    
}


void ASEPSimulation::compute_lower_bound_prob() {
    // Process probability files for l_0 in {1,2,...,9}
    std::vector<int> y = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    for(size_t fileIndex = 0; fileIndex < y.size(); fileIndex++) {
        std::string filename = "prob" + std::to_string(N)+ "_L" + std::to_string(y[fileIndex]) + ".dat";
        std::string outfilename = "rel.probs" + std::to_string(N)+ "_L" + std::to_string(y[fileIndex]) + ".dat";
        
        std::ifstream infile(filename);
        std::ofstream outfile(outfilename);
        
        if(!infile.is_open()) {
            std::cout << "Error: Cannot open" << filename << std::endl;
            continue;
        }
        // Read all simulation results
        std::vector<std::string> lines;
        std::string line;
        while(std::getline(infile, line)) {
            lines.push_back(line);
        }
        infile.close();
        
        if(lines.empty()) continue;
        
        int num_c_values = lines[0].length();
        int num_trials = lines.size();
        
        // Calculate empirical probabilities for each c value
        for(int col = 0; col < num_c_values; col++) {
            int count_ones = 0;
            for(int row = 0; row < num_trials; row++) {
                if(lines[row][col] == '1') {
                    count_ones++;
                }
            }
            double prob = (double)count_ones / num_trials;
            outfile << std::fixed << std::setprecision(6) << prob << std::endl;
        }
        
        outfile.close();
    }
}
