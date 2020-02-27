#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include "Frontend.h"
#include "diver.hpp"


const int nPar = 1;                            // Dimensionality of the parameter space
const double lowerbounds[nPar] = {0.2};                 // Lower boundaries of parameter space
const double upperbounds[nPar] = {0.22};                 // Upper boundaries of parameter space
const char path[] = "/mnt/d/Uni/Lectures/thesis/diver_output/example"; // Path to save samples, resume files, etc
const int nDerived = 0;                            // Number of derived quantities to output
const int nDiscrete = 0;                            // Number of parameters that are to be treated as discrete
const int discrete[] = {
        none};                       // Indices of discrete parameters, Fortran style, i.e. starting at 1!!
const bool partitionDiscrete = false;                        // Split the population evenly amongst discrete parameters and evolve separately
const int maxciv = 1;                          // Maximum number of civilisations
const int maxgen = 100;                          // Maximum number of generations per civilisation
const int NP = 100;                          // Population size (individuals per generation)
const int nF = 1;                            // Size of the array indicating scale factors
const double F[nF] = {0.6};                        // Scale factor(s).  Note that this must be entered as an array.
const double Cr = 0.9;                          // Crossover factor
const double lambda = 0.8;                          // Mixing factor between best and rand/current
const bool current = false;                        // Use current vector for mutation
const bool expon = false;                        // Use exponential crossover
const int bndry = 3;                            // Boundary constraint: 1=brick wall, 2=random re-initialization, 3=reflection
const bool jDE = true;                         // Use self-adaptive choices for rand/1/bin parameters as per Brest et al 2006
const bool lambdajDE = true;                         // Use self-adaptive rand-to-best/1/bin parameters; based on Brest et al 2006
const double convthresh = 1.e-3;                        // Threshold for gen-level convergence: smoothed fractional improvement in the mean population value
const int convsteps = 10;                           // Number of steps to smooth over when checking convergence
const bool removeDuplicates = true;                         // Weed out duplicate vectors within a single generation
const bool doBayesian = true;                         // Calculate approximate log evidence and posterior weightings
const double maxNodePop = 1.9;                          // Population at which node is partitioned in binary space partitioning for posterior
const double Ztolerance = 0.1;                          // Input tolerance in log-evidence
const int savecount = 1;                            // Save progress every savecount generations
const bool resume = false;                        // Restart from a previous run
const bool outputSamples = true;                         // Write output .raw and .sam (if nDerived != 0) files
const int init_pop_strategy = 2;                            // Initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
const bool discard_unfit_points = false;                        // Recalculate any trial vector whose fitness is above max_acceptable_value
const int max_init_attempts = 10000;                        // Maximum number of times to try to find a valid vector for each slot in the initial population.
const double max_acceptable_val = 1e6;                          // Maximum fitness to accept for the initial generation if init_population_strategy > 0, or any generation if discard_unfit_points = true.
const int seed = -1;                           // base seed for random number generation; non-positive or absent means seed from the system clock
const int verbose = 1;                            // Output verbosity: 0=only error messages, 1=basic info, 2=civ-level info, 3+=population info

Frontend frontend;

void construct_setup() {
    frontend.input_settings("dbar", false, true, 7.);
   // frontend.max_simulations(45);
    frontend.set_H_data();
    frontend.set_H_sim();
    frontend.load_pythia_data();
}

const double Pi = 3.14159265359;                                            // Tasty
typedef double (*likelihood)(double[], const int, int &, bool &,
                             const bool); // This example's internal standard likelihood function signature

double my_minus_lnlike(double params[], const int param_dim, int &fcall,
                       bool &quit, const bool validvector) {
    double result = 0;
    result = frontend.chi_sqaured(params);
    // What to do if the vector of parameters were not
    // within our parameter bounds? Just return the largest
    // possible number
    if (not validvector) { result = std::numeric_limits<double>::max(); }
    // Increment the function call counter that Diver provided;
    fcall += 1;
    // Tell Diver that everything went fine, so no need to quit.
    quit = false;
    return result;
}


//Function to be minimized.  Corresponds to -ln(Likelihood).  Redirects to the target of context pointer.
double
objective(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void *&context) {
    likelihood *like_ptr = static_cast<likelihood *>(context);
    likelihood minus_lnlike = *like_ptr;
    return minus_lnlike(params, param_dim, fcall, quit, validvector);
}

//Plain Gaussian likelihood centred at the origin, good for any number of dimensions.
double gauss(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector) {
    double result = 0.0;
    for (int i = 0; i < param_dim; i++) result += params[i] * params[i];
    result += 0.5 * nPar * log(Pi);
    if (not validvector) result = std::numeric_limits<double>::max();
    fcall += 1;
    quit = false;
    return result;
}


//Gaussian shell likelihood, good for any number of dimensions (just remember to expand the subarrays in c).
double gauss_shell(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector) {
    double result, temp, dist, loclike;
    int i, j;
    double *greater, lesser;
    const int nRings = 1;                                      // Number of rings
    const double w[nRings] = {0.1};                        // Gaussian widths of the shells
    const double r[nRings] = {2.0};                        // Widths of the rings
    const double c[nRings][nPar] = {{0}}; // Positions of ring centres

    result = -std::numeric_limits<double>::max() * 1e-5;
    for (i = 0; i < nRings; i++) {
        temp = 0.0;
        for (j = 0; j < nPar; j++) temp += (params[j] - c[i][j]) * (params[j] - c[i][j]);
        dist = pow(pow(temp, 0.5) - r[i], 2);
        loclike = -dist / (2.0 * w[i] * w[i]) - 0.5 * log(2.0 * Pi * w[i] * w[i]);
        if (result > loclike) result = result + log(1.0 + exp(loclike - result));
        else result = loclike + log(1.0 + exp(result - loclike));
    }
    if (not validvector) result = std::numeric_limits<double>::max();
    fcall += 1;
    quit = false;
    return -result;
}

//Flat prior function
double flatprior(const double real_params[], const int real_param_dim, void *&context) {
    int result = 1.0;
    for (int i = 0; i < real_param_dim; i++) result *= upperbounds[i] - lowerbounds[i];
    return 1.0 / result;
}

//Log prior function.  Remember it won't work if any of lowerbounds are <= 0!
double logprior(const double real_params[], const int real_param_dim, void *&context) {
    int result = 1.0;
    for (int i = 0; i < real_param_dim; i++) result /= real_params[i] * log(upperbounds[i] / lowerbounds[i]);
    return result;
}


int main(int argc, char **argv) {


    //Scan the shell likelihood if 'shell' is given as the first command-line argument, gauss if not (illustrates use of the context pointer).
    //likelihood minus_lnlike;
    //if (argc > 1 and strcmp(argv[1], "shell") == 0) minus_lnlike = gauss_shell; else minus_lnlike = &gauss;
    //void *context = &minus_lnlike;

    // Save the address of the my_minus_lnlike function as the void pointer 'context'.
    // In the 'objective' function, this pointer will later be cast back
    // to the correct function pointer type (see typedef 'likelihood') and
    // then called.
    construct_setup();
    likelihood minus_lnlike;
    minus_lnlike = my_minus_lnlike;
    void *context = &minus_lnlike;

    cdiver(objective, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, partitionDiscrete,
           maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE, convthresh, convsteps,
           removeDuplicates, doBayesian, flatprior, maxNodePop, Ztolerance, savecount, resume, outputSamples,
           init_pop_strategy, discard_unfit_points, max_init_attempts, max_acceptable_val, seed, context, verbose);

    return 0;
}