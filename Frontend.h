//
// Created by rwagn on 19.02.2020.
//

#ifndef MYSIMULATION_FRONTEND_H
#define MYSIMULATION_FRONTEND_H

#include "mySimulation.h"
#include "histogram.h"
#include "ALICE_input.h"

class Frontend {
private:

    Input input;
    histogram H_sim;
    histogram H_data;
    mySimulation* S;

    static void read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error);

public:

    Frontend();

    void input_settings(string particle_string, bool particle_type, bool reduced_data, double cms);

    void max_simulations(int N_max);

    void set_H_data();

    void set_H_sim();

    void load_pythia_data();

    double chi_sqaured(double params[]);

};


#endif //MYSIMULATION_FRONTEND_H
