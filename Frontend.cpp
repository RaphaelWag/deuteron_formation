//
// Created by rwagn on 19.02.2020.
//

#include "Frontend.h"

#include <utility>

Frontend::Frontend() {
    S = new mySimulation[0];
}

void Frontend::input_settings(string particle_string, bool particle_type, bool reduced_data, double cms) {
    if (reduced_data) { input.reduce_data(); }
    input.set_particle(std::move(particle_string), particle_type);
    input.set_cms(cms);
}

void Frontend::max_simulations(int N_max) {
    input.max_simulations(N_max);
}

void Frontend::set_H_data() {

    vector<double> x_axis;
    vector<double> data;
    vector<double> error;
    read_alice_data(input.ALICE_data, x_axis, data, error);

    H_data.set_xaxis(x_axis);
    //set data from ALICE
    H_data.not_scale();
    H_data.set_data(data);
    H_data.set_error(error);

}

void
Frontend::read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error) {

    string line;
    string token;
    vector<double> temporary;

    ifstream myfile(file);

    if (!myfile.is_open()) {
        cout << "Error in: main -> read_alice_data" << endl;
        cout << "Unable to open file" << endl;
    } else {
        while (getline(myfile, line)) {
            stringstream ss(line);
            while (getline(ss, token, ' ')) {

                if (token.empty()) { continue; }
                temporary.push_back(stod(token));
            }
        }
    }
    myfile.close();

    unsigned i;
    for (unsigned j = 0; j < temporary.size() / 3; ++j) {
        i = j * 3;
        x_axis.push_back(temporary[i]);
        data.push_back(temporary[i + 1]);
        error.push_back(temporary[i + 2]);
    }

}

void Frontend::set_H_sim() {
    vector<double> x_axis;
    vector<double> data;
    vector<double> error;
    read_alice_data(input.ALICE_data, x_axis, data, error);

    H_sim.set_xaxis(x_axis);
    H_sim.rescale_data(input.ff / (2 * M_PI * input.N_event_total)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    H_sim.rescale_data(x_axis); //bin individual rescale factors
    H_sim.normalize(); //Normalize per bin width

}

void Frontend::load_pythia_data() {
    delete[]S;
    S = new mySimulation[input.N_simulations];

    //load weights
    for (int k = 0; k < input.N_simulations; ++k) {
        S[k].set_ALICE_weights_discrete(input.cms_string, input.particle_type);
        S[k].set_ALICE_weights_lt(input.cms_string, input.particle_type);
    }
    //read in data from pythia and set weights for the events
    for (int l = 0; l < input.N_simulations; ++l) {
        S[l].load_txt(input.dataset_folder + input.files[l], input.N_events[l], true);
        S[l].rescale_spectrum();
    }
}

double Frontend::chi_sqaured(double *params) {

    double chi_squared;
    for (int m = 0; m < input.N_simulations; ++m) {
        S[m].set_cutoff_momentum(params[0]);
        S[m].coalescence();

        //fill data in Histograms
        for (auto &dbar:S[m].deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_sim.fill(dbar.pT(), dbar.wLT());
            }
        }
    }
    H_sim.rescale();
    chi_squared = H_sim.chi_squared_raw(H_data)/(H_data.N_bins()-1);
    H_sim.reset();
    return chi_squared;
}
