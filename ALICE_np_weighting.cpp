//
// Created by rwagn on 27.08.2019.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
#include "ALICE_input.h"
#include <vector>
#include <cmath>

using namespace std;

void read_alice_data_x(const string &file, vector<double> &x_axis);

int main() {

    vector<double> x_axis;
    double cms = 2.76;

    Input input;
    input.full_data();
    input.set_cms(cms);
    int particle_id = -2212;
    string particle_string = "pbar";

    double N_events_total = 0;

    for (int j = 0; j < input.N_simulations; ++j) {
        N_events_total += input.N_events[j];
    }

    histogram H_pbar;
    read_alice_data_x(
            "/mnt/d/Uni/Lectures/thesis/ALICE_data/ALICE_" + particle_string + "_data_" + input.cms_string + ".dat",
            x_axis);
    H_pbar.set_xaxis(x_axis);
    H_pbar.rescale_data(1. / (2 * M_PI * N_events_total)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    H_pbar.rescale_data(x_axis); //bin individual rescale factors
    H_pbar.normalize(); //Normalize per bin width

    auto *S = new mySimulation[input.N_simulations];

    //analyze events

    //read in data from pythia
    for (int l = 0; l < input.N_simulations; ++l) {
        S[l].load_txt(input.dataset_folder + input.files[l], input.N_events[l], particle_id);
    }

    //fill Histogram
    for (int i = 0; i < input.N_simulations; ++i) {
        for (int k = 0; k < input.N_events[i]; ++k) {
            for (auto &pbar:S[i].Particles[k][0]) {
                if (abs(pbar.y()) <= 0.5) {
                    H_pbar.fill(pbar.pT());
                }
            }
        }
    }

    //print histogram
    H_pbar.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/" +
                 particle_string + "_data_" + input.cms_string + ".txt");

    return 0;
}

void read_alice_data_x(const string &file, vector<double> &x_axis) {

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
                x_axis.push_back(stod(token));
                break;
            }
        }
    }

    myfile.close();
}