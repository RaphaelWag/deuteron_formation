//
// Created by rwagn on 07.11.2019.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
//#include <fstream>
//#include <sstream>
#include <vector>
#include <cmath>
#include "ALICE_input.h"


using namespace std;

void read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error);

int main() {

    vector<double> x_axis;
    vector<double> data;
    vector<double> error;
    histogram H_dw;
    histogram H_nw;
    histogram H_cw;


    double cms = 7;
    double cutoff_momentum_nw = 0.203;
    double cutoff_momentum_dw = 0.211;
    double cutoff_momentum_cw = 0.211;
    string particle = "dbar";
    Input input;
    input.set_particle(particle, false);
    input.reduce_data();
    input.set_cms(cms);
    x_axis.clear();
    error.clear();
    data.clear();

    cout << input.cms_string << endl;

    read_alice_data(input.ALICE_data, x_axis, data, error);

    H_dw.set_xaxis(x_axis);
    H_dw.rescale_data(input.ff / (2 * M_PI * input.N_event_total)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    H_dw.rescale_data(x_axis); //bin individual rescale factors
    H_dw.normalize(); //Normalize per bin width
    H_nw = H_cw = H_dw; //normalize and set x_axis equal to all histograms

    auto *S = new mySimulation[input.N_simulations];

    for (int k = 0; k < input.N_simulations; ++k) {
        S[k].set_ALICE_weights_discrete(input.cms_string, input.particle_type);
        S[k].set_ALICE_weights_lt(input.cms_string, input.particle_type);
    }

    //analyze events

    //read in data from pythia and set weights for the events
    for (int l = 0; l < input.N_simulations; ++l) {
        S[l].load_txt(input.dataset_folder + input.files[l], input.N_events[l]);
        S[l].rescale_spectrum();

        S[l].set_cutoff_momentum(cutoff_momentum_nw);
        S[l].coalescence();
        for (auto &dbar:S[l].deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_nw.fill(dbar.pT());
            }
        }

        S[l].set_cutoff_momentum(cutoff_momentum_dw);
        S[l].coalescence();
        for (auto &dbar:S[l].deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_dw.fill(dbar.pT(),dbar.wA());
            }
        }

        S[l].set_cutoff_momentum(cutoff_momentum_cw);
        S[l].coalescence();
        for (auto &dbar:S[l].deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_cw.fill(dbar.pT(),dbar.wLT());
            }
        }

    }
    H_dw.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_dw.txt");
    H_cw.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_cw.txt");
    H_nw.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_nw.txt");

    return 0;
}

void read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error) {

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
