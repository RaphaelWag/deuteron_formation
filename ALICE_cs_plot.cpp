//
// Created by rwagn on 01.02.2020.
//

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
    histogram H_nw;
    histogram H_cw;
    histogram H_comb;
    histogram H_comb_w;


    double cms = 13;
    double simga_inv_nw = 3.3 * 1e-6; // 5.3, 4.0, 3.6, 3.1
    double sigma_inv_cw = 3.9 * 1e-6; // 5.1, 3.4, 3.1, 3.9
    double sigma_inv_comb = 3.4 * 1e-6;
    double sigma_inv_comb_w = 3.4 * 1e-6;
    string particle = "dbar";
    Input input;
    input.set_particle(particle, false);
    input.cs_data();
    input.set_cms(cms);
    x_axis.clear();
    error.clear();
    data.clear();
    int mc_cycles = 10;
    cout << input.cms_string << endl;

    read_alice_data(input.ALICE_data, x_axis, data, error);

    H_cw.set_xaxis(x_axis);
    H_cw.rescale_data(
            input.ff * input.nsdtoinel / (/*2 * M_PI */ input.N_event_total * mc_cycles)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    //H_cw.rescale_data(x_axis); //bin individual rescale factors
    H_cw.normalize(); //Normalize per bin width
    H_nw = H_comb = H_comb_w = H_cw; //normalize and set x_axis equal to all histograms

    mySimulation S;

    for (int k = 0; k < input.N_simulations; ++k) {
        S.set_ALICE_weights_lt(input.cms_string, input.particle_type);
    }

    //analyze events

    //read in data from pythia and set weights for the events
    for (int l = 0; l < input.N_simulations; ++l) {
        S.load_txt(input.dataset_folder + input.files[l], input.N_events[l], true);
        S.rescale_spectrum();

        for (int i = 0; i < mc_cycles; ++i) {
            S.set_sigma0(simga_inv_nw);
            S.cs_model_formation();
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    H_nw.fill(dbar.pT());
                }
            }
        }
        for (int i = 0; i < mc_cycles; ++i) {
            S.set_sigma0(sigma_inv_cw);
            S.cs_model_formation();
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    H_cw.fill(dbar.pT(), dbar.wLT());
                }
            }
        }
        for (int i = 0; i < mc_cycles; ++i) {
            S.set_sigma0(sigma_inv_comb);
            S.cs_model_formation();
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    H_comb.fill(dbar.pT());
                }
            }
        }
        for (int i = 0; i < mc_cycles; ++i) {
            S.set_sigma0(sigma_inv_comb_w);
            S.cs_model_formation();
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    H_comb_w.fill(dbar.pT(), dbar.wLT());
                }
            }
        }

    }

    H_cw.print(
            "/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/cs_model/" + particle + "_data_" + input.cms_string + "_cw.txt");
    H_nw.print(
            "/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/cs_model/" + particle + "_data_" + input.cms_string + "_nw.txt");
    H_comb.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/cs_model/" + particle + "_data_" + input.cms_string +
                 "_comb.txt");
    H_comb_w.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/cs_model/" + particle + "_data_" + input.cms_string +
                   "_comb_w.txt");

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
