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


    double cms = 7;
    double cutoff_momentum_nw = 0.204; // 205, 197, 204, 191
    double cutoff_momentum_cw = 0.196; // 222 , 197 , 196. 216
    double cutoff_momentum_comb = 0.196;
    double cutoff_momentum_comb_w = 0.203;
    string particle = "dbar";
    Input input;
    input.set_particle(particle, false);
    input.reduce_data();
    input.set_cms(cms);

    //vector<double> class_bin_limits_weights = {23255, 23.255814, 17.9069767, 14.8837209, 13.0232558,11.3953488, 8.8372093,
    //                                          6.97674419, 5.34883721, 3.48837209, 0};
    //string class_string_weights[10] = {"i", "ii", "iii", "iv","v", "vi", "vii", "iix", "ix", "x"};


    x_axis.clear();
    error.clear();
    data.clear();

    cout << input.cms_string << endl;

    read_alice_data(input.ALICE_data, x_axis, data, error);

    H_cw.set_xaxis(x_axis);
    H_cw.rescale_data(input.ff*input.nsdtoinel / (2* M_PI *input.N_event_total)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    H_cw.rescale_data(x_axis); //bin individual rescale factors
    H_cw.normalize(); //Normalize per bin width
    H_nw = H_comb = H_comb_w = H_cw; //normalize and set x_axis equal to all histograms

    mySimulation S;


    S.set_ALICE_weights_lt(input.cms_string, input.particle_type);
    //S.set_ALICE_weights_multi(input.cms_string, input.particle_type, class_string_weights, 10);


    //analyze events

    //read in data from pythia and set weights for the events
    for (int l = 0; l < input.N_simulations; ++l) {
        S.load_txt(input.dataset_folder + input.files[l], input.N_events[l],true);
        S.rescale_spectrum();
        //S.rescale_spectrum_multi(class_bin_limits_weights);

        S.set_cutoff_momentum(cutoff_momentum_nw);
        S.coalescence();
        for (auto &dbar:S.deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_nw.fill(dbar.pT());
            }
        }

        S.set_cutoff_momentum(cutoff_momentum_cw);
        S.coalescence();
        for (auto &dbar:S.deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_cw.fill(dbar.pT(),dbar.wLT());
            }
        }

        S.set_cutoff_momentum(cutoff_momentum_comb);
        S.coalescence();
        for (auto &dbar:S.deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_comb.fill(dbar.pT());
            }
        }

        S.set_cutoff_momentum(cutoff_momentum_comb_w);
        S.coalescence();
        for (auto &dbar:S.deuteron) {
            if (abs(dbar.y()) <= 0.5) {
                H_comb_w.fill(dbar.pT(),dbar.wLT());
            }
        }

    }

    H_cw.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_cw.txt");
    H_nw.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_nw.txt");
    H_comb.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_comb.txt");
    H_comb_w.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/"+particle+"_data_"+input.cms_string+"_comb_w.txt");

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
