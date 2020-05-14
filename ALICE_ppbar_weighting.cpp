//
// Created by rwagn on 20.04.2020.
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
    double cms = 7;

    Input input;
    input.full_data();
    input.set_cms(cms);
    input.max_simulations(5);
    int particle_id = 2212;
    string particle_string = "p";

    histogram H_pbar;
    read_alice_data_x("/mnt/d/Uni/Lectures/thesis/ALICE_data_new/7TeVppbar/ppbar_7.txt", x_axis);

    H_pbar.set_xaxis(x_axis);
    H_pbar.rescale_data(1. / ( input.N_event_total)); //Overall normalization

    //need inverse x axis values for rescaling
    H_pbar.normalize(); //Normalize per bin width

    auto *S = new mySimulation[input.N_simulations];

    //analyze events

    //read in data from pythia
    for (int l = 0; l < input.N_simulations; ++l) {
        S[l].load_txt(input.dataset_folder + input.files[l], input.N_events[l], particle_id);
    }

    //fill Histogram
    for (int i = 0; i < input.N_simulations; ++i) {
        for (auto& event:S[i].Event) {
            for (auto &pbar:event.protons) {
                if (abs(pbar.y()) <= 0.5) {
                    H_pbar.fill(pbar.pT());
                }
            }
        }
    }

    //print histogram
    H_pbar.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_new/7TeVppbar/" +
                 particle_string + "_pythia_" + input.cms_string + ".txt");

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