//
// Created by rwagn on 27.02.2020.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
#include <vector>
#include <cmath>
#include "ALICE_input.h"


void read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error);


int main() {

    double cms = 7;
    vector<double> x_axis;
    vector<double> data;
    vector<double> error;
    histogram H;
    histogram H_data;
    histogram H_test;

    Input input;
    input.reduce_data();
    input.set_particle("dbar", false);

    input.set_cms(cms);
    string cms_string = input.cms_string;
    x_axis.clear();
    error.clear();
    data.clear();

    cout << input.cms_string << endl;

    read_alice_data(input.ALICE_data, x_axis, data, error);

    H.set_xaxis(x_axis);
    H.rescale_data(input.ff*input.nsdtoinel/ (2 * M_PI * input.N_event_total)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    H.rescale_data(x_axis); //bin individual rescale factors
    H.normalize(); //Normalize per bin width
    H_test = H_data = H; //normalize and set x_axis equal to all histograms

    double threshold = 200;
    int Np_values = 101;
    int Nx_values = 101;
    double x_min = 1e-14;
    double x_max = 5e-12;
    double p_min = 0.180;
    double p_max = 1;
    vector<double> x_cutoff;
    vector<double> p_cutoff;

    for (int i = 0; i < Nx_values; ++i) {
        x_cutoff.push_back(x_min + i / (double(Np_values) - 1.) * (x_max - x_min));
    }

    for (int i = 0; i < Np_values; ++i) {
        p_cutoff.push_back(p_min + i / (double(Np_values) - 1.) * (p_max - p_min));
    }

    auto **Hist = new histogram *[Np_values];
    for (int i = 0; i < Np_values; ++i) {
        Hist[i] = new histogram[Nx_values];

    }
    auto **use = new bool *[Np_values];
    for (int i = 0; i < Np_values; ++i) {
        use[i] = new bool[Nx_values];

    }
    for (int j = 0; j < Np_values; ++j) {
        for (int i = 0; i < Nx_values; ++i) {
            Hist[j][i] = H;
        }
    }

    //set data from ALICE
    H_data.not_scale();
    H_data.set_data(data);
    H_data.set_error(error);

    mySimulation S;
    S.set_ALICE_weights_lt(input.cms_string, input.particle_type);

    //use one data file to see which combination of values is usefull
    S.load_txt(input.dataset_folder + input.files[0], input.N_events[0], true);
    S.rescale_spectrum();

    for (int i = 0; i < Np_values; ++i) {
        S.set_cutoff_momentum(p_cutoff[i]);
        for (int j = 0; j < Nx_values; ++j) {
            S.set_cutoff_position(x_cutoff[j]);
            S.coalescence_px();
            //fill data in Histograms
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    Hist[i][j].fill(dbar.pT(), dbar.wLT());
                }
            }
        }
    }
    double chi_squared;
    //check chi squared and throw out all useless combinations
    for (int j = 0; j < Np_values; ++j) {
        for (int i = 0; i < Nx_values; ++i) {
            H_test = Hist[j][i];
            H_test.rescale_data(input.ff / (2 * M_PI * input.N_events[0])); //Overall normalization
            H_test.rescale();
            chi_squared = H_test.chi_squared_raw_norm_err(H_data, input.norm_err_above,input.norm_err_below);
            chi_squared = chi_squared / (H_data.N_bins() - 2.);
            use[j][i] = chi_squared <= threshold;
        }
    }
    //do complete analysis with all usefull combinations
    for (int k = 1; k < input.N_simulations; ++k) {

        S.load_txt(input.dataset_folder + input.files[k], input.N_events[k], true);
        S.rescale_spectrum();

        for (int i = 0; i < Np_values; ++i) {
            S.set_cutoff_momentum(p_cutoff[i]);
            for (int j = 0; j < Nx_values; ++j) {
                S.set_cutoff_position(x_cutoff[j]);
                if (use[i][j]) {
                    S.coalescence_px();
                    //fill data in Histograms
                    for (auto &dbar:S.deuteron) {
                        if (abs(dbar.y()) <= 0.5) {
                            Hist[i][j].fill(dbar.pT(), dbar.wLT());
                        }
                    }
                }
            }
        }
    }

    string file = "/mnt/d/Uni/Lectures/thesis/ALICE_results/p_x_corr/p_x_coalescence/";
    file = file + "data_" + cms_string + ".txt";
    //print results
    ofstream myfile(file);
    if (!myfile.is_open()) {
        cout << "unable to open file in write_txt" << "\n";
    } else {
        for (int i = 0; i < Np_values; ++i) {
            for (int j = 0; j < Nx_values; ++j) {
                if (use[i][j]) {
                    Hist[i][j].rescale();

                    chi_squared = Hist[i][j].chi_squared_raw_norm_err(H_data,input.norm_err_above, input.norm_err_below) / (H_data.N_bins() - 2.);
                    myfile << p_cutoff[i] << " " << x_cutoff[j] << " " << chi_squared << endl;
                    cout << p_cutoff[i] << " " << x_cutoff[j] << " " << chi_squared << endl;
                } else {
                    myfile << p_cutoff[i] << " " << x_cutoff[j] << " " << threshold << endl;
                    cout << p_cutoff[i] << " " << x_cutoff[j] << " " << threshold << endl;
                }

            }
        }


    }
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