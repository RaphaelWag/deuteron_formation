//
// Created by rwagn on 28.01.2020.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
#include <vector>
#include <cmath>
#include "ALICE_input.h"


using namespace std;

void read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error);

void
print_results_txt(double *chi_squared, double *cutoff_momentum, int N, const string &cms_string, const string &wtype);

int main() {

    vector<double> x_axis;
    vector<double> data;
    vector<double> error;
    histogram H_rescaled_LT;
    histogram H_nscaled;
    histogram H_data;

    vector<double> cms_list = {2.76};
    Input input;
    input.cs_data();
    input.set_particle("dbar", false);


    int mc_cycles = 10;
    int N_values = 50;
    double sigma_0_inv[N_values];
    double chi_squared[N_values];
    double chi_squared_wLT[N_values];
    for (int n = 0; n < N_values; ++n) {
        sigma_0_inv[n] = (3.0 + n / 10.) * 1e-6;
    }

    for (auto &cms:cms_list) {
        input.set_cms(cms);
        x_axis.clear();
        error.clear();
        data.clear();

        cout << input.cms_string << endl;

        read_alice_data(input.ALICE_data, x_axis, data, error);

        H_rescaled_LT.set_xaxis(x_axis);
        H_rescaled_LT.rescale_data(input.ff / (2 * M_PI * input.N_event_total * mc_cycles)); //Overall normalization
        //need inverse x axis values for rescaling
        for (auto &values:x_axis) {
            values = 1. / values;
        }

        H_rescaled_LT.rescale_data(x_axis); //bin individual rescale factors
        H_rescaled_LT.normalize(); //Normalize per bin width
        H_nscaled = H_data = H_rescaled_LT; //normalize and set x_axis equal to all histograms

        //set data from ALICE
        H_data.not_scale();
        H_data.set_data(data);
        H_data.set_error(error);

        auto *S = new mySimulation[input.N_simulations];

        for (int k = 0; k < input.N_simulations; ++k) {
            S[k].set_ALICE_weights_discrete(input.cms_string, input.particle_type);
            S[k].set_ALICE_weights_lt(input.cms_string, input.particle_type);
        }

        //analyze events

        //read in data from pythia and set weights for the events
        for (int l = 0; l < input.N_simulations; ++l) {
            S[l].load_txt(input.dataset_folder + input.files[l], input.N_events[l], true);
            S[l].rescale_spectrum();
        }

        //loop over cutoff values
        for (int i = 0; i < N_values; ++i) {
            for (int m = 0; m < input.N_simulations; ++m) {
                for (int j = 0; j < mc_cycles; ++j) {
                    S[m].set_sigma0(sigma_0_inv[i]);
                    S[m].cs_model_formation();

                    //fill data in Histograms
                    for (auto &dbar:S[m].deuteron) {
                        if (abs(dbar.y()) <= 0.5) {
                            H_rescaled_LT.fill(dbar.pT(), dbar.wLT());
                            H_nscaled.fill(dbar.pT());
                        }
                    }
                }
            }
            //rescale
            H_rescaled_LT.rescale();
            H_nscaled.rescale();


            //calculate chi squared and safe in array
            chi_squared[i] = H_nscaled.chi_squared_raw(H_data);
            chi_squared_wLT[i] = H_rescaled_LT.chi_squared_raw(H_data);

            //reset Histograms
            H_nscaled.reset();
            H_rescaled_LT.reset();

            cout << sigma_0_inv[i] << " " << chi_squared[i] / (H_data.N_bins() - 1.) << " "
                 << chi_squared_wLT[i] / (H_data.N_bins() - 1.)
                 << endl;
            //end loop
        }

        //print chi squared values to plot in python
        print_results_txt(chi_squared, sigma_0_inv, N_values, input.cms_string, "nw"); //nw = no weights
        print_results_txt(chi_squared_wLT, sigma_0_inv, N_values, input.cms_string, "cw"); //dw = discrete weights

        delete[] S;
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

void
print_results_txt(double *chi_squared, double *cutoff_momentum, int N, const string &cms_string, const string &wtype) {

    string path;
    path = "/mnt/d/Uni/Lectures/thesis/ALICE_results/chi_squared_cs/chi_squared_" + cms_string + "_" + wtype + ".txt";

    ofstream myfile(path);
    if (!myfile.is_open())
        cout << "Unable to open file";
    else {
        for (int k = 0; k < N; k++) {
            myfile << chi_squared[k] << " " << cutoff_momentum[k] << endl;
        }
    }
    myfile.close();
}
