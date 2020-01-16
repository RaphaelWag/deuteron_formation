//
// Created by rwagn on 03.12.2019.
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
print_results_txt(double **chi_squared, double *cutoff_momentum, int N, string *class_string,
                  const string &particle_type);

int main() {

    const int N_classes = 5;
    vector<double> x_axis;
    vector<double> data;
    vector<double> error;
    vector<double> counter(N_classes,0);
    vector<histogram> alice(N_classes);
    vector<histogram> pythia(N_classes);

    double cms = 7;
    Input input;
    input.reduce_data();
    input.set_particle("dbar", false);
    input.set_cms(cms);

    int N_values = 30;
    auto **chi_squared = new double *[N_classes];
    for (int i = 0; i < N_classes; ++i) { chi_squared[i] = new double[N_values]; }

    double cutoff_momentum[N_values];
    for (int n = 0; n < N_values; ++n) {
        cutoff_momentum[n] = (190. + n) / 1000.;
    }
    string class_string[N_classes] = {"_0_5_", "_5_10_", "_10_20_", "_20_40_", "_40_100_"};
    double rescale_factors[N_classes] = {16. / input.N_event_total/0.05, 8. / input.N_event_total/0.05, 4. / input.N_event_total/0.1,
                                         2. / input.N_event_total/0.2, 1. / input.N_event_total/0.6};

    double class_bin_limits[N_classes + 1] = {10000000., 14.88372093, 12.55813953, 9.53488372, 6.04651163, 0};

    string error_string = "2.txt";
    string path = "/mnt/d/Uni/Lectures/thesis/ALICE_multiplicity_data/V0M_";


    x_axis.clear();
    error.clear();
    data.clear();

    for (int j = 0; j < N_classes; ++j) {
        //read in alice data
        read_alice_data(path + input.particle + class_string[j] + error_string, x_axis, data, error);

        //set alice data in histograms
        alice[j].set_xaxis(x_axis);
        alice[j].set_data(data);
        alice[j].set_error(error);
        alice[j].not_scale();


        //set rescale factors for pythia data
        pythia[j].rescale_data(rescale_factors[j]);
        //set x axis for pythia
        pythia[j].set_xaxis(x_axis);
        pythia[j].normalize(); //differential in pT

        x_axis.clear();
        error.clear();
        data.clear();
    }

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
    }

    //loop over cutoff values
    for (int i = 0; i < N_values; ++i) {
        for (int m = 0; m < input.N_simulations; ++m) {
            S[m].set_cutoff_momentum(cutoff_momentum[i]);
            S[m].coalescence();

            //fill data in Histograms
            for (auto &dbar:S[m].deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    for (int j = 0; j < N_classes; ++j) {
                        if (dbar.Nch_eta() > class_bin_limits[j]) { continue; }
                        if (dbar.Nch_eta() < class_bin_limits[j + 1]) { continue; }
                        pythia[j].fill(dbar.pT(), dbar.wLT());
                        counter[j]++;
                    }
                }
            }
        }

        //rescale
        for (auto &hist:pythia) { hist.rescale(); }

        cout << cutoff_momentum[i] << " ";
        //calculate chi squared and safe in array
        for (int k = 0; k < N_classes; ++k) {
            chi_squared[k][i] = pythia[k].chi_squared_raw(alice[k]);
            cout << chi_squared[k][i] / (pythia[k].N_bins() - 1.) << " " ;//<< counter[k] << " ";
            counter[k] =0;
        }
        //reset Histograms
        for (auto &hist:pythia) { hist.reset(); }
        cout << endl;

    }

    //print chi squared values to plot in python
    print_results_txt(chi_squared, cutoff_momentum, N_values, class_string, input.particle);
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
print_results_txt(double **chi_squared, double *cutoff_momentum, int N, string *class_string,
                  const string &particle_type) {

    string path;
    int N_class = 5;
    for (int i = 0; i < N_class; ++i) {

        path = "/mnt/d/Uni/Lectures/thesis/ALICE_results/multiplicity/chi_squared/chi_squared_" + class_string[i] +
               "_" +
               particle_type + ".txt";

        ofstream myfile(path);
        if (!myfile.is_open())
            cout << "Unable to open file";
        else {
            for (int k = 0; k < N; k++) {
                myfile << cutoff_momentum[k] << " " << chi_squared[i][k] << endl;
            }
        }
        myfile.close();

    }
}