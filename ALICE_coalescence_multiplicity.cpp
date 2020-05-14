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
                  const string &particle_type, const string &cms_string);

int main() {

    const int N_classes = 9;
    vector<double> x_axis;
    vector<double> data;
    vector<double> error;

    double cms = 13;
    Input input;
    input.reduce_data();
    input.set_particle("dbar", false);
    input.set_cms(cms);

    int N_values = 150;
    vector<histogram> alice(N_classes);
    auto **chi_squared = new double *[N_classes];
    auto **pythia = new histogram *[N_classes];
    for (int i = 0; i < N_classes; ++i) {
        chi_squared[i] = new double[N_values];
        pythia[i] = new histogram[N_values];
    }


    double cutoff_momentum[N_values];
    for (int n = 0; n < N_values; ++n) {
        cutoff_momentum[n] = (150. + n) / 1000.;
    }
    //string class_string[N_classes] = {"i+ii", "iii", "iv+v", "vi+vii","iix+ix+x"};
    //double cs_frac[N_classes] = {0.047, 0.048, 0.095, 0.19, 0.62};
    double cs_frac[N_classes] = {0.0092, 0.0368, 0.046, 0.092, 0.092, 0.092, 0.092, 0.185, 0.355};
    string class_string[N_classes] = {"i", "ii", "iii", "iv+v", "vi", "vii", "iix", "ix", "x"};
    string class_string_weights[10] = {"i", "ii", "iii", "iv","v", "vi", "vii", "iix", "ix", "x"};

    double rescale_factors[N_classes];

    for (int i = 0; i < N_classes; ++i) {
        rescale_factors[i] = input.ff * input.nsdtoinel / (input.N_event_total * cs_frac[i] * 0.75);
    }

    vector<double> class_bin_limits = {23255, 23.255814, 17.9069767, 14.8837209, 11.3953488, 8.8372093,
                                              6.97674419, 5.34883721, 3.48837209, 0};
    vector<double> class_bin_limits_weights = {23255, 23.255814, 17.9069767, 14.8837209, 13.0232558,11.3953488, 8.8372093,
                                       6.97674419, 5.34883721, 3.48837209, 0};

    string path = "/mnt/d/Uni/Lectures/thesis/ALICE_multiplicity_data/" + input.cms_string + "TeV/V0M_dbar_";

    x_axis.clear();
    error.clear();
    data.clear();

    for (int j = 0; j < N_classes; ++j) {
        //read in alice data
        read_alice_data(path + class_string[j] + ".txt", x_axis, data, error);

        //set alice data in histograms
        alice[j].set_xaxis(x_axis);
        alice[j].set_data(data);
        alice[j].set_error(error);
        alice[j].not_scale();


        //set rescale factors for pythia data
        for (int k = 0; k < N_values; ++k) {
            pythia[j][k].rescale_data(rescale_factors[j]);
            //set x axis for pythia
            pythia[j][k].set_xaxis(x_axis);
            pythia[j][k].normalize(); //differential in pT
        }
        x_axis.clear();
        error.clear();
        data.clear();
    }

    mySimulation S;

    //S.set_ALICE_weights_lt(input.cms_string, input.particle_type);
    S.set_ALICE_weights_multi(input.cms_string, input.particle_type, class_string_weights, 10);

    //analyze events

    //read in data from pythia and set weights for the events
    for (int l = 0; l < input.N_simulations; ++l) {
        S.load_txt(input.dataset_folder + input.files[l], input.N_events[l], true);
        //S.rescale_spectrum();
        S.rescale_spectrum_multi(class_bin_limits_weights);

        //loop over cutoff values
        for (int i = 0; i < N_values; ++i) {
            S.set_cutoff_momentum(cutoff_momentum[i]);
            S.coalescence();
            //fill data in Histograms
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    for (int j = 0; j < N_classes; ++j) {
                        if (dbar.Nch_eta() > class_bin_limits[j]) { continue; }
                        if (dbar.Nch_eta() < class_bin_limits[j + 1]) { continue; }
                        pythia[j][i].fill(dbar.pT(), dbar.wLT());
                        break;
                    }
                }
            }
        }
    }

    //rescale
    for (int m = 0; m < N_classes; ++m) {
        for (int i = 0; i < N_values; ++i) {
            pythia[m][i].rescale();
            chi_squared[m][i] = pythia[m][i].chi_squared_raw_norm_err(alice[m], input.norm_err_above,
                                                                      input.norm_err_below);
        }
    }

    //print chi squared values to plot in python
    print_results_txt(chi_squared, cutoff_momentum, N_values, class_string, input.particle, input.cms_string);
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
                  const string &particle_type, const string &cms_string) {

    string path;
    int N_class = 9;
    for (int i = 0; i < N_class; ++i) {

        path = "/mnt/d/Uni/Lectures/thesis/ALICE_results/multiplicity/chi_squared/chi_squared_" + class_string[i] +
               "_" + cms_string + "_" + particle_type + ".txt";

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