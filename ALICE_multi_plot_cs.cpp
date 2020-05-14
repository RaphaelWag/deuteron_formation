//
// Created by rwagn on 05.05.2020.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
#include <vector>
#include <cmath>
#include "ALICE_input.h"


using namespace std;

void read_alice_data(const string &file, vector<double> &x_axis, vector<double> &data, vector<double> &error);

int main() {


    const int N_classes = 9;
    vector<double> x_axis;
    vector<double> data;
    vector<double> error;

    double cms = 13;
    Input input;
    //input.cs_data();
    input.reduce_data();
    input.set_particle("dbar", false);
    input.set_cms(cms);

    //vector<double> sigma_0_inv = {2.2, 2.55, 2.6, 3.75, 4.8};
    //vector<double> p_0 = {184,190,188,207,209};
    vector<double> sigma_0_inv = {2.3, 3.1, 3.1, 3.4, 4.4, 3.7, 5.7, 3.9, 6.5};
    vector<double> p_0 = {190, 206, 202, 206, 220, 202, 229, 198, 221};


    int mc_cycles = 10;

    vector<histogram> pythia(N_classes);

    //string class_string[N_classes] = {"i+ii", "iii", "iv+v", "vi+vii", "iix+ix+x"};
    //double cs_frac[N_classes] = {0.047, 0.048, 0.095, 0.19, 0.62};
    //double class_bin_limits[N_classes + 1] = {23255, 14.4186047, 12.0930233, 9.30232558, 5.81395349, 0};

    double cs_frac[N_classes] = {0.0092, 0.0368, 0.046, 0.092, 0.092, 0.092, 0.092, 0.185, 0.355};
    string class_string[N_classes] = {"i", "ii", "iii", "iv+v", "vi", "vii", "iix", "ix", "x"};
    string class_string_weights[10] = {"i", "ii", "iii", "iv", "v", "vi", "vii", "iix", "ix", "x"};
    vector<double> class_bin_limits = {23255, 23.255814, 17.9069767, 14.8837209, 11.3953488, 8.8372093,
                                       6.97674419, 5.34883721, 3.48837209, 0};
    vector<double> class_bin_limits_weights = {23255, 23.255814, 17.9069767, 14.8837209, 13.0232558, 11.3953488,
                                               8.8372093,
                                               6.97674419, 5.34883721, 3.48837209, 0};


    double rescale_factors[N_classes];
    for (int i = 0; i < N_classes; ++i) {
        rescale_factors[i] = input.ff * input.nsdtoinel / (input.N_event_total * cs_frac[i] * 0.75 * mc_cycles);
    }

    string path = "/mnt/d/Uni/Lectures/thesis/ALICE_multiplicity_data/" + input.cms_string + "TeV/V0M_dbar_";

    x_axis.clear();
    error.clear();
    data.clear();


    for (int j = 0; j < N_classes; ++j) {
        //read in alice data
        read_alice_data(path + class_string[j] + ".txt", x_axis, data, error);

        pythia[j].rescale_data(rescale_factors[j]);
        pythia[j].set_xaxis(x_axis);
        pythia[j].normalize(); //differential in pT

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
        for (int i = 0; i < N_classes; ++i) {
            //loop over cutoff values
            for (int k = 0; k < mc_cycles; ++k) {

                //S.set_sigma0(sigma_0_inv[i] * 1e-6);
                //S.cs_model_formation();

                S.set_cutoff_momentum(p_0[i] * 1e-3);
                S.coalescence();

                //fill data in Histograms
                for (auto &dbar:S.deuteron) {
                    if (abs(dbar.y()) <= 0.5) {
                        if (dbar.Nch_eta() > class_bin_limits[i]) { continue; }
                        if (dbar.Nch_eta() < class_bin_limits[i + 1]) { continue; }
                        pythia[i].fill(dbar.pT(), dbar.wLT());
                    }
                }
            }

        }
    }


    for (int m = 0; m < N_classes; ++m) {
        //rescale
        pythia[m].rescale();
        //print histograms to plot in python
        //pythia[m].print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/cs_model/dbar_data_" + input.cms_string +
         //             "_" + class_string[m] + ".txt");
        pythia[m].print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/dbar_data_" + input.cms_string +
                        "_" + class_string[m] + ".txt");
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

