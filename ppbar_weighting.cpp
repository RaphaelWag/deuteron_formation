//
// Created by rwagn on 03.05.2020.
//

#include <iostream>
#include "histogram.h"
#include "mySimulation.h"
#include "ALICE_input.h"
#include <vector>
#include <cmath>

using namespace std;

double levy_tsallis(double pT, vector<double> params);

int main() {

    vector<double> x_axis = {0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13,
                             0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22,
                             0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.325,
                             0.350, 0.375, 0.400, 0.425, 0.475, 0.525, 0.575, 0.625,
                             0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025,
                             1.075, 1.125, 1.175, 1.225, 1.275, 1.35, 1.45, 1.55,
                             1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45,
                             2.55, 2.65, 2.75, 2.85, 2.95};
    double cms = 13;
    vector<double> params_i = {1.358, 8.62, 0.3454, 0.9382720813}; // N,n,C,m0
    vector<double> params_ii = {1.061, 8.35, 0.3146, 0.9382720813}; // N,n,C,m0
    vector<double> params_iii = {0.871, 8.01, 0.2900, 0.9382720813}; // N,n,C,m0
    vector<double> params_iv = {0.750, 7.84, 0.2732, 0.9382720813}; // N,n,C,m0
    vector<double> params_v = {0.660, 7.70, 0.2600, 0.9382720813}; // N,n,C,m0
    vector<double> params_vi = {0.555, 7.49, 0.2428, 0.9382720813}; // N,n,C,m0
    vector<double> params_vii = {0.443, 7.36, 0.2239, 0.9382720813}; // N,n,C,m0
    vector<double> params_viii = {0.3539, 7.17,0.2054, 0.9382720813}; // N,n,C,m0
    vector<double> params_ix = {0.2484, 7.016, 0.1837, 0.9382720813}; // N,n,C,m0
    vector<double> params_x = {0.1296, 6.81, 0.1437, 0.9382720813}; // N,n,C,m0

    vector<vector<double>> params = {params_i, params_ii, params_iii, params_iv, params_v, params_vi, params_vii,
                                     params_viii, params_ix, params_x};

    Input input;
    input.full_data();
    input.set_cms(cms);

    int particle_id = -2212;
    string particle_string = "pbar";

    string class_string[10] = {"i", "ii", "iii", "iv", "v", "vi", "vii", "iix", "ix", "x"};
    double cs_frac[10] = {0.0092, 0.0368, 0.046, 0.046, 0.046, 0.092, 0.092, 0.092, 0.185, 0.355};
    double rescale_factors[10];

    for (int i = 0; i < 10; ++i) {
        rescale_factors[i] = input.ff * input.nsdtoinel / (input.N_event_total * cs_frac[i] * 0.75);
    }

    double class_bin_limits[11] = {23255, 23.255814, 17.9069767, 14.8837209, 13.0232558, 11.3953488, 8.8372093,
                                   6.97674419, 5.34883721, 3.48837209, 0};

    auto H_pbar = new histogram[10];

    for (int j = 0; j < 10; ++j) {
        H_pbar[j].set_xaxis(x_axis);
        H_pbar[j].rescale_data(rescale_factors[j]); //Overall normalization 13TeV data INEL>0
        H_pbar[j].normalize(); //Normalize per bin width
    }
    mySimulation S;
    //analyze events

    //read in data from pythia
    for (int l = 0; l < input.N_simulations; ++l) {
        S.load_txt(input.dataset_folder + input.files[l], input.N_events[l], particle_id);

        //fill Histogram
        for (auto &event:S.Event) {
            for (auto &pbar:event.protons) {
                if (abs(pbar.y()) <= 0.5) {
                    for (int j = 0; j < 10; ++j) {
                        if (pbar.Nch_eta() > class_bin_limits[j]) { continue; }
                        if (pbar.Nch_eta() < class_bin_limits[j + 1]) { continue; }
                        H_pbar[j].fill(pbar.pT(), pbar.wLT());
                        break;
                    }
                }
            }
        }
    }


    //print histogram
    for (int k = 0; k < 10; ++k) {
        H_pbar[k].rescale();
        H_pbar[k].print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/" +
                        particle_string + "_data_" + input.cms_string + "_" + class_string[k] + ".txt");

        //write weights in file
        string path = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/levy_tsallis_weights/";
        string file = "weights_" + particle_string + "_" + input.cms_string + "_" + class_string[k] + ".txt";
        string filename = path + file;

        ofstream myfile(filename);
        if (!myfile.is_open()) {
            cout << "write weightss in file" << endl;
            cout << "Unable to open file" << endl;
        } else {
            myfile << x_axis.size() << endl;
            for (int i = 0; i < H_pbar[k].N_bins(); ++i) {
                myfile << H_pbar[k].x_axis_[i] << " "
                       << levy_tsallis(H_pbar[k].x_axis_[i], params[k]) / H_pbar[k].data_[i] << endl;
            }
            myfile.close();
        }
    }
    return 0;
}


double levy_tsallis(double pT, vector<double> params) {
    double N = params[0]/2.;
    double n = params[1];
    double C = params[2];
    double m0 = params[3];
    double mT = sqrt(pT * pT + m0 * m0);
    double nC = n * C;

    double term1 = pT*(n - 1.) * (n - 2.);
    double term2 =  nC * (nC + m0 * (n - 2.));
    double term3 = 1. + (mT - m0) / nC;
    return N * term1 / term2 * pow(term3, -n);
}
