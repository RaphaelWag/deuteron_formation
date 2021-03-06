//
// Created by rwagn on 12.11.2019.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
#include "ALICE_input.h"
#include <vector>
#include <cmath>

using namespace std;

double levy_tsallis(double pT, vector<double> params);

double levy_tsallis_new(double pT, vector<double> params);

int main() {

    vector<double> x_axis = {0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13,
                             0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22,
                             0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.325,
                             0.350, 0.375, 0.400, 0.425, 0.475, 0.525, 0.575, 0.625,
                             0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025,
                             1.075, 1.125, 1.175, 1.225, 1.275, 1.35, 1.45, 1.55,
                             1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45,
                             2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25, 3.35,
                             3.45, 3.55, 3.65, 3.75, 3.85, 3.95};
    double cms = 13;
    vector<double> params = {0.3691/2., 7.23, 0.2442, 0.9382720813}; // N,n,C,m0

    // anti protons
    // 7TeV
    // [0.16840874, 6.64441328, 0.22235897]
    // [0.16840874, 6.64441328, 0.22235897]
    // 2.76 TeV
    // [0.1352278,  6.81995388, 0.19849261]
    // [0.1352278,  6.81995388, 0.19849261]
    // 0.9 TeV
    // [0.09917695, 8.44442985, 0.19216617]
    // [0.09917695, 8.44442985, 0.19216617]

    //13 TeV different norm
    // [0.3691/2., 7.23, 0.2442]
    // [0.3691/2., 7.23, 0.2442]


    // new 0.9 TeV data with different normalization
    // [0.07854844 8.17249344 0.18937513]
    // [0.07854844 8.17249344 0.18937513]

    Input input;
    input.full_data();
    input.set_cms(cms);

    int particle_id = -2212;
    string particle_string = "pbar";

    histogram H_pbar;

    H_pbar.set_xaxis(x_axis);
    //H_pbar.rescale_data(input.ff / (2 * M_PI * input.N_event_total)); //Overall normalization to NSD
    H_pbar.rescale_data(input.ff*input.nsdtoinel / (input.N_event_total*0.75)); //Overall normalization 13TeV data

    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    //H_pbar.rescale_data(x_axis); //bin individual rescale factors, not used in 13TeV data

    H_pbar.normalize(); //Normalize per bin width

    mySimulation S;

    //analyze events

    //read in data from pythia
    for (int l = 0; l < input.N_simulations; ++l) {
        S.load_txt(input.dataset_folder + input.files[l], input.N_events[l], particle_id);

        //fill Histogram
        for (auto &event:S.Event) {
            for (auto &pbar:event.protons) {
                if (abs(pbar.y()) <= 0.5) {
                    H_pbar.fill(pbar.pT());
                }
            }
        }
    }


    //print histogram
    H_pbar.rescale();

    H_pbar.print("/mnt/d/Uni/Lectures/thesis/ALICE_data_mc/" +
                 particle_string + "_data_" + input.cms_string + ".txt");

    //write weights in file
    string path = "/mnt/d/Uni/Lectures/thesis/ALICE_weighting/levy_tsallis_weights/";
    string file = "weights_" + particle_string + "_" + input.cms_string + ".txt";
    string filename = path + file;

    ofstream myfile(filename);
    if (!myfile.is_open()) {
        cout << "write weightss in file" << endl;
        cout << "Unable to open file" << endl;
    } else {
        myfile << x_axis.size() << endl;
        for (int i = 0; i < H_pbar.N_bins(); ++i) {
            myfile << H_pbar.x_axis_[i] << " " << levy_tsallis_new(H_pbar.x_axis_[i], params) / H_pbar.data_[i] << endl;
        }
        myfile.close();
    }
    return 0;
}


double levy_tsallis(double pT, vector<double> params) {
    double N = params[0];
    double n = params[1];
    double C = params[2];
    double m0 = params[3];
    double mT = sqrt(pT * pT + m0 * m0);
    double nC = n * C;

    double term1 = (n - 1.) * (n - 2.);
    double term2 = 2. * M_PI * nC * (nC + m0 * (n - 2.));
    double term3 = 1. + (mT - m0) / nC;
    return N * term1 / term2 * pow(term3, -n);
}

double levy_tsallis_new(double pT, vector<double> params) {
    double N = params[0];
    double n = params[1];
    double C = params[2];
    double m0 = params[3];
    double mT = sqrt(pT * pT + m0 * m0);
    double nC = n * C;

    double term1 = pT * (n - 1.) * (n - 2.);
    double term2 = nC * (nC + m0 * (n - 2.));
    double term3 = 1. + (mT - m0) / nC;
    return N * term1 / term2 * pow(term3, -n);
}