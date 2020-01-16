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
    histogram H_rescaled_A;
    histogram H_rescaled_LT;
    histogram H_nscaled;
    histogram H_data;

    vector<double> cms_list = {7};
    Input input;
    input.reduce_data();
    input.set_particle("dbar", false);

    int N_values = 65;
    double cutoff_momentum[N_values];
    double chi_squared[N_values];
    double chi_squared_wA[N_values];
    double chi_squared_wLT[N_values];
    for (int n = 0; n < N_values; ++n) {
        cutoff_momentum[n] = (190. + n) / 1000.;
    }
    double N_deuterons = 0;
    double N_LT = 0;
    double N_A = 0;

    for (auto &cms:cms_list) {
        input.set_cms(cms);
        x_axis.clear();
        error.clear();
        data.clear();

        cout << input.cms_string << endl;

        read_alice_data(input.ALICE_data, x_axis, data, error);

        H_rescaled_A.set_xaxis(x_axis);
        H_rescaled_A.rescale_data(input.ff / (2 * M_PI * input.N_event_total)); //Overall normalization
        //need inverse x axis values for rescaling
        for (auto &values:x_axis) {
            values = 1. / values;
        }

        H_rescaled_A.rescale_data(x_axis); //bin individual rescale factors
        H_rescaled_A.normalize(); //Normalize per bin width
        H_nscaled = H_rescaled_LT = H_data = H_rescaled_A; //normalize and set x_axis equal to all histograms

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
                        H_rescaled_A.fill(dbar.pT(), dbar.wA());
                        H_rescaled_LT.fill(dbar.pT(), dbar.wLT());
                        H_nscaled.fill(dbar.pT());
                        N_deuterons++;
                        if (dbar.wLT() == 1.) { N_LT++; }
                        if (dbar.wA() == 1.) { N_A++; }
                    }
                }
            }
            //rescale
            H_rescaled_A.rescale();
            H_rescaled_LT.rescale();
            H_nscaled.rescale();


            //calculate chi squared and safe in array
            chi_squared[i] = H_nscaled.chi_squared_raw(H_data);
            chi_squared_wA[i] = H_rescaled_A.chi_squared_raw(H_data);
            chi_squared_wLT[i] = H_rescaled_LT.chi_squared_raw(H_data);

            //reset Histograms
            H_nscaled.reset();
            H_rescaled_A.reset();
            H_rescaled_LT.reset();

            cout << cutoff_momentum[i] << " " << chi_squared[i] / (H_data.N_bins() - 1.) << " "
                 << chi_squared_wA[i] / (H_data.N_bins() - 1.) << " " << chi_squared_wLT[i] / (H_data.N_bins() - 1.)
                 << " " << N_A / N_deuterons << " " << N_LT / N_deuterons << endl;

            N_deuterons = N_A = N_LT = 0;
            //end loop
        }

        //print chi squared values to plot in python
        print_results_txt(chi_squared, cutoff_momentum, N_values, input.cms_string, "nw"); //nw = no weights
        print_results_txt(chi_squared_wA, cutoff_momentum, N_values, input.cms_string, "dw"); //dw = discrete weights
        print_results_txt(chi_squared_wLT, cutoff_momentum, N_values, input.cms_string, "cw"); //dw = discrete weights

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
    path = "/mnt/d/Uni/Lectures/thesis/ALICE_results/chi_squared/chi_squared_" + cms_string + "_" + wtype + ".txt";

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
