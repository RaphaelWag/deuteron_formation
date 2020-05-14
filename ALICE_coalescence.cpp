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
    histogram H;
    histogram H_data;

    double cms = 13;
    Input input;
    input.reduce_data();
    input.set_particle("dbar", false);

    int N_values = 50;
    double cutoff_momentum[N_values];
    double chi_squared[N_values];
    double chi_squared_wLT[N_values];
    histogram H_w[N_values];
    histogram H_nw[N_values];
    for (int n = 0; n < N_values; ++n) {
        cutoff_momentum[n] = (185. + n) / 1000.;
    }

    vector<double> class_bin_limits_weights = {23255, 23.255814, 17.9069767, 14.8837209, 13.0232558,11.3953488, 8.8372093,
                                               6.97674419, 5.34883721, 3.48837209, 0};
    string class_string_weights[10] = {"i", "ii", "iii", "iv","v", "vi", "vii", "iix", "ix", "x"};

    input.set_cms(cms);


    x_axis.clear();
    error.clear();
    data.clear();

    cout << input.cms_string << endl;

    read_alice_data(input.ALICE_data, x_axis, data, error);

    H.set_xaxis(x_axis);
    //H.rescale_data(input.ff*input.nsdtoinel / (2 * M_PI * input.N_event_total)); //Overall normalization
    //13 TeV
    H.rescale_data(input.ff*input.nsdtoinel / (input.N_event_total)); //Overall normalization
    //need inverse x axis values for rescaling
    for (auto &values:x_axis) {
        values = 1. / values;
    }

    //H.rescale_data(x_axis); //bin individual rescale factors
    H.normalize(); //Normalize per bin width
    H_data = H; //normalize and set x_axis equal to all histograms

    for (auto &element:H_w) { element = H; }
    for (auto &element:H_nw) { element = H; }

    //set data from ALICE
    H_data.not_scale();
    H_data.set_data(data);
    H_data.set_error(error);

    mySimulation S;
    S.set_ALICE_weights_lt(input.cms_string, input.particle_type);
    //S.set_ALICE_weights_multi(input.cms_string, input.particle_type, class_string_weights, 10);

    //analyze events

    for (int j = 0; j < input.N_simulations; ++j) {
        //read in data from pythia and set weights for the events
        S.load_txt(input.dataset_folder + input.files[j], input.N_events[j], true);
        S.rescale_spectrum();
        //S.rescale_spectrum_multi(class_bin_limits_weights);

        //loop over cutoff values
        for (int i = 0; i < N_values; ++i) {
            S.set_cutoff_momentum(cutoff_momentum[i]);
            S.coalescence();

            //fill data in Histograms
            for (auto &dbar:S.deuteron) {
                if (abs(dbar.y()) <= 0.5) {
                    H_w[i].fill(dbar.pT(), dbar.wLT());
                    H_nw[i].fill(dbar.pT());

                }
            }
        }
    }
    //rescale
    for (auto &element:H_w) { element.rescale(); }
    for (auto &element:H_nw) { element.rescale(); }


    //calculate chi squared and safe in array
    for (int i = 0; i < N_values; ++i) {
        chi_squared[i] = H_nw[i].chi_squared_raw_norm_err(H_data,input.norm_err_above, input.norm_err_below);
        chi_squared_wLT[i] = H_w[i].chi_squared_raw_norm_err(H_data,input.norm_err_above, input.norm_err_below);

        cout << cutoff_momentum[i] << " " << chi_squared[i] / (H_data.N_bins() - 1.) << " "
             << chi_squared_wLT[i] / (H_data.N_bins() - 1.) << endl;
    }




//print chi squared values to plot in python
    print_results_txt(chi_squared, cutoff_momentum, N_values, input.cms_string, "nw"); //nw = no weights
    print_results_txt(chi_squared_wLT, cutoff_momentum, N_values, input.cms_string, "cw"); //dw = discrete weights


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
