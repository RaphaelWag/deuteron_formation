//
// Created by rwagn on 06.08.2019.
//

#include <iostream>
#include "mySimulation.h"
#include <vector>
#include <complex>
#include "ALICE_input.h"

using namespace std;

struct com_data {
    vector<double> n_x = {0, 0, 0};
    vector<double> p_x = {0, 0, 0};
    vector<double> n_p = {0, 0, 0};
    vector<double> p_p = {0, 0, 0};
    double dx=0;
    double dp=0;
};

void print_data(const string &file, vector<com_data> &Data);



int main() {

    double cms = 7;

    Input input;
    input.full_data();
    input.set_cms(cms);
    input.max_simulations(1);
    int pbar_id = -2212;
    int nbar_id = -2112;
    vector<com_data> Data;
    int max_events = 10000;
    int used_events = 0;

    auto *S = new mySimulation[2];

    //read in data from pythia

    S[0].load_txt(input.dataset_folder + input.files[0], input.N_events[0], pbar_id);
    S[1].load_txt(input.dataset_folder + input.files[0], input.N_events[0], nbar_id);

    // loop over all particles
    for (int i = 0; i < input.N_events[0]; ++i) {
        for (auto &pbar:S[0].Particles[i][0]) {
            for (auto &nbar:S[1].Particles[i][0]) {
                if(used_events==max_events){ continue;}
                myParticle neutron = nbar;
                myParticle proton = pbar;
                myParticle boost_particle = nbar + pbar;
                // boost each pair in com frame
                neutron.bstback(boost_particle);
                proton.bstback(boost_particle);


                //check if it works
                double dpx = pow(neutron.p(1) - proton.p(1), 2);
                double dpy = pow(neutron.p(2) - proton.p(2), 2);
                double dpz = pow(neutron.p(3) - proton.p(3), 2);

                double dx = pow(neutron.x(1) - proton.x(1), 2);
                double dy = pow(neutron.x(2) - proton.x(2), 2);
                double dz = pow(neutron.x(3) - proton.x(3), 2);

                dx = sqrt(dx+dy+dz);
                double dp = sqrt(dpx+dpy+dpz);

                //get data
                com_data temp;
                temp.n_x = {neutron.x(1),neutron.x(2),neutron.x(3)};
                temp.n_p = {neutron.p(1),neutron.p(2),neutron.p(3)};
                temp.p_x = {proton.x(1),proton.x(2),proton.x(3)};
                temp.p_p = {proton.p(1),proton.p(2),proton.p(3)};
                temp.dx = dx;
                temp.dp = dp;

                Data.push_back(temp);
                ++used_events;
            }
        }
    }


    //compute correlation factors in python

    //print data in txt file to plot in python
    print_data("/mnt/d/Uni/Lectures/thesis/ALICE_python/p_x_corr_data/"+input.cms_string+".txt",Data);

    return 0;
}

void print_data(const string &file, vector<com_data> &Data) {
    ofstream myfile(file);
    if (!myfile.is_open()) {
        cout << "Error in: print_results" << endl;
        cout << "Unable to open file" << endl;
    } else {
        for (auto &element:Data) {
            myfile << element.n_x[0] << " " <<element.n_x[1] << " " << element.n_x[2] << " ";
            myfile << element.n_p[0] << " " <<element.n_p[1] << " " << element.n_p[2] << " ";
            myfile << element.p_x[0] << " " <<element.p_x[1] << " " << element.p_x[2] << " ";
            myfile << element.p_p[0] << " " <<element.p_p[1] << " " << element.p_p[2] << " ";
            myfile << element.dx << " " << element.dp ;
            myfile << endl;
        }
        myfile.close();
    }
}