//
// Created by rwagn on 07.11.2019.
//

#include <iostream>
#include "mySimulation.h"
#include "histogram.h"
#include <complex>
#include <vector>
#include <cmath>
#include "ALICE_input.h"


using namespace std;

struct Data {
    double dx;
    double dp;
    double w;
};

bool coming_closer(myParticle &proton, myParticle &neutron);

void print_results_txt(const vector<Data> &data, const string &cms_string);

int main() {

    double cms = 7;

    string particle = "dbar";
    Input input;
    input.set_particle(particle, false);
    input.reduce_data();
    input.set_cms(cms);
    input.max_simulations(17);
    vector<Data> data;
    histogram H;
    vector<double> x_axis;
    double min =0;
    double max = 2e-11;
    double bins = 200;
    for (int i = 0; i < bins; ++i) {
        x_axis.push_back(min+i/bins*max);
    }
    H.set_boundaries(x_axis);
    H.probability();

    cout << input.cms_string << endl;

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

        //loop over pairs in each event and get dx dp
        for (auto& event:S[l].Event) {
            Data temp{};
            double dp;
            double dx;
            for (auto &proton:event.protons) {
                //if(data.size()>100000){break;}
                for (auto &neutron:event.neutrons) {
                    //if(data.size()>100000){break;}
                    //check which particle was created first
                    double t = neutron.x(0) - proton.x(0);
                    if (t > 0) {
                        proton.move(t);
                    } else {
                        neutron.move(-t);
                    }
                    //boost to CoM frame
                    myParticle boostparticle = neutron + proton;
                    neutron.bstback(boostparticle);
                    proton.bstback(boostparticle);
                    // in CoM frame get closest distance
                    //check if they are coming closer or not
                    if (coming_closer(proton, neutron)) {
                        //get closest distance
                        // cross product of direction vectors (= momentum vectors)
                        double nx = proton.p(2) * neutron.p(3) - proton.p(3) * neutron.p(2);
                        double ny = proton.p(3) * neutron.p(1) - proton.p(1) * neutron.p(3);
                        double nz = proton.p(1) * neutron.p(2) - proton.p(2) * neutron.p(1);

                        // r1 - r2 = r
                        double rx = proton.x(1) - neutron.x(1);
                        double ry = proton.x(2) - neutron.x(2);
                        double rz = proton.x(3) - neutron.x(3);

                        // d =n*r/n
                        double d = abs((nx * rx + ny * ry + nz * rz) / sqrt(nx * nx + ny * ny + nz * nz));

                        double kx = proton.p(1) - neutron.p(1);
                        double ky = proton.p(2) - neutron.p(2);
                        double kz = proton.p(3) - neutron.p(3);

                        dp = sqrt(kx * kx + ky * ky + kz * kz);


                        //save dx and dp
                        temp.dx = d;
                        temp.dp = dp;
                        temp.w = neutron.wLT() * proton.wLT();
                        //data.push_back(temp);
                        H.fill(temp.dx,temp.w);

                    } else {
                        //save dx and dp
                        // r1 - r2 = r
                        double rx = proton.x(1) - neutron.x(1);
                        double ry = proton.x(2) - neutron.x(2);
                        double rz = proton.x(3) - neutron.x(3);

                        double kx = proton.p(1) - neutron.p(1);
                        double ky = proton.p(2) - neutron.p(2);
                        double kz = proton.p(3) - neutron.p(3);
                        dp = sqrt(kx * kx + ky * ky + kz * kz);
                        dx = sqrt(rx * rx + ry * ry + rz * rz);
                        //save dx and dp
                        temp.dx = dx;
                        temp.dp = dp;
                        temp.w = neutron.wLT() * proton.wLT();
                        //data.push_back(temp);
                        H.fill(temp.dx,temp.w);
                    }

                }

            }
        }


    }
    //write results to plot in python
    //print_results_txt(data, input.cms_string);
    H.print("/mnt/d/Uni/Lectures/thesis/ALICE_results/p_x_corr/dx_dist_"+input.cms_string+".txt");
    return 0;
}

bool coming_closer(myParticle &proton, myParticle &neutron) {
    double x, y, z;
    x = neutron.x(1) - proton.x(1);
    y = neutron.x(2) - proton.x(2);
    z = neutron.x(3) - proton.x(3);

    double sp = x * proton.p(1) + y * proton.p(2) + z * proton.p(3);

    return sp > 0;
}

void
print_results_txt(const vector<Data> &data, const string &cms_string) {

    string path;
    path = "/mnt/d/Uni/Lectures/thesis/ALICE_results/p_x_corr/p_x_corr_" + cms_string + ".txt";

    ofstream myfile(path);
    if (!myfile.is_open())
        cout << "Unable to open file";
    else {
        for (auto &element:data) {
            myfile << element.dx << " " << element.dp << " " << element.w << endl;
        }
    }
    myfile.close();
}
