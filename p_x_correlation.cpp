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

bool coming_closer(myParticle &proton, myParticle &neutron);

int main() {

    double cms = 7;

    string particle = "dbar";
    Input input;
    input.set_particle(particle, false);
    input.reduce_data();
    input.set_cms(cms);
    vector<complex<double>> data; // real part = dx, imaginary part = dp

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
        for (int k = 0; k < S[l].N_events(); ++k) {
            for (auto &proton:S[l].Particles[k][0]) {
                for (auto &neutron:S[l].Particles[k][1]) {
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

                        //save dx and dp

                    } else {
                        //save dx and dp

                    }
                }
            }
        }


    }
    //write results to plot in python

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