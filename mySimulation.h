//
// Created by rwagn on 07.03.2019.
//

#ifndef MYSIMULATION_MYSIMULATION_H
#define MYSIMULATION_MYSIMULATION_H

#include "myParticle.h"
#include "histogram.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

struct weights {
    vector<double> weights_;
    vector<double> weight_boundaries_;
    vector<double> weight_centre_;
    bool rescaled = false;
};

class mySimulation {

private:

    int n_id_ = 2112; // neutron
    int p_id_ = 2212; // proton
    vector<int> id_array_ = {p_id_, n_id_};
    int N_events_ = 0;
    double cutoff_momentum_2_; //coalescence cutoff for momentum
    double cutoff_position_2_; //coalescence cutoff for position

    weights weights_A;
    weights weights_LT;

    bool cascade_momentum(myParticle &neutron, myParticle &proton);

    static bool cascade_position(myParticle &neutron, myParticle &proton);

    void form_deuterons(vector<myParticle> Protons, vector<myParticle> Neutrons, vector<myParticle> &Deuterons);

    void load_weights(const string &file, weights &weight_struct);

public:
    vector<myParticle> **Particles; //Particle[Event][Type][i-th particle of this event and this typ]
    //Type: 0 proton, 1 neutron
    vector<myParticle> deuteron;


    mySimulation();

    void load_txt(const string &file, int N_Events);

    void load_txt(const string &file, int N_Events, int particle_id);

    void set_cutoff_momentum(double cutoff);

    void set_cutoff_position(double cutoff);

    void coalescence();

    void set_ALICE_weights_discrete(const string &cms, bool particle_type);

    void set_ALICE_weights_lt(const string &cms, bool particle_type);

    void rescale_spectrum();
};

#endif //MY&SIMULATION_MYSIMULATION_H
