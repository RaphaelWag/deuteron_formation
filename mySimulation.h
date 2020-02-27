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

struct temp_d {
    int type; // same as probability functions
    int n1; // number of first nucleus in array
    int n2; // number of second nucleus in array
    double probability;
};

struct event {
    vector<myParticle> neutrons;
    vector<myParticle> protons;
};

class mySimulation {

private:

    int n_id_ = 2112; // neutron
    int p_id_ = 2212; // proton
    vector<int> id_array_ = {p_id_, n_id_};
    vector<temp_d> temporary_deuteron;
    int N_events_ = 0;
    double cutoff_momentum_2_; //coalescence cutoff for momentum
    double cutoff_position_2_; //coalescence cutoff for position
    double sigma_0_inv_; //normalization for cross section model

    weights weights_A;
    weights weights_LT;

    bool cascade_momentum(myParticle &neutron, myParticle &proton);

    bool cascade_position(myParticle &neutron, myParticle &proton);

    void form_deuterons(vector<myParticle> Protons, vector<myParticle> Neutrons);

    void form_deuterons_cs(vector<myParticle> Protons, vector<myParticle> Neutrons);

    void load_weights(const string &file, weights &weight_struct);

    double prob_1(double k); //pn-> d gamma

    double prob_2(double q, double inverse_pi_mass); // pn -> d pi0

    double prob_3(double k); // pn -> d pi0 pi0

    double prob_4(double k); // pn -> d pi+ pi-

    double prob_5(double q, double inverse_pi_mass); // pp -> d pi-

    double prob_6(double k); // pp -> d pi- pi0

    double prob_7(double q, double inverse_pi_mass); // nn -> d pi+

    double prob_8(double k); // nn -> pi+ pi0

    double prob(double k, double a, double b, double c, double d, double e);

public:
    //vector<myParticle> **Particles; //Particle[Event][Type][i-th particle of this event and this typ]
    //Type: 0 proton, 1 neutron
    vector<event> Event;
    vector<myParticle> deuteron;

    mySimulation();

    int N_events();

    void load_txt(const string &file, int N_Events, bool reduced_data = false);

    void load_txt(const string &file, int N_Events, int particle_id);

    void set_cutoff_momentum(double cutoff);

    void set_cutoff_position(double cutoff);

    void set_sigma0(double sigma0);

    void coalescence();

    void cs_model_formation();

    void set_ALICE_weights_discrete(const string &cms, bool particle_type);

    void set_ALICE_weights_lt(const string &cms, bool particle_type);

    void rescale_spectrum();
};

#endif //MY&SIMULATION_MYSIMULATION_H
